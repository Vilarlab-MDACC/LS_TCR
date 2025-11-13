#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(furrr)

options(future.globals.maxSize = 3 * 1024^3)  # 3 GiB
plan(multisession, workers = availableCores() - 1)
library(progressr)
mem.maxVSize(vsize = 160000) # Set limit to 16 GB (example`)

fast_fisher <- function(training_set, training_set_meta, positive_cutoff) {
    data.cnt <- training_set[, .(n = .N), by = .(Name, Contrast)]
    n_Contrast_TRUE <- training_set_meta[Contrast == TRUE, .N]
    n_Contrast_FALSE <- training_set_meta[Contrast == FALSE, .N]



    # Complete missing combinations
    all_combinations <- CJ(
        Name = unique(data.cnt$Name),
        Contrast = unique(data.cnt$Contrast)
    )
    data.cnt <- merge(all_combinations, data.cnt, all.x = TRUE)
    data.cnt[is.na(n), n := 0]

    # Pivot wider using dcast
    data.cnt <- dcast(data.cnt, Name ~ Contrast, value.var = "n", fill = 0)
    setnames(data.cnt, c("TRUE", "FALSE"), c("Contrast_TRUE_Present_TRUE", "Contrast_FALSE_Present_TRUE"))




    # Add absent counts
    data.cnt[, `:=`(
        Contrast_TRUE_Present_FALSE = n_Contrast_TRUE - Contrast_TRUE_Present_TRUE,
        Contrast_FALSE_Present_FALSE = n_Contrast_FALSE - Contrast_FALSE_Present_TRUE
    )]

    data.cnt <- data.cnt[Contrast_TRUE_Present_TRUE / (Contrast_TRUE_Present_TRUE + Contrast_TRUE_Present_FALSE) > positive_cutoff &
        Contrast_TRUE_Present_TRUE / (Contrast_TRUE_Present_TRUE + Contrast_TRUE_Present_FALSE) > Contrast_FALSE_Present_TRUE / (Contrast_FALSE_Present_TRUE + Contrast_FALSE_Present_FALSE)]
    if (nrow(data.cnt) == 0) {
        return(data.cnt)
    }

    fisher.test.result <- unique(data.cnt[, .(
        Contrast_TRUE_Present_TRUE, Contrast_TRUE_Present_FALSE,
        Contrast_FALSE_Present_TRUE, Contrast_FALSE_Present_FALSE
    )])

    fisher.test.result[, fisher := {
        fisher.test(matrix(c(
            Contrast_TRUE_Present_TRUE, Contrast_TRUE_Present_FALSE,
            Contrast_FALSE_Present_TRUE, Contrast_FALSE_Present_FALSE
        ), nrow = 2))$p.value
    }, by = 1:nrow(fisher.test.result)]


    fisher.result <- merge(data.cnt, fisher.test.result,
        by = c(
            "Contrast_TRUE_Present_TRUE", "Contrast_TRUE_Present_FALSE",
            "Contrast_FALSE_Present_TRUE", "Contrast_FALSE_Present_FALSE"
        )
    )

    return(fisher.result)
}


calculate_AUC <- function(TCR_sum) {
    TCR_sum <- TCR_sum %>%
        ungroup() %>%
        arrange(desc(metric))
    0:nrow(TCR_sum) %>%
        map_dfr(
            ~ TCR_sum |>
                mutate(
                    prediction =
                        case_when(
                            row_number() > .x ~ F,
                            T ~ T
                        )
                ) %>%
                summarise(
                    TPR = sum(Contrast &
                        prediction) / (sum(Contrast & prediction) + sum(Contrast & !prediction)),
                    FPR = sum(!Contrast &
                        prediction) / (sum(!Contrast & prediction) + sum(!Contrast & !prediction))
                )
        ) |>
        mutate(diff_FPR = c(0, diff(FPR))) %>%
        mutate(AUC = sum(TPR * diff_FPR))
}


wilcox_test <- function(training_set, training_set_meta, target_tcr) {
    all_combinations <- CJ(
        Name = target_tcr,
        Sample = training_set_meta$Sample
    )
    training_expanded <- merge(all_combinations, training_set[, .(Name, Sample, Proportion)], all.x = TRUE)
    training_expanded[is.na(Proportion), Proportion := 0]
    
    training_expanded <- merge(training_expanded, training_set_meta[, .(Sample, Contrast)], by = "Sample", all.x = TRUE)

    wrong_mean = training_expanded[,.(Proportion = mean(Proportion)), by = .(Name, Contrast)]%>%
      dcast(Name ~ Contrast, value.var = "Proportion", fill = 0)%>%
    `[`(`FALSE` > `TRUE`)
    
    
    wilcox_result <- training_expanded[!Name%in%wrong_mean][, .(wilcox = wilcox.test(Proportion ~ Contrast)$p.value), by = Name]
    return(wilcox_result)
}

args <- commandArgs(trailingOnly = TRUE)



calculater_logp <- function(training_set, validate_set, Signature) {
    training_set <- setDT(training_set)
    validate_set <- setDT(validate_set)

    validate_proportion <- asin(sqrt(sum(validate_set[Name %in% Signature]$Proportion)))

    training_pos_proportion <- asin(sqrt(training_set[Name %in% Signature & Contrast == TRUE][, .(Proportion = sum(Proportion)), by = .(Sample)]$Proportion))
    training_neg_proportion <- asin(sqrt(training_set[Name %in% Signature & Contrast == FALSE][, .(Proportion = sum(Proportion)), by = .(Sample)]$Proportion))
    # validate_proportion <- sum(validate_set[Name %in% Signature]$Proportion)
    # training_pos_proportion <- training_set[Name %in% Signature & Contrast == TRUE][, .(Proportion = sum(Proportion)), by = .(Sample)]$Proportion
    # training_neg_proportion <- training_set[Name %in% Signature & Contrast == FALSE][, .(Proportion = sum(Proportion)), by = .(Sample)]$Proportion

    pos <- dnorm(validate_proportion,
        mean = mean(training_pos_proportion),
        sd = sd(training_pos_proportion)
    )
    neg <- dnorm(validate_proportion,
        mean = mean(training_neg_proportion),
        sd = sd(training_neg_proportion)
    )
    return(log10(pos / neg))
}
# 0.2,0.001,0.1,"aa_j_v",0.172384510869565,NA,NA,66666,580,"c_vs_ps"

    positive_cutoff <- as.numeric(args[1])
    fisher_cutoff <- as.numeric(args[2])
    wilcox_cutoff <- as.numeric(args[3])
    target <- as.character(args[4])
    compare <- as.character(args[5])

print(positive_cutoff)
print(fisher_cutoff)
print(wilcox_cutoff)
print(target)
print(compare)



meta_data <- fread("meta.txt")


samples <- fread("./final_samples_for_classification.csv")[Sample%in%meta_data$Sample,]



print("read files done")
samples <- switch(target,
    aa_j = samples[, .(Name = paste(CDR3.aa, J.name, sep = "_"), Sample, Proportion, Status_MMRd_Cancer, Institution)],
    aa = samples[, .(Name = CDR3.aa, Sample, Proportion, Status_MMRd_Cancer, Institution)],
    aa_v = samples[, .(Name = paste(CDR3.aa, V.name, sep = "_"), Sample, Proportion, Status_MMRd_Cancer, Institution)],
    aa_j_v = samples[, .(Name = paste(CDR3.aa, J.name, V.name, sep = "_"), Sample, Proportion, Status_MMRd_Cancer, Institution)],
    stop("Invalid target value. Must be one of: aa_j, aa, aa_v, aa_j_v")
)

pop <- switch(compare,
    c_vs_ps = list(positive_pop = c("Previvor", "Survivor"), negative_pop = c("Control")),
    c_vs_p = list(negative_pop = c("Control"), positive_pop = c("Previvor")),
    p_vs_s = list(positive_pop = c("Survivor"), negative_pop = c("Previvor")),
    stop("Invalid compare value. Must be one of: c_vs_ps, c_vs_p, p_vs_s")
)

print("prepare samples done")


sample_status <- samples[, .(
    Status_MMRd_Cancer = data.table::first(Status_MMRd_Cancer),
    Institution = data.table::first(Institution)
), by = Sample][, Sample := factor(Sample, levels = meta_data$Sample)]

setorder(sample_status, Sample)

set.seed(1234)
sample_status[, validate := sample(c(rep(TRUE, ceiling(.N * 0.20)), rep(FALSE, .N - ceiling(.N * 0.2)))),
    by = .(Status_MMRd_Cancer, Institution)
]
print("selected sample done")


filename <- paste0("sample_split_", positive_cutoff, "_", fisher_cutoff, "_", wilcox_cutoff,  "_", target, "_", compare, ".csv")
write.csv(sample_status, filename)


samples <- samples[Status_MMRd_Cancer %in% c(pop$positive_pop, pop$negative_pop)][sample_status[, .(Sample, validate)], on = "Sample", nomatch = NULL][, .(Proportion = sum(Proportion, na.rm = TRUE)),
    by = .(Name, Sample, Status_MMRd_Cancer, validate)
]

cut_n <- max(2,nrow(sample_status[!validate & Status_MMRd_Cancer %in% pop$positive_pop]) * positive_cutoff)
samples <- samples[, if (.N >= cut_n) .SD, by = Name]


print("filtered samples done")


testing_meta <- sample_status[!validate & Status_MMRd_Cancer %in% c(pop$positive_pop, pop$negative_pop)][, Contrast := Status_MMRd_Cancer %in% pop$positive_pop]
testing_samples <- samples[Status_MMRd_Cancer %in% c(pop$positive_pop, pop$negative_pop)][, Contrast := Status_MMRd_Cancer %in% pop$positive_pop][validate == FALSE]



with_progress({
    p <- progressor(steps = nrow(testing_meta))
    loovc_HL <- future_map_dfr(seq_len(nrow(testing_meta)), function(i) {
        p(sprintf("Processing fold %d", i))

        training_fold <- testing_samples[Sample %in% testing_meta[-i, Sample]]
        validation_fold <- testing_samples[Sample %in% testing_meta[i, Sample]]
        training_meta_fold <- testing_meta[-i, ]
        # Run Fisher's exact test
        fisher_result <- fast_fisher(training_fold, training_meta_fold, positive_cutoff)
        if (nrow(fisher_result) == 0) {
            return(data.table(LH = NA, Sample = testing_meta[i, Sample]))
        }

        fisher_filtered <- fisher_result[fisher <= fisher_cutoff]$Name

        if (length(fisher_filtered) == 0) {
            return(data.table(LH = NA, Sample = testing_meta[i, Sample]))
        }
        # Run Wilcoxon test
        wilcox_result <- wilcox_test(training_fold, training_meta_fold, fisher_filtered)

        Signature <- wilcox_result[wilcox <= wilcox_cutoff]$Name

        if (length(Signature) == 0) {
            return(data.table(LH = NA, Sample = testing_meta[i, Sample]))
        }
        # Calculate AUC on validation fold
        LH <- calculater_logp(training_fold, validation_fold, Signature)
        return(data.table(LH = LH, Sample = testing_meta[i, Sample]))
    })
})


loocv_auc <- loovc_HL[testing_meta, on = "Sample", nomatch = NULL] %>%
    rename(metric = LH) %>%
    arrange(desc(metric)) %>%
    calculate_AUC()
filename <- paste0("loocv_AUC_", positive_cutoff, "_", fisher_cutoff, "_", wilcox_cutoff, "_", target, "_", compare, ".csv")
write.csv(loocv_auc, filename, row.names = FALSE)


print("loocv done")


fisher_result <- fast_fisher(testing_samples, testing_meta, positive_cutoff)



if (nrow(fisher_result) == 0) {
    validate_AUC <- NA
    test_AUC <- NA
    size <- NA
    Signature <- NULL
} else {
    fisher_filtered <- fisher_result[fisher <= fisher_cutoff]$Name

    if (length(fisher_filtered) == 0) {
        validate_AUC <- NA
        test_AUC <- NA
        size <- NA
        Signature <- NULL
    } else {
        wilcox_result <- wilcox_test(testing_samples, testing_meta, fisher_filtered)

        Signature <- wilcox_result[wilcox <= wilcox_cutoff]$Name

        if (length(Signature) == 0) {
            validate_AUC <- NA
            test_AUC <- NA
            size <- NA
            Signature <- NULL
        }
    }
}




filename <- paste0("Signature_", positive_cutoff, "_", fisher_cutoff, "_", wilcox_cutoff, "_",  target, "_", compare, ".RDS")
saveRDS(Signature, filename)



print("signature saved")

if (length(Signature) == 0) {
    results <- data.frame(
        positive_cutoff = positive_cutoff,
        fisher_cutoff = fisher_cutoff,
        wilcox_cutoff = wilcox_cutoff,
        target = target,
        loocv_auc = loocv_auc$AUC[1],
        training_auc = NA,
        validation_auc = NA,
        size = NA,
        compare = compare
    )
} else {
    validation_meta <- sample_status[validate & Status_MMRd_Cancer %in% c(pop$positive_pop, pop$negative_pop)][, Contrast := Status_MMRd_Cancer %in% pop$positive_pop]
    validation_samples <- samples[Status_MMRd_Cancer %in% c(pop$positive_pop, pop$negative_pop)][, Contrast := Status_MMRd_Cancer %in% pop$positive_pop][validate == TRUE]





    testing_comb <- CJ(
        Name = Signature,
        Sample = unique(testing_samples$Sample)
    )

    testing_samples_full <- testing_samples[testing_comb, on = c("Sample", "Name")][is.na(Proportion), Proportion := 0][, .(Name, Sample, Proportion)][testing_meta, on = "Sample"]


    validation_comb <- CJ(
        Name = Signature,
        Sample = unique(validation_samples$Sample)
    )

    validation_LH <- validation_samples[validation_comb, on = c("Sample", "Name")][is.na(Proportion), Proportion := 0][, .(Name, Sample, Proportion)][validation_meta, on = "Sample"] %>%
        group_by(Sample) %>%
        group_map(
            \(df, idx){
                data.table(
                    LH = calculater_logp(testing_samples_full, df, Signature),
                    Sample = idx$Sample
                )
            }
        ) %>%
        bind_rows() %>%
        arrange(desc(LH))




    validation_auc <- calculate_AUC(validation_LH[validation_meta, on = "Sample", nomatch = NULL] %>% rename(metric = LH) %>% arrange(desc(metric)))
    filename <- paste0("validation_AUC_", positive_cutoff, "_", fisher_cutoff, "_", wilcox_cutoff, "_", target, "_", compare, ".csv")
    write.csv(validation_auc, filename, row.names = FALSE)
    print("valdation auc saved")



    training_LH <- testing_samples_full %>%
        group_by(Sample) %>%
        group_map(
            \(df, idx){
                data.table(
                    LH = calculater_logp(testing_samples_full, df, Signature),
                    Sample = idx$Sample
                )
            }
        ) %>%
        bind_rows() %>%
        arrange(desc(LH))

    training_auc <- calculate_AUC(training_LH[testing_meta, on = "Sample", nomatch = NULL] %>% rename(metric = LH) %>% arrange(desc(metric)))
    filename <- paste0("training_AUC_", positive_cutoff, "_", fisher_cutoff, "_", wilcox_cutoff, "_",  target, "_", compare, ".csv")
    write.csv(training_auc, filename, row.names = FALSE)
    print("traing auc saved")


    results <- data.frame(
        positive_cutoff = positive_cutoff,
        fisher_cutoff = fisher_cutoff,
        wilcox_cutoff = wilcox_cutoff,
        target = target,
        loocv_auc = loocv_auc$AUC[1],
        training_auc = training_auc$AUC[1],
        validation_auc = validation_auc$AUC[1],
        size = length(Signature),
        compare = compare
    )
}



filename <- paste0("summarize_", positive_cutoff, "_", fisher_cutoff, "_", wilcox_cutoff, "_",  target, "_", compare, ".csv")
write.csv(results, filename, row.names = FALSE)


