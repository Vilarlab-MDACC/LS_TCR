library(tidyverse)
source("functions.R")



file = commandArgs(trailingOnly=TRUE)[1]

message(file)

#file = "test_result/test.p.vs.s.Sample.rds"
#file = "test_result/test.c.vs.ps.P_ABDL25.rds"

type = file|>dirname()|>(\(x) str_split_1(x[1],"\\.")[-1]|>str_c(collapse="."))()
LOO = file|>str_split_1("\\.")|>(\(x) tail(x,2)[1])()
compare = file|>str_split_1("\\.")|>(\(x) tail(x,5)[1:3]|>str_c(collapse="."))()
meta = readRDS(str_glue("meta.{compare}.rds"))

message(type)
message(LOO)
message(compare)


RESULTFOLDER=str_glue("summarized_result.{type}")
if (!file.exists(RESULTFOLDER)) dir.create(RESULTFOLDER)

tcr = readRDS(str_glue("public.3.{type}.rds"))|>right_join(meta)
mx <- readRDS(file)



expand_grid(
  fisher_cutoff = 10^-(0:5),
  wilcox_cutoff = 10^-(-0:5),
  positive_cutoff = (0:6)/20
)|>
  filter(fisher_cutoff!=1 | fisher_cutoff!=1)|>
  pmap(
    \(fisher_cutoff,wilcox_cutoff,positive_cutoff) {
      para = as_tibble(
        list(fisher_cutoff = fisher_cutoff,
             wilcox_cutoff = wilcox_cutoff,
             positive_cutoff = positive_cutoff)
      )
      
      result_list <- list()
      mx |> filter(
        Contrast_TRUE_Present_TRUE/(Contrast_TRUE_Present_TRUE+Contrast_TRUE_Present_FALSE)>=positive_cutoff,
        fisher <=fisher_cutoff,
        wilcox <= wilcox_cutoff)|>
        bind_cols(para)-> result_list[["signature"]]
      
      if (nrow(result_list[["signature"]])>0){
        result_list[["data"]] <- tcr|>
          filter(Name%in%result_list[["signature"]]$Name)|>
          bind_cols(para)
        
        if(LOO!="Sample"){
          train = result_list[["data"]]|>filter(Sample!=LOO)
          cv = result_list[["data"]]|>filter(Sample==LOO)
        }else{
          train = result_list[["data"]]
          cv = result_list[["data"]]
        }
        
        result_list[["training_summary"]] <- train|>
          group_by(Sample,Contrast)|>
          summarise(mean_freq = mean(frequency))|>
          group_by(Contrast)|>
          summarise(
            mean = mean(mean_freq),
            sd = sd(mean_freq)
          )|>
          bind_cols(para)
        
        result_list[["cv_summary"]] <- 
          cv|>
          group_by(Sample,Contrast)|>
          summarise(frequency = mean(frequency))|>
          ungroup()|>
          mutate(prob.pos = dnorm(frequency,
                                  mean = result_list[["training_summary"]]|>filter(Contrast)|>pull(mean),
                                  sd = result_list[["training_summary"]]|>filter(Contrast)|>pull(sd)),
                 prob.neg = dnorm(frequency,
                                  mean = result_list[["training_summary"]]|>filter(!Contrast)|>pull(mean),
                                  sd = result_list[["training_summary"]]|>filter(!Contrast)|>pull(sd)),
                 log10_prob_ratio = log10(prob.pos/prob.neg)
          )|>
          bind_cols(para)
        
        result_list
      }}
  )|>
  list_transpose()|>
  imap(~bind_rows(.x))|>
  saveRDS(str_glue("{RESULTFOLDER}/{compare}.{LOO}.rds"))