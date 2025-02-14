#TCRseq manuscript code

library(immunarch)
library(readxl)
library(stats)
library(tidyverse)
library(ggplot2)
library(ggdist)
library(DHARMa)
library(gap)
library(sjPlot)
library(lme4)
library(glmmTMB)

#Diversity analysis for PBMCs/Blood samples
#Data loading: all samples should be in a single folder with a metadata file containing a column named "Sample" with the exact same names of samples without extensions
data_final <- repLoad("/home/ccp/vilarsanchez/tcrseq/Data_final/")

#Effective number of clonotypes and Simpson Clonality Index values were obtained from Adaptive's ImmunoSEQ analyzer

#Homeostasis analysis
imm_hom_all_final <- repClonality(data_final$data,
                                  .method = "homeo",
                                  .clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1))

#Exported values and used Prism to create graphs
write.csv(imm_hom_all_final, file = "Homeostasis_final.csv") 

#Gamma Generalized Linear Model analysis of Simpson Clonality Index
#Density plot of Simpson Clonality Index
plot(density(metadata_for_manova$Simpson_index), main = "Density Plot of Simpson Index")

#GLM
model <- glm(Simpson_index ~ Status2 + Institution2 + Gender +  MMRgene +  + Input_material + Age_Range3, 
              family = Gamma(link = "log"), data = metadata_for_manova)
summary_model <- summary(model)
simulationOutput <- simulateResiduals(fittedModel = model)
plot(simulationOutput)
tab_model(model)

#All other diversity graphs were created with graphpad Prism



#Public Repertoire
#A public TCR is defined as a TCR shared by at least 2 individuals by looking for exact matches of the CDR3, V and J segments at aminoacid level
public_repertoire_proportion_final_2 <- pubRep(data_final$data, "aa+v+j", .verbose = F, .quant = "prop", .coding = TRUE, .min.samples = 2)

#Public TCRs filtered for Survivors
public_survivors_proportion_final_2 <- pubRepFilter(public_repertoire_proportion_final_2, data_final$meta, c(Status_MMRd_Cancer = "Survivor"))

#Public TCRs filtered for Previvors
public_previvors_proportion_final_2 <- pubRepFilter(public_repertoire_proportion_final_2, data_final$meta, c(Status_MMRd_Cancer = "Previvor"))

#Public TCRs filtered for Controls
public_controls_proportion_final_2 <- pubRepFilter(public_repertoire_proportion_final_2, data_final$meta, c(Status_MMRd_Cancer = "Control"))

#Join CDR3, V and J columns into a single column
public_survivors_proportion_final_2$CDR3VJ <- paste(public_survivors_proportion_final_2$CDR3.aa, ",", public_survivors_proportion_final_2$V.name, ",", public_survivors_proportion_final_2$J.name)
public_previvors_proportion_final_2$CDR3VJ <- paste(public_previvors_proportion_final_2$CDR3.aa, ",", public_previvors_proportion_final_2$V.name, ",", public_previvors_proportion_final_2$J.name)
public_controls_proportion_final_2$CDR3VJ <- paste(public_controls_proportion_final_2$CDR3.aa, ",", public_controls_proportion_final_2$V.name, ",", public_controls_proportion_final_2$J.name)


#export public repertoires for each group
write.csv(public_survivors_proportion_final_2, file= "public_survivors_proportion_final_2.csv")
write.csv(public_previvors_proportion_final_2, file= "public_previvors_proportion_final_2.csv")
write.csv(public_controls_proportion_final_2, file= "public_controls_proportion_final_2.csv")

#Nan helped mi with upset plot

#Percentage of individuals by which each TCR is shared graph
#With these values, it was transformed to percetanges and grpahed in Prism
Controls_distribution_public_TCRs <- as.data.frame(table(public_controls_proportion_final_2$Samples))
Previvors_distribution_public_TCRs <- as.data.frame(table(public_previvors_proportion_final_2$Samples))
Survivors_distribution_public_TCRs <- as.data.frame(table(public_survivors_proportion_final_2$Samples))


#Nan obtained Jaccard and Morisita indexes

#Tissue and PBMCs concordance
P_1008 <- repLoad("/home/ccp/vilarsanchez/tcrseq/Tissue_PBMC_pairs/P_1008")
Track_1008 <- trackClonotypes(P_1008$data, list(1,100), .col = "aa+v+j")
graph_Track_1008 <- vis(Track_1008, .order = c(1, 2)) +
  theme(axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        text = element_text(size = 4),
        axis.text.y = element_text(size = 10)) 
graph_Track_1008

P_2468<- repLoad("/home/ccp/vilarsanchez/tcrseq/Tissue_PBMC_pairs/P_2468")
Track_2468 <- trackClonotypes(P_2468$data, list(1,100), .col = "aa+v+j")
graph_Track_2468 <- vis(Track_2468, .order = c(1, 2)) +
  theme(axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        text = element_text(size = 5),
        axis.text.y = element_text(size = 10))
graph_Track_2468

P_Lyn16 <- repLoad("/home/ccp/vilarsanchez/tcrseq/Tissue_PBMC_pairs/P_Lyn16")
Track_Lyn16 <- trackClonotypes(P_Lyn16$data, list(1,100), .col = "aa+v+j")
graph_Track_Lyn16 <- vis(Track_Lyn16, .order = c(1, 2)) +
  theme(axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        text = element_text(size = 4),
        axis.text.y = element_text(size = 10)) 
graph_Track_Lyn16

P_Lyn22 <- repLoad("/home/ccp/vilarsanchez/tcrseq/Tissue_PBMC_pairs/P_Lyn22")
Track_Lyn22 <- trackClonotypes(P_Lyn22$data, list(1,100), .col = "aa+v+j")
graph_Track_Lyn22 <- vis(Track_Lyn22, .order = c(1, 2)) +
  theme(axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        text = element_text(size = 5),
        axis.text.y = element_text(size = 10)) 
graph_Track_Lyn22 

S_1305 <- repLoad("/home/ccp/vilarsanchez/tcrseq/Tissue_PBMC_pairs/S_1305")
Track_1305_T1 <- trackClonotypes(S_1305$data, list(1,100), .col = "aa+v+j")
Track_1305_T2 <- trackClonotypes(S_1305$data, list(2,100), .col = "aa+v+j")
Track_1305_TA <- trackClonotypes(S_1305$data, list(3,100), .col = "aa+v+j")
Track_1305_TA1 <- trackClonotypes(S_1305$data, list(4,100), .col = "aa+v+j")
Track_1305_TA2 <- trackClonotypes(S_1305$data, list(5,100), .col = "aa+v+j")
graph_Track_1305_T1 <- vis(Track_1305_T1, .order = c(3,2,4,1,5,6)) +
  theme(axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        text = element_text(size = 4),
        axis.text.y = element_text(size = 10)) 
graph_Track_1305_T1

graph_Track_1305_T2 <- vis(Track_1305_T2, .order = c(3,2,4,1,5,6)) +
  theme(axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        text = element_text(size = 4),
        axis.text.y = element_text(size = 10)) 
graph_Track_1305_T2

graph_Track_1305_TA <- vis(Track_1305_TA, .order = c(3,2,4,1,5,6)) +
  theme(axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        text = element_text(size = 4),
        axis.text.y = element_text(size = 10)) 
graph_Track_1305_TA

graph_Track_1305_TA1 <- vis(Track_1305_TA1, .order = c(3,2,4,1,5,6)) +
  theme(axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        text = element_text(size = 4),
        axis.text.y = element_text(size = 10)) 
graph_Track_1305_TA1

graph_Track_1305_TA2 <- vis(Track_1305_TA2, .order = c(3,2,4,1,5,6)) +
  theme(axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        text = element_text(size = 4),
        axis.text.y = element_text(size = 10)) 
graph_Track_1305_TA2

S_8421 <- repLoad("/home/ccp/vilarsanchez/tcrseq/Tissue_PBMC_pairs/S_8421")
Track_8421 <- trackClonotypes(S_8421$data, list(1,100), .col = "aa+v+j")
graph_Track_8421 <- vis(Track_8421, .order = c(1,2)) +
  theme(axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        text = element_text(size = 4),
        axis.text.y = element_text(size = 10)) 
graph_Track_8421

S_4670 <- repLoad("/home/ccp/vilarsanchez/tcrseq/Tissue_PBMC_pairs/S_4670")
Track_4670 <- trackClonotypes(S_4670$data, list(1,100), .col = "aa+v+j")
graph_Track_4670 <- vis(Track_4670, .order = c(1,2)) +
  theme(axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        text = element_text(size = 4),
        axis.text.y = element_text(size = 10)) 
graph_Track_4670


S_5642 <- repLoad("/home/ccp/vilarsanchez/tcrseq/Tissue_PBMC_pairs/S_5642")
Track_5642 <- trackClonotypes(S_5642$data, list(1,100), .col = "aa+v+j")
graph_Track_5642 <- vis(Track_5642, .order = c(1, 2)) +
  theme(axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        text = element_text(size = 5),
        axis.text.y = element_text(size = 10)) 
graph_Track_5642

S_8228 <- repLoad("/home/ccp/vilarsanchez/tcrseq/Tissue_PBMC_pairs/S_8228")
Track_8228 <- trackClonotypes(S_8228$data, list(1,100), .col = "aa+v+j")
graph_Track_8421 <- vis(Track_8421, .order = c(1,2)) +
  theme(axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        text = element_text(size = 4),
        axis.text.y = element_text(size = 10))  
graph_Track_8421

S_lyn11 <- repLoad("/home/ccp/vilarsanchez/tcrseq/Tissue_PBMC_pairs/S_Lyn11")
Track_Lyn11 <- trackClonotypes(S_lyn11$data, list(1,100), .col = "aa+v+j")
graph_Track_8421 <- vis(Track_8421, .order = c(1,2)) +
  theme(axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        text = element_text(size = 4),
        axis.text.y = element_text(size = 10))
graph_Track_8421 

#The numbers in the grpah were added manually based on the number of overlapping TCRs
#Overlap with public repertoire
Track_8421_merged <- merge(Track_8421, public_repertoire_proportion_final_2, 
                           by.x = "CDR3VJ", by.j = "CDR3VJ", all.x = TRUE, all.y = FALSE)
write.csv(Track_8421_merged, "8421_blood.csv")

Track_1008$CDR3VJ <- paste(Track_1008$CDR3.aa, ",", Track_1008$V.name, ",", Track_1008$J.name)
Track_1008_merged <- merge(Track_1008, public_repertoire_proportion_final_2, 
                           by = "CDR3VJ")
write.csv(Track_1008_merged, "1008_blood.csv")

Track_2468$CDR3VJ <- paste(Track_2468$CDR3.aa, ",", Track_2468$V.name, ",", Track_2468$J.name)
Track_2468_merged <- merge(Track_2468, public_repertoire_proportion_final_2, 
                           by = "CDR3VJ")
write.csv(Track_2468_merged, "2468_blood.csv")

Track_Lyn16$CDR3VJ <- paste(Track_Lyn16$CDR3.aa, ",", Track_Lyn16$V.name, ",", Track_Lyn16$J.name)
Track_Lyn16_merged <- merge(Track_Lyn16, public_repertoire_proportion_final_2, 
                            by = "CDR3VJ")
write.csv(Track_Lyn16_merged, "Lyn16_blood.csv")

Track_Lyn22$CDR3VJ <- paste(Track_Lyn22$CDR3.aa, ",", Track_Lyn22$V.name, ",", Track_Lyn22$J.name)
Track_Lyn22_merged <- merge(Track_Lyn22, public_repertoire_proportion_final_2, 
                            by = "CDR3VJ")
write.csv(Track_Lyn22_merged, "Lyn22_blood.csv")


Track_1305_T1$CDR3VJ <- paste(Track_1305_T1$CDR3.aa, ",", Track_1305_T1$V.name, ",", Track_1305_T1$J.name)
Track_1305_T1_merged <- merge(Track_1305_T1, public_repertoire_proportion_final_2, 
                              by = "CDR3VJ")
write.csv(Track_1305_T1_merged, "1305_T1_blood.csv")

Track_1305_T2$CDR3VJ <- paste(Track_1305_T2$CDR3.aa, ",", Track_1305_T2$V.name, ",", Track_1305_T2$J.name)
Track_1305_T2_merged <- merge(Track_1305_T2, public_repertoire_proportion_final_2, 
                              by = "CDR3VJ")
write.csv(Track_1305_T2_merged, "1305_T2_blood.csv")

Track_1305_TA$CDR3VJ <- paste(Track_1305_TA$CDR3.aa, ",", Track_1305_TA$V.name, ",", Track_1305_TA$J.name)
Track_1305_TA_merged <- merge(Track_1305_TA, public_repertoire_proportion_final_2, 
                              by = "CDR3VJ")
write.csv(Track_1305_TA_merged, "1305_TA_blood.csv")


Track_1305_TA1$CDR3VJ <- paste(Track_1305_TA1$CDR3.aa, ",", Track_1305_TA1$V.name, ",", Track_1305_TA1$J.name)
Track_1305_TA1_merged <- merge(Track_1305_TA1, public_repertoire_proportion_final_2, 
                               by = "CDR3VJ")
write.csv(Track_1305_TA1_merged, "1305_TA1_blood.csv")


Track_1305_TA2$CDR3VJ <- paste(Track_1305_TA2$CDR3.aa, ",", Track_1305_TA2$V.name, ",", Track_1305_TA2$J.name)
Track_1305_TA2_merged <- merge(Track_1305_TA2, public_repertoire_proportion_final_2, 
                               by = "CDR3VJ")
write.csv(Track_1305_TA2_merged, "1305_TA2_blood.csv")

Track_4670$CDR3VJ <- paste(Track_4670$CDR3.aa, ",", Track_4670$V.name, ",", Track_4670$J.name)
Track_4670_merged <- merge(Track_4670, public_repertoire_proportion_final_2, 
                           by = "CDR3VJ")
write.csv(Track_4670_merged, "4670_blood.csv")



Track_5642$CDR3VJ <- paste(Track_5642$CDR3.aa, ",", Track_5642$V.name, ",", Track_5642$J.name)
Track_5642_merged <- merge(Track_5642, public_repertoire_proportion_final_2, 
                           by = "CDR3VJ")
write.csv(Track_5642_merged, "5642_blood.csv")

Track_8228$CDR3VJ <- paste(Track_8228$CDR3.aa, ",", Track_8228$V.name, ",", Track_8228$J.name)
Track_8228_merged <- merge(Track_8228, public_repertoire_proportion_final_2, 
                           by = "CDR3VJ")
write.csv(Track_8228_merged, "8228_blood.csv")

Track_Lyn11$CDR3VJ <- paste(Track_Lyn11$CDR3.aa, ",", Track_Lyn11$V.name, ",", Track_Lyn11$J.name)
Track_Lyn11_merged <- merge(Track_Lyn11, public_repertoire_proportion_final_2, 
                            by = "CDR3VJ")
write.csv(Track_Lyn11_merged, "Lyn11_blood.csv")



#Annotation against McPAS TCR database and dMMR CRC

#I downloaded the database only for humans from the website and loaded it here as follows
mcpas = dbLoad("/home/ccp/vilarsanchez/tcrseq/Annotation/MCPAS_database/McPas_Database.csv", "mcpas")
mcpas$CDR3VJ <- paste(mcpas$CDR3.beta.aa, ",", mcpas$TRBV, ",", mcpas$TRBJ)

#Annotation against McPAS. I then manually checked for the numbers in Pathogens category and in Automimmune conditions
annotatated_survivors_final_mcpas <- inner_join(public_survivors_proportion_final_2, mcpas, by= c("CDR3.aa" = "CDR3.beta.aa", "V.name" = "TRBV",  "J.name" = "TRBJ"))
annotatated_previvors_final_mcpas <- inner_join(public_previvors_proportion_final_2, mcpas, by= c("CDR3.aa" = "CDR3.beta.aa", "V.name" = "TRBV",  "J.name" = "TRBJ"))
annotatated_controls_final_mcpas <- inner_join(public_controls_proportion_final_2, mcpas, by= c("CDR3.aa" = "CDR3.beta.aa", "V.name" = "TRBV",  "J.name" = "TRBJ"))

#Annotation against dMMR tumors

#data loading
tumors_for_annotation <- repLoad("/ccp_projects/Laboratory-Vilar/Ana/TCRseq_Project/Analysis/Annotation/Tumors_for_annotation")

#Pool of dMMR resident TCRs
repertoires_dMMR_tumors_for_annotation <- pubRep(tumors_for_annotation$data, "aa+v+j", .verbose = F, .quant = "prop", .coding = TRUE, .min.samples = 1)

#Annotation
annotatated_controls_final_dMMR_tumors <- inner_join(public_controls_proportion_final_2, repertoires_dMMR_tumors_for_annotation, by= c("CDR3.aa", "V.name",  "J.name"))
annotatated_previvors_final_dMMR_tumors <- inner_join(public_previvors_proportion_final_2, repertoires_dMMR_tumors_for_annotation, by= c("CDR3.aa", "V.name",  "J.name"))
annotatated_survivors_final_dMMR_tumors <- inner_join(public_survivors_proportion_final_2, repertoires_dMMR_tumors_for_annotation, by= c("CDR3.aa", "V.name",  "J.name"))

#Nan did all the classifier analyses and graphs

#Single cell validation analyses
#Viral peptide for EBV
BMLF1_EBV_results <- read.csv("/ccp_projects/Laboratory-Vilar/Ana/TCRseq_Project/Viral_peptides/Clonotypes_sc_BMLF1_EBV.csv", header = TRUE)
BMLF1_EBV_results$CDR3VJ <- paste(BMLF1_EBV_results$CDR3.aa, ",", BMLF1_EBV_results$trb_v_genes, ",", BMLF1_EBV_results$trb_j_genes)
BMLF1_public_merge <- merge(BMLF1_EBV_results, public_repertoire_proportion_final_2, 
                            by = "CDR3VJ")

write.csv(BMLF1_public_merge, file = "BMLF1_public_merge.csv")

#graphs were made on graphpad

#RNF43 validation of TCR signatures

#Signatures load
format_gene <- function(.col) {
  .col%>%
    str_replace(",",", ")%>%
    str_replace("-([0])([0-9])", "-\\2")%>%
    str_replace("([VDJ])([0])([0-9])", "\\1\\3")%>%
    str_replace("TCR", "TR")
}


read_tsv("/ccp_projects/Laboratory-Vilar/Ana/TCRseq_Project/Analysis/ML_Classifier/classifier/c.vs.ps_0.001_0.001_0.15_aa.j.tsv")%>%
  separate(Name,into=c("CDR3.aa","J.name","V.name"),sep="_")%>%
  mutate(across(ends_with("name"),format_gene)) -> target_TCR_lynch

read_tsv("/ccp_projects/Laboratory-Vilar/Ana/TCRseq_Project/Analysis/ML_Classifier/classifier/p.vs.s_0.1_0.01_0.15_aa.j.tsv")%>%
  separate(Name,into=c("CDR3.aa","J.name","V.name"),sep="_")%>%
  mutate(across(ends_with("name"),format_gene)) -> target_TCR_survivors

filtered_data = fulldata

for (i in names(fulldata$data)){
  filtered_data$data[[i]] = fulldata$data[[i]]%>%inner_join(target_TCR,na_matches ="na" )
}

target_TCR_lynch$CDR3VJ <- paste(target_TCR_lynch$CDR3.aa, ",", target_TCR_lynch$V.name, ",", target_TCR_lynch$J.name)

#Overlap of signatures against dMMR tumors pool
overlap_Lynch_signature_with_dMMR_tumors <- inner_join(target_TCR_lynch, repertoires_dMMR_tumors_for_annotation, by= c("CDR3.aa", "J.name"))
write.csv(overlap_Lynch_signature_with_dMMR_tumors, "overlap_Lynch_signature_with_dMMR_tumors.csv")
overlap_survivors_signature_with_dMMR_tumors <- inner_join(target_TCR_survivors, repertoires_dMMR_tumors_for_annotation, by= c("CDR3.aa", "J.name"))

#overlap of signatures against RNF43 TCRs
RNF43_results <- read.csv("/ccp_projects/Laboratory-Vilar/Ana/TCRseq_Project/Analysis/In_vitro_validation/Clonotypes_RNF43.csv", header = TRUE)
overlap_Lynch_signature_with_RNF43_results <- inner_join(target_TCR_lynch, RNF43_results, by= c("CDR3.aa", "J.name"))
write.csv(overlap_Lynch_signature_with_RNF43_results, "overlap_Lynch_signature_with_RNF43_results.csv")
overlap_survivors_signature_with_RNF43_results <- inner_join(target_TCR_survivors,RNF43_results, by= c("CDR3.aa", "J.name"))


