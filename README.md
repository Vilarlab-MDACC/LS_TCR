- `read_data.R`:  R code to read TCRseq data and save to RData for following process

- `summarize.R`: summarize result from `process.R` and generate ROC curve
- `Manuscript_code.R`: code to generate figures in the manuscript

- renv.lock: renv file for R enviroment used in the paper.
- main_LOOCV.nf: nextflow script for processing the classifer parapmeter selection, training and validation AUC estimate.
- select_LOOCV.R: R script called by main_LOOCV.nf
