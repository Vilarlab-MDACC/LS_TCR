library(tidyverse)
source("functions.R")

compare = commandArgs(trailingOnly=TRUE)[1]
LOO = commandArgs(trailingOnly=TRUE)[2]
input = commandArgs(trailingOnly=TRUE)[3]


RESULTFOLDER = str_glue("result.{input}")


if (!file.exists(RESULTFOLDER)) dir.create(RESULTFOLDER)



result = readRDS(str_glue("public.3.{input}.rds"))|>filter(Sample!=LOO)|>right_join(readRDS(str_glue("meta.{compare}.rds")))


result|>fast_fisher() ->tmp_fisher
result|>filter(Name%in%tmp_fisher$Name) -> result

result|>group_by(Name,Contrast)|>summarise(mean_freq = mean(frequency))|>ungroup()|>pivot_wider(id_cols = Name, names_from =Contrast,names_prefix="mean_freq_",values_from =mean_freq)|>mutate(mean_freq_diff = mean_freq_TRUE - mean_freq_FALSE)|>filter(mean_freq_diff>0) -> tmp_mean

result|>filter(Name%in%tmp_mean$Name)->result

result|>group_by(Name)|>summarise(wilcox = wilcox.test(frequency~Contrast)$p.value)->tmp_wilcox

right_join(tmp_fisher,tmp_mean)|>right_join(tmp_wilcox)|>saveRDS(str_glue("{RESULTFOLDER}/test.{compare}.{LOO}.rds"))