library(immunarch)
library(tidyverse)

repLoad("../TCR_all/ALL_SAMPLES/")$data[[1]]%>%
  pubRep(immundata$data[[1]], "aa+v+j")%>%
  saveRDS(pr.aav,"public.3.aa.j.v.rds")


repLoad("../TCR_all/ALL_SAMPLES/")$data[[1]]%>%
  pubRep(immundata$data[[1]], "aa+v")%>%
  saveRDS(pr.aav,"public.3.aa.v.rds")


repLoad("../TCR_all/ALL_SAMPLES/")$data[[1]]%>%
  pubRep(immundata$data[[1]], "aa+j")%>%
  saveRDS(pr.aav,"public.3.aa.j.rds")


repLoad("../TCR_all/ALL_SAMPLES/")$data[[1]]%>%
  pubRep(immundata$data[[1]], "aa")%>%
  saveRDS(pr.aav,"public.3.aa.rds")