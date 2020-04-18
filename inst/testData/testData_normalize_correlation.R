library(spqn)
library(spqnData)
library(SummarizedExperiment)
library(RUnit)

dir_data = "/users/ywang/packages/data_generation/" 
data(gtex.4k)
cor_m <- cor(t(assay(gtex.4k)))

ave_logrpkm <- rowData(gtex.4k)$ave_logrpkm

cor_m_spqn <- normalize_correlation(cor_m, ave_exp=ave_logrpkm, ngrp=20, size_grp=300, ref_grp=18)

cor_m_spqn_sample_linux = cor_m_spqn[1:500,1:500]
cor_m_spqn_sample_linux = round(cor_m_spqn_sample_linux,4)
save(cor_m_spqn_sample_linux,file=paste0(dir_data,"cor_m_spqn_sample_linux.rda"))


cor_m_spqn_sample = cor_m_spqn_sample_linux
save(cor_m_spqn_sample,file=paste0(dir_data,"cor_m_spqn_sample.rda"))

sessionInfo()
