library(spqn)
library(spqnData)
library(SummarizedExperiment)

data(gtex.4k)
cor_m <- cor(t(assay(gtex.4k)))
ave_logrpkm <- rowData(gtex.4k)$ave_logrpkm

cor_m_spqn <- normalize_correlation(cor_m, ave_exp=ave_logrpkm, ngrp=20, size_grp=300, ref_grp=18)

cor_m_spqn_sample = cor_m_spqn[1:500,1:500]
save(cor_m_spqn_sample, file=file.path("..", "unitTests", "cor_m_spqn_sample.rda"))

sessionInfo()
