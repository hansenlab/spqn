test_normalize_correlation <- function(){
data(gtex.4k)
gtex.4k=gtex.4k[1:100,]
cor_m_sample <- cor(t(assay(gtex.4k)))
ave_logrpkm <- rowData(gtex.4k)$ave_logrpkm
cor_m_spqn_sample <- normalize_correlation(cor_m, ave_exp=ave_logrpkm, ngrp=10, size_grp=15, ref_grp=9)
}
