test_normalize_correlation <- function(){
  load(file.path(path.package("spqn"), "unitTests", "cor_m_spqn_sample.rda"))
  data(gtex.4k)
  cor_m <- cor(t(assay(gtex.4k)))
  ave_logrpkm <- rowData(gtex.4k)$ave_logrpkm
  cor_m_spqn_test <- normalize_correlation(cor_m, ave_exp=ave_logrpkm, ngrp=20, size_grp=300, ref_grp=18)
  idx <- seq_len(500)
  cor_m_spqn_test_sample <- cor_m_spqn_test[idx, idx]
  checkEquals(cor_m_spqn_sample, cor_m_spqn_test_sample)
}

 
