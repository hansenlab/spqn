test_normalize_correlation <- function(){
  data(gtex.4k)
  load(file.path(path.package("spqn"), "useful", "grp_loc_target.rda"))
  cor_m <- cor(t(assay(gtex.4k)))
  ave_logrpkm <- rowData(gtex.4k)$ave_logrpkm
  cor_m_spqn_test <- normalize_correlation(cor_m, ave_exp=ave_logrpkm, ngrp=20, size_grp=300, ref_grp=18)
  checkEquals(cor_m_spqn_hash, .digest.matrix(cor_m_spqn_test))
}
