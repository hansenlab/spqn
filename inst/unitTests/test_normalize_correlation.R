test_normalize_correlation <- function(){
  load(file.path(path.package("spqn"), "unitTests", "ave_logrpkm.rda"))
  load(file.path(path.package("spqn"), "unitTests", "exp_mat.rda"))
  load(file.path(path.package("spqn"), "unitTests", "cor_m_spqn_hash.rda"))
  
  cor_m <- cor(t(exp_mat))

  cor_m_spqn_test <- normalize_correlation(cor_m, ave_exp=ave_logrpkm, ngrp=20, size_grp=300, ref_grp=18)
  checkEquals(cor_m_spqn_hash, spqn:::.digestMatrix(cor_m_spqn_test))
}

 
