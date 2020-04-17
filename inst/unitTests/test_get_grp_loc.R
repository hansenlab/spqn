test_get_grps <- function(){
  data(grp_loc_target)
  load(file.path(path.package("spqn"), "unitTests", "grp_loc_target.rda"))

  grp_loc_assigend = spqn:::get_grp_loc(cor_mat=array(dim=c(100,100)),ngrp=10)
  grp_loc_assigend =  matrix(unlist(grp_loc_assigend), ncol = 10, byrow = TRUE)
    
  checkEquals(grp_loc_target, grp_loc_assigend)
}

