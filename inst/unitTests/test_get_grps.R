test_get_grps <- function(){
  data(grp_loc_target)
  load(file.path(path.package("spqn"), "useful", "grp_loc_target.rda"))

  grp_loc_assigend = get_grp_loc(cor_mat=array(dim=c(100,100)),ngrp=10)

  all(compare.list(grp_loc_target, grp_loc_assigend))
}


