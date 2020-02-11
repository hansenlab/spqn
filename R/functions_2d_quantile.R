get_grps<-function(avelog2cpm,ngrp=20,size_grp=400){
  ngene=nrow(avelog2cpm) 
  grp_label=cut(1:(ngene-size_grp+1),ngrp-1)
  grp_loc0=split(1:(ngene-size_grp+1),grp_label)
  grp_loc=list()

  if(size_grp-length(grp_loc0[[1]])<5){
    grp_label=cut(1:ngene,ngrp)
    grp_loc0=split(1:ngene,grp_label)
    for(i in 1:(ngrp)){
      grp_loc[[i]]=grp_loc0[[i]]
    }
  }else{
    grp_label=cut(1:(ngene-size_grp+1),ngrp-1)
    grp_loc0=split(1:(ngene-size_grp+1),grp_label)
    for(i in 1:(ngrp-1)){
      grp_loc[[i]]=c(grp_loc0[[i]][1]:(grp_loc0[[i]][1]+size_grp-1))
    }
    grp_loc[[ngrp]]=c((ngene-size_grp+1):ngene)
  }
  grp_loc
}
