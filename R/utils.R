.get_grps <- function(cor_mat, ngrp=20, size_grp=400){
    ngene <- nrow(cor_mat)
    grp_label <- cut(seq_len(ngene-size_grp+1), ngrp-1)
    grp_loc0 <- split(seq_len(ngene-size_grp+1), grp_label)
    grp_loc <- list()

    if(size_grp-length(grp_loc0[[1]])<5){
        grp_label <- cut(seq_len(ngene), ngrp)
        grp_loc0 <- split(seq_len(ngene), grp_label)
        for(i in seq_len(ngrp)){
            grp_loc[[i]] <- grp_loc0[[i]]
        }
    }else{
        grp_label <- cut(seq_len(ngene-size_grp+1), ngrp-1)
        grp_loc0 <- split(seq_len(ngene-size_grp+1), grp_label)
        for(i in seq_len(ngrp-1)){
            grp_loc[[i]] <- c(grp_loc0[[i]][1]:(grp_loc0[[i]][1]+size_grp-1))
        }
        grp_loc[[ngrp]] <- c((ngene-size_grp+1):ngene)
    }
    grp_loc
}
