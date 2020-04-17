.get_grps <- function(cor_mat, ngrp=20, size_grp=400){
    ngene <- nrow(cor_mat)
    idx <- seq_len(ngene-size_grp+1)
    grp_label <- cut(idx, ngrp-1)
    grp_loc0 <- split(idx, grp_label)
    
    if(size_grp-length(grp_loc0[[1]])<5){
        idx <- seq_len(ngene)
        grp_label <- cut(idx, ngrp)
        grp_loc <- split(idx, grp_label)
    }else{
        grp_loc <- vector(mode = "list", length = ngrp)
        idx <- seq_len(ngene-size_grp+1)
        grp_label <- cut(idx, ngrp-1)
        grp_loc0 <- split(idx, grp_label)
        for(i in seq_len(ngrp-1)){
            id1 <- grp_loc0[[i]][1]
            id2 <- grp_loc0[[i]][1]+size_grp-1
            grp_loc[[i]] <- c(id1:id2)
        }
        id1 <- ngene-size_grp+1
        grp_loc[[ngrp]] <- c(id1:ngene)
    }
    grp_loc
}


##### Asssign inner bins
## In each running group, get a inner group
.get_grps_inner <- function(grp_loc){
    grp_loc_inner <- list()
    ngrp <- length(grp_loc)
    size_bin <- length(grp_loc[[1]])
    ngene <- max(grp_loc[[ngrp]])
    
    width_tmp <- grp_loc[[2]][1] - grp_loc[[1]][1]
    length_inner_1 <- round(size_bin/2 + width_tmp/2)
    grp_loc_inner[[1]] <- seq_len(length_inner_1)
    
    for(i in 2:(ngrp-1)){
        id1 <- grp_loc[[i+1]][1]
        id2 <- grp_loc[[i]][1]
        width_tmp <- id1 - id2
        tail <- tail(grp_loc_inner[[i-1]], 1)
        grp_loc_inner[[i]] <- (tail + 1):(tail + width_tmp)
    }
    tail <- tail(grp_loc_inner[[ngrp-1]],1)
    grp_loc_inner[[ngrp]] <- c((tail+1) : ngene)
    
    grp_loc_inner
}


##### Get rank for each running bin
.get_bin_rank <- function(cor_obs, grp_loc, grp_loc_inner, cor_ref){
    rank_bin <- rank_bin_pre <- array(dim=dim(cor_obs))
    ngrp <- length(grp_loc)
    size_bin <- length(grp_loc[[1]])
    
    l_cor_tmp_ref <- length(cor_ref[upper.tri(cor_ref)])
    
    for(i in seq_len(ngrp)){
        for(j in i:ngrp){
            cor_bin_tmp <- cor_obs[grp_loc[[i]],grp_loc[[j]]]
            rank_bin_tmp <- array(dim=dim(cor_bin_tmp))
            l_cor_tmp <- length(cor_bin_tmp)
            
            idx <- seq_len(l_cor_tmp)
            rank_bin_tmp[idx] <- rank(cor_bin_tmp[idx])
            
            ## Number of diagonals(of full correlation matrix) contained in the bin
            n_diag <- sum(grp_loc[[i]] %in% grp_loc[[j]])
            
            ## Scale the rank of each bin to same scale as rank(cor_ref),
            ## After scaling, rank_bin could contain non-integars
            rank_bin_pre[grp_loc[[i]], grp_loc[[j]]] <- rank_bin_tmp
            rank_replace <- rank_bin_pre[grp_loc_inner[[i]], grp_loc_inner[[j]]]
            # scale the rank to range [0,1]
            nreplace <- l_cor_tmp-n_diag-1
            rank_replace <- (rank_replace-1)/ nreplace
            # scale the rank to the same scale as the referent bin
            rank_replace <- 1 + (rank_replace * (l_cor_tmp_ref-1))
            # change the ranks that larger than the size of referent bin, to be the size of referent bin
            rank_replace[which(rank_replace>l_cor_tmp_ref)] <- l_cor_tmp_ref
            
            # rank_bin_pre[grp_loc_inner[[i]],grp_loc_inner[[j]]] <-  rank_replace
            # rank_bin_pre[grp_loc_inner[[i]],grp_loc_inner[[j]]][which(rank_replace>l_cor_tmp_ref)] <- l_cor_tmp_ref
       
            rank_bin[grp_loc_inner[[i]],grp_loc_inner[[j]]] <- rank_replace
            # rank_bin[grp_loc_inner[[i]],grp_loc_inner[[j]]] <- rank_bin_pre[grp_loc_inner[[i]],grp_loc_inner[[j]]]
        }
    }
    
    rank_bin
}

##### Transform rank to cor_est
.est_cor <- function(rank_bin, cor_ref){
    cor_adj <- array(dim=dim(rank_bin))
    
    cor_ref_upper_tri <- cor_ref[upper.tri(cor_ref)]
    cor_ref_sorted <- sort(cor_ref_upper_tri)
    
    up_tri <- upper.tri(cor_adj)
    
    ## Find to nearest integars to rank_bin, since rank_bin could contain non-integars
    rank_bin <- rank_bin[up_tri]
    rank_bin1 <- rank_bin %/% 1

    ## Assign weights to each rank according to the distance to rank_bin
    rank_bin_w2 <- (rank_bin - rank_bin1)

    rank_bin2 <- rank_bin1+1
    rank_bin_w1 <- 1-rank_bin_w2

    ## Find the correlations in the cor_ref corresponding to the two nearest ranks to rank_bin
    ## Estimate the correlation using weighted average based on the distance to rank_bin
    length_cor_ref <- length(cor_ref_sorted)
    rank_bin2[which(rank_bin2>length_cor_ref)] <- length_cor_ref
    cor_adj[up_tri] <- rank_bin_w1*cor_ref_sorted[rank_bin1]+ rank_bin_w2*cor_ref_sorted[rank_bin2]

    low_tri <- lower.tri(cor_adj)
    cor_adj[low_tri] <- t(cor_adj)[low_tri]
    diag(cor_adj) <- 1
    cor_adj
}
normalize_correlation <- function(cor_mat, ave_exp, ngrp, size_grp, ref_grp){
  ## stopifnot(is.integer(ngrp), is.integer(size_grp), is.integer(ref_grp))
  stopifnot(length(ngrp) == 1, length(size_grp) == 1, length(ref_grp) == 1)
  stopifnot(ref_grp <= ngrp)
  stopifnot(0 < ngrp, 0 < size_grp, 0 < ref_grp)
  stopifnot(is.matrix(cor_mat))
  stopifnot(nrow(cor_mat) == ncol(cor_mat), nrow(cor_mat) == length(ave_exp))
  
  idx <- order(ave_exp)
  ngene <- length(ave_exp)
  ref_vec <- seq_len(ngene)
  
  cor_mat <- cor_mat[idx, idx]
  ref_vec_rearranged <- ref_vec[idx]
  
  group_loc <- .get_grps(cor_mat, ngrp, size_grp)
  
  ## Asssign inner bins
  group_loc_adj <- .get_grps_inner(group_loc)
  
  ## Get rank for each running bin
  cor_ref <- cor_mat[group_loc[[ref_grp]], group_loc[[ref_grp]]]
  rank_bin <- .get_bin_rank(cor_mat, group_loc, group_loc_adj, cor_ref)
  
  ## Transform rank to cor_adj
  cor_est <- .est_cor(rank_bin, cor_ref)
  
  idx <- order(ref_vec_rearranged)
  cor_est <- cor_est[idx, idx]
  cor_est
}

