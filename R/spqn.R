.get_grps <- function(cor_mat, ngrp=20, size_grp=400){
    ngene <- nrow(cor_mat)
    #idx <- seq_len(ngene-size_grp+1)
    #grp_label <- cut(idx, ngrp-1)
    #grp_loc0 <- split(idx, grp_label)
    #grp_loc <- list()
    
    if(size_grp-length(grp_loc0[[1]])<5){
        idx <- seq_len(ngene)
        grp_label <- cut(idx, ngrp)
        grp_loc <- split(idx, grp_label)
    }else{
        grp_loc <- vector(mode = "list", length = ngrp)
        idx = seq_len(ngene-size_grp+1)
        grp_label <- cut(idx, ngrp-1)
        grp_loc0 <- split(idx, grp_label)
        for(i in seq_len(ngrp-1)){
            grp_loc[[i]] <- c(grp_loc0[[i]][1]:(grp_loc0[[i]][1]+size_grp-1))
        }
        grp_loc[[ngrp]] <- c((ngene-size_grp+1):ngene)
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
    grp_loc_inner[[1]] <- c(seq_len(round(size_bin/2+width_tmp/2)))
    
    for(i in 2:(ngrp-1)){
        width_tmp <- grp_loc[[i+1]][1] - grp_loc[[i]][1]
        # grp_loc_inner[[i]] <- c( (tail(grp_loc_inner[[i-1]],1)+1) :(tail(grp_loc_inner[[i-1]],1) + width_tmp))
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
            rank_bin_pre[grp_loc_inner[[i]],grp_loc_inner[[j]]] <-  1+((rank_bin_pre[grp_loc_inner[[i]], grp_loc_inner[[j]]]-1)/(l_cor_tmp-n_diag-1) *(l_cor_tmp_ref-1))
            rank_bin_pre[grp_loc_inner[[i]],grp_loc_inner[[j]]][which(rank_bin_pre[grp_loc_inner[[i]],grp_loc_inner[[j]]]>l_cor_tmp_ref)] <- l_cor_tmp_ref
            
            rank_bin[grp_loc_inner[[i]],grp_loc_inner[[j]]] <- rank_bin_pre[grp_loc_inner[[i]],grp_loc_inner[[j]]]
        }
    }
    
    rank_bin
}

##### Transform rank to cor_est
.est_cor <- function(rank_bin, cor_ref){
    cor_adj <- array(dim=dim(rank_bin))
    cor_ref_sorted <- sort(cor_ref[upper.tri(cor_ref)])
    
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
    rank_bin2[which(rank_bin2>length(cor_ref_sorted))] <- length(cor_ref_sorted)
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

    rownames(cor_mat) <- colnames(cor_mat) <- seq_len(length(ave_exp))
    cor_mat <- cor_mat[order(ave_exp), order(ave_exp)]
    
    group_loc <- .get_grps(cor_mat, ngrp, size_grp)
    
    ## Asssign inner bins
    group_loc_adj <- .get_grps_inner(group_loc)
    
    ## Get rank for each running bin
    cor_ref <- cor_mat[group_loc[[ref_grp]], group_loc[[ref_grp]]]
    rank_bin <- .get_bin_rank(cor_mat, group_loc, group_loc_adj, cor_ref)
    
    ## Transform rank to cor_adj
    cor_est <- .est_cor(rank_bin, cor_ref)
    
    cor_est <- cor_est[order(as.numeric(rownames(cor_mat))),
                       order(as.numeric(colnames(cor_mat)))]
    cor_est
}

