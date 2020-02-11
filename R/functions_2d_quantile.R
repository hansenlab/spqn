cpm<-function(counts){
  counts_cpm=t( (t(counts)/colSums(counts)) * (1e+6)   )
  counts_cpm
}


rpm_filt_v2 <- function(counts_raw){
  counts_log2cpm=log2(cpm(counts_raw)+0.5)
  keep=which(rowMedians(data.matrix(counts_log2cpm))>0)
  counts_log2cpm_keep=counts_log2cpm[keep,]
  counts_raw_keep=counts_raw[keep,]
  return(list(counts_raw_keep=counts_raw_keep,counts_log2cpm_keep=counts_log2cpm_keep,keep=keep))
}

removePCs <- function(counts_logrpm, num_PCs = 4){
  if(num_PCs > 0) {
    counts_removePCs= t(removePrincipalComponents( scale(t(counts_logrpm)), num_PCs))
    return(counts_removePCs)} else{return(counts_logrpm)}
}




## step 1, asssign running bins
## group the sorted genes into equal-size running groups, each group containing 400 genes 

create_bins <- function(n_genes, n_bins = 20, size_bins = 400) {

    bin_label <- cut(1:(n_genes - size_bins + 1), n_bins-1)
    grp_loc0 <- split(1:(n_genes - size_bins+1), bin_label)

    outer_bins <- list()
    for(i in 1:(n_bins - 1)){
        outer_bins[[i]] <- c(grp_loc0[[i]][1]:(grp_loc0[[i]][1] + size_bins - 1))
    }
    outer_bins[[n_bins]] <- c((n_genes-size_bins+1):n_genes)
    ## outer bins done
    
    inner_bins <- list()
    size_bins <- length(outer_bins[[1]])
        
    width_tmp <- outer_bins[[2]][1] - outer_bins[[1]][1] 
    inner_bins[[1]] <- c(1:round(size_bins/2+width_tmp/2))
    
    for(i in 2:(n_bins-1)){
        width_tmp <- outer_bins[[i+1]][1] - outer_bins[[i]][1] 
        inner_bins[[i]] <- c( (tail(inner_bins[[i-1]],1)+1) :(tail(inner_bins[[i-1]],1) + width_tmp))
    }  
    inner_bins[[n_bins]] <- c( (tail(inner_bins[[n_bins-1]],1)+1) : n_genes)
    ## inner bins done

    
    list(inner = inner_bins, outer = outer_bins)
}






get_grps<-function(avelog2cpm,ngrp=20,size_grp=400){
  ngene=length(avelog2cpm) 
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


##### step 2, asssign inner bins
# in each running group, get a inner group
get_grps_inner <- function(grp_loc){
    grp_loc_inner <- list()
    ngrp <- length(grp_loc)
    size_bin <- length(grp_loc[[1]])
    ngene <- max(grp_loc[[ngrp]])
    
    width_tmp <- grp_loc[[2]][1] - grp_loc[[1]][1] 
    grp_loc_inner[[1]] <- c(1:round(size_bin/2+width_tmp/2))
    
    for(i in 2:(ngrp-1)){
        width_tmp <- grp_loc[[i+1]][1] - grp_loc[[i]][1] 
        grp_loc_inner[[i]] <- c( (tail(grp_loc_inner[[i-1]],1)+1) :(tail(grp_loc_inner[[i-1]],1) + width_tmp))
    }  
    
    grp_loc_inner[[ngrp]] <- c( (tail(grp_loc_inner[[ngrp-1]],1)+1) : ngene)
    
    grp_loc_inner
}


##### step 3, get rank for each running bin
get_bin_rank <- function(cor_obs, grp_loc, grp_loc_inner, cor_ref) {
    rank_bin <- rank_bin_pre <- array(dim=dim(cor_obs))
    ngrp <- length(grp_loc)
    size_bin <- length(grp_loc[[1]])
  
    l_cor_tmp_ref <- length(cor_ref[upper.tri(cor_ref)])
  
    for(i in 1:ngrp) {
        for(j in i:ngrp) {
            cor_bin_tmp <- cor_obs[grp_loc[[i]],grp_loc[[j]]]
            rank_bin_tmp <- array(dim=dim(cor_bin_tmp))
            l_cor_tmp <- length(cor_bin_tmp)
      
            rank_bin_tmp[1:l_cor_tmp] <- rank(cor_bin_tmp[1:l_cor_tmp])
      
            ## number of diagonals(of full correlation matrix) contained in the bin
            n_diag <- sum(grp_loc[[i]] %in% grp_loc[[j]])
      
            ## scale the rank of each bin to same scale as rank(cor_ref), 
            ## after scaling, rank_bin could contain non-integars 
            rank_bin_pre[grp_loc[[i]],grp_loc[[j]]] <- rank_bin_tmp
            rank_bin_pre[grp_loc_inner[[i]],grp_loc_inner[[j]]] <- 1 + ((rank_bin_pre[grp_loc_inner[[i]],grp_loc_inner[[j]]]-1)/(l_cor_tmp-n_diag) *(l_cor_tmp_ref-1))
      
            rank_bin[grp_loc_inner[[i]],grp_loc_inner[[j]]] <- rank_bin_pre[grp_loc_inner[[i]],grp_loc_inner[[j]]]
      
        }
    }
    
    rank_bin
}

##### step 4, transform rank to cor_est 
est_cor <- function(rank_bin, cor_ref){
  
    cor_adj <- array(dim=dim(rank_bin))
    cor_ref_sorted <- sort(cor_ref[upper.tri(cor_ref)])
  
    up_tri <- upper.tri(cor_adj)#15s for all genes
  
    ## find to nearest integars to rank_bin, since rank_bin could contain non-integars
    rank_bin <- rank_bin[up_tri]
    rank_bin1 <- rank_bin %/% 1
  
    ## assign weights to each rank according to the distance to rank_bin
    rank_bin_w2  <-  (rank_bin - rank_bin1)
  
    rank_bin2 <- rank_bin1+1
    rank_bin_w1 <- 1-rank_bin_w2
    
    ## find the correlations in the cor_ref corresponding to the two nearest ranks to rank_bin
    ## estimate the correlation using weighted average based on the distance to rank_bin
    cor_adj[up_tri] <- rank_bin_w1*cor_ref_sorted[rank_bin1]+ rank_bin_w2*cor_ref_sorted[rank_bin2]
    
    remove(rank_bin1)
    remove(rank_bin2)
    remove(rank_bin_w1)
    remove(rank_bin_w2)
    remove(up_tri)
    remove(rank_bin)
    
    low_tri <- lower.tri(cor_adj)
    cor_adj[low_tri] <- t(cor_adj)[low_tri]
    diag(cor_adj) <- 1
    cor_adj
}

quantile_norm <- function(cor_mat, ngrp,size_grp, ref_grp){
    group_loc <- get_grps(cor_mat, ngrp, size_grp)
    
    ## step 2, asssign inner bins
    group_loc_adj <- get_grps_inner(group_loc)

    ## step 3, get rank for each running bin
    cor_ref <- cor_mat[group_loc[[ref_grp]], group_loc[[ref_grp]]]
    rank_bin <- get_bin_rank(cor_mat, group_loc, group_loc_adj, cor_ref)
    ##6s for 20*20 bins, 4000 genes; 2min for all genes,ngrp=20,size_grp=1000

    ## step 4, transform rank to cor_adj
    cor_est <- est_cor(rank_bin, cor_ref)
    ##6s for 20*20 bins, 4000 genes; 3.5min for all genes, 20*20 bins,siez_grp=100

    cor_est
}

