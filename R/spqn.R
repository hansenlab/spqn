.get_grps <- function(cor_mat, ngrp=20, size_grp=400){
    ngene <- nrow(cor_mat)
    grp_label <- cut(1:(ngene-size_grp+1), ngrp-1)
    grp_loc0 <- split(1:(ngene-size_grp+1), grp_label)
    grp_loc <- list()

    if(size_grp-length(grp_loc0[[1]])<5){
        grp_label <- cut(1:ngene, ngrp)
        grp_loc0 <- split(1:ngene, grp_label)
        for(i in 1:(ngrp)){
            grp_loc[[i]] <- grp_loc0[[i]]
        }
    }else{
        grp_label <- cut(1:(ngene-size_grp+1), ngrp-1)
        grp_loc0 <- split(1:(ngene-size_grp+1), grp_label)
        for(i in 1:(ngrp-1)){
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
    grp_loc_inner[[1]] <- c(1:round(size_bin/2+width_tmp/2))
    
    for(i in 2:(ngrp-1)){
        width_tmp <- grp_loc[[i+1]][1] - grp_loc[[i]][1]
        grp_loc_inner[[i]] <- c( (tail(grp_loc_inner[[i-1]],1)+1) :(tail(grp_loc_inner[[i-1]],1) + width_tmp))
    }
    
    grp_loc_inner[[ngrp]] <- c( (tail(grp_loc_inner[[ngrp-1]],1)+1) : ngene)
    
    grp_loc_inner
}


##### Get rank for each running bin
.get_bin_rank <- function(cor_obs,grp_loc,grp_loc_inner,cor_ref){
    rank_bin <- rank_bin_pre <- array(dim=dim(cor_obs))
    ngrp <- length(grp_loc)
    size_bin <- length(grp_loc[[1]])
    
    l_cor_tmp_ref <- length(cor_ref[upper.tri(cor_ref)])
    
    for(i in 1:ngrp){
        for(j in i:ngrp){
            cor_bin_tmp <- cor_obs[grp_loc[[i]],grp_loc[[j]]]
            rank_bin_tmp <- array(dim=dim(cor_bin_tmp))
            l_cor_tmp <- length(cor_bin_tmp)
            
            rank_bin_tmp[1:l_cor_tmp] <- rank(cor_bin_tmp[1:l_cor_tmp])
            
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
.est_cor<-function(rank_bin, cor_ref){
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
    cor_adj[up_tri] <- rank_bin_w1*cor_ref_sorted[rank_bin1]+ rank_bin_w2*cor_ref_sorted[rank_bin2]#1.3min for all genes
    
    remove(rank_bin1)
    remove(rank_bin2)
    remove(rank_bin_w1)
    remove(rank_bin_w2)
    remove(up_tri)
    remove(rank_bin)

    low_tri <- upper.tri(cor_adj)
    cor_adj[low_tri] <- t(cor_adj)[low_tri]
    diag(cor_adj) <- 1
    cor_adj
}


normalize_correlation <- function(cor_mat, ngrp, size_grp, ref_grp){
    group_loc <- .get_grps(cor_mat, ngrp, size_grp)
    
    ## Asssign inner bins
    group_loc_adj <- .get_grps_inner(group_loc)
    
    ## Get rank for each running bin
    cor_ref <- cor_mat[group_loc[[ref_grp]], group_loc[[ref_grp]]]
    rank_bin <- .get_bin_rank(cor_mat, group_loc, group_loc_adj, cor_ref)
    
    ## Transform rank to cor_adj
    cor_est <- .est_cor(rank_bin, cor_ref)
    
    cor_est
}

boxplot<-function(group_mat,grp_medians){
  group_mat2=round(group_mat,3)
  group_mat=group_mat/max(group_mat)
  max_sd=max(group_mat)
  sd_grps_offset_y=t(matrix(rep(max_sd,100),nrow=10))
  for(i in 2:10){
    sd_grps_offset_y[i,]=rep(max_sd,10)+sd_grps_offset_y[i-1,]+min(group_mat)/10
  }
  
  sd_grps_offset_x=t(matrix(rep(max_sd,100),nrow=10))
  for(i in 2:10){
    sd_grps_offset_x[,i]=t(rep(max_sd,10)+sd_grps_offset_y[i-1,])+min(group_mat)/10
  }
  sd=as.numeric(group_mat)
  x1=as.numeric((-group_mat/2)+sd_grps_offset_x)#+c(0:9)*min(group_mat)/10
  x2=as.numeric((group_mat/2)+sd_grps_offset_x)
  y1=as.numeric((-group_mat/2)+sd_grps_offset_y)
  y2=as.numeric((group_mat/2)+sd_grps_offset_y)
  
  group_mat=round(group_mat,3)
  d=data.frame(x1,x2,y1,y2,sd)
  ggplot() + 
    geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="black", alpha=0) +
    xlab("average(log2RPKM)")+ylab("average(log2RPKM)")   +
    scale_y_discrete(limits=c(1:10)+c(0:9)*min(group_mat)/10,
                     labels=grp_medians)+
    scale_x_discrete(limits=c(1:10)+c(0:9)*min(group_mat)/10,
                     labels=grp_medians) + theme_bw()+
    theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text=element_text(size=20))
}
