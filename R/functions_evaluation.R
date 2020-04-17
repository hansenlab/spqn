get_grp_loc <- function(cor_matrix, ngrp=10){
    ngene <- nrow(cor_matrix)
    idx <- seq_len(ngene)
    grp_label <- cut(idx, ngrp)
    grp_loc <- split(idx, grp_label)
    return(grp_loc)
}

get_IQR_condition_exp <- function(cor_mat, ave_exp){
    idx <- order(ave_exp)
    cor_mat <- cor_mat[idx, idx]
    ave_exp <- ave_exp[idx]
    grp_loc <- get_grp_loc(cor_mat)
    IQR_cor_mat <- array(dim=c(10,10))
    grp_mean <- array(dim=10)
    for(i in seq_len(10)) {
        cor_tmp <- cor_mat[grp_loc[[i]],grp_loc[[i]]]
        cor_tmp_up_tri <- cor_tmp[upper.tri(cor_tmp)]
        IQR_cor_mat[i,i] <- IQR(cor_tmp_up_tri)
        if(i < 10){ 
            for(j in (i+1):10) {
                cor_tmp <- cor_mat[grp_loc[[i]],grp_loc[[j]]]
                IQR_cor_mat[i,j] <- IQR(cor_tmp)
            }
        }
        ave_exp_grp <- ave_exp[grp_loc[[i]]]
        grp_mean[i] <-  mean(ave_exp_grp)
    }
    idx <- lower.tri(IQR_cor_mat)
    IQR_cor_mat[idx] <- t(IQR_cor_mat)[idx]
    return(list(IQR_cor_mat=IQR_cor_mat,
                grp_mean=grp_mean))
}


qqplot_condition_exp <- function(cor_mat, ave_exp, i, j){
    idx <- order(ave_exp)
    cor_mat <- cor_mat[idx, idx]
    group_loc <- get_grp_loc(cor_mat)
    cor_ref <- cor_mat[group_loc[[9]], group_loc[[9]]]
    cor_ref_vec <- cor_ref[upper.tri(cor_ref)]
    cor_tmp <- cor_mat[group_loc[[i]], group_loc[[j]]]
    if(i==j) {
        cor_tmp <- cor_tmp[upper.tri(cor_tmp)]
    }
    qqplot(cor_ref_vec, as.vector(cor_tmp),
           xlab="reference correlations - (9,9)",
           ylab=paste0("correlations in (",i,",",j,")"),
           cex.lab=1, cex.axis=1, cex=0.2)
    abline(0, 1, col="blue")
}

plot_signal_condition_exp <- function(cor_mat, ave_exp, signal) { 
  idx <- order(ave_exp)
  cor_mat <- cor_mat[idx, idx]
  if(isTRUE(all.equal(signal, 0))){
    ngrp <- 10
    ngene <- ncol(cor_mat)
    idx <- seq_len(ngene)
    grp_label <- cut(idx, ngrp)
    grp_loc <- split(idx, grp_label) 
    
    length_group <- lapply(grp_loc, length)
    length_group <- unlist(length_group)
    length_group <- length_group*(length_group-1)/2
    length_group_cumulate <- cumsum(length_group)
    
    cor_vec_all <- data.frame(matrix(ncol = 2, 
                                     nrow = length_group_cumulate[10]))
    colnames(cor_vec_all) <- c("correlation",  "group")
    
    for(i in seq_len(10)) {
      cor_tmp <- cor_mat[grp_loc[[i]], grp_loc[[i]]]
      cor_tmp <- cor_tmp[upper.tri(cor_tmp)]
      cor_tmp <- as.numeric(cor_tmp)
      if(i==1){
        idx1 <- 1
      }else{
        idx1 <- length_group_cumulate[i-1]+1
      }
      idx2 <- length_group_cumulate[i]
      idx <- c(idx1 : idx2)
      cor_vec_all$correlation[idx] <- cor_tmp
      length_tmp <- length_group[i]
      cor_vec_all$group[idx] <- rep(i,length_tmp)
    }
    cor_vec_all <- data.frame(cor_vec_all)
    names(cor_vec_all) <- c("correlation", "group")
    cor_vec_all$group <- as.factor(cor_vec_all$group)
    ggplot(cor_vec_all, aes_string(x = "correlation", y = "group")) +
      geom_density_ridges2(fill="blue") +
      theme_ridges( grid = TRUE) + theme(axis.title.x = element_blank()) + 
      geom_vline(xintercept = 0, linetype="dotted", color = "black", size=0.3)
  }else{
    ncor <- length(which(upper.tri(cor_mat)))
    nsignal <- ncor*signal
    nsignal_grp <- round(nsignal/100)
    
    # list_cor_sig <- list_cor_back <- grp_back <- grp_sig <- c()
    grp_locs <- get_grp_loc(cor_mat)
    list_cor_sig <- grp_sig <- numeric(10*nsignal_grp)
    nbackground_group <- lapply(grp_locs, length)
    nbackground_group <- unlist(nbackground_group)
    nbackground_group <- nbackground_group*(nbackground_group-1)/2
    list_cor_back_cumulate <- cumsum(nbackground_group)
    list_cor_back <- grp_back <- numeric(list_cor_back_cumulate[10])
    for(ngrp in seq_len(10)) {
      loc_back <- grp_locs[[ngrp]]
      cor_back <- cor_mat[loc_back,loc_back]
      cor_back <- cor_back[upper.tri(cor_back)]
      order_cor_ori_grp <- order(abs(cor_back), decreasing=TRUE)
      idx <- seq_len(nsignal_grp)
      id_signal <- order_cor_ori_grp[idx]
      cor_signal_ori <- cor_back[id_signal]
      
      idx1 <- (ngrp-1)*nsignal_grp+1
      idx2 <- ngrp*nsignal_grp
      idx <- c(idx1 : idx2)
      list_cor_sig[idx] <- cor_signal_ori
      grp_sig[idx] <- ngrp
      
      if(ngrp==1){
        idx <- seq_len(list_cor_back_cumulate[ngrp])
        list_cor_back[idx] <- cor_back
        grp_back[idx] <- ngrp
      }else{
        idx1 <- list_cor_back_cumulate[ngrp-1]+1
        idx2 <- list_cor_back_cumulate[ngrp]
        idx <- c(idx1:idx2)
        list_cor_back[idx] <- cor_back
        grp_back[idx] <- ngrp }                                                                                                       
    }
    df_sig <- data.frame(correlation=list_cor_sig, bin=grp_sig, group="signal")
    df_back <- data.frame(correlation=list_cor_back, 
                          bin=grp_back, group="background")
    df_merge <- rbind(df_sig,df_back)
    df_merge$bin_group  <-  paste(df_merge$bin, df_merge$group)
    ggplot(df_merge, aes_string(y = "bin")) +
      geom_density_ridges(
        aes_string(x = "correlation", fill = "bin_group"),
        alpha = .8, color = "white") +
      labs(x = "correlation", y = "bin")  +
      scale_y_discrete(limits=seq_len(10)) +
      scale_fill_cyclical(
        labels = c("background","signal"),
        values = c("#0000ff", "#ff0000"),
        name = "group", guide = "legend") +
      theme_ridges(grid = FALSE)
  }
}


plot_IQR_condition_exp <- function(IQR_list){
  IQR_cor_mat <- IQR_list$IQR_cor_mat
  grp_mean <- IQR_list$grp_mean
  IQR_cor_mat2 <- round(IQR_cor_mat, 3)
  IQR_cor_mat <- IQR_cor_mat/max(IQR_cor_mat)
  max_IQR <- max(IQR_cor_mat)
  margin <- min(IQR_cor_mat)/10
  sd_grps_offset_y <- t(matrix(rep(max_IQR,100), nrow=10))
  for(i in 2:10){
    sd_grps_offset_y[i,] <- max_IQR + sd_grps_offset_y[i-1,] + margin
  }
  sd_grps_offset_x=t(matrix(rep(max_IQR,100),nrow=10))
  for(i in 2:10){
    sd_grps_offset_x[,i]=max_IQR+sd_grps_offset_y[i-1,]+margin
  }
  sd <- as.numeric(IQR_cor_mat)
  x1 <- as.numeric((-IQR_cor_mat/2)+sd_grps_offset_x)
  x2 <- as.numeric((IQR_cor_mat/2)+sd_grps_offset_x)
  y1 <- as.numeric((-IQR_cor_mat/2)+sd_grps_offset_y)
  y2 <- as.numeric((IQR_cor_mat/2)+sd_grps_offset_y)
  IQR_cor_mat <- round(IQR_cor_mat,3)
  d <- data.frame(x1,x2,y1,y2,sd)
  
  ggplot() + 
    geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), 
              color="black", alpha=0) +
    xlab("average(log2RPKM)")+ylab("average(log2RPKM)") +
    scale_y_discrete(limits=c(seq_len(10))+c(0:9)*min(IQR_cor_mat)/10,
                     labels=round(grp_mean,1)) +
    scale_x_discrete(limits=c(seq_len(10))+c(0:9)*min(IQR_cor_mat)/10,
                     labels=round(grp_mean,1)) + theme_bw() +
    theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),
          axis.text=element_text(size=20))
}
