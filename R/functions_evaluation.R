get_grp_loc <- function(cor_matrix, ngrp=10){
    ngene <- nrow(cor_matrix)
    grp_label <- cut(1:ngene, ngrp)
    grp_loc <- split(1:ngene, grp_label)
    return(grp_loc)
}

get_IQR_condition_exp <- function(cor_matrix, ave_logcpm){
    grp_loc <- get_grp_loc(cor_matrix)
  
    IQR_cor_mat= array(dim=c(10,10))
    grp_mean <- c()
    for(i in 1:10) {
        cor_tmp <- cor_matrix[grp_loc[[i]],grp_loc[[i]]]
        IQR_cor_mat[i,i] <- IQR(cor_tmp[upper.tri(cor_tmp)])
        if(i < 10){ 
            for(j in (i+1):10) {
                cor_tmp <- cor_matrix[grp_loc[[i]],grp_loc[[j]]]
                IQR_cor_mat[i,j] <- IQR(cor_tmp)
            }
        }
        grp_mean <- c(grp_mean, mean(ave_logcpm[grp_loc[[i]]]))
    }
    IQR_cor_mat[lower.tri(IQR_cor_mat)]= t(IQR_cor_mat)[lower.tri(IQR_cor_mat)]
    return(list(IQR_cor_mat=IQR_cor_mat,
                grp_mean=grp_mean))
}


qqplot_condition_exp <- function(cor_matrix, i, j){
    group_loc <- get_grp_loc(cor_matrix)
  
    cor_ref <- cor_matrix[group_loc[[9]], group_loc[[9]]]
    cor_ref_vec <- cor_ref[upper.tri(cor_ref)]
  
    cor_tmp <- cor_matrix[group_loc[[i]], group_loc[[j]]]
    if(i==j) {
        cor_tmp <- cor_tmp[upper.tri(cor_tmp)]
    }
  
    qqplot(cor_ref_vec, as.vector(cor_tmp),
           xlab="reference correlations - (9,9)",
           ylab=paste0("correlations in (",i,",",j,")"),
           cex.lab=1, cex.axis=1, cex=0.2)
    abline(0, 1, col="blue")
}

plot_density_condition_exp <- function(cor_matrix, ngrp=10){
    ngene <- ncol(cor_matrix)
    grp_label <- cut(1:ngene, ngrp)
    grp_loc <- split(1:ngene, grp_label) 
    for(i in 1:10) {
        cor_tmp <- cor_matrix[grp_loc[[i]], grp_loc[[i]]]
        cor_tmp <- cor_tmp[upper.tri(cor_tmp)]
        if(i==1) {
            cor_vec_all <- cbind(as.numeric(cor_tmp), rep(i,length(as.numeric(cor_tmp))))
            names(cor_vec_all) <- c("correlation","group")
        } else {
            cor_vec_all <- rbind(cor_vec_all, cbind(as.numeric(cor_tmp), rep(i,length(as.numeric(cor_tmp)))))
        }
    }
    cor_vec_all <- data.frame(cor_vec_all)
    names(cor_vec_all) <- c("correlation","group")
    cor_vec_all$group <- as.factor(cor_vec_all$group)
  
    ggplot(cor_vec_all, aes(x = `correlation`, y = `group`)) +
        geom_density_ridges2(fill="blue") +
        theme_ridges( grid = TRUE) + theme(axis.title.x = element_blank()) + 
        geom_vline(xintercept = 0, linetype="dotted", color = "black", size=0.3)
}


plot_signal_condition_exp <- function(cor_mat, percent_sig) { 
    ncor <- length(which(upper.tri(cor_mat)))
    nsignal <- ncor*percent_sig/100
    nsignal_grp <- nsignal/100
  
    list_cor_sig <- list_cor_back <- grp_back <- grp_sig <- c()
    grp_locs <- get_grp_loc(cor_mat)
  
    for(ngrp in 1:10) {
        loc_back <- grp_locs[[ngrp]]
        cor_back <- cor_mat[loc_back,loc_back]
        cor_back <- cor_back[upper.tri(cor_back)]
        order_cor_ori_grp <- order(abs(cor_back), decreasing=TRUE)
        cor_signal_ori <- cor_back[order_cor_ori_grp[1:nsignal_grp]]
    
        list_cor_sig <- c(list_cor_sig,cor_signal_ori)
        list_cor_back <- c(list_cor_back, cor_back)
        grp_back <- c(grp_back, rep(ngrp, length(cor_back)))
        grp_sig <- c(grp_sig, rep(ngrp, length(cor_signal_ori)))
    }
  
    df_sig <- data.frame(correlation=list_cor_sig, bin=grp_sig, group="signal")
    df_back <- data.frame(correlation=list_cor_back, bin=grp_back, group="background")
    
    df_merge <- rbind(df_sig,df_back)
  
    ggplot( df_merge,aes(y = bin)) +
        geom_density_ridges(
            aes(x = correlation, fill = paste(bin, group)),
            alpha = .8, color = "white") +
        labs(x = "correlation",
             y = "bin" )  +
        scale_y_discrete(limits=c(1:10)) +
        scale_fill_cyclical(
            labels = c("background","signal"),
            values = c("#0000ff", "#ff0000"),
            name = "group", guide = "legend") +
        theme_ridges(grid = FALSE)
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
