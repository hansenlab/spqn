\name{plot_IQR_condition_exp}
\alias{plot_IQR_condition_exp}
\title{Plot the IQRs of bins across the correlation matrix}
\usage{plot_IQR_condition_exp(IQR_list)}
\arguments{
\item{IQR_list}{List, output of get_IQR_condition_exp}
}
\description{This function allows users to compare the IQR of correlations conditional on expression matrix.
}
\value{A plot with boxes that shows the IQR of each expression level}
\examples{
a=rnorm(10000)
b=a/max(abs(a))
cor_mat=array(b,dim=c(100,100))
ave_logrpkm=c(1:100)
IQR_list=get_IQR_condition_exp(cor_mat,ave_logrpkm)
plot_IQR_condition_exp(IQR_list)  
}
