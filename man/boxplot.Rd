\name{boxplot}
\alias{boxplot}
\title{Plot the IQRs of bins across the correlation matrix}
\usage{boxplot(IQR_cor_mat)}
\arguments{
\item{IQR_cor_mat}{Matrix, IQR of bins across the correlation matrix}
}
\value{}
\description{This function allows users to compare the IQR of correlations conditional on expression matrix.
}
\examples{
a=rnorm(10000)
b=a/max(abs(a))
cor_ori=array(b,dim=c(100,100))
IQR_list=get_IQR_condition_exp(cor_m, c(1:100))
boxplot(IQR_list$IQR_cor_mat)
}
