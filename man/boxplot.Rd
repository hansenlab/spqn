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
a=rnorm(100)
IQR_cor_mat=array(a,dim=c(10,10))
grp_mean=c(1:10)
boxplot(IQR_cor_mat,grp_mean)
}
