\name{plot_signal_condition_exp}
\alias{plot_signal_condition_exp}
\title{Plot the signal and background distribution}
\usage{
plot_signal_condition_exp(cor_mat, percent_sig) 
}
\arguments{
\item{cor_mat}{Matrix, correlation matrix, generated by gene expression matrix, with genes sorted by average expression levels}

\item{percent_sig}{Numeric, percentage of signals in each group of the correlation matrix} 
}
\value{A ridge plot visualizing the distribution of background and signals in each group of correlation matrix}
\description{
This function allows users to visualize the separation between signals and background, conditional on expression levels.
}
\examples{
a=rnorm(10000)
b=a/max(abs(a))
cor_ori=array(b,dim=c(100,100))
plot_signal_condition_exp(cor_ori,0.1)
}