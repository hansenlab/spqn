\name{plot_IQR_condition_exp}
\alias{plot_IQR_condition_exp}
\alias{get_IQR_condition_exp}
\title{
  Get and plot the IQRs of submatrices of the correlation matrix.
}
\description{
  The \code{get_IQR_condition_exp} function computes the IQRs of a set
  of 10 by 10 same-size bins that partition the correlation matrix,
  ordered according to expression level.
  
  The \code{plot_IQR_condition_exp} function plots the IQR for each bin
  among a set of 10 by 10 same-size bins that partition the correlation
  matrix, with IQR denoted by the width of boxes in the plot.
}
\usage{
get_IQR_condition_exp(cor_mat, ave_exp)
plot_IQR_condition_exp(IQR_list)
}
\arguments{
\item{cor_mat}{correlation matrix, generated by gene expression matrix,
  with genes sorted by average expression levels.}
\item{ave_exp}{vector, average expression level of each gene for the
  normalized gene expression matrix.}
\item{IQR_list}{List, output of \code{get_IQR_condition_exp}.}
}
\note{
  The mnemonic for \code{condition_exp} is \sQuote{conditional on
    expression}.
}
\value{
  A plot with boxes that shows the IQR of each bin
}
\examples{
if(require(spqnData)) {
    data(gtex.4k)
    cor_mat <- cor(t(assay(gtex.4k)))
    ave_logrpkm <- rowData(gtex.4k)$ave_logrpkm
    IQR_list <- get_IQR_condition_exp(cor_mat, ave_exp = ave_logrpkm)
    plot_IQR_condition_exp(IQR_list)
    }
}
