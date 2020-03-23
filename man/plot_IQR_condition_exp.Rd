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
if(require(spqnData)) {
  data(gtex.4k)
  cor_mat <- cor(t(assay(gtex.4k)))
  ave_logrpkm <- rowData(gtex.4k)$ave_logrpkm
  IQR_list=get_IQR_condition_exp(cor_mat, ave_logrpkm)
  plot_IQR_condition_exp(IQR_list)}
}
