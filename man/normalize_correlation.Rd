\name{normalize_correlation}
\alias{normalize_correlation}
\title{Spatial quantile normalization (SpQN)}
\description{
  This method was developed to remove a mean-correlation relationship in
  correlation matrices built from gene expression data. It can serve as
  pre-processing step prior to a co-expression analysis.
}
\usage{
normalize_correlation(cor_mat, ave_exp, ngrp, size_grp, ref_grp)
}
\arguments{
  \item{cor_mat}{A (square and symmetrix) correlation matrix.}
  \item{ave_exp}{A vector of expression levels, same length as the
    number of rows of the correlation matrix in \code{cor_mat}. For
    other types of data, \code{ave_exp} can be the vector corresponding to the
    row/column of the correlation matrix, whose dependency with the
    distribution of correlations need to be removed.}
  \item{ngrp}{Number of bins in each row/column to be used to
    partition the correlation matrix, integer.}
  \item{size_grp}{Size of the outer bins to be used to appriximate the
    distribution of the inner bins, in order to smooth the
    normalization. Note that the product of size_grp and ngrp must be
    equal or larger than than the row/column number of cor_mat, and
    there is no smoothness in the normalization when they are equal.}
  \item{ref_grp}{Location of the reference bin on the diagonal, whose
    distribution will be used as target distribution in the
    normalization, an integer.}
}
\value{
  A normalized correlation matrix.
}
\examples{
if(require(spqnData)){
  data(gtex.4k)
  cor_ori <- cor(t(assay(gtex.4k)))
  ave_logrpkm <- rowData(gtex.4k)$ave_logrpkm
  normalize_correlation(cor_ori, ave_exp = ave_logrpkm,
                        ngrp=10, size_grp=15, ref_grp=9)}
}
