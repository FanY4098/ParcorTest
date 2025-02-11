% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TestCCA.R
\name{TestCCA}
\alias{TestCCA}
\title{Spectrum Based Test}
\usage{
TestCCA(res1, res2)
}
\arguments{
\item{res1, res2}{Two data matrice with the same sample size. The goal is to test
if the covariance between them are 0 or not.}
}
\value{
Returns an list with \item{tstats}{ The test statistic.}
\item{pval}{The p-value computed based on the test statistic. The reference
distribution is weighted Chi-square. }
\item{weights}{The estimated weights for the weighted Chi-square of reference.}
}
\description{
This function employs the Spectrum based test statistic to infer
if covariance between two random vectors are zero or not.
}
