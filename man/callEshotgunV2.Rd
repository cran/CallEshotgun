% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CallEshotgun.R
\name{callEshotgunV2}
\alias{callEshotgunV2}
\title{Call the e-shotgun version 2 with random selection function in python}
\usage{
callEshotgunV2(Xtr, Ytr, f_lb, f_ub, q = 10L, epsilon = 0.1, pf = FALSE)
}
\arguments{
\item{Xtr}{a matrix containing the initial points}

\item{Ytr}{a matrix containing the evaluation of Xtr with a given function}

\item{f_lb}{a vector with the values of the lower bounds}

\item{f_ub}{a vector with the values of the upper bounds}

\item{q}{the amount if points that the e-shotgun should evaluate}

\item{epsilon}{the epsilon value of the e-shotgun}

\item{pf}{boolean that decides if pareto front is used}
}
\value{
a matrix or a vector
}
\description{
The function checks the passed parameter and than calls the e-shotgun Python
implementation and returns a matrix with the evaluated points.
}
\details{
George De Ath, Richard M. Everson, Jonathan E. Fieldsend, and Alma A. M. Rahat. 2020.
e-shotgun : e-greedy Batch Bayesian Optimisation. In Genetic and Evolutionary Computation Conference (GECCO ’20), July 8–12, 2020, Cancún, Mexico.
ACM, New York, NY, USA, 9 pages.
https://doi.org/10.1145/3377930.3390154
https://github.com/georgedeath/eshotgun
}
