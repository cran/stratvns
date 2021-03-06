% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/STRATENUM.R
\name{STRATENUM}
\alias{STRATENUM}
\title{Enumeration Algorithm}
\usage{
STRATENUM(X, L, cvt = 0.1, nhmin = 2)
}
\arguments{
\item{X}{Stratification Variable}

\item{L}{Number of strata}

\item{cvt}{Target cv}

\item{nhmin}{Mininum sample size by stratum}
}
\value{
\item{n}{Sample size}

\item{nh}{Sample size by strata}

\item{cv}{coefficient of variation}

\item{Nh}{Strata sizes}

\item{Vh}{Strata variances}

\item{totoptg}{Total global optimal solutions}

\item{tfeasible}{Total feasible solutions}

\item{cputime}{Runtime in seconds}
}
\description{
This function enumerates all feasible solutions to
the stratification problem and produces the global optimum,
applying an integer formulation proposed by Brito et al (2015).
}
\details{
STRATENUM
}
\examples{
\dontrun{
Example1:
s<-STRATENUM(U21,L=3,cvt=0.05)
Example2:
s<-STRATENUM(U15,L=4)
Example3:
s<-STRATENUM(U1,L=3,nhmin=4)
}
}
\references{
1. Brito, J.A.M., Silva, P.L.N., Semaan, G.S., Maculan, N.,
2015. Integer programming formulations applied to optimal allocation
in stratified sampling. Survey Methodology 41, 2, 427–442.
}
\author{
Leonardo de Lima, Jose Brito, Pedro Gonzalez and Breno Oliveira
}
