\name{Data.COST}
\alias{Data.COST}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Data Generation
%%  ~~function to do ... ~~
}
\description{Generating data from COST DGP, assuming Markov process of order one
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Data.COST(n,n.total,seed1,coord,par.t)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{number of time points for parameter estimation
%%     ~~Describe \code{x} here~~
}
\item{n.total}{number of total time points, with a burning sequence}
\item{seed1}{random seed to generate a data set, for reproducibility}
\item{coord}{coordinates of the locations}
\item{par.t}{the true copula parameters}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{Y.all}{data from all locations and time points, may include data at time point n+1, or data from new locations}
\item{mean.true}{true conditional mean of observed locations at time point n+1}
%% ...
}
\references{
Yanlin Tang, Huixia Judy Wang, Ying Sun, Amanda Hering. Copula-based semiparametric models for spatio-temporal data.
}
\author{
Yanlin Tang, Huixia Judy Wang
}

\examples{
library(COST)
n = 500
n.total = 1001
seed1 = 22222
coord = cbind(rep(c(1,3,5)/6,each=3),rep(c(1,3,5)/6,3))
par.t = c(0,1,1,0.5,1.5,100)
dat = Data.COST(n,n.total,seed1,coord,par.t)
#it returns a data set with dimension 501*9
}
