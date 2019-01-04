\name{example.forecast}
\alias{example.forecast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{example for one-step ahead forecast
%%  ~~function to do ... ~~
}
\description{Example for one-step ahead forecast for Gaussian Process and our COST method with Gaussian and t copulas, where the data are generated from COST DGP, where the parameters are assumed to be known; the parameters can be obtained by the ``optim" function. Assuming that data are observed at d=9 locations, and n+1 time points, where the last time point is for validation.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
example.forecast(n,n.total,seed1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{n}{number of time points for parameter estimation}
\item{n.total}{number of total time points, with a burning sequence}
\item{seed1}{random seed to generate a data set, for reproducibility}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{COST.t.fore.ECP}{a vector of length d, with value 1 or 0, 1 means the verifying value from the corresponding location lies in the 95\% forecast interval, 0 means not}
  \item{COST.t.fore.ML}{a vector of length d, each element is the length of forecast interval of the corresponding location}
  \item{COST.t.fore.rank}{multivariate rank of the verifying vector by t copula}
  \item{COST.G.fore.ECP}{same as COST.t.fore.ECP}
  \item{COST.G.fore.ML}{same as COST.t.fore.ML}
  \item{COST.G.fore.rank}{multivariate rank of the verifying vector by Gaussian copula}
  \item{GP.fore.ECP}{same as COST.t.fore.ECP}
  \item{GP.fore.ML}{same as COST.t.fore.ML}
  \item{GP.fore.rank}{multivariate rank of the verifying vector by Gaussian process method}
%% ...
}
\references{
Yanlin Tang, Huixia Judy Wang, Ying Sun, Amanda Hering. Copula-based semiparametric models for spatio-temporal data.
}
\author{Yanlin Tang and Huixia Judy Wang
%%  ~~who you are~~
}

\examples{
library(COST)
#settings
seed1 = 2222222
n.total = 101 #number of total time points, including the burning sequence
n = 50 #number of time points we observed
example.forecast(n,n.total,seed1)
#OUTPUTS

# $COST.t.fore.ECP #whether the forecast interval includes the true value at n+1
# [1] 1 1 1 1 1 1 1 1 1
#
# $COST.t.fore.ML #length of the forecast interval
# [1] 0.7036 4.1318 4.8749 2.7615 3.7398 5.8186 4.4532 4.9251 6.3757
#
# $COST.t.fore.rank #multivariate rank
# [1] 162
#
#
# $COST.G.fore.ECP #whether the forecast interval includes the true value at n+1
# [1] 1 1 1 1 1 1 1 1 1
#
# $COST.G.fore.ML #length of the forecast interval
# [1]  0.7035 4.1316 4.8656 2.7611 3.7388 5.7913 4.4458 4.9036 6.3727
#
# $COST.G.fore.rank #multivariate rank
# [1] 186
#

# $GP.fore.ECP #whether the forecast interval includes the true value at n+1
# [1] 1 0 0 1 1 1 1 1 1
#
# $GP.fore.ML #length of the forecast interval
# [1] 0.4879 2.0449 3.4436 2.2107 2.9170 4.4537 4.2169 5.5789 7.3689
#
# $GP.fore.rank #multivariate rank
# [1] 17
}
