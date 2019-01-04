\name{example.prediction}
\alias{example.prediction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{example for new location prediction
%%  ~~function to do ... ~~
}
\description{Example for new location prediction, Gaussian process method, and our COST method with Gaussian and t copulas, where the parameters are assumed to be known; the parameters can be obtained by the ``optim" function. Data are generated at 13 locations and n time points, and assume that 9 locations are observed, and 4 new locations need prediction at time n, conditional on 9 locations at time points n-1 and n.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
example.prediction(n,n.total,seed1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{number of time points for parameter estimation
%%     ~~Describe \code{x} here~~
}
\item{n.total}{number of total time points, with a burning sequence}
\item{seed1}{random seed to generate a data set, for reproducibility}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{COST.t.pre.ECP}{a vector of length K=4 (number of new locations), with value 1 or 0, 1 means the verifying value from the corresponding location lies in the 95\% prediction interval, 0 means not}
\item{COST.t.pre.ML}{a vector of length K=4, each element is the length of prediction interval of the corresponding location}
\item{COST.t.pre.med.error}{prediction error based on conditional median}
\item{COST.G.pre.ECP}{same as COST.t.pre.ECP}
\item{COST.G.pre.ML}{same as COST.t.pre.ML}
\item{COST.G.pre.med.error}{same as COST.t.pre.med.error}
\item{GP.pre.ECP}{same as COST.t.pre.ECP}
\item{GP.pre.ML}{same as COST.t.pre.ML}
\item{GP.pre.med.error}{same as COST.t.pre.med.error}
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
n.total = 101 #number of total time points, including the burning sequence
n = 50 #number of time points we observed
seed1 = 22222
example.prediction(n,n.total,seed1)

#OUTPUTS

# $COST.t.pre.ECP #whether the prediction interval includes the true value, time point n
# [1] 1 1 1 1
#
# $COST.t.pre.ML #length of the prediction interval
# [1] 1.445576 2.146452 2.260688 2.706681
#
# $COST.t.pre.med.error #point prediction error, using conditional median
# [1]  0.01127162 -0.03222058 -0.22081051  0.57831480
#
# $COST.G.pre.ECP #whether the prediction interval includes the true value, time point n
# [1] 1 1 1 1
#
# $COST.G.pre.ML #length of the prediction interval
# [1] 1.445576 2.432646 2.260688 2.914887
#
# $COST.G.pre.med.error #point prediction error, using conditional median
# [1] 0.01127162 -0.03222058 -0.22081051  0.57831480
#
# $GP.pre.ECP #whether the prediction interval includes the true value, time point n
# [1] 1 1 1 1
#
# $GP.pre.ML #length of the prediction interval
# [1] 0.8345359 1.4096642 1.5948724 2.3419428
#
# $GP.pre.med.error #point prediction error, using conditional median
# [1] 0.09447685 -0.05889409 -0.08923935  0.58494684
}

