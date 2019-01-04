\name{Forecasts.CF}
\alias{Forecasts.CF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{one-step ahead forecast by separate time series analysis
%%  ~~function to do ... ~~
}
\description{one-step ahead forecast, analyzing the time series at each location separately with a t copula, including: (i) point forecast, either conditional median or mean; (ii) 95\% forecast intervals, which can also be adjusted by the users; (iii) m (m=500 by default) random draws from the conditional distribution for each location, can be used for multivariate rank after combining all the locations together
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Forecasts.CF(par,Y,seed1,m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{parameters in the copula function
%%     ~~Describe \code{x} here~~
}
\item{Y}{observed data}
\item{seed1}{random seed used to generate random draws from the conditional distribution, for reproducibility}
\item{m}{number of random draws to approximate the conditional distribution}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{y.qq}{0.025-, 0.975- and 0.5-th conditional quantiles of the conditional distribution for each location}
\item{mean.est}{conditional mean estimate for each location}
\item{y.draw.random}{m random draws from the conditional distribution}
%% ...
}
\references{
Yanlin Tang, Huixia Judy Wang, Ying Sun, Amanda Hering. Copula-based semiparametric models for spatio-temporal data.
}
\author{Yanlin Tang and Huixia Judy Wang
%%  ~~who you are~~
}

