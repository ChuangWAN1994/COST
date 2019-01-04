\name{Forecasts.GP}
\alias{Forecasts.GP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{one-step ahead forecast by Gaussian process fitting
%%  ~~function to do ... ~~
}
\description{one-step ahead forecast by Gaussian process fitting, including: (i) point forecast, either conditional mean; (ii) 95\% forecast intervals, which can also be adjusted by the users; (iii) m (m=500 by default) random draws from the conditional distribution, can be used for multivariate rank
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Forecasts.GP(par,Y,s.ob,seed1,m,isotropic)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{parameters in the copula function
%%     ~~Describe \code{x} here~~
}
\item{Y}{observed data}
\item{s.ob}{coordinates of observed locations}
\item{seed1}{random seed used to generate random draws from the conditional distribution, for reproducibility}
\item{m}{number of random draws to approximate the conditional distribution}
\item{isotropic}{indicator, True for isotropic correlation matrix, False for anisotropic correlation matrix, and we usually choose False for flexibility}
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
