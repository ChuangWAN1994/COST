\name{Predictions.COST.G}
\alias{Predictions.COST.G}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{new location prediction by Gaussian copula
%%  ~~function to do ... ~~
}
\description{new location prediction by Gaussian copula, where the copula dimension is extended, and the marginal CDF of the new location is estimated by neighboring information; it gives 0.025-, 0.975- and 0.5-th conditional quantiles of the conditional distribution for each new location, at time n, conditional on observed locations at time n-1 and n; both point and interval predictions are provided
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Predictions.COST.G(par,Y,s.ob,s.new,isotropic)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{parameters in the copula function
%%     ~~Describe \code{x} here~~
}
\item{Y}{observed data}
\item{s.ob}{coordinates of observed locations}
\item{s.new}{coordinates of new locations}
\item{isotropic}{indicator, True for isotropic correlation matrix, False for anisotropic correlation matrix, and we usually choose False for flexibility}
}

\value{0.025-, 0.975- and 0.5-th conditional quantiles of the conditional distribution for each new location, at time n
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Yanlin Tang, Huixia Judy Wang, Ying Sun, Amanda Hering. Copula-based semiparametric models for spatio-temporal data.
%% ~put references to the literature/web site here ~
}
\author{Yanlin Tang and Huixia Judy Wang
%%  ~~who you are~~
}
