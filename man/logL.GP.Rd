\name{logL.GP}
\alias{logL.GP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{negtive log-likelihood of Gaussian process
%%  ~~function to do ... ~~
}
\description{negtive log-likelihood of Gaussian process, with mean vector and variance vector obtained by the empirical version, and it is for parameter estimation in the correlation matrix
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
logL.GP(par,Y,s.ob)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{parameters in the copula function, will be obtained by minimizing the negtive log-likelihood
%%     ~~Describe \code{x} here~~
}
\item{Y}{the data set from observed locations, used for parameter estimation}
\item{s.ob}{coordinates of observed locations}
}

\value{the negative log-likelihood
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
