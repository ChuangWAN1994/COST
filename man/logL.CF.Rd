\name{logL.CF}
\alias{logL.CF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{negtive log-likelihood for separate time series analysis
%%  ~~function to do ... ~~
}
\description{negtive log-likelihood for separate time series analysis, copula-based semiparametric method from Chen and Fan (2006), assuming t copula for each time series and Markov process of order one, with marginal distribution estimated by espirical CDF, and it is for correlation parameter estimation
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
logL.CF(par,Yk,dfs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{correlation parameter in the t copula function, will be obtained by minimizing the negtive log-likelihood
%%     ~~Describe \code{x} here~~
}
\item{Yk}{observed data from k-th location}
\item{dfs}{degrees of freedom for the t copula, obtained from COST method with t copula}
}

\value{the negative log-likelihood
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{1.Chen, X. and Fan, Y. (2006). Estimation of copula-based semiparametric time series models. Journal of Econometrics 130, 307–335.\\
2.Yanlin Tang, Huixia Judy Wang, Ying Sun, Amanda Hering. Copula-based semiparametric models for spatio-temporal data.
%% ~put references to the literature/web site here ~
}
\author{Yanlin Tang and Huixia Judy Wang
%%  ~~who you are~~
}
