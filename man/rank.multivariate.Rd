\name{rank.multivariate}
\alias{rank.multivariate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{multivariate rank of a vector
%%  ~~function to do ... ~~
}
\description{calculating the multivariate rank of a vector among a set of vectors, used to evaluate the performance of conditional distribution, and the rank would be uniform when the conditional distribution is estimated well
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rank.multivariate(y.test,y.random,seed1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y.test}{the observed (verifying) vector at time n+1
%%     ~~Describe \code{x} here~~
}
\item{y.random}{m random draws from the conditional distribution}
\item{seed1}{random seed to solve tie at random}
}

\value{the multivariate rank of the observed (verifying) vector at time n+1
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
