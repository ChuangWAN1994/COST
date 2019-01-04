\name{Wind6month}
\alias{Wind6month}
\docType{data}
\title{Wind speed data from 10 sites
%%   ~~ data name/kind ... ~~
}
\description{The data set is a subset of the data we used in the paper, with 10 sites and 6-month long time series.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{
data(Wind6month)
}
\format{
  A 4320*10 matrix from 10 locations, date ranges from Sep 22, 2014 to Dec 20, 2014, 180 days
  \describe{
    \item{\code{BiddleButte}}{wind speed from site BiddleButte}
    \item{\code{ForestGrove}}{wind speed from site ForestGrove}
    \item{\code{HoodRiver}}{wind speed from site HoodRiver}
    \item{\code{HorseHeaven}}{wind speed from site HorseHeaven}
    \item{\code{Megler}}{wind speed from site Megler}
    \item{\code{NaselleRidge}}{wind speed from site NaselleRidge}
    \item{\code{Roosevelt}}{wind speed from site Roosevelt}
    \item{\code{Shaniko}}{wind speed from site Shaniko}
    \item{\code{Sunnyside}}{wind speed from site Sunnyside}
    \item{\code{Tillamook}}{wind speed from site Tillamook}
  }
  }
\source{\url{https://transmission.bpa.gov/business/operations/wind/MetData/default.aspx}
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{Yanlin Tang, Huixia Judy Wang, Ying Sun, Amanda Hering. Copula-based semiparametric models for spatio-temporal data.
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(Wind6month)
Y.ob = Wind6month[,-3]
Y.newloc = Wind6month[,3]
dim(Y.ob) #4320*9, data at 9 locations, with length 4320 (hours)
length(Y.newloc) #4320, time series at the new location
}
\keyword{datasets}
