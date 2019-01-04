\name{location}
\alias{location}
\docType{data}
\title{Locations of 10 sites
%%   ~~ data name/kind ... ~~
}
\description{Locations of 10 sites.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{
data(location)
}
\format{
  Locations of 10 sites, 10*2 matrix in Cartesian coordinate system
  }
\source{\url{https://transmission.bpa.gov/business/operations/wind/MetData/default.aspx}
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{Yanlin Tang, Huixia Judy Wang, Ying Sun, Amanda Hering. Copula-based semiparametric models for spatio-temporal data.
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
s.ob = location[-3,2:3]
s.new = location[3,2:3]
}
\keyword{datasets}
