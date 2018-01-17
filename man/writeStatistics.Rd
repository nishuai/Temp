\name{writeStatistics}
\alias{writeMainReport}
\alias{writeStatisticsWholeExperiment}
\alias{writeStatistics}
\alias{writeBindingCategory}
\title{Write RBDmap report}
\description{Write RBDmap report.}
\usage{
writeMainReport                (dirname = "results")
writeStatisticsWholeExperiment (E,
                                filename = "indexOverall.html", 
                                dirname = "results")
writeStatistics                (E, efit, cond, label = cond,
                                filename = sprintf("index\%s.html", label),
                                dirname = "results", maxFDR = 0.01)
writeBindingCategory           (E, efit,
                                filename = "indexBindingCategory.html",
                                dirname = "results", maxFDR = 0.01)
}
\arguments{
  \item{E}{}
  \item{efit}{}
  \item{cond}{}
  \item{label}{}
  \item{filename}{}
  \item{dirname}{}
  \item{maxFDR}{}
}
\details{}
\value{}
\author{Bernd Fischer}
\seealso{
\code{\link{getBindingCategory}},\code{\link{RBDmap-package}}
}
\examples{
# writeMainReport()
}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
