\name{getComplexity}
\alias{getComplexity}
\title{Local sequence complexity}
\description{Computes the local sequence complexity of a protein sequence.}
\usage{
getComplexity(ProtSeq, w = 10)
}
\arguments{
  \item{ProtSeq}{An AAStringSet object.}
  \item{w}{window size to compute the local complexity}
}
\value{A list of real-valued vectors containing the sequence entropy.}
\author{Bernd Fischer}
\examples{
   library(Biostrings)
   proteins = AAStringSet(c(paste(LETTERS[1:26],collapse=""),
                            paste(rep("A",26),collapse="")))
   getComplexity(proteins)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
