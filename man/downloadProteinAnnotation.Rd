\name{downloadProteinAnnotation}
\alias{downloadProteinAnnotation}
\title{Download protein annotation}
\description{Download a whole organism protein sequence database and sequence features from Interpro, Uniprot, and compute disorderdness and sequence complexity.}
\usage{
downloadProteinAnnotation(organism = "HUMAN",
                          maindatadir = ".",
                          ProtFeaturesFile = "ProtFeatures.rda",
                          RBDmapFTfile = "RBDmapFT.rda",
                          UniprotVersionFile = "UniprotVersion.rda",
                          idmappingfile,
                          fastafile)
}
\arguments{
  \item{organism}{The name of the organism. One out of ARATH, CAEEL, CHICK, DANRE, DROME, HUMAN, RAT, YEAST}
  \item{maindatadir}{The directory for output files.}
  \item{ProtFeaturesFile}{Output file for protein features used for mRNA-interactome capture.}
  \item{RBDmapFTfile}{Output file for protein features with protein position used by RBDmap.}
  \item{UniprotVersionFile}{Outputfile containing version information of Uniprot and Interpro.}
  \item{idmappingfile}{Unused. Will be filled according to organism.}
  \item{fastafile}{Unused. Will be filled according to organism.}
}
\details{Protein sequence database is downloaded from Uniprot as well as Uniprot sequence features and Interpro domains. Disorderdned is computed by iupred and sequence complexity by the function \code{\link{getComplexity}}.}
\section{Warning}{The download can take up to several hours to be completed.}
\value{Nothing is returned, but the function as a site effect the function saves the three files \code{ProtFeaturesFile}, \code{RBDmapFTfile} and \code{UniprotVersionFile}.}
\author{Bernd Fischer}
\examples{
## downloadProteinAnnotation("YEAST")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
