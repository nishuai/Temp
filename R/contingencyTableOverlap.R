
contingencyTableOverlap <- function(Peptides, RBDmapFTsel, 
                                    selected = Peptides$category %in% c("CandidateRBDpep","RBDpep"),
                                    columns="ShortName") {
  Pep2PTM = mapPeptide2Domains(PeptideSet=Peptides, RBDmapFTsel,win=0,columns=columns,verbose=FALSE)

  overall = table(factor(listLen(Pep2PTM) > 0,levels=c(FALSE,TRUE)),
                  factor(selected,levels=c(FALSE,TRUE)))
  dimnames(overall) = list(c("released", "bound"), c("hasnotFT","hasFT"))

  testPTMs = unique(unlist(Pep2PTM))
  X = t(sapply(Pep2PTM, function(x) {
    testPTMs %in% x
  }))
  colnames(X) = testPTMs
  FT = t(apply(X, 2, function(x) {
    table(factor(selected,levels=c(FALSE,TRUE)),factor(x,levels=c(FALSE,TRUE)))
  }))
  dim(FT) = c(dim(FT)[1], 2, 2)
  dimnames(FT) = list(testPTMs, c("released", "bound"), c("hasnotFT","hasFT"))

  PID = apply(X[selected,],2,function(x) { unique(Peptides$ProtID[selected][x]) })
  
  res = list(overall = overall, FT=FT, PID=PID)
  res
}

convertName <- function(names, prefix="") {
  for (i in seq_along(names)) {
    name = names[i]
    name=gsub(" ","",name,fixed=TRUE)
    name=gsub("-","",name,fixed=TRUE)
    name=gsub(",","",name,fixed=TRUE)
    name=gsub("/","x",name,fixed=TRUE)
    name=gsub(".","x",name,fixed=TRUE)
    name=gsub("(","",name,fixed=TRUE)
    name=gsub(")","",name,fixed=TRUE)
    name = paste(prefix,name,sep="")
    names[i] = name
  }
  names
}

