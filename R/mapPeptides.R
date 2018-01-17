
toKmers <- function(proteins, verbose=TRUE) {
  if (verbose) cat("Indexing. Step 1\n")
  SP1 = strsplit(as.character(subseq(proteins,1,nchar(proteins)-3)),split="")
  SP2 = strsplit(as.character(subseq(proteins,2,nchar(proteins)-2)),split="")
  SP3 = strsplit(as.character(subseq(proteins,3,nchar(proteins)-1)),split="")
  SP4 = strsplit(as.character(subseq(proteins,4,nchar(proteins))),split="")
  if (verbose) cat("Indexing. Step 2\n")
  SP = mapply(function(a,b,c,d) {
    unique(paste(a,b,c,d,sep=""))
  }, SP1, SP2, SP3, SP4)
  n = listLen(SP)
  n = rep(1:length(proteins), n)
  SP = unlist(SP)
  if (verbose) cat("Indexing. Step 3\n")
  L = tapply(X=n, INDEX=SP, function(x) { x } )
  if (verbose) cat("Done.\n")
  L
}

mapPeptides <- function(PeptideSet, ProteinSet, proteinNames=TRUE, Index, verbose=TRUE) {
  if (length(PeptideSet) > 300) {
    useIndex = TRUE
    if (missing(Index)) {
      if (verbose) cat("Indexing protein database.\n")
      h = toKmers(ProteinSet, verbose=verbose)
    } else {
      h = Index
    }
  } else {
    useIndex = FALSE
  }

  if (verbose) cat("\n")
  Pep2Protein = list()
  N = length(PeptideSet)
  for (i in seq_along(PeptideSet)) {
    if (verbose) cat("\rpeptide ",i, " from ",N)
    if (useIndex) {
      K1 = h[[ substr(PeptideSet[i],1,4) ]]
      n = nchar(PeptideSet[i])
      K2 = h[[ substr(PeptideSet[i],n-3,n) ]]
      K = K1[K1 %in% K2]
      I = K[grep(PeptideSet[i],ProteinSet[K])]
    } else {
      I = grep(PeptideSet[i],ProteinSet)
    }
    if (length(I) > 0) {
      Pep2Protein[[i]] = I
    } else {
      Pep2Protein[i] = list(NULL)
    }
  }
  if (verbose) cat("\n")
  
  if (proteinNames) {
    if (verbose) cat("Convert indices to protein names.\n")
    Pep2Protein = lapply(Pep2Protein, function(x) { names(ProteinSet)[x] } )
  }
  
  ## if (useIndex) {
  ##   if (missing(Index)) {
  ##     if (verbose) cat("Release index.\n")
  ##     rm(h)
  ##   }
  ## }
  Pep2Protein
}

