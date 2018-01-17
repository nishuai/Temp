
mapPeptide2Domains <- function(PeptideSet, RBDmapFT, win=0, columns=c("Type","ShortName"),verbose=TRUE) {
  PID = unique(PeptideSet$ProtID[!is.na(PeptideSet$ProtID)])
  FT = split(RBDmapFT,factor(RBDmapFT$UniProtID, levels=PID))

  Pep2domain = list()
  if (verbose) cat("\n")
  for (i in seq_len(nrow(PeptideSet))) {
    if (verbose) cat("\rcollect domains on peptide ",i," from ",nrow(PeptideSet))
    id = PeptideSet$ProtID[i]
    if (!is.na(id)) {
      a = PeptideSet$Start[i]
      b = PeptideSet$Stop[i]
      J = which((FT[[id]]$Start <= (PeptideSet$Stop[i]+win)) &
                  ((PeptideSet$Start[i]-win) <= FT[[id]]$Stop))
      Pep2domain[i] = list(apply(as.matrix(FT[[id]][,columns])[J,,drop=FALSE],1,paste,collapse="::"))
    } else {
      Pep2domain[i] = list(NULL)
    }
  }
  if (verbose) cat("\n")
  Pep2domain
}

getKnownRBD <- function() {
  data(RNA_binding_domains, package="RBDmap")
  classicalDomains = RNA_binding_domains$pFamName[ RNA_binding_domains$Type == "Classical"]
  classicalDomains = c(classicalDomains,"KH")
  classicalDomains = c(classicalDomains, "Ribosomal_S6e", "Ribosomal_S16", "Ribosomal_S2", "Ribosomal_S4e", 
                       "Ribosomal_L23", "Ribosomal_L37ae", "Ribosomal_L16", "Ribosomal_L30_N", 
                       "Ribosomal_L30", "Ribosomal_L23eN", "Ribosomal_S4", "Ribosomal_L7Ae", 
                       "Ribosomal_L27", "Ribosomal_L12", "Ribosomal_S10", "Ribosomal_S7e", 
                       "Ribosomal_L3", "Ribosomal_L50", "Ribosomal_S21e", "Ribosomal_L31e", 
                       "Ribosomal_L20", "Ribosomal_L14", "Ribosomal_L2", "Ribosomal_L2_C", 
                       "Ribosomal_L29e", "Ribosomal_S27e", "Ribosomal_L24e", "Ribosomal_L22e", 
                       "Ribosomal_L35Ae", "Ribosomal_L32e", "Ribosomal_L35p", "Ribosomal_S3Ae", 
                       "Ribosomal_L37e", "Ribosomal_L6", "Ribosomal_S18", "Ribosomal_L13", 
                       "Ribosomal_S11", "Ribosomal_L18e", "Ribosomal_L15e", "Ribosomal_L14e", 
                       "Ribosomal_L22", "Ribosomal_S24e", "Ribosomal_S17", "Ribosomal_S3_C", 
                       "Ribosomal_L17", "Ribosomal_S5", "Ribosomal_S5_C", "Ribosomal_S30", 
                       "Ribosomal_L11_N", "Ribosomal_S13_N", "Ribosomal_L29", "Ribosomal_L21p", 
                       "Ribosomal_L10", "Ribosomal_L6e_N", "Ribosomal_60s", "Ribosomal_L21e", 
                       "Ribosomal_L44", "Ribosomal_L1", "Ribosomal_L28e", "Ribosomal_L11", 
                       "Ribosomal_S17e", "Ribosomal_L4", "Ribosomal_S8", "Ribosomal_L13e", 
                       "Ribosomal_S13", "Ribosomal_S15", "Ribosomal_S7", "Ribosomal_L38e", 
                       "Ribosomal_L36e", "Ribosomal_L19e", "Ribosomal_S19", "Ribosomal_L27e", 
                       "Ribosomal_S9", "Ribosomal_S19e", "Ribosomal_L18ae", "Ribosomal_L40e", 
                       "Ribosomal_L34", "Ribosomal_S14", "Ribosomal_L33", "Ribosomal_S8e", 
                       "Ribosomal_L18p", "Ribosomal_L18_c", "Ribosomal_L34e", "Ribosomal_L19", 
                       "Ribosomal_S25", "Ribosomal_S26e", "Ribosomal_S28e", "Ribosomal_L39", 
                       "Ribosomal_L5", "Ribosomal_L5_C", "Ribosomal_L41", "Ribosomal_S27", 
                       "Ribosomal_S21", "Ribosomal_S6", "Ribosomal_L6e", "Ribosomal_L9_N", 
                       "Ribosomal_L37", "Ribosomal_L32p", "Ribosomal_L36", "SAM", "Ribosomal")
  nonclassicalDomains = RNA_binding_domains$pFamName[ RNA_binding_domains$Type == "Non-Canonical"]
  knownRBD = c(classicalDomains,nonclassicalDomains)
  knownRBD
}

mapPeptide2DomainClass <- function(Pep2domain) {
  knownRBD = getKnownRBD()
  domain = rep("unknownRBD", length(Pep2domain))
  domain[sapply(Pep2domain,function(x) { any(x %in% paste("Pfam",knownRBD,sep="::")) })] = "knownRBD"
  domain = factor(domain, levels = c("unknownRBD", "knownRBD"))
  domain
}

