## library(hwriter)
## library(Biostrings)

# proteolytic enzymes

CleavagePattern = c(
  ArgC =    ".{1}.{1}.{1}[R]{1}.{1}.{1}",
  AspNendopeptidase  = ".{1}.{1}.{1}.{1}[D]{1}.{1}",
  BNPSSkatole = ".{1}.{1}.{1}[W]{1}.{1}.{1}",
  #   Caspase1 = "[FWYL]{1}.{1}[HAT]{1}[D]{1}[^PEDQKR]{1}.{1}",
  #   Caspase2 = "[D]{1}[V]{1}[A]{1}[D]{1}[^PEDQKR]{1}.{1}",
  #   Caspase3 = "[D]{1}[M]{1}[Q]{1}[D]{1}[^PEDQKR]{1}.{1}",
  #   Caspase4 = "[L]{1}[E]{1}[V]{1}[D]{1}[^PEDQKR]{1}.{1}",
  #   Caspase5 = "[LW]{1}[E]{1}[H]{1}[D]{1}.{1}.{1}",
  #   Caspase6 = "[V]{1}[E]{1}[HI]{1}[D]{1}[^PEDQKR]{1}.{1}",
  #   Caspase7 = "[D]{1}[E]{1}[V]{1}[D]{1}[^PEDQKR]{1}.{1}",
  #   Caspase8 = "[IL]{1}[E]{1}[T]{1}[D]{1}[^PEDQKR]{1}.{1}",
  #   Caspase9 = "[L]{1}[E]{1}[H]{1}[D]{1}.{1}.{1}",
  #   Caspase10 = "[I]{1}[E]{1}[A]{1}[D]{1}.{1}.{1}",
  ChymotrypsinHighSpec = ".{1}.{1}.{1}[FY]{1}[^P]{1}.{1}|.{1}.{1}.{1}[W]{1}[^MP]{1}.{1}",
  ChymotrypsinLowSpec = ".{1}.{1}.{1}[FLY]{1}[^P]{1}.{1}|.{1}.{1}.{1}[W]{1}[^MP]{1}.{1}|.{1}.{1}.{1}[M]{1}[^PY]{1}.{1}|.{1}.{1}.{1}[H]{1}[^DMPW]{1}.{1}",
  #   Clostripain = ".{1}.{1}.{1}[R]{1}.{1}.{1}",
  #   CNBr = ".{1}.{1}.{1}[M]{1}.{1}.{1}",
  #   Enterokinase = "[DN]{1}[DN]{1}[DN]{1}[K]{1}.{1}.{1}",
  #   FactorXa = "[AFGILTVM]{1}[DE]{1}[G]{1}[R]{1}.{1}.{1}",
  #   FormicAcid = ".{1}.{1}.{1}[D]{1}.{1}.{1}",
  GlutamylEndopeptidase = ".{1}.{1}.{1}[E]{1}.{1}.{1}",
  #   GranzymeB = "[I]{1}[E]{1}[P]{1}[D]{1}.{1}.{1}",
  Hydroxylamine = ".{1}.{1}.{1}[N]{1}[G]{1}.{1}",
  #   IodosobenzoicAcid = ".{1}.{1}.{1}[W]{1}.{1}.{1}",
  Trypsin = ".{1}.{1}.{1}[KR]{1}.{1}.{1}",
  LysC = ".{1}.{1}.{1}[K]{1}.{1}.{1}",
  #   NTCB = ".{1}.{1}.{1}.{1}[C]{1}.{1}",
  #   PepsinpH13 = ".{1}[^HKR]{1}[^P]{1}[^R]{1}[FLWY]{1}[^P]{1}|.{1}[^HKR]{1}[^P]{1}[FLWY]{1}.{1}[^P]{1}",
  #   PepsinpH20 = ".{1}[^HKR]{1}[^P]{1}[^R]{1}[FL]{1}[^P]{1}|.{1}[^HKR]{1}[^P]{1}[FL]{1}.{1}[^P]{1}",
  ProlineEndopeptidase = ".{1}.{1}[HKR]{1}[P]{1}[^P]{1}.{1}",
  #   ProteinaseK = ".{1}.{1}.{1}[AEFILTVWY]{1}.{1}.{1}",
  #   StaphylococcalPeptidaseI = ".{1}.{1}[^E]{1}[E]{1}.{1}.{1}",
  Thermolysin = ".{1}.{1}.{1}[^DE]{1}[AFILMV]{1}.{1}"
  #   Thrombin = ".{1}.{1}[G]{1}[R]{1}[G]{1}.{1}|[AFGILTVM]{1}[AFGILTVWA]{1}[P]{1}[R]{1}[^DE]{1}[^DE]{1}"
)

CleavageSites = function(ProtSeq, cleavagePattern, offsetCleavageSite=4) {
  L <- gregexpr(pattern=cleavagePattern,text=as.character(ProtSeq))
  M <- sapply(L, function(x) { x[1] } )
  I <- which(M > 0)
  J <- which(M <= 0)
  if (length(I) > 0) {
    L[I] <- lapply(L[I],function(x) {
      x=c(1,x+offsetCleavageSite)
      })
    n = nchar(ProtSeq[I])+1
    K = which(sapply(L[I], function(x) { max(x) }) != n)
    if (length(K) > 0) {
      L[I[K]] = mapply(function(x,y) {
          x = c(x,y)
        }, L[I[K]], n[K])
    }
  }
  if (length(J) > 0) {
    n = nchar(ProtSeq[J])+1
    L[J] = lapply(n, function(x) { c(1,x) })
  }
  names(L) = names(ProtSeq)
  L
}

mapToFragments = function(Peptides, category, cleavageSites, UniprotSeq) {
  Peptides$category = category
  I = which(!is.na(Peptides$ProtID))
  tmp = t(mapply(function(pid, a,b) {
    rs = cleavageSites[[pid]]
    rs1 = rs[1:(length(rs)-1)]
    rs2 = rs[2:length(rs)]-1
    A = min(which((a >= rs1) & (a <= rs2)))
    B = max(which((b >= rs1) & (b <= rs2)))
    res = c(A,B)
    res
  }, Peptides$ProtID[I], Peptides$Start[I], Peptides$Stop[I]))

  Peptides$id = NA
  Peptides$id[I] = sprintf("%s.%d.%d",Peptides$ProtID[I], tmp[,1], tmp[,2])
  Peptides$Start[I] = mapply(function(pid, a) {
    cleavageSites[[pid]][a]
  }, Peptides$ProtID[I], tmp[,1])
  Peptides$Stop[I] = mapply(function(pid, a) {
    cleavageSites[[pid]][a+1]-1
  }, Peptides$ProtID[I], tmp[,2])
  PID = unique(Peptides$ProtID[I])
  CHR = as.character(UniprotSeq[PID])
  CHR = CHR[Peptides$ProtID[I]]
  Peptides$proteolyticFragment = NA
  Peptides$fragmentStart = NA
  Peptides$fragmentStop = NA
  Peptides$distToPeptide = NA
  Peptides$proteolyticFragment[I] = substr(CHR, Peptides$Start[I], Peptides$Stop[I])
  Peptides$fragmentStart[I] = tmp[,1]
  Peptides$fragmentStop[I] = tmp[,2]
  Peptides$distToPeptide[I] = 0

  Peptides
}

splitSpanningFragments = function(Peptides, cleavageSites, UniprotSeq) {
  Peptides = Peptides[which(!is.na(Peptides$ProtID)),]
  I = rep(seq_len(nrow(Peptides)),Peptides$fragmentStop-Peptides$fragmentStart+1)
  tmp2 = unlist(apply(cbind(Peptides$fragmentStart,Peptides$fragmentStop),1,function(x) { x[1]:x[2] } ))
  Peptides = Peptides[I,]
  Peptides$id = sprintf("%s.%d",Peptides$ProtID, tmp2)
  Peptides$Start = mapply(function(pid, a) {
    cleavageSites[[pid]][a]
  }, Peptides$ProtID, tmp2)
  Peptides$Stop = mapply(function(pid, a) {
    cleavageSites[[pid]][a+1]-1
  }, Peptides$ProtID, tmp2)
  Peptides$fragment = tmp2
  Peptides$fragmentStart = Peptides$fragmentStop = NULL
  PID = unique(Peptides$ProtID)
  CHR = as.character(UniprotSeq[PID])
  CHR = CHR[Peptides$ProtID]
  Peptides$proteolyticFragment = substr(CHR, Peptides$Start, Peptides$Stop)
  Peptides
}

# This function maps all peptides to its proteolytics fragment, if the peptide is uniquely mapping to a gene.
# The resulting ExpressionSet contains one row per fragment that is spanned by these peptides. In the case, 
# a peptide spans two or more fragments, all fragments are reported.
## test if pid in names(cleavageSites)
## test if a,b outside of protein length of cleavage sites
# mapToFragments = function(Peptides, catgory, cleavageSites, UniprotSeq) {
#   Peptides$category = catgory
#   Peptides = Peptides[which(!is.na(Peptides$ProtID)),]
#   tmp = t(mapply(function(pid, a,b) {
#     rs = cleavageSites[[pid]]
#     rs1 = rs[1:(length(rs)-1)]
#     rs2 = rs[2:length(rs)]-1
#     A = min(which((a >= rs1) & (a <= rs2)))
#     B = max(which((b >= rs1) & (b <= rs2)))
#     res = c(A,B)
#     res
#   }, Peptides$ProtID, Peptides$Start, Peptides$Stop))
# 
#   I = rep(seq_len(nrow(Peptides)),tmp[,2]-tmp[,1]+1)
#   tmp2 = unlist(apply(tmp,1,function(x) { x[1]:x[2] } ))
#   Peptides = Peptides[I,]
#   Peptides$id = sprintf("%s.%d",Peptides$ProtID, tmp2)
#   Peptides$Start = mapply(function(pid, a) {
#     cleavageSites[[pid]][a]
#   }, Peptides$ProtID, tmp2)
#   Peptides$Stop = mapply(function(pid, a) {
#     cleavageSites[[pid]][a+1]-1
#   }, Peptides$ProtID, tmp2)
#   PID = unique(Peptides$ProtID)
#   CHR = as.character(UniprotSeq[PID])
#   CHR = CHR[Peptides$ProtID]
#   Peptides$proteolyticFragment = substr(CHR, Peptides$Start, Peptides$Stop)
#   Peptides$fragment = tmp2
#   Peptides$distToPeptide = 0
# 
#   Peptides
# }
# 
addNeighbors = function(Fragments, cleavageSites, UniprotSeq, neighbors=1) {
  FragmentsNew = NULL
  N = seq_len(neighbors)
  N = c(-N,N)
  for (n in N) {
    FragmentsTmp = Fragments
    Fr = Fragments$fragment+n
    Nfr = sapply(cleavageSites[Fragments$ProtID], length)-1
    FragmentsTmp = FragmentsTmp[(Fr > 0) & (Fr <= Nfr),]
    Fr = FragmentsTmp$fragment+n
    FragmentsTmp$id = sprintf("%s.%d",FragmentsTmp$ProtID, Fr)
    FragmentsTmp$Start = mapply(function(pid, a) {
      cleavageSites[[pid]][a]
    }, FragmentsTmp$ProtID, Fr)
    FragmentsTmp$Stop = mapply(function(pid, a) {
      cleavageSites[[pid]][a+1]-1
    }, FragmentsTmp$ProtID, Fr)
    
    PID = unique(FragmentsTmp$ProtID)
    CHR = as.character(UniprotSeq[PID])
    CHR = CHR[FragmentsTmp$ProtID]
    FragmentsTmp$proteolyticFragment = substr(CHR, FragmentsTmp$Start, FragmentsTmp$Stop)

    FragmentsTmp$fragment = Fr
    FragmentsTmp$distToPeptide = abs(n)
    FragmentsTmp = FragmentsTmp[ which(!(FragmentsTmp$id %in% Fragments$id)), ]
    FragmentsNew = rbind(FragmentsNew, FragmentsTmp)
  }
  FragmentsNew = rbind(Fragments,FragmentsNew)
  FragmentsNew
}

summarizeFragments <- function(Fragments, useCategories = c("Input", "CandidateRBDpep", "RBDpep")) {

  if (!("category" %in% colnames(Fragments))) {
    error("Column category not found in Fragments.")
  }
  I = which((Fragments$category %in% useCategories) & !is.na(Fragments$ProtID))
  Fragments = Fragments[I,]
  Fragments = Fragments[order( Fragments$category ), ]

  SP = split(Fragments, Fragments$ProtID)
  SP = lapply(SP, function(F) {
    C = as.integer(F$category)
    Cat = rep(-1, max(F$Stop))
    for (i in seq_len(nrow(F))) {
      Cat[F$Start[i]:F$Stop[i]] = C[i]
    }
    RLE = rle(Cat)
    Stop = cumsum(RLE$length)
    Start = c(1,Stop[-length(Stop)]+1)
    category = RLE$values
    df = data.frame(ENSG = F$ENSG[1], ProtID = F$ProtID[1], Symbol = F$Symbol[1], Start=Start, Stop=Stop, category = category, distToPeptide = 0, stringsAsFactors=FALSE)
    df = df[df$category >0,]
    df
  })
  F = do.call(rbind, SP)
  F$category = factor(F$category, levels = seq_len(nlevels(Fragments$category)))
  levels(F$category) = levels(Fragments$category)
  F
}

# 
# summarizeFragments <- function(Fragments) {
#   
#   #  key = paste(Fragments$ProtID, Fragments$fragment, sep="__")
#   key = Fragments$id
#   Key = unique(key)
#   
#   I = match(Key, Fragments$id)
#   FragmentsSummary = Fragments[I,]
#   FragmentsSummary$category = tapply(as.integer(Fragments$category),
#                                      factor(Fragments$id, levels=Key),                                      
#                                      max, na.rm=TRUE)
#   FragmentsSummary$category = factor(FragmentsSummary$category,levels=1:4)
#   levels(FragmentsSummary$category) = levels(Fragments$category)
#   FragmentsSummary
# }
# 
FragmentFeatures <- function(ProtSeq, cleavagePattern, offsetCleavageSite=4) {

  cleavageSites = CleavageSites(ProtSeq=ProtSeq,
                                      cleavagePattern=cleavagePattern, 
                                      offsetCleavageSite=offsetCleavageSite)
  Start = lapply(cleavageSites,function(x) { x[1:(length(x)-1)] })
  Stop = lapply(cleavageSites,function(x) { x[2:length(x)]-1 })
  Fragments = mapply(function(seq, A, B) {
    f = mapply(function(a,b) {
      substr(seq,a,b)
    },A,B)
  }, as.character(ProtSeq), Start, Stop )
  UniProtID = rep(names(cleavageSites), sapply(cleavageSites, length)-1)
  pf = data.frame(UniProtID = UniProtID,
                  Source = "Cleavage",
                  Type = paste(names(cleavagePattern),"Fragment",sep=""),
#                  Name = unlist(Fragments),
                  Name = paste(names(cleavagePattern),"Fragment",sep=""),
                  Start = unlist(Start),
                  Stop = unlist(Stop),
                  ShortName = paste(names(cleavagePattern),"Fragment",sep=""),
                  Info = "",
                  stringsAsFactors = FALSE)
  pf
}

