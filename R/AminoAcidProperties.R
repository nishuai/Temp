
AA <- function() {
  res = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", 
         "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  res
}
  
AAproperties <- function() {
  res = list(
   polar = c("R","N","D","C","E","Q","H","K","S","T","Y"),
   positive = c("R","K"),
   negative = c("D","E"),
   hydrophobic = c("V","I","L","M","F","W","C"),
   aromatic = c("F","W","Y","H"),
   aliphatic = c("I","L","V"),
   tiny = c("A","G","S","T","C"),
   small = c("V","P","N","D"))
  res
}
  
AAcomposition <- function(Seq, type="AA") {
  t = table(factor(unlist(strsplit(Seq, split="")), levels=AA()))
  if (type == "prop") {
    t = sapply(AAproperties(), function(aa) { sum(t[aa]) })
  }
  t
}

AAcompositionDifferential <- function(SeqFG, SeqBG, type = "AA") {
  aaFG = AAcomposition(SeqFG)
  aaBG = AAcomposition(SeqBG)
  nFG = sum(aaFG)
  nBG = sum(aaBG)
  if (type == "prop") {
    aaFG = AAcomposition(SeqFG, type="prop")
    aaBG = AAcomposition(SeqBG, type="prop")
  }
  
  pval = oddsRatio = rep(NA, length(aaFG))
  names(pval) = names(oddsRatio) = names(aaFG)
  for (aa in names(aaFG)) {
    M = matrix(NA, nr=2, nc=2)
    M[1,2] = nFG - aaFG[aa]
    M[1,1] = nBG - aaBG[aa]
    M[2,2] = aaFG[aa]
    M[2,1] = aaBG[aa]
    t = fisher.test(M)
    pval[aa]=t$p.value
    oddsRatio[aa]=t$estimate
  }
  padj = p.adjust(pval, method="BH")
  library(RColorBrewer)
  sig = rep("", length(pval))
  names(sig) = names(pval)
  sig[which(padj <= 0.1)] = "*"
  sig[which(padj <= 0.01)] = "**"
  col = rep("gray90", length(pval))
  col[(sig != "") & (oddsRatio > 1.0)] = brewer.pal(3,"Pastel1")[1]
  col[(sig != "") & (oddsRatio < 1.0)] = brewer.pal(3,"Pastel1")[2]
  aaFG = as.vector(aaFG)
  aaBG = as.vector(aaBG)
  df = data.frame(aaFG, aaBG, oddsRatio, pval, padj, sig, col, stringsAsFactors=FALSE)
  df
}
