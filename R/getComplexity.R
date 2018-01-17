
getComplexity <- function(ProtSeq, w=10) {
  AA = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", 
         "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  PS = strsplit(as.character(ProtSeq),split="")
  PS2 = lapply(PS, function(x) { c(rep("X",w),x,rep("X",w)) } )

  f = match(unlist(PS),AA,nomatch=0)
  q = table(factor(f, levels=1:length(AA)))
  q = q / sum(q)
  
  f = match(unlist(PS2),AA,nomatch=0)
  nfactor = filter(f > 0,rep(1.0,2*w+1))+1
  nfactor[is.na(nfactor)] = 1
  h = rep(0.0, length(f))
  for (aa in 1:length(AA)) {
    p = (filter(f==aa,rep(1.0,2*w+1)) + q[aa]) / nfactor
    h = h - p*log2(p)
  }
  H = relist(flesh=h, skeleton=PS2)
  H = lapply(H, function(x) { x[(w+1):(length(x)-w)] } )
  H
}
