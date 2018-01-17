
parseTree <- function(T, CLroot, height=0.0,level=0) {
  if (missing(CLroot)) {
    CLroot = unique(T$edge[which( !(T$edge[,1] %in% T$edge[,2]) ),1])
    if (length(CLroot) > 1) {
      stop("Multiple root nodes detected!")
    }
  }
#  cat(rep("-",level),"CLroot = ",CLroot, " h=",height,"\n",collapse="")
  
  if (CLroot <= (T$Nnode+1)) {
    res = CLroot
#     if (T$edge.length[CLroot] < 0) {
#       h = height
#     } else {
      h = height + T$edge.length[CLroot]
#     }
    attributes(res) <- list(label = T$tip.label[CLroot],
                            members = 1L,
                            height = h,
                            leaf = TRUE)
    class(res) = "dendrogram"
  } else {
#     if (T$edge.length[CLroot-1] < 0) {
#       h = height
#     } else {
      h = height + T$edge.length[CLroot-1]
#     }

R = list()
    CL = T$edge[which(T$edge[,1] == CLroot),2]
    for (i in seq_along(CL)) {
      R[[i]] = parseTree(T, CLroot=CL[i], height=h,level=level+1)
    }
    res = do.call(merge, R)
  }
#  cat(rep("-",level),"CLroot = ",CLroot, " n=",attr(res,"members"),"\n",collapse="")
  res
}

