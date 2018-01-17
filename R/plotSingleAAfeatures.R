
makeTransparent <- function(col, alpha=0.5) {
  colnew = col
  for (i in seq_along(col)) {
    RGB = col2rgb(col[i])
    RGB = (1-alpha)*col2rgb("white") + alpha*RGB
    colnew[i] = rgb(t(RGB)/255)
  }
  colnew
}

plotSingleAAfeatures <- function(ProtID, Peptides, ProtFeatures, RBDmapFT,
                                 Type = "MOD_RES", column = "ShortName",
                                 main = "Post-translational modifications",
                                 cex = 0.5, xoffset = 0.13, N=1000, alpha=0.35,
                                 size = c(Dom = 1.5, Disorder = 0.5,
                                          Prot = 1.5, PTM = 1.0) * 1.3) {
  RBDmapFTdom = RBDmapFT[ which(RBDmapFT$Type == "Pfam"), ]
  RBDmapFTdisorder = RBDmapFT[ which(RBDmapFT$Type == "disorder"), ]
  RBDmapFT = RBDmapFT[ which(RBDmapFT$Type == Type), ]

  pageWidth = convertY(unit(1,"npc"),"lines",valueOnly=TRUE)/cex
  VPoffset = 5.0
  page=1

  partmp = par(xpd=NA)
  pushViewport(viewport())
  grid.text(label=main,
            y = unit(1,"npc")-unit(2,"lines"),just=c("center","center"),gp=gpar(fontface="bold"))
  grid.text(label=sprintf("page %d",page),
            y = unit(4,"lines"),just=c("center","center"),gp=gpar(cex=0.5))
  popViewport()
  for (protID in ProtID) {
    n = nchar(ProtFeatures$ProtSeq[protID])
    name = ProtFeatures$Symbol[ protID ]
    I = which(RBDmapFT$UniProtID == protID)
    ptms = names(sort(table(RBDmapFT[[column]][I]),decreasing=TRUE))
    S = sum(size)+length(ptms)*size["PTM"]

    if (VPoffset + S + 7 >= pageWidth) {
      grid.newpage()
      pushViewport(viewport())
      page = page+1
      grid.text(label=main,
                y = unit(1,"npc")-unit(2,"lines"),just=c("center","center"),gp=gpar(fontface="bold"))
      grid.text(label=sprintf("page %d",page),
                y = unit(4,"lines"),just=c("center","center"),gp=gpar(cex=0.5))
      popViewport()

      VPoffset = 5.0
    }
    vp = viewport(y = unit(1.0,"npc")-unit(1,"lines")-unit(VPoffset*cex,"lines"),
                  height = unit(S*cex,"lines"),
                  x=unit(xoffset, "npc"),width=unit(1-xoffset,"npc"),
                  xscale=c(0.5,N+0.5),yscale=c(0,1),
                  clip=FALSE, just=c("left","top"))
    pushViewport(vp)
  
    ######################
    ## domains
    ######################
    I = which(RBDmapFTdom$UniProtID == protID)
    for (i in I) {
      grid.rect(x=unit(RBDmapFTdom$Start[i],"native"),
                width=unit(RBDmapFTdom$Stop[i]-RBDmapFTdom$Start[i]+1,"native"),
                y=unit(size["Dom"]/2,"lines"),
                height=unit(size["Dom"],"lines"),
                just=c("left","center"), gp=gpar(cex=cex))
      grid.text(label=RBDmapFTdom$ShortName[i],
                x=unit(RBDmapFTdom$Start[i]+(RBDmapFTdom$Stop[i]-RBDmapFTdom$Start[i])/2,"native"),
                y=unit(size["Dom"]/2,"lines"),just=c("center","center"),gp=gpar(cex=cex))
    }
    offset = size["Dom"]

    ######################
    ## disorders
    ######################
    I = which(RBDmapFTdisorder$UniProtID == protID)
    for (i in I) {
      grid.lines(x=unit(c(RBDmapFTdisorder$Start[i],RBDmapFTdisorder$Stop[i]),"native"),
                 y=unit(rep(size["Disorder"]/2+offset,2),"lines"),
                 gp=gpar(col="blue",lwd=3,cex=cex))
    }
    offset = offset + size["Disorder"]

    ######################
    ## protein
    ######################
    grid.rect(x=unit(0.0, "native"),
              width=unit(n,"native"),
              y=unit(size["Prot"]/2+offset,"lines"),
              height=unit(size["Prot"],"lines"),
              gp=gpar(fill="gray95",cex=cex),just=c("left","center"))
    I = which(Peptides$ProtID == protID)
    for (i in I) {
      col = colRBDmap[as.character(Peptides$category[i])]
      if (Peptides$distToPeptide[i] > 0) {
        col = makeTransparent(col, alpha=alpha)
      }
      if ((Peptides$category[i] == "onlyRBpeptide") | (Peptides$category[i] == "RBpeptideAndInput")) {
        col = ifelse(Peptides$distToPeptide[i] == 0, "tomato4", "tomato1")
      }
      grid.rect(x=unit(Peptides$Start[i],"native"),
                width=unit(Peptides$Stop[i]-Peptides$Start[i]+1,"native"),
                y=unit(size["Prot"]/2+offset,"lines"),
                height=unit(size["Prot"],"lines"),
                just=c("left","center"),
                gp=gpar(fill=col,cex=cex))
    }
    grid.text(label=name,
              x = unit(-0.6,"lines"),
              y = unit((size["Prot"]/2+offset)*cex,"lines"),
              just=c("right","center"))
    offset = offset + size["Prot"]

    ######################
    ## PTMS
    ######################
    I = which(RBDmapFT$UniProtID == protID)
    ptms = names(sort(table(RBDmapFT[[column]][I]),decreasing=TRUE))
    for (i in rev(seq_along(ptms))) {
      J = I[which(RBDmapFT[[column]][I] == ptms[i])]
      grid.polyline(x=unit(rep(RBDmapFT$Start[J],2),"native"),
                    y=unit(rep(c(i*size["PTM"]-1,i*size["PTM"]),each=length(J))+offset,"lines"),
                    id=rep(seq_along(J),2),gp=gpar(cex=cex))
      r = c(min(RBDmapFT$Start[I])-4,max(RBDmapFT$Start[J]))
      grid.lines(x=unit(r,"native"),
                 y=unit(c(i,i)*size["PTM"]+offset,"lines"),gp=gpar(cex=cex))
      grid.text(label=ptms[i],
                x=unit(r[1]-4,"native"),
                y=unit(i*size["PTM"]+offset,"lines"),just=c("right","center"),gp=gpar(cex=cex))
    }
    VPoffset = VPoffset + offset + (2+length(ptms))*size["PTM"]

    popViewport()
  }
  par(partmp)
}
