#   vp = list()
#   for (i in seq_len(C)) {
#     vp[[i]] = viewport(name=sprintf("pos.%d",i), layout.pos.row=1, layout.pos.col=i, clip=FALSE)
#   }

makeSequenceLogoVP <- function(PFM, PFMneg, cex, ylim) {
  C = length(PFM)

#  print("A")
  
  vp = list()
  for (i in seq_along(PFM)) {
    R = length(PFM[[i]])
    if(R > 0) {
      vp2 = list()
      for (j in seq_len(R)) {
        vp2[[j]] = viewport(name=sprintf("aa.%d.%d",i,j), layout.pos.row=j, layout.pos.col=1, clip=FALSE)
      }
      vplStack = do.call(vpList, vp2)
      vp[[length(vp)+1]] = vpTree(viewport(name=sprintf("pos.%d",i), layout.pos.row=1, layout.pos.col=i, clip=FALSE,yscale=c(0.0,ylim[2])),
                       vpList(vpTree(viewport(name=sprintf("aaLayout.%d",i), layout = grid.layout(R,1,
                                                                           widths = unit(rep(1,R),rep("npc",R)),
                                                                           heights = unit(PFM[[i]]/ylim[2],rep("native",R)),
                                                                           just=c("center","bottom"),
                                                                           respect=FALSE)
                       ),vplStack)))
    }
  }

  vpl = do.call(vpList, vp)

#  print("B")
  
  vp = list()
  for (i in seq_along(PFMneg)) {
    R = length(PFMneg[[i]])
    if(R > 0) {
      vp2 = list()
      for (j in seq_len(R)) {
        vp2[[j]] = viewport(name=sprintf("aa.%d.%d",i,j), layout.pos.row=j, layout.pos.col=1, clip=FALSE)
      }
      vplStack = do.call(vpList, vp2)
      vp[[length(vp)+1]] = vpTree(viewport(name=sprintf("pos.%d",i), layout.pos.row=1, layout.pos.col=i, clip=FALSE,yscale=c(ylim[1],0)),
                       vpList(vpTree(viewport(name=sprintf("aaLayout.%d",i), layout = grid.layout(R,1,
                                                                                                  widths = unit(rep(1,R),rep("npc",R)),
                                                                                                  heights = unit(PFMneg[[i]]/(-ylim[1]),rep("native",R)),
                                                                                                  just=c("center","top"),
                                                                                                  respect=FALSE)
                       ),vplStack)))
    }
  }
  
  vplneg = do.call(vpList, vp)
  
#  print("C")
  
  vplogoPos = vpTree(viewport(name="posLayoutPos", layout = grid.layout(1,C,
                                                        widths = unit(rep(1,C),rep("null",C)),
                                                        heights = unit(rep(1,C),rep("npc",C)),
                                                        respect=FALSE),
                           xscale=c(0.5,C+0.5),
                           yscale=ylim,
                           layout.pos.row=1,layout.pos.col=1
                  ),vpl)
  
#  print("D")

  vplogoNeg = vpTree(viewport(name="posLayoutNeg", layout = grid.layout(1,C,
                                                                  widths = unit(rep(1,C),rep("null",C)),
                                                                  heights = unit(rep(1,C),rep("npc",C)),
                                                                  respect=FALSE),
                           xscale=c(0.5,C+0.5),
                           yscale=ylim,
                           layout.pos.row=2,layout.pos.col=1
  ),vplneg)
  
  h = rev(ylim)
  h[2] = -h[2]
  h = h / sum(h)
  vplogoList = vpList(vplogoNeg, vplogoPos)
  if (length(vplneg) == 0) {
    vplogoList = vpList(vplogoPos)
  }
  if (length(vpl) == 0) {
    vplogoList = vpList(vplogoNeg)
  }
  vplogo = vpTree(viewport(name="posLayout", layout = grid.layout(2,1,
                                                                  widths = unit(1,"npc"),
                                                                  heights = unit(h,rep("npc",2)),
                                                                  respect=FALSE),
                            xscale=c(0.5,C+0.5),
                            yscale=ylim,
                            layout.pos.row=2,layout.pos.col=2
  ),vplogoList)
  
  
#  print("E")
  
  vpall = vpTree(viewport(name="seqLogo", layout = grid.layout(3,3,
                                                        widths = unit(c(4*cex,1,1),c("lines","null","lines")),
                                                        heights = unit(c(1,1,4*cex),c("lines","null","lines")),
                                                        respect=FALSE)
                          ),vpList(vplogo))
  vpall
}

makeSequenceLogo <- function(PFM, PFMneg, Plet,cex, ylim, ylab="probability", xlab="position") {
  gl = gList()
#   for (i in seq_len(length(PFM))) {
#     col = rep(c("red","blue","yellow"),length.out=length(PFM))[i]
#     gl = gList(gl, rectGrob(gp=gpar(fill=col),vp=vpPath("posLayout",sprintf("pos.%d",i))) )
#   }
  g=linesGrob(x=c(0,1), 
              y=c(0,0),
              default.units="npc",
              name = "GRID.zeroline",
              vp=vpPath("seqLogo","posLayout","posLayoutPos") )
  gl = gList(gl, g)
  col = c("#000000", "#00811B", "#D00001", "#D00001", "#000000", 
          "#00811B", "#2000C7", "#000000", "#2000C7", "#000000",
          "#000000", "#800080", "#000000", "#800080", "#2000C7",
          "#00811B", "#00811B", "#000000", "#000000", "#00811B")
  names(col) = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                 "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  for (i in seq_len(length(PFM))) {
    for (j in seq_len(length(PFM[[i]]))) {
      if (PFM[[i]][[j]] >= 0.01) {
      g=polygonGrob(x=Plet[[names(PFM[[i]])[j]]][,1], 
                    y=Plet[[names(PFM[[i]])[j]]][,2],
                    default.units="npc",
                    gp=gpar(fill=col[names(PFM[[i]])[j]],col=NA),
                    name = sprintf("GRID.aapicture.%d.%d",i,j),
                    vp=vpPath("seqLogo","posLayout","posLayoutPos",sprintf("pos.%d",i),
                              sprintf("aaLayout.%d",i),sprintf("aa.%d.%d",i,j))
                    )
      gl = gList(gl, g)
      }
    }
  }
  for (i in seq_len(length(PFMneg))) {
    for (j in seq_len(length(PFMneg[[i]]))) {
      if (PFMneg[[i]][[j]] >= 0.01) {
        g=polygonGrob(x=Plet[[names(PFMneg[[i]])[j]]][,1], 
                      y=1-Plet[[names(PFMneg[[i]])[j]]][,2],
                      default.units="npc",
                      gp=gpar(fill=col[names(PFMneg[[i]])[j]],col=NA),
                      name = sprintf("GRID.aapictureneg.%d.%d",i,j),
                      vp=vpPath("seqLogo","posLayout","posLayoutNeg",sprintf("pos.%d",i),
                                sprintf("aaLayout.%d",i),sprintf("aa.%d.%d",i,j))
        )
        gl = gList(gl, g)
      }
    }
  }
  gl = gList(gl, xaxisGrob(at=1:ncol(pfm),label=colnames(pfm),
                           vp=vpPath("seqLogo","posLayout"),name="position",
                           gp=gpar(cex=cex)),
              textGrob(xlab,y=unit(-3,"lines"),
                       vp=vpPath("seqLogo","posLayout"),gp=gpar(cex=cex)),
             yaxisGrob(at=pretty(ylim),label=pretty(ylim),
                       vp=vpPath("seqLogo","posLayout"),name="probability",
                       gp=gpar(cex=cex)),
             textGrob(ylab,x=unit(-3,"lines"),rot=90,
                      vp=vpPath("seqLogo","posLayout"),gp=gpar(cex=cex))
  )
  gl
}

sequenceLogoGrob <- function (pfm, Plet,cex=1.0, ylim, ylab="probability", xlab="position", ...) {
  PFM = list()
  for (i in seq_len(ncol(pfm))) {
    PFM[[i]] = pfm[,i][which(pfm[,i] > 0)]
    PFM[[i]] = sort(PFM[[i]], decreasing=TRUE)
  }
  PFMneg = list()
  for (i in seq_len(ncol(pfm))) {
    PFMneg[[i]] = -pfm[,i][which(pfm[,i] < 0)]
    PFMneg[[i]] = sort(PFMneg[[i]], decreasing=TRUE)
  }

  if (missing(ylim)) {
    ylim = c(-max(sapply(PFMneg, sum, na.rm=TRUE)),max(sapply(PFM, sum, na.rm=TRUE)))
    ylim = range(pretty(ylim))
  }
  
#  name=name, gp=gp, vp=vp,

  igt <- gTree(childrenvp=makeSequenceLogoVP(PFM,PFMneg,cex,ylim,...),
               children=makeSequenceLogo(PFM,PFMneg,Plet,cex, ylim,ylab=ylab,xlab=xlab),
               cl="sequenceLogo")
  igt
}

grid.sequenceLogo <- function(..., draw = TRUE) {
  igt <- sequenceLogoGrob(...)
  if (draw) {
    grid.draw(igt)
  }
  igt
}

