
## library(grid)
## library(gridSVG)
## library(RColorBrewer)

# write a table to SVG: copySubstitute
# writing objects from R to javascript: RJSONIO
plotProtein <- function(pid, Peptides, PF, UniprotSeq, dir=".",prefix=paste("prot_",pid,sep=""), main, profiles, lettercode) {
  t1=Sys.time()
  Peptides$col = colRBDmap[as.character(Peptides$category)]

  maxlr = attr(Peptides, "range")["NA"]
  maxlr2 = attr(Peptides, "range")["max"]
  yline = 0.1
  ydata = 1.5

  I = which((Peptides$ProtID == pid) & (Peptides$Start > 0))
  N = nchar(UniprotSeq[pid])
  if (missing(main)) {
    main = Peptides$Symbol[I[1]]
  }
  if (missing(profiles)) {
     profiles = list()
  }
  profiles_arg = list()
  for (i in seq_along(profiles)) {
    a = list(y = 1, yrange=c(0,1),
             gp = gpar(fill=brewer.pal(3,"Pastel1")[2]),
             gp_bg = gpar(fill="gray95",col=NA))
#     I = names(a)[names(a) %in% names(profiles[[i]])]
#     a[I] = profiles[[i]][I]
    profiles_arg[[i]] = a
  }

  P = which(names(UniprotSeq) == pid)

  pf = PF[ which(PF$UniProtID == pid),]

  if (nrow(pf) > 1) {
    newline = c(TRUE,(pf$Type[2:nrow(pf)] != pf$Type[2:nrow(pf)-1]) | 
      (pf$Start[2:nrow(pf)] <= pf$Stop[2:nrow(pf)-1]))
    Line = cumsum(newline)
    L = Line[length(Line)]
  } else {
    L = nrow(pf)
    Line = rep(1, L)
  }

  pdf(file=file.path(dir,paste(prefix,"pdf",sep=".")),width=N/72+2,height=ydata+L*yline+4*yline)
  yp = unlist(sapply(profiles_arg, function(x) { x$y }))
  vp = viewport(width = (N/72)/(N/72+2),layout=grid.layout(nrow=3+length(profiles)+1,ncol=1,
                                                           heights=c(ydata,yp,4*yline,(length(lettercode)+2)*yline,L*yline)),name = "panelvpall")
  pushViewport(vp)
  vp = viewport(xscale=c(1,N),yscale=c(maxlr-0.1,maxlr2+0.1),layout.pos.row=1,layout.pos.col=1,name = "panelvpdata")
  pushViewport(vp)
  if ("col" %in% colnames(Peptides)) {
    col = Peptides$col    
  } else {
    col = RBDmapCol[as.character(Peptides$category)]
  }
  grid.yaxis()
  grid.lines(c(1,N),rep(0.0,2),default.units="native",gp=gpar(col="black",lwd=1))
  for (i in I) {
    nr = Peptides$lrF3F2[i]
    peptt = sprintf("%s: %s",Peptides$Enzyme[i],Peptides$proteolyticFragment[i])
    grid.lines(c(Peptides$Start[i],Peptides$Stop[i]),
                rep(nr,2),
                default.units="native",
                name=paste("pep",i,sep=""),
                gp=gpar(col=col[i],lwd=3))
    grid.garnish(paste("pep", i, sep=""),
                 onmouseover=paste("pephighlight('", i,".1','",peptt,"')", sep=""),
                 onmouseout=paste('pepdim(', i, '.1)', sep=""))
  }

  upViewport()

  for (i in seq_along(profiles)) {
    vp = viewport(xscale=c(1,N),yscale=profiles_arg[[i]]$yrange,layout.pos.row=1+i,layout.pos.col=1,name = sprintf("panelprofile%d",i))
    pushViewport(vp)
    grid.polygon(c(0,0,N,N),c(0.0,1.0,1.0,0.0),gp = profiles_arg[[i]]$gp_bg, default.units="native")
    grid.polygon(c(0,1:N,N),c(0.0,profiles[[i]],0.0),gp = profiles_arg[[i]]$gp, default.units="native")
    upViewport()
  }

  ############################
  ## domains
  ############################
  vp = viewport(xscale=c(1,N),yscale=c(0,1),layout.pos.row=2+length(profiles),layout.pos.col=1,name = "paneldom")
  pushViewport(vp)
  I = which(pf$Type == "Pfam")
  if (length(I) > 0) {
    for (i in I) {
      grid.polygon(c(pf$Start[i],pf$Start[i],pf$Stop[i],pf$Stop[i]),
                   c(0.1,0.9,0.9,0.1),default.units="native")
      grid.text(label=pf$ShortName[i],(pf$Start[i]+pf$Stop[i])/2,0.5,default.units="native")
    }
  }
  upViewport()

  ##############################
  if (length(lettercode) > 0) {
    vp = viewport(xscale=c(1,N),yscale=c(length(lettercode)+0.5,0.5),layout.pos.row=3+length(profiles),
                                         layout.pos.col=1,name = "panelvplettercode")
    pushViewport(vp)
    for (i in seq_len(length(lettercode))) {
      SP = strsplit(lettercode[[i]],split="")[[1]]
      grid.polyline(x=rep(1:length(SP),each=2),
                    y=rep(rep(i,length(SP)),each=2),
                    id = rep(1:length(SP),each=2),
                    default.units="native",gp=gpar(col=attr(lettercode,"col")[SP]))
    }
    upViewport()
  }

  ################################
  ## interval features
  ################################
  vp = viewport(xscale=c(1,N),yscale=c(L+0.5,0.5),layout.pos.row=4+length(profiles),layout.pos.col=1,name = "panelvpfeatures")
  pushViewport(vp)
  if (length(Line) > 0) {
    grid.polyline(x=rep(c(1,N),each=L),
                  y=rep(seq_len(L),times=2),
                  id = rep(seq_len(L),times=2),
                  default.units="native",gp=gpar(col="gray80"))
    for (i in seq_len(nrow(pf))) {
      peptt = sprintf("[\'%d - %d\',\'%s\',\'%s\',\'%s\']",pf$Start[i],pf$Stop[i],
                      pf$Source[i],pf$Type[i], pf$Name[i])
      grid.lines(x=c(pf$Start[i],pf$Stop[i]),y=c(Line[i],Line[i]),
                 gp=gpar(col=pf$col[i],lwd=3),
                 name=paste("gf",P,i,sep="."),
                 default.units="native")
      grid.garnish(paste("gf",P,i,sep="."),
                   onmouseover=paste("gfhighlight('", P,".",i,".1',",peptt,")", sep=""),
                   onmouseout=paste("gfdim('", P,".",i,".1')", sep=""))
    }

  }
  upViewport()

  upViewport()
  vp2 = viewport(name = "panelvp2")
  pushViewport(vp2)
  grid.text(label=main,x=unit(1,"npc")-unit(0.2,"lines"),y=unit(1,"npc")-unit(0.2,"lines"),just = c("right","top"))
  upViewport()

  ## grid.script(file="highlight.js")
  grid.script(file=system.file("highlight","highlight.js",package="RBDmap"),inline=TRUE)
  grid.export(file.path(dir,paste(prefix,"svg",sep=".")))

  dev.off()
  t2=Sys.time()
#   print(t2-t1)
  invisible(NULL)
}

printProteinRegions <- function(pid, Peptides, PF, UniprotSeq, dir=".",prefix=paste("prot_",pid,sep="")) {
  I = which((Peptides$ProtID == pid) & (Peptides$Start > 0))
  N = nchar(UniprotSeq[pid])
  P = which(names(UniprotSeq) == pid)

  pf = PF[ PF$UniProtID == pid,]
  
  I = I[order(Peptides$Start[I])]
  res = paste("proteolyticFragment",paste(Peptides$Enzyme[I], Peptides$category[I],sep=":"),
              Peptides$Start[I],Peptides$Stop[I],
              Peptides$proteolyticFragment[I],sep="\t")
  res = c(res,"")
  if (nrow(pf) > 0) {
    seq = substr(rep(UniprotSeq[pid],nrow(pf)),pf$Start,pf$Stop)
    res = c(res,
            paste(pf$Type,pf$ShortName,pf$Start,pf$Stop,seq,sep="\t"))
    res = c(res,"")
  }

  writeLines(res,file.path(dir,paste(prefix,"txt",sep=".")))
  writeXStringSet(UniprotSeq[pid],filepath=file.path(dir,paste(prefix,"txt",sep=".")),append=TRUE)
}

orderProteins <- function(ProtID, E, category, ProtRanges, names=TRUE) {
  x = sapply(levels(E$Enzyme), function(enzyme) {
    rowMeans(exprs(E[, (E$Cond == "lrF3F2") & (E$Enzyme == enzyme)]),na.rm=TRUE)
  })
  x[is.na(x)] = -5.0
  x[x < -5] = -5.0
  x = rowMax(x)
  if (missing(ProtRanges)) {
    M = cbind(fData(E)[,c("Sequence", "ProtID")], x=x, category=category)
  } else {  
    M = merge(data.frame(Sequence = fData(E)$Sequence, x=x, category = category, stringsAsFactors=FALSE),
              ProtRanges[,c("Sequence","ProtID")],all.y=TRUE)
  }
  sortKey = tapply(as.integer(M$category)*2*diff(range(x,finite=TRUE))+M$x, M$ProtID, max, na.rm=TRUE)
  I = order(sortKey[ProtID], decreasing=TRUE)
  if (names) {
    I = ProtID[I]
  }
  I
}

writeIndexPage <- function(ProtID, ProtFeatures, fileprefix, title, description,
                           titleSummary = "Protein profiles", K=50,
                           dirname = file.path("result","profileIndex"),
                           dirprofiles = file.path("..","profiles")) {
  if (!is.list(ProtID)) {
    stop("ProtID is expected to be a list of Protein IDs.")
  }
  if (is.null(names(ProtID))) {
    stop("ProtID is expected to have a names attribute.")
  }
  if (length(ProtID) != length(title)) {
    stop("length of ProtID must be same as length of title.")
  }
  if (missing(fileprefix)) {
    stop("fileprefix is missing.")
  }
  if (missing(description)) {
    description = c()
  }
  dir.create(dirname, recursive=TRUE,showWarnings=FALSE)
  if (!file.exists(file.path(dirname,"hwriter.css"))) {
    file.copy(system.file("images","hwriter.css",package="hwriter"),
              file.path(dirname,"hwriter.css"))
  }
  if (!file.exists(file.path(dirname,"highlight.js"))) {
    file.copy(system.file("highlight","highlight.js",package="RBDmap"),
              file.path(dirname,"highlight.js"))
  }
  if (!missing(fileprefix)) {
    filename = sprintf("%sIndex.html",fileprefix)
    page = openPage(filename=filename,dirname=dirname,link.css="hwriter.css")
    hwrite(titleSummary,heading=1,page=page)
    for (i in seq_along(ProtID)) {
      hwrite(title[i], link=sprintf("%sIndex_%s.html",fileprefix,names(ProtID)[i]),page=page,br=TRUE)
    }
    closePage(page, splash=FALSE)
  }

  for (k in seq_along(ProtID)) {
    label = names(ProtID)[k]
    PID = ProtID[[k]]

    m = match(PID, ProtFeatures$names)
    names(m) = PID
    name = ProtFeatures$Symbol[m]
    ensg = ProtFeatures$GeneName[m]
    L = sprintf("<object data=\"%s/prot_%s.svg\" type=\"image/svg+xml\"></object>",dirprofiles,PID)
    names(L) = PID
    H = rep("",length(PID))
    I = which(PID %in% names(description))
    H[I] = description[PID[I]]
    LL = sprintf("<tr><td>%s<br><a href=%s/prot_%s.txt target=text>%s</a><br>%s</td><td><br>%s<br></td></tr>",
                 name,dirprofiles,PID,ensg,H,L)

    N = ceiling(length(LL) / K)
    for (i in 1:N) {
      I = (1:K) + (i-1)*K
      I = I[I %in% seq_along(LL)]
      if (i > 1) {
        prev = sprintf("<a href=%sPage_%s_%d.html>[previous page]</a> ",fileprefix, label, i-1)
      } else {
        prev = "[previous page]"
      }
      if (i < N) {
        np = sprintf(" <a href=%sPage_%s_%d.html>[next page]</a> ",fileprefix, label, i+1)
      } else {
        np = "[next page]"
      }
      overview = sprintf(" <a href=%sIndex_%s.html>[overview]</a> ",fileprefix, label)
      pp = sprintf("%s - %s - %s - %s",prev,
                   paste(sprintf(" <a href=%sPage_%s_%d.html>[%d]</a> ",
                                 fileprefix, label,seq_len(N),seq_len(N)),collapse=" "),
                   overview,np)
      LL2 = c("<html><body>",pp,"<br><br><h3>",title[k]," page ",i,"</h3>",
              "<table>",LL[I],"</table><br><br>",pp,"</body></html>")
      fn = sprintf("%sPage_%s_%d.html",fileprefix, label,i)
      writeLines(LL2,file.path(dirname,fn))
    }

    filename = sprintf("%sIndex_%s.html",fileprefix, label)
    page = openPage(filename=filename,dirname=dirname,link.css="hwriter.css")
    hwrite(title[k], heading=1, page=page)
    hwrite("Protein profiles",heading=3,page=page)
    hwrite(matrix(sprintf("page %d",seq_len(N)),ncol=1),
           link=sprintf("%sPage_%s_%d.html",fileprefix,label,seq_len(N)),page=page)
    hwrite("Protein profiles by gene",heading=3,page=page)
    hn = name
    l = sprintf("%sPage_%s_%d.html",fileprefix,label,ceiling(seq_along(PID) / K))
    I = order(hn)
    R = 8
    M = L = matrix("",nrow=ceiling(length(I)/R),ncol=R)
    M[seq_along(I)] = hn[I]
    L[seq_along(I)] = l[I]
    hwrite(M,link=L,page=page)
    closePage(page,splash=FALSE)
  }
  invisible(NULL)
}
