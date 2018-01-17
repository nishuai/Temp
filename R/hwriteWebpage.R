hwriteStatisticsWholeExperiment <- function(E, filename="indexOverall.html",
                                           dirname="result") {
  dir.create(dirname, showWarnings=FALSE)
  if (!file.exists(file.path(dirname,"hwriter.css"))) {
    file.copy(system.file("images","hwriter.css",package="hwriter"),
              file.path(dirname,"hwriter.css"))
  }
  page = openPage(filename=filename,dirname=dirname,link.css="hwriter.css")
  hwrite("Peptides identified in the whole experiment (F1, F1-noCL, F2, F3, F3-noCL, all replicates, all enzymes) ",heading=1,page=page)
  hwrite(sprintf("In total there are %d peptides identified. Peptides are mapped to the Uniprot proteins and the corresponding genes. Some peptides do not map to a unique gene or the gene id is unknown. For each gene the highest-covered protein is chosen. Some peptides map uniquely to a gene, but not to the selected representative protein.",
                 length(fData(E)$Uniqueness)),page=page, br=TRUE)
  hwrite(as.matrix(table(fData(E)$Uniqueness)), page=page, br=TRUE)
  hwrite(paste("The unique peptides cover ", 
               length(unique(fData(E)$ENSG[!is.na(fData(E)$ENSG) & (fData(E)$Uniqueness == "UniqueGene")])),
               " different genes.",sep=""),
         page=page,br=TRUE)
  hwrite("The following histograms show the uniqueness of the peptide mapping. The x-axis shows the number of genes (proteins) a peptide maps to. The y-axis shows the number of peptides that map to as many genes (proteins).",page=page,br=TRUE)
  hwriteImage(sprintf("PeptideMappingBarplotNumGenes-%d.png",1:2),
              link=sprintf("PeptideMappingBarplotNumGenes-%d.pdf",1:2),
              page=page)
  closePage(page,splash=FALSE)
}

hwriteNormalization <- function(EF3F2,dirname="result") {
  dir.create(dirname, showWarnings=FALSE)
  if (!file.exists(file.path(dirname,"hwriter.css"))) {
    file.copy(system.file("images","hwriter.css",package="hwriter"),
              file.path(dirname,"hwriter.css"))
  }
  n = ncol(EF3F2)
  page = openPage(filename="index.html",dirname=dirname,link.css="hwriter.css")
  hwrite("Normalization",heading=1,page=page)
  hwriteImage(matrix(sprintf("NormalizationShorth-%d.png",seq_len(n)),ncol=1),
              link=matrix(sprintf("NormalizationShorth-%d.pdf",seq_len(n)),ncol=1),
              page=page)
  closePage(page,splash=FALSE)
}

hwriteStatistics <- function(E, efit, label="DPE",name=label,
                            dirname="result",maxFDR=0.01) {
  dir.create(dirname, showWarnings=FALSE)
  if (!file.exists(file.path(dirname,"hwriter.css"))) {
    file.copy(system.file("images","hwriter.css",package="hwriter"),
              file.path(dirname,"hwriter.css"))
  }
  filename=sprintf("index%s.html",label)
  page = openPage(filename=filename,dirname=dirname,link.css="hwriter.css")
  hwrite(sprintf("Peptides with enriched intensity ratio in %s",name),heading=1,page=page)
  
  if (missing(efit)) {
    col = "black"
    warning("No test given for ",name)
  } else {
    col = ifelse((efit$padj <= maxFDR) & (efit$coefficients > 0.0),"red","black")
  }
  N = dim(E)[2]
  if (N == 0) {
    warning("No expression data for ",name," found.")
  } else {  # N > 0
    r = range(exprs(E),finite=TRUE)
    if (N == 1) {
      warning("Only one replicate found for ",name,".")
    } else {  # N > 1
      if (N > 2) {
        hwriteImage(matrix(sprintf("StatisticsScatter%s-%d.png",label,seq_len(N*(N-1)/2)),ncol=1),
                    link=matrix(sprintf("StatisticsScatter%s-%d.pdf",label,seq_len(N*(N-1)/2)),ncol=1),
                    page=page,br=TRUE)
      } else {
        hwriteImage(matrix(sprintf("StatisticsScatter%s.png",label),ncol=1),
                    link=matrix(sprintf("StatisticsScatter%s.pdf",label),ncol=1),
                    page=page,br=TRUE)
      }
      t = table(apply(!is.na(exprs(E)),1,sum))
      t = t[names(t) != "0"]
      hwrite(paste("There are ",paste(t,collapse=",")," peptides with exactly ",
                   paste(names(t),collapse=",")," ratios available.",sep=""),
             page=page,br=TRUE)

      if (!missing(efit)) {
        hwrite(paste("There are ",sum((efit$padj <= maxFDR) &
                                        (efit$coefficients > 0.0),na.rm=TRUE),
                     " peptides enriched in ",name," at an FDR of ",sprintf("%0.2f",maxFDR),"%",sep=""),
               page=page,br=TRUE)
        hwrite(paste(sum((efit$padj  <= maxFDR) &
                           (efit$coefficients > 0.0) &
                           (fData(E)$numENSG == 1),na.rm=TRUE),
                     " of these peptides map uniquely to one gene id.",sep=""),
               page=page,br=TRUE)
        hwrite(paste("These peptides cover ",
                     length(unique(fData(E)$ENSG[which((efit$padj  <= maxFDR) &
                                                         (efit$coefficients > 0.0) &
                                                         (fData(E)$numENSG == 1))])),
                     " different genes.",sep=""),
               page=page,br=TRUE)
        hwrite("<br><br>",page=page)
        Peptides = fData(E)
        Peptides = Peptides[,c("ENSG", "ProtID", "Symbol",
            "Sequence", "Start", "Stop", "Uniqueness")]
        e = exprs(E)
        colnames(e) = paste(E$Cond, E$Experiment, E$Enzyme, sep="_")
        Peptides = cbind(Peptides, e)
        Peptides$log2FC = efit$coefficients[,1]
        Peptides$p.value = efit$p.value
        Peptides$p.adj = efit$padj
        filenameTab=sprintf("resTable%s.txt",label)
        write.table(Peptides, file=file.path(resultdir, filenameTab), sep="\t", quote=FALSE)
        hwrite("[Download table statistics]", link=filenameTab, br=TRUE, page=page)
      }
    }
  }
  closePage(page,splash=FALSE)
  invisible(NULL)
}


hwriteBindingCategory <- function(E, category,
                                 filename="indexBindingCategory.html",
                                 dirname="result") {
  dir.create(dirname, showWarnings=FALSE)
  if (!file.exists(file.path(dirname,"hwriter.css"))) {
    file.copy(system.file("images","hwriter.css",package="hwriter"),
              file.path(dirname,"hwriter.css"))
  }
  page = openPage(filename=filename,dirname=dirname,link.css="hwriter.css")
  hwrite("Binding category of peptides",heading=1,page=page)
  hwrite("The peptides are devided in four groups. The table shows the number of unique peptides.",page=page,br=TRUE)
  I = which(fData(E)$Uniqueness == "UniqueGene")
  M = sapply(category[I,], table)
  hwrite(M,page=page,br=TRUE)
  
  hwrite("The unique mRNA-binding peptides cover the following number of genes.",page=page,br=TRUE)
  M = sapply(category[I,], function(x) {
    length(unique(fData(E)$ENSG[I][(x == "RBDpep") & 
                                     !is.na(fData(E)$ENSG[I])])) })
  hwrite(M,page=page,br=TRUE)

  d = dir(dirname,pattern="StatisticsCategoryPlotsF3input")
  dpng = sort(d[grep("png", d)])
  dpdf = sort(d[grep("pdf", d)])
  hwriteImage(dpng,link=dpdf,page=page,br=TRUE)

  d = dir(dirname,pattern="StatisticsCategoryPlotsEnzyme")

  if (length(d) >= 1) {
    dpng = sort(d[grep("png", d)])
    dpdf = sort(d[grep("pdf", d)])
    ## hwrite("The log2-ratios of F3 / F2 are compared between experimental conditions.",page=page,br=TRUE,heading=1)
    ## M = as.matrix(table(category[ fData(E)$Uniqueness == "UniqueGene",levels(E$Enzyme)[1]],
    ##                     category[ fData(E)$Uniqueness == "UniqueGene",levels(E$Enzyme)[2]])[])
    ## A = hwrite(M)
    ## M = matrix(A, nrow=1, ncol=1)
    ## row.names(M) = levels(E$Enzyme)[1]
    ## colnames(M) = levels(E$Enzyme)[2]
    ## hwrite(M, page=page, br=TRUE)
    M = matrix(hwriteImage(dpng, table=FALSE), nrow=1)
    hwrite(M,link=dpdf,page=page,br=TRUE)
  }

  closePage(page,splash=FALSE)
}


hwriteKnownDomains <- function(E, filename="indexDomainAnalysis.html",
                                dirname="result") {
  dir.create(dirname, showWarnings=FALSE)
  if (!file.exists(file.path(dirname,"hwriter.css"))) {
    file.copy(system.file("images","hwriter.css",package="hwriter"),
              file.path(dirname,"hwriter.css"))
  }
  page = openPage(filename=filename,dirname=dirname,link.css="hwriter.css")
  hwrite("Intensity ratio of peptides for known RNA-binding peptides",heading=1,page=page)

  dD = dir(dirname,pattern="KnownDomainsValidationSetDensity")
  dDpng = sort(dD[grep("png", dD)])
  dDpdf = sort(dD[grep("pdf", dD)])

  dP = dir(dirname,pattern="KnownDomainsValidationSetProbability")
  dPpng = sort(dP[grep("png", dP)])
  dPpdf = sort(dP[grep("pdf", dP)])

  M = matrix(c(dDpng, dPpng), ncol=2)
  L = matrix(c(dDpdf, dPpdf), ncol=2)
  hwriteImage(M, link=L, page=page)
  closePage(page, splash=FALSE)
  invisible(NULL)
}


hwriteDomainEnrichment <- function(resultdir, Fragments, T, FTtest, ProtFeatures) {
  fnpng = sprintf("DomainEnrichmentBarplot-%d.png",seq_along(names(Fragments)))
  fnpdf = sprintf("DomainEnrichmentBarplot-%d.pdf",seq_along(names(Fragments)))
  names(fnpng) = names(fnpdf) = names(Fragments)
  for (enzyme in names(Fragments)) {
    PID = T[[enzyme]]$PID[row.names(FTtest[[enzyme]])]
    names(PID) = convertName(names(PID), prefix="Pfam")
    writeIndexPage(PID, ProtFeatures=ProtFeatures,
                   fileprefix=sprintf("indexPfam%s",enzyme),
                   dirname=resultdir,title=sprintf("Pfam %s",rownames(FTtest[[enzyme]])))
    
    M = as.matrix(FTtest[[enzyme]])
    row.names(M) = 1:nrow(M)
    L = matrix(NA, nrow=nrow(M), ncol=ncol(M))
    L[,1] = sprintf("indexPfam%sIndex_%s.html",enzyme,names(PID))
    page2 = openPage(sprintf("indexPfam%s.html",enzyme),dirname=resultdir,
                     link.css="hwriter.css")
    hwrite("Pfam domain enrichment by Fisher's exact test",heading=1,page=page2)
    #   hwrite("Overall enrichment of PTMs",
    #          heading=3,page=page2)
    #   hwrite(paste(c("<PRE>",capture.output(overAllTest[[enzyme]]),"</PRE>"),
    #                collapse="\n"),page=page2, pre=TRUE,br=TRUE)
    hwrite("In the following barcharts, all under- or overrepresented Pfam domains at an FDR of 0.1 are shown (excluding Pfam-B domains).", page=page2, br=TRUE)
    hwriteImage(fnpng[enzyme], link=fnpdf[enzyme], page=page2, center=TRUE, br=TRUE)
    hwrite(M, link=L,page=page2,row.bgcolor=ifelse(FTtest[[enzyme]]$padj <= 0.1,
                                                   "#FFFFB7",NA))
    closePage(page2,splash=FALSE)
  }
}


hwritePTM <- function(resultdir, Fragments, T, FTtest, overAllTest, ProtFeatures) {
  fnpng = sprintf("PTMbarplot-%d.png",seq_along(names(Fragments)))
  fnpdf = sprintf("PTMbarplot-%d.pdf",seq_along(names(Fragments)))
  names(fnpng) = names(fnpdf) = names(Fragments)
  for (enzyme in names(Fragments)) {
    PID = T[[enzyme]]$PID[row.names(FTtest[[enzyme]])]
    names(PID) = convertName(names(PID), prefix="PTM")
    writeIndexPage(PID, ProtFeatures=ProtFeatures,
                   fileprefix=sprintf("indexPTM%s",enzyme),
                   dirname=resultdir,title=sprintf("PTM %s",rownames(FTtest[[enzyme]])))
    
    M = as.matrix(FTtest[[enzyme]])
    row.names(M) = 1:nrow(M)
    L = matrix(NA, nrow=nrow(M), ncol=ncol(M))
    L[,1] = sprintf("indexPTM%sIndex_%s.html",enzyme,names(PID))
    page2 = openPage(sprintf("indexPTM%s.html",enzyme),dirname=resultdir,
                     link.css="hwriter.css")
    hwrite("PTM enrichment by Fisher's exact test",heading=1,page=page2)
    hwrite("[PTM enrichment profiles]", link=sprintf("PTMprofiles%s.pdf",enzyme), 
           page=page2,br=TRUE)
    hwrite("Overall enrichment of PTMs",
           heading=3,page=page2)
    hwrite(paste(c("<PRE>",capture.output(overAllTest[[enzyme]]),"</PRE>"),
                 collapse="\n"),page=page2, pre=TRUE,br=TRUE)
    hwrite("In the following barcharts, all under- or overrepresented PTMs at an FDR of 0.2 are shown.", page=page2, br=TRUE)
    hwriteImage(fnpng[enzyme], link=fnpdf[enzyme], page=page2, center=TRUE, br=TRUE)
    hwrite(M, link=L,page=page2,row.bgcolor=ifelse(FTtest[[enzyme]]$padj <= 0.1,
                                                   "#FFFFB7",NA))
    closePage(page2,splash=FALSE)
  }
}

hwriteDisease <- function(resultdir, Fragments, T, Tdisease, Tsink, overAllTest, FTtest, FTtestSink, ProtFeatures) {
  for (enzyme in names(Fragments)) {
    page2 = openPage(sprintf("indexDisease%s.html",enzyme),dirname=resultdir,
                     link.css="hwriter.css")
    hwrite("Enrichment of disease associated variants by Fisher's exact test",heading=1,page=page2)
    hwrite("[MIM enrichment profiles] ", link=sprintf("MIMprofiles%s.pdf",enzyme), 
           page=page2,br=FALSE)
    hwrite("[MIM enrichment profiles simple]", link=sprintf("MIMprofilesSimple%s.pdf",enzyme), 
           page=page2,br=TRUE)
    hwrite("Overall enrichment of disease associated variants",
           heading=3,page=page2)
    hwrite(paste(c("<PRE>",capture.output(overAllTest[[enzyme]]),"</PRE>"),
                 collapse="\n"),page=page2, pre=TRUE,br=TRUE)
    
    # Number of peptides covering disease associated variants
    PID = Tdisease[[enzyme]]$PID
    names(PID) = convertName(names(PID), prefix="MIM")
    writeIndexPage(PID, ProtFeatures=ProtFeatures,
                   fileprefix=sprintf("indexDisease%s",enzyme),
                   dirname=resultdir,title=names(Tdisease[[enzyme]]$PID))
    M = cbind(row.names(Tdisease[[enzyme]]$Tab),Tdisease[[enzyme]]$Tab)
    row.names(M) = 1:nrow(M)
    L = matrix(NA, nrow=nrow(M), ncol=ncol(M))
    L[,1] = sprintf("indexDisease%sIndex_%s.html",enzyme,convertName(row.names(Tdisease[[enzyme]]$Tab), prefix="MIM"))
    hwrite("Number of peptides covering disease associated variants",
           heading=3,page=page2)
    hwrite(M, link=L,page=page2)
    
    # Number of peptides covering amino acids as source
    PID = T[[enzyme]]$PID
    names(PID) = convertName(names(PID), prefix="MIMSource")
    writeIndexPage(PID, ProtFeatures=ProtFeatures,
                   fileprefix=sprintf("indexMIMSource%s",enzyme),
                   dirname=resultdir,title=names(T[[enzyme]]$PID))
    M = as.matrix(FTtest[[enzyme]])
    row.names(M) = 1:nrow(M)
    L = matrix(NA, nrow=nrow(M), ncol=ncol(M))
    L[,1] = sprintf("indexMIMSource%sIndex_MIMSource%s.html",enzyme,FTtest[[enzyme]]$name)
    hwrite("Number of peptides covering single amino acids used as source of variation",
           heading=3,page=page2)
    hwrite(M, link=L,page=page2)
    
    # Number of peptides covering amino acids as sink
    PID = Tsink[[enzyme]]$PID
    names(PID) = convertName(names(PID), prefix="MIMSink")
    writeIndexPage(PID, ProtFeatures=ProtFeatures,
                   fileprefix=sprintf("indexMIMSink%s",enzyme),
                   dirname=resultdir,title=names(Tsink[[enzyme]]$PID))
    M = as.matrix(FTtestSink[[enzyme]])
    row.names(M) = 1:nrow(M)
    L = matrix(NA, nrow=nrow(M), ncol=ncol(M))
    L[,1] = sprintf("indexMIMSink%sIndex_MIMSink%s.html",enzyme,FTtestSink[[enzyme]]$name)
    hwrite("Number of peptides covering single amino acids used as sink of variation",
           heading=3,page=page2)
    hwrite(M, link=L,page=page2)
    
    closePage(page2,splash=FALSE)
  }  
}

hwriteKMers <- function(resultdir, Fragments, X, ProtFeatures, binomTest) {
  file.copy(system.file("images","hwriter.css", package="hwriter"),
            file.path(resultdir, "hwriter.css"))
  for (enzyme in names(Fragments)) {
    ProtID = Fragments[[enzyme]]$ProtID[Fragments[[enzyme]]$category != "Input"]
    PID = apply(X[[enzyme]][,Fragments[[enzyme]]$category != "Input"],1,
                function(x) {
                  ProtID[which(x > 0)]
                } )
    writeIndexPage(PID, ProtFeatures=ProtFeatures,
                   fileprefix=sprintf("indexKMers%s",enzyme),
                   dirname=resultdir,title=sprintf("K-mer %s",names(PID)))
    
    M = as.matrix(binomTest[[enzyme]])
    row.names(M) = 1:nrow(binomTest[[enzyme]])
    L = matrix(NA, nrow=nrow(M), ncol=ncol(M))
    L[,1] = sprintf("indexKMers%sIndex_%s.html",enzyme,binomTest[[enzyme]]$name)
    I = which(binomTest[[enzyme]]$p.adj <= 0.1)
    I = I[order(binomTest[[enzyme]]$p.adj[I])]
    M = M[I,]
    L = L[I,]
    page2 = openPage(file.path(resultdir,sprintf("index%s.html", enzyme)),
                     link.css="hwriter.css")
    hwrite("k-mer enrichment by binomial test",heading=1,page=page2)
    #   hwriteImage(matrix(sprintf("hwriteKMersScatter%d.pdf",seq_along(Fragments)),nrow=1),
    #               page=page2)
    hwriteImage(matrix(sprintf("KMersScatterText-%d.png",match(enzyme, names(Fragments))),nrow=1),
                link=matrix(sprintf("KMersScatterText-%d.pdf",match(enzyme, names(Fragments))),nrow=1),
                page=page2)
    #   hwriteImage(matrix(sprintf("hwriteKMersVolcano%d.pdf",seq_along(Fragments)),nrow=1),
    #               page=page2)
    #   hwriteImage(matrix(sprintf("hwriteKMersVolcanoText%d.pdf",seq_along(Fragments)),nrow=1),
    #               page=page2)
    hwrite(M, link=L,page=page2)
    closePage(page2,splash=FALSE)
  }
}

hwriteCV <- function(resultdir) {
  page = openPage(file.path(resultdir, "indexPrediction.html"), link.css="hwriter.css")
  hwrite("Cross validation of binding site prediction", heading=1,page=page)
  L = matrix("", nrow=3, ncol=1)
  L[1,1] = "PredictionPrecisionRecallCurve.png"
  L[2,1] = "PredictionMultiDensity.png"
  L[3,1] = "PredictionVsMeasurement.png"
  M = hwriteImage(L, table=FALSE)
  row.names(M) = c("Precision-recall<br>curve",
                   "density of<br>posterior probability",
                   "prediction compared<br>to<br>mass spec intensities")
  L = matrix("", nrow=3, ncol=1)
  L[1,1] = "PredictionPrecisionRecallCurve.pdf"
  L[2,1] = "PredictionMultiDensity.pdf"
  L[3,1] = "PredictionVsMeasurement.pdf"
  hwrite(M, link=L, page=page, br=TRUE)
  closePage(page, splash=FALSE)  
}
