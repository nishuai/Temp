
BindingCategory <- function(E, efit, enzyme = "LysC", FDR.RBDpep=0.01,FDR.candidateRBDpep=0.1, FDR.Input = FDR.RBDpep) {
  category = factor(rep("Background", nrow(E)),
                  levels=names(colRBDmap))
  if ("lrF1" %in% names(efit[[enzyme]])) {
    category[which((efit[[enzyme]]$lrF1$padj  <= FDR.Input) &
                 (efit[[enzyme]]$lrF1$coefficients > 0))] = "Input"
  }
  if ("lrF3F2" %in% names(efit[[enzyme]])) {
    category[which((efit[[enzyme]]$lrF3F2$padj  <= FDR.candidateRBDpep) & 
                 (efit[[enzyme]]$lrF3F2$coefficients > 0))] = "CandidateRBDpep"
    category[which((efit[[enzyme]]$lrF3F2$padj  <= FDR.RBDpep) & 
                 (efit[[enzyme]]$lrF3F2$coefficients > 0))] = "RBDpep"
  }

  category
}

