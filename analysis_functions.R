ObtainLowNacaSampleNames = function() {
>>>>>>> c4a81a66080eca2e76d5e9d4e3ab303fc93a4d5e
  # Return sample IDs with missing exon 6 in NACA
  
  vobj = readRDS("expression_data/matrix.exon.RDS")
  info = readRDS("expression_data/info.exon.RDS")
  fit = lmFit(vobj, model.matrix(~ as.factor(Batch), info))
  R = residuals(fit, vobj)
  datExpr.exon = as.matrix(t(R))
  colnames(datExpr.exon) = rownames(R)
  rownames(datExpr.exon) = colnames(R)
  samples.low.naca = rownames(datExpr.exon[datExpr.exon[, "NACA.6"] < -4, ])
  return(samples.low.naca)
}
