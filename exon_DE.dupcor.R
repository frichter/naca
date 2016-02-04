

options(stringsAsFactors=FALSE)
setwd("naca/")

# Might have to take blocking into account in voom
# > nf <- calcNormFactors(Counts)
# > design <- model.matrix(~Group)
# > y <- voom(Counts,design,lib.size=colSums(Counts)*nf)
# > corfit <- duplicateCorrelation(y,design,block=Subject)
# > y <- voom(Counts,design,plot=TRUE,lib.size=colSums(Counts)*nf,block=
#               Subject,correlation=corfit$consensus)
# Voom is run twice - the first time to obtain log-counts-per-million
# values and observational level weights to feed into
# duplicateCorrelation(); and a second time with the blocking variable
# and estimated correlation

# explanation on why to use duplicateCorrelation twice
# https://support.bioconductor.org/p/59700/

dupcor.out = duplicateCorrelation(vobj, design, block = info$ID)
fit <- lmFit(vobj, design, block = info$ID, 
             correlation = dupcor.out$consensus.correlation)

#cont.matrix = makeContrasts(HLHSvsOther = L2HLHS, levels = design)
#fit2 = contrasts.fit(fit, cont.matrix)

fit <- eBayes(fit)