

############################
# obtain subset sample names
############################
library("WGCNA")
library("limma")
options(stringsAsFactors=FALSE)
setwd("~/Documents/phd/variancepartitionpcgc")

vobj = readRDS("matrix.exon.RDS")
info = readRDS("info.exon.RDS")
info$L2 <- as.character(info$L2)
info[info$L2 %in% c("Fallot"), ]$L2 = "TOF"
info$L2 <- droplevels(factor(info$L2))
fit = lmFit(vobj, model.matrix(~ as.factor(Batch), info))
R = residuals(fit, vobj)
datExpr.exon = as.matrix(t(R))
colnames(datExpr.exon) = rownames(R)
rownames(datExpr.exon) = colnames(R)
exon.list = grep("NACA\\>", colnames(datExpr.exon))
data = as.data.frame(datExpr.exon[, exon.list])
data$sample = rownames(data)
info$sample = rownames(info)
info$ID = as.character(info$ID)
data.all = merge(info, data)
samples.low.naca = unique(data.all[data.all$NACA.6 < -4, "sample"])
sample.file = "analysisExon/DE_Exon_NACA/naca/samples_low_naca.txt"
write.table(samples.low.naca, file = sample.file, sep = ", ", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)

#########
# Gene DE
#########

# Generate DE
# NACA down in these HLHS: 1-00425, 1-01026
# also this other: 1-05824
# run DE between naca.all vs all and naca.hlhs vs all

setwd("~/Documents/phd/variancepartitionpcgc")
options(stringsAsFactors=FALSE)
vobj <- readRDS("matrix.gene.RDS")
info <- readRDS("info.gene.RDS")

naca.col = rep(0, nrow(info))
naca.col[rownames(info) %in% samples.low.naca[1:5]] = 1
info$naca = naca.col

# shorten Chr 
vobj$genes$ChrShort <- sapply(strsplit(vobj$genes$Chr, ";"), function(x) x[1])

# fit linear model
design <- model.matrix( ~ Batch + naca, data = info)
colnames(design)

# fit limmar linear model
fit <- lmFit(vobj, design)
fit <- eBayes(fit)
summary(decideTests(fit))
topSet <- topTable(fit, coef = "naca", p.value = 0.05, 
                   number = 17000)
# filters by adj.P.Val
topSet[nrow(topSet),]
colnames(topSet)
de.up.hlhs = topSet[topSet$logFC > 1, c(1,7:14)]
de.down.hlhs = topSet[topSet$logFC < -2, c(1,7:14)]

# de.up = topSet[topSet$logFC > 1, c(1,7:14)]
# de.down = topSet[topSet$logFC < 1, c(1,7:14)]

#########################################
# enrichment of genes that NACA regulates
#########################################

# don't use kolmogorov-smirnov tot test for enrichment because it assumes ordering
# and just looking for absolute increase in 

# down:
# Irx4 (IRXA3). grep "IRX3"
# Myoglobin (MB, PVALB). grep "\\<MB\\>"
# Periostin (POSTN, PN, OSF-2). grep "POSTN"
# Cxcl12 (Sdf-1). grep "CXCL12"
# kcne1 (JLNS2, LQT5, JLNS, LQT2/5, MINK, ISK). grep "KCNE1\\>"
# Myh7 (MYHCB). grep "MYH7"
# Myh2 (MYH2A). grep "MYH2"

downGenes = c("IRX3", "\\<MB\\>", "POSTN", "CXCL12", "KCNE1\\>", "MYH7\\>", "MYH2", "NACA\\>")

down.interest = length(grep(paste(downGenes, collapse = "|"), de.down.hlhs$GeneSymbol))
total.genes = nrow(vobj$genes)
total.de = nrow(de.down.hlhs)
p = total.de/total.genes
1-pbinom(down.interest-1,8,p)

upGenes = c("KRT8\\>", "SHOX2", "MYH4")
normalGenes = c("SHOX2", "MYH4")

#######################
# Plot downstream genes
#######################
setwd("~/Documents/phd/variancepartitionpcgc")
options(stringsAsFactors=FALSE)
vobj <- readRDS("matrix.gene.RDS")
info <- readRDS("info.gene.RDS")

sample.file = "analysisExon/DE_Exon_NACA/naca/samples_low_naca.txt"
samples.low = read.table(file = sample.file)

# [1:5]
naca.col = rep(0, nrow(info))
naca.col[rownames(info) %in% samples.low$V1] = 1
info$naca = naca.col

fit = lmFit(vobj, model.matrix(~ as.factor(Batch), info))
R = residuals(fit, vobj)
datExpr = as.matrix(t(R))
colnames(datExpr) = rownames(R)
rownames(datExpr) = colnames(R)

downGenes = c("IRX3", "\\<MB\\>", "POSTN", "CXCL12", "KCNE1\\>", "MYH7\\>", "MYH2", "NACA\\>")
upGenes = c("KRT8\\>", "SHOX2", "MYH4")
normalGenes = c("SHOX2", "MYH4")

down.indices = grep(paste(c(downGenes, upGenes, normalGenes), collapse = "|"), colnames(datExpr))
data = as.data.frame(datExpr[, down.indices])
data$sample = rownames(data)
info$sample = rownames(info)
info$ID = as.character(info$ID)
data.all = merge(info, data)

library(reshape)
data.long = melt(data)
data.long.all = merge(data.long, info)
# plot
library(ggplot2)

#data.long.all[data.long.all$sample %in% samples.low$V1, ]
ggplot(data.long.all, aes(x = variable, y = value, col = as.factor(naca))) + 
  geom_point()

# provide boxplot + error bars
# plot all genes not just down, as well as just down

  xlab ("NACA exon") +
  ylab("Expression (normalized)") + 
  scale_x_discrete(breaks = c("NACA", "NACA.1", "NACA.2", "NACA.3", "NACA.4", 
                              "NACA.5", "NACA.6", "NACA.7", "NACA.8"),
                   labels = c("9", "8", "7", "6", "5", "4", "3", "2", "1")) + 
  theme_bw() +
  theme(legend.title = element_blank()) 



