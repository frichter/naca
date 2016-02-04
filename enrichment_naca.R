

############################

############################
library("WGCNA")
library("limma")
options(stringsAsFactors=FALSE)
<<<<<<< HEAD
setwd("/media/Storage/PhD/naca")
source("analysis_functions.R")

# obtain subset sample names (and save to file)
samples.low.naca = ObtainLowNacaSampleNames(sample.file)

sample.file = "id_info/samples_low_naca.txt"
write.table(samples.low.naca, file = sample.file, sep = ", ", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)
setwd("naca/")
#setwd("/media/Storage/PhD/naca")
source("analysis_functions.R")

# obtain subset sample names (and save to file)
# samples.low.naca = ObtainLowNacaSampleNames()
# sample.file = "id_info/samples_low_naca.txt"
# write.table(samples.low.naca, file = sample.file, sep = ",", quote = FALSE, 
#             row.names = FALSE, col.names = FALSE)


#########
# Gene DE
#########

# Generate DE
# NACA down in these HLHS: 1-00425, 1-01026
# also this other: 1-05824
# run DE between naca.all vs all and naca.hlhs vs all

options(stringsAsFactors=FALSE)
samples.low.naca = read.table(file = "id_info/samples_low_naca.txt")$V1

>>>>>>> c4a81a66080eca2e76d5e9d4e3ab303fc93a4d5e
vobj <- readRDS("matrix.gene.RDS")
info <- readRDS("info.gene.RDS")

naca.col = rep(0, nrow(info))
naca.col[rownames(info) %in% samples.low.naca, ] = 1
info$naca = naca.col

naca.hlhs.col = rep(0, nrow(info))
naca.hlhs.col[rownames(info) %in% samples.low.naca, ] = 1
info$naca.hlhs = naca.hlhs.col

# shorten Chr 
vobj$genes$ChrShort <- sapply(strsplit(vobj$genes$Chr, ";"), function(x) x[1])

# fit linear model
design <- model.matrix( ~ Batch + naca.hlhs, data = info)
colnames(design)

# fit limmar linear model
fit <- lmFit(vobj, design)
fit <- eBayes(fit)
summary(decideTests(fit))
topSet <- topTable(fit, coef = "naca.hlhs", p.value = 0.05, 
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

setwd("/media/Storage/PhD/naca")
library("limma")

# setwd("/media/Storage/PhD/naca")
setwd("D:/PhD/naca")
>>>>>>> c4a81a66080eca2e76d5e9d4e3ab303fc93a4d5e
options(stringsAsFactors=FALSE)
vobj <- readRDS("expression_data/matrix.gene.RDS")
info <- readRDS("expression_data/info.gene.RDS")

sample.file = "id_info/samples_low_naca.txt"
samples.low = read.table(file = sample.file)$V1

# [1:5]
naca.col = rep(0, nrow(info))
naca.col[rownames(info) %in% samples.low] = 1
info$naca = naca.col

fit = lmFit(vobj, model.matrix(~ as.factor(Batch), info))
R = residuals(fit, vobj)
datExpr = as.matrix(t(R))
colnames(datExpr) = rownames(R)
rownames(datExpr) = colnames(R)

downGenes = c("IRX3", "\\<MB\\>", "POSTN", "CXCL12", "KCNE1\\>", "MYH7\\>", "MYH2", "NACA\\>")
# upGenes = c("KRT8\\>", "SHOX2", "MYH4")
# normalGenes = c("SHOX2", "MYH4")
# c(downGenes, upGenes, normalGenes)
down.indices = grep(paste(downGenes, collapse = "|"), colnames(datExpr))
data = as.data.frame(datExpr[, down.indices])
data$sample = rownames(data)
info$sample = rownames(info)
info$ID = as.character(info$ID)

library(reshape)
data.long = melt(data)
data.long.all = merge(data.long, info)
# plot

library(ggplot2)
#data.long.all[data.long.all$sample %in% samples.low$V1, ]
<<<<<<< HEAD
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



plot.down = ggplot(data=subset(data.long.all, naca == 0), aes(x = variable, y = value)) + 
  geom_point(color = "red") +
  geom_point(data=subset(data.long.all, naca == 1), aes(x = variable, y = value), color = "blue") +
  xlab ("NACA downstream target") +
  ylab("Expression (normalized)") + 
  theme_bw() +
  theme(legend.title = element_blank())
f = "figures/expr.downstream_targets.png"
ggsave(filename = f, plot = plot.down, width = 5, height = 5, dpi = 300)

# provide boxplot + error bars
# plot all genes not just down, as well as just down


  

