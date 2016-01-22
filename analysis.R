# L2HLHS

###########################################
# graph NACA exon expression per individual
###########################################

library("WGCNA")
library("limma")
options(stringsAsFactors=FALSE)
setwd("~/Documents/phd/variancepartitionpcgc")

vobj = readRDS("matrix.exon.RDS")
info = readRDS("info.exon.RDS")
# 
# vobj = readRDS("vobj.tech.exon.RDS")
# info = readRDS("info.tech.exon.RDS")

info$L2 <- as.character(info$L2)
info[info$L2 %in% c("Fallot"), ]$L2 = "TOF"
info$L2 <- droplevels(factor(info$L2))

#fit = lmFit(vobj, model.matrix(~ as.factor(Batch), info))
fit = lmFit(vobj, model.matrix(~ as.factor(Batch) + Storage, info))
R = residuals(fit, vobj)
datExpr.exon = as.matrix(t(R))
colnames(datExpr.exon) = rownames(R)
rownames(datExpr.exon) = colnames(R)

# paste(grep("NACA\\.", colnames(datExpr.exon), value = T), collapse = "\", \"")
# grep("NACA", colnames(datExpr.exon), value = T)
# exon.list = grep("NACA\\.", colnames(datExpr.exon))[-8]
# exon.list = grep("NACA\\.", colnames(datExpr.exon))
grep("NACA\\>", colnames(datExpr.exon), value = T)
exon.list = grep("NACA\\>", colnames(datExpr.exon))
data = as.data.frame(datExpr.exon[, exon.list])
data$sample = rownames(data)
info$sample = rownames(info)

# i = sapply(info, is.factor)
# info[i] = lapply(info[i], as.character)
info$ID = as.character(info$ID)
data.all = merge(info, data)
library(reshape)
data.long = melt(data)
data.long.all = merge(data.long, info)
unique(data.long.all$variable)
unique(data.long.all$L2)
# plot
library(ggplot2)
# L3 and TissueMeta only work if they are factors, ID only works as character

paste(grep("NACA\\>", colnames(datExpr.exon), value = T), collapse = "\", \"")
paste(rev(seq(1:9)), collapse = "\", \"")
# opened in new window and replaced all occurrences of slash

naca.all.plot = ggplot(data.long.all, aes(x = variable, y = value, col = L2)) + 
  geom_point() + 
  xlab ("NACA exon") +
  ylab("Expression (normalized)") + 
  scale_x_discrete(breaks = c("NACA", "NACA.1", "NACA.2", "NACA.3", "NACA.4", 
                              "NACA.5", "NACA.6", "NACA.7", "NACA.8"),
                   labels = c("9", "8", "7", "6", "5", "4", "3", "2", "1")) + 
  theme_bw() +
  theme(legend.title = element_blank()) 
  
ggsave("analysisExon/DE_Exon_F30/naca/expr_by_exon.all.L2.png", 
       plot = naca.all.plot, width = 4, height = 3, dpi = 300)


data.all[data.all$NACA.6 < -2.5, "sample"]
mu = mean(data.all$NACA.6)
se = sd(data.all$NACA.6)/(sqrt(30))
xbar = mean(data.all[data.all$NACA.6 < -4, "NACA.6"])
(xbar - mu) / se

# use -2.5 to incorporate a bad sample, -4 for only positives
data.all.naca.low.samples = unique(data.all[data.all$NACA.6 < -4, "sample"])
data.long.all.naca.low = data.long.all[data.long.all$sample %in% data.all.naca.low.samples, ]

naca.low.ID.plot = ggplot(data.long.all.naca.low, aes(x = variable, y = value, col = ID)) + 
  geom_point() +
  geom_point(data.long.all, aes(x = variable, y = value))
ggsave("analysisExon/DE_Exon_F30/naca/expr_by_exon.low.1badsample.png", 
       plot = naca.low.ID.plot, width = 7, height = 6, dpi = 300)

###############
# Exon graphing
###############

i = "NACA"
plotSplice(ex, geneid=i)

ex[grep("NACA\\.", rownames(ex)), ]
fit = ex
i = "NACA"
j <- fit$gene.firstexon[i]:fit$gene.lastexon[i]
coef = ncol(fit)
exoncolname <- fit$exoncolname
exon.start = fit$genes[j, exoncolname]
exon.end = fit$genes[j, exoncolname] + 112
exon.end = fit$genes[j, exoncolname] + 5000
exon.pos = rbind(exon.start, exon.end)
apply(exon.pos, 2, function(x) paste("chr12:",x[1],"-",x[2],sep=""))
# list of exon locations for UCSC

exon.id <- fit$genes[j, exoncolname]
geneid = NULL
xlab <- paste("Exon", exoncolname, sep = " ")
fit$coefficients[j, coef] + 113
# log(FC) relative to other exons: -0.87318295
plot(fit$coefficients[j, coef], xlab = "", ylab = "logFC (this exon vs rest)", 
     main = geneid, type = "b", xaxt = "n")
axis(1, at = 1:length(j), labels = exon.id, las = 2, 
     cex.axis = 0.5)

info[info$ID %in% c("1-00425", "1-01026", "1-05824"),]

topSpliceSimes[grep("NACA", topSpliceSimes$GeneSymbol),]
