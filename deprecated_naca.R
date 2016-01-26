# package installation commands

install.packages("ggplot2")
source("https://bioconductor.org/biocLite.R")
biocLite("limma")

source("http://bioconductor.org/biocLite.R")
biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore"))
install.packages("WGCNA") 

install.packages("reshape")

# failed statistics for downstream genes


up.interest = length(grep(paste(upGenes, collapse = "|"), de.up.hlhs$GeneSymbol))
normal.interest = length(normalGenes) - 
  (length(grep(paste(normalGenes, collapse = "|"), de.down.hlhs$GeneSymbol)) +
     length(grep(paste(upGenes, collapse = "|"), de.up.hlhs$GeneSymbol)) )

grep(paste(downGenes, collapse = "|"), de.down.hlhs$GeneSymbol, value = TRUE)

total.genes = nrow(vobj$genes)
total.de = nrow(de.down.hlhs)
total.interest = length(downGenes)
interests.de = length(grep(paste(downGenes, collapse = "|"), de.down.hlhs$GeneSymbol))
# binomial test of proportions to see if significant

obs = interests.de/total.de
pi = total.interest/total.genes
num = (obs) - pi
den = sqrt ( (pi*(1-pi))/7)
1 - pnorm(num/den)

downGenes = c("IRX3", "\\<MB\\>", "POSTN", "CXCL12", "KCNE1\\>", "MYH7\\>", "MYH2", "NACA\\>")
fit$p.value[grep(paste(downGenes, collapse = "|"), fit$genes$GeneSymbol), ]
topSet <- topTable(fit, coef = snowflake_colname, p.value = 1, 
                   number = 17032)
topSet[grep(paste(downGenes, collapse = "|"), topSet$GeneSymbol), c(7:14)]
dim(topSet[topSet$logFC < 0, ])
# binomial test of proportions to see if significant, divide pi by 2 to only look for down genes
# x/n = 3/7
pi = (6658/17032)/2 #p = 0.13
pi = (4730/17032)/2 # p = 0.05
pi = (3739/17032)/2 # p = 0.025
pi = 8339/17032  # number of trending down genes, vs 8693 up
num = (2/7) - pi
den = sqrt ( (pi*(1-pi))/7)
1 - pnorm(num/den)

# only down genes
pi = (3311/8339) #p = 0.13
pi = (2374/8339) # p = 0.05
pi = (1877/8339) # p = 0.025
num = (2/4) - pi
den = sqrt ( (pi*(1-pi))/7)
1 - pnorm(num/den)
# looking at all down trending genes is disengenuous bc a lot of down could be up and vice versa

# conclusion: with p-value < 0.05 MYH2 and MB are down, no change observed in the others
# 13% probability this occured by chance
# from a biological perspective

qt(0.975, 1000)


# up:
# Krt2.8
# Shox2
# Myh4
upGenes = c("KRT8\\>", "SHOX2", "MYH4")
fit$p.value[grep(paste(upGenes, collapse = "|"), fit$genes$GeneSymbol), ]
topSet[grep(paste(upGenes, collapse = "|"), topSet$GeneSymbol), c(7:14)]
summary(topSet$AveExpr)
# are up, but not significantly
# none of the expected up genes are up :( although fold change is positive)

# Normal:
# Smyd1
# Hand2
normalGenes = c("SHOX2", "MYH4")
fit$p.value[grep(paste(normalGenes, collapse = "|"), fit$genes$GeneSymbol), ]
topSet[grep(paste(normalGenes, collapse = "|"), topSet$GeneSymbol), c(7:14)]
# genes that are neither up or down
# don't know which samples are alternatively spliced

# binomial for all findings
pi = (4730/17032) # p = 0.05
num = (4/13) - pi
den = sqrt ( (pi*(1-pi))/7)
1 - pnorm(num/den)


