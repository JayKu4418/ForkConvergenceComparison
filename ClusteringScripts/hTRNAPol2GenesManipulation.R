source('~/Documents/SeqData/DSBioinformaticsSpyder/CompareOkazakiFragDist/ForkConvergenceComparison/ClusteringScripts/clusterTroughs.R')
rna.pol2.ht.genes <- read.table('~/Documents/SeqData/DSBioinformaticsSpyder/CompareOkazakiFragDist/ForkConvergenceComparison/DataToUse/HighlyTxnRNAPol2GenesSorted.txt',header=FALSE)
colnames(rna.pol2.ht.genes) <- c('Name','Chr','Start','End','Strand')
rownames(rna.pol2.ht.genes) <- rna.pol2.ht.genes$Name
setwd("~/Documents/SeqData/DSBioinformaticsSpyder/CompareOkazakiFragDist/TestData")

yeastsize <- list(230218,813184,316620,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066,12157105)

# Based on window
rna.pol2.ht.genes <- rna.pol2.ht.genes[rna.pol2.ht.genes$Start-5000 >= 0,]
rna.pol2.ht.genes.valid <- data.frame()
for (i in 1:nrow(rna.pol2.ht.genes)){
  if ((rna.pol2.ht.genes[i,]$Start + 5000) <= yeastsize[[rna.pol2.ht.genes[i,]$Chr]]){
    rna.pol2.ht.genes.valid <- rbind(rna.pol2.ht.genes.valid,rna.pol2.ht.genes[i,])
  }
  else{
    print(rna.pol2.ht.genes[i,]$Start + 5000)
    print(yeastsize[[rna.pol2.ht.genes[i,]$Chr]])
  }
}

writeConvDivFilesForStrain('cdc9deg','HighlyTransRNAPol2Genes',rna.pol2.ht.genes.valid,c(6,7,16,19,26,30,42),c())
writeConvDivFilesForStrain('rrm3d','HighlyTransRNAPol2Genes',rna.pol2.ht.genes.valid,c(6,7,16,19,26,30,42),c())
writeConvDivFilesForStrain('pif1m2','HighlyTransRNAPol2Genes',rna.pol2.ht.genes.valid,c(6,7,16,19,26,30,42),c())
writeConvDivFilesForStrain('rrm3d_pif1m2','HighlyTransRNAPol2Genes',rna.pol2.ht.genes.valid,c(6,7,16,19,26,30,42),c())

