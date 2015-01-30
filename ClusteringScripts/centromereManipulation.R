source('~/Documents/SeqData/DSBioinformaticsSpyder/CompareOkazakiFragDist/ForkConvergenceComparison/ClusteringScripts/clusterTroughs.R')
centromeres <-  read.table('~/Documents/SeqData/DSBioinformaticsSpyder/CompareOkazakiFragDist/ForkConvergenceComparison/DataToUse/allcent.txt',header=F)
colnames(centromeres) <- c('Chr','Start','End','Strand')
centromeres <- centromeres[order(centromeres$Chr),]
centnames <- vector()
for (i in 1:16){
  namefori <- paste('CEN',i,sep='')
  centnames[i] <- namefori
}
centwatson <- centromeres[centromeres$Strand=='W',]
centcrick <- centromeres[centromeres$Strand=='C',]
centcrick <- centcrick[,c(1,2,4,3,5)]
colnames(centcrick) <- c('Name',Chr','Start','End','Strand')
centromeres <- rbind(centwatson,centcrick)
rownames(centromeres) <- centnames
centromeres <- cbind(Name=centnames,centromeres)
writeConvDivFilesForStrain('cdc9deg','Centromeres',centromeres,c(),c())
writeConvDivFilesForStrain('rrm3d','Centromeres',centromeres,c(),c())
writeConvDivFilesForStrain('pif1m2','Centromeres',centromeres,c(),c())
writeConvDivFilesForStrain('rrm3d_pif1m2','Centromeres',centromeres,c(),c())
