source('~/Documents/SeqData/DSBioinformaticsSpyder/CompareOkazakiFragDist/ForkConvergenceComparison/ClusteringScripts/clusterTroughs.R')
setwd("~/Documents/SeqData/DSBioinformaticsSpyder/CompareOkazakiFragDist/TestData")
gettRNAs <- function(trnaswatsoncrickratiofile,strand,posorneg,window){
  trnasfile <- 'trnasOrdered2Valid.txt'
  strain.trnas <- read.table(trnaswatsoncrickratiofile,header=F)
  trnas <- read.table(trnasfile,header=T)
  strain.trnasJustRatios <- strain.trnas[,5:(2*window+4)]
  rownames(strain.trnas) <- trnas$Name
  rownames(strain.trnasJustRatios) <- trnas$Name
  rownames(trnas) <- trnas$Name
  if (strand=='W'){
      trnas.names.strand <- rownames(trnas[trnas$Strand=='W',])
  }
  else{
    trnas.names.strand <- rownames(trnas[trnas$Strand=='C',])
  }
  
  strain.trnasJustRatios.strand <- strain.trnasJustRatios[trnas.names.strand,]
  rowmean.WCstrain.strand <- apply(strain.trnasJustRatios.strand,1,mean)
  
  if (posorneg=='pos'){
    trnas.avg.strand.names <- names(rowmean.WCstrain.strand[rowmean.WCstrain.strand>0])
  }
  else{
    trnas.avg.strand.names <- names(rowmean.WCstrain.strand[rowmean.WCstrain.strand<0])
  }
  return(trnas.avg.strand.names)
}

watson.pos.names.cdc9deg <- gettRNAs('cdc9degtRNAGenesValidRawWCRatio-5000.txt','W','pos',5000)
watson.neg.names.cdc9deg <- gettRNAs('cdc9degtRNAGenesValidRawWCRatio-5000.txt','W','neg',5000)
watson.pos.names.rrm3d <- gettRNAs('rrm3dtRNAGenesValidRawWCRatio-5000.txt','W','pos',5000)
watson.neg.names.rrm3d <- gettRNAs('rrm3dtRNAGenesValidRawWCRatio-5000.txt','W','neg',5000)
watson.pos.names.pif1m2 <- gettRNAs('pif1m2tRNAGenesValidRawWCRatio-5000.txt','W','pos',5000)
watson.neg.names.pif1m2 <- gettRNAs('pif1m2tRNAGenesValidRawWCRatio-5000.txt','W','neg',5000)
watson.pos.names.rrm3dpif1m2 <- gettRNAs('rrm3d_pif1m2tRNAGenesValidRawWCRatio-5000.txt','W','pos',5000)
watson.neg.names.rrm3dpif1m2 <- gettRNAs('rrm3d_pif1m2tRNAGenesValidRawWCRatio-5000.txt','W','neg',5000)

trnas <- read.table('trnasOrdered2Valid.txt',header=TRUE)
trnas.watsoncrick.ratio.cdc9deg <- read.table('cdc9degtRNAGenesValidRawWCRatio-5000.txt',header=FALSE)
rownames(trnas.watsoncrick.ratio.cdc9deg) <- trnas$Name
trnas.watsoncrick.ratio.rrm3d <- read.table('rrm3dtRNAGenesValidRawWCRatio-5000.txt',header=FALSE)
rownames(trnas.watsoncrick.ratio.rrm3d) <- trnas$Name
trnas.watsoncrick.ratio.pif1m2 <- read.table('pif1m2tRNAGenesValidRawWCRatio-5000.txt',header=FALSE)
rownames(trnas.watsoncrick.ratio.pif1m2) <- trnas$Name
trnas.watsoncrick.ratio.rrm3dpif1m2 <- read.table('rrm3d_pif1m2tRNAGenesValidRawWCRatio-5000.txt',header=FALSE)
rownames(trnas.watsoncrick.ratio.rrm3dpif1m2) <- trnas$Name
write.table(trnas.watsoncrick.ratio.cdc9deg[watson.neg.names.cdc9deg,],'cdc9degConvergentWatsontRNAGenesRawWCRatio-5000.txt',quote = FALSE,sep = '\t',row.names = FALSE,col.names = FALSE)
write.table(trnas.watsoncrick.ratio.cdc9deg[watson.pos.names.cdc9deg,],'cdc9degDivergentWatsontRNAGenesRawWCRatio-5000.txt',quote = FALSE,sep = '\t',row.names = FALSE,col.names = FALSE)
write.table(trnas.watsoncrick.ratio.rrm3d[watson.neg.names.rrm3d,],'rrm3dConvergentWatsontRNAGenesRawWCRatio-5000.txt',quote = FALSE,sep = '\t',row.names = FALSE,col.names = FALSE)
write.table(trnas.watsoncrick.ratio.rrm3d[watson.pos.names.rrm3d,],'rrm3dDivergentWatsontRNAGenesRawWCRatio-5000.txt',quote = FALSE,sep = '\t',row.names = FALSE,col.names = FALSE)
write.table(trnas.watsoncrick.ratio.pif1m2[watson.neg.names.pif1m2,],'pif1m2ConvergentWatsontRNAGenesRawWCRatio-5000.txt',quote = FALSE,sep = '\t',row.names = FALSE,col.names = FALSE)
write.table(trnas.watsoncrick.ratio.pif1m2[watson.pos.names.pif1m2,],'pif1m2DivergentWatsontRNAGenesRawWCRatio-5000.txt',quote = FALSE,sep = '\t',row.names = FALSE,col.names = FALSE)
write.table(trnas.watsoncrick.ratio.rrm3dpif1m2[watson.neg.names.rrm3dpif1m2,],'rrm3d_pif1m2ConvergentWatsontRNAGenesRawWCRatio-5000.txt',quote = FALSE,sep = '\t',row.names = FALSE,col.names = FALSE)
write.table(trnas.watsoncrick.ratio.rrm3dpif1m2[watson.pos.names.rrm3dpif1m2,],'rrm3d_pif1m2DivergentWatsontRNAGenesRawWCRatio-5000.txt',quote = FALSE,sep = '\t',row.names = FALSE,col.names = FALSE)

crick.pos.names.cdc9deg <- gettRNAs('cdc9degtRNAGenesValidRawWCRatio-5000.txt','C','pos',5000)
crick.neg.names.cdc9deg <- gettRNAs('cdc9degtRNAGenesValidRawWCRatio-5000.txt','C','neg',5000)
crick.pos.names.rrm3d <- gettRNAs('rrm3dtRNAGenesValidRawWCRatio-5000.txt','C','pos',5000)
crick.neg.names.rrm3d <- gettRNAs('rrm3dtRNAGenesValidRawWCRatio-5000.txt','C','neg',5000)
crick.pos.names.pif1m2 <- gettRNAs('pif1m2tRNAGenesValidRawWCRatio-5000.txt','C','pos',5000)
crick.neg.names.pif1m2 <- gettRNAs('pif1m2tRNAGenesValidRawWCRatio-5000.txt','C','neg',5000)
crick.pos.names.rrm3dpif1m2 <- gettRNAs('rrm3d_pif1m2tRNAGenesValidRawWCRatio-5000.txt','C','pos',5000)
crick.neg.names.rrm3dpif1m2 <- gettRNAs('rrm3d_pif1m2tRNAGenesValidRawWCRatio-5000.txt','C','neg',5000)

write.table(trnas.watsoncrick.ratio.cdc9deg[crick.pos.names.cdc9deg,],'cdc9degConvergentCricktRNAGenesRawWCRatio-5000.txt',quote = FALSE,sep = '\t',row.names = FALSE,col.names = FALSE)
write.table(trnas.watsoncrick.ratio.cdc9deg[crick.neg.names.cdc9deg,],'cdc9degDivergentCricktRNAGenesRawWCRatio-5000.txt',quote = FALSE,sep = '\t',row.names = FALSE,col.names = FALSE)
write.table(trnas.watsoncrick.ratio.rrm3d[crick.pos.names.rrm3d,],'rrm3dConvergentCricktRNAGenesRawWCRatio-5000.txt',quote = FALSE,sep = '\t',row.names = FALSE,col.names = FALSE)
write.table(trnas.watsoncrick.ratio.rrm3d[crick.neg.names.rrm3d,],'rrm3dDivergentCricktRNAGenesRawWCRatio-5000.txt',quote = FALSE,sep = '\t',row.names = FALSE,col.names = FALSE)
write.table(trnas.watsoncrick.ratio.pif1m2[crick.pos.names.pif1m2,],'pif1m2ConvergentCricktRNAGenesRawWCRatio-5000.txt',quote = FALSE,sep = '\t',row.names = FALSE,col.names = FALSE)
write.table(trnas.watsoncrick.ratio.pif1m2[crick.neg.names.pif1m2,],'pif1m2DivergentCricktRNAGenesRawWCRatio-5000.txt',quote = FALSE,sep = '\t',row.names = FALSE,col.names = FALSE)
write.table(trnas.watsoncrick.ratio.rrm3dpif1m2[crick.pos.names.rrm3dpif1m2,],'rrm3d_pif1m2ConvergentCricktRNAGenesRawWCRatio-5000.txt',quote = FALSE,sep = '\t',row.names = FALSE,col.names = FALSE)
write.table(trnas.watsoncrick.ratio.rrm3dpif1m2[crick.neg.names.rrm3dpif1m2,],'rrm3d_pif1m2DivergentCricktRNAGenesRawWCRatio-5000.txt',quote = FALSE,sep = '\t',row.names = FALSE,col.names = FALSE)


strain.trnasJustRatios.crick <- strain.trnasJustRatios[cricktrnas.names,]
rowmean.WCstrain.Crick <- apply(strain.trnasJustRatios.crick,1,mean)

watson.pos.names.cdc9deg <- gettRNAs('cdc9degtRNAGenesRawWCRatio-1000.txt','W','pos')
watson.pos.names.rrm3d <- gettRNAs('rrm3dtRNAGenesRawWCRatio-1000.txt','W','pos')
watson.pos.names.rrm3d.pif1m2 <- gettRNAs('rrm3d_pif1m2tRNAGenesRawWCRatio-1000.txt','W','pos')
watson.pos.names.pif1m2 <- gettRNAs('pif1m2tRNAGenesRawWCRatio-1000.txt','W','pos')
x <- intersect(watson.pos.names.cdc9deg,watson.pos.names.pif1m2)
y <- intersect(watson.pos.names.cdc9deg,watson.pos.names.rrm3d)
z <- intersect(watson.pos.names.cdc9deg,watson.pos.names.rrm3d.pif1m2)
diff(x,watson.pos.names.cdc9deg)


trnas.pos.Watson <- trnas[trnas.avgpos.Watson,]
trnas.neg.Watson <- trnas[trnas.avgneg.Watson,]
rowstd.WCstrain.Watson.pos <- apply(trnas.strain.pos.Watson,1,sd)
boxplot(rowstd.WCstrain.Watson.pos)
clusterTroughsUsingHierarchical(strain.trnasJustRatios.watson)
clusts <- kmeans.cluster.grps(strain.trnasJustRatios.watson,2)
plot(apply(strain.trnasJustRatios.watson[clusts[[1]],],2,mean))
plot(apply(strain.trnasJustRatios.watson[clusts[[2]],],2,mean))
x <- clusterTroughsUsingHierarchical(strain.trnasJustRatios[clusts[[1]],])
siggroups = cutree(x, k=2)
table(siggroups)
silhouetteWidthKmeans(strain.trnasJustRatios.watson,10)
#}



trnasavgposWatson <- names(meanWCcdc9degWatson[meanWCcdc9degWatson>0])

trnasavgnegWatson <- names(meanWCcdc9degWatson[meanWCcdc9degWatson><0])
trnasavgnegWatson <- names(meanWCcdc9degWatson[meanWCcdc9degWatson<0])


plot(colmeantrnascdc9degposWatson)
plot(colmeantrnascdc9degnegWatson)
colmeantrnascdc9degWatson <- apply(cdc9degtrnasWatson,2,mean)
plot(colmeantrnascdc9degWatson)


cl1.further <- kmeans.cluster.grps(strain.trnasJustRatios.watson[cls[[1]],],3)
cl1.further.1 <- apply(strain.trnasJustRatios.watson[cl1.further[[1]],],2,mean)
cl1.further.2 <- apply(strain.trnasJustRatios.watson[cl1.further[[2]],],2,mean)
cl1.further.3 <- apply(strain.trnasJustRatios.watson[cl1.further[[3]],],2,mean)
boxplot(cl1.further.1,cl1.further.2,cl1.further.3)


cl1.further <- kmeans.cluster.grps(strain.trnasJustRatios.watson[cls[[1]],],3)
cl1.further.1 <- apply(strain.trnasJustRatios.watson[cl1.further[[1]],],2,mean)
cl1.further.2 <- apply(strain.trnasJustRatios.watson[cl1.further[[2]],],2,mean)
cl1.further.3 <- apply(strain.trnasJustRatios.watson[cl1.further[[3]],],2,mean)
boxplot(cl1.further.1,cl1.further.2,cl1.further.3)

x <- clusterTroughsUsingHierarchical(trnas.strain.pos.Watson)
siggroups = cutree(x, k=2)
table(siggroups)
grp1 <- names(siggroups[siggroups==1])
grp2 <- names(siggroups[siggroups==2])
grp2.colmean <- apply(strain.trnasJustRatios.watson[grp2,],2,mean)
grp1.colmean <- apply(strain.trnasJustRatios.watson[grp1,],2,mean)
plot(grp2.colmean,col='red',ylim=c(min(grp1.colmean,grp2.colmean),max(grp1.colmean,grp2.colmean)))
lines(grp1.colmean,col='blue')

trnas.watsoncrick.ratio.rrm3dpif1m2.justwatson <- trnas.watsoncrick.ratio.rrm3dpif1m2[trnas.watsoncrick.ratio.rrm3dpif1m2$V4=='W',]
trnas.watsoncrick.ratio.rrm3dpif1m2.justcrick <- trnas.watsoncrick.ratio.rrm3dpif1m2[trnas.watsoncrick.ratio.rrm3dpif1m2$V4=='C',]
trnas.watsoncrick.ratio.rrm3dpif1m2.valid.watson <- trnas.watsoncrick.ratio.rrm3dpif1m2.justwatson[-c(1,2,4,13,15,17,66,70,114,115,119,120),]
trnas.watsoncrick.ratio.rrm3dpif1m2.valid.crick <- trnas.watsoncrick.ratio.rrm3dpif1m2.justcrick[-c(5,7,8,31,33,35,39,40,57,59,64,111,113,114,124),]
trnas.watsoncrick.ratio.rrm3dpif1m2.valid <- rbind(trnas.watsoncrick.ratio.rrm3dpif1m2.valid.watson,trnas.watsoncrick.ratio.rrm3dpif1m2.valid.crick)
trnas.watsoncrick.ratio.rrm3dpif1m2.valid <- trnas.watsoncrick.ratio.rrm3dpif1m2.valid[order(trnas.watsoncrick.ratio.rrm3dpif1m2.valid$V1,trnas.watsoncrick.ratio.rrm3dpif1m2.valid$V2),]
write.table(trnas.watsoncrick.ratio.rrm3dpif1m2.valid,'rrm3d_pif1m2tRNAGenesValidRawWCRatio-5000.txt',quote = F,sep = '\t',row.names = F,col.names = F)

