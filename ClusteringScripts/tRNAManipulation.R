source('~/Documents/SeqData/DSBioinformaticsSpyder/CompareOkazakiFragDist/ForkConvergenceComparison/ClusteringScripts/clusterTroughs.R')
setwd("~/Documents/SeqData/DSBioinformaticsSpyder/CompareOkazakiFragDist/TestData")
gettRNAs <- function(trnaswatsoncrickratiofile,strand,posorneg){
  trnasfile <- 'trnasOrdered2.txt'
  strain.trnas <- read.table(trnaswatsoncrickratiofile,header=F)
  trnas <- read.table(trnasfile,header=T)
  strain.trnasJustRatios <- strain.trnas[,5:2004]
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

watson.pos.names.cdc9deg <- gettRNAs('cdc9degtRNAGenesRawWCRatio-1000.txt','W','pos')
watson.pos.names.rrm3d <- gettRNAs('rrm3dtRNAGenesRawWCRatio-1000.txt','W','pos')

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

