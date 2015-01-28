############FUNCTIONS#########################

# This function reads in the troughcomparisons file which contains all the paired troughs and different features
# describing it. It assigns Trough Names to the file and writes a file that has troughNames for each row to the 
# writefilename if writefile is TRUE. It returns a data frame with rownames that are trough names. 
readTroughComparisonsAndAddTroughNames <- function(readfilename,writefile,writefilename){
  troughComparisons <- read.table(readfilename,header = TRUE,sep = "\t")
  tn <- rep(NA,nrow(troughComparisons))
  for (i in 1:nrow(troughComparisons)){
    tn[i] <- paste('Trough',i,sep="")
  }
  row.names(troughComparisons) <- tn
  if (writefile){
    write.table(troughComparisons,file = writefilename,sep = "\t",col.names = FALSE,row.names=TRUE)
  }
  return(troughComparisons)
}

# This function returns a data frame that contains just the features and the row names 
# which can be used for clustering purposes
minimizeColsInTroughs <- function(troughComparisons,valsabs){
  tc <- troughComparisons[,c(6:ncol(troughComparisons))]
  if (valsabs){
    WTMUTdiff <- abs(tc$WT.Term.Ratio - tc$MUT.Term.Ratio)
  }
  else{
    WTMUTdiff <- tc$WT.Term.Ratio - tc$MUT.Term.Ratio
  }
  tc <- cbind(WTMUTdiff,tc[,3:ncol(tc)])
  row.names(tc) <- row.names(troughComparisons)
  return(tc)
}

# This function clusters the troughs using hierarchical clustering and then plots the clusters in a pdf.
# It also looks to see the silhoutte width of the clusters as you cutoff the tree at various distances
# to create a different number of clusters
clusterTroughsUsingHierarchical <- function(troughs){
  tc.t <- t(troughs)
  colnames(tc.t) <- rownames(troughs)
  tc.cor <- cor(tc.t,method='pearson')
  tc.dist <- as.dist(1-tc.cor)
  tc.clust <- hclust(tc.dist,method='ave')
  pdf('troughComparisonCluster.pdf')
  plot(tc.clust,cex=0.6)
  dev.off()
  library(cluster)
  avg.sil.values=numeric()
  avg.sil.values[1]=0
  for (i in 2:20) {
    temp.clusters = cutree(tc.clust, k=i)
    silhouette(temp.clusters, dist=tc.dist)-> temp.cluster.sil
    avg.sil.values[i]=mean(temp.cluster.sil[,"sil_width"])
  }
  print(avg.sil.values)
  pdf("sigsilplot.pdf")
  plot(avg.sil.values)
  dev.off()
  return(tc.clust)
}

#siggroups = cutree(tc.clust, k=6)
#table(siggroups)
silhouetteWidthKmeans <- function(dat,k){
  sil=rep(0,k)
  for (i in 2:k){
    km=kmeans(dat, centers=i)
    mem=km$cluster
    aa=silhouette(mem,dist(dat))
    sil[i]=mean(aa[,3])
  }
plot(1:k, sil)
lines(1:k, sil)
}


# This function uses kmeans clustering to cluster troughs. Need to provide the number of clusters you want.
kmeans.cluster.grps <- function(troughs,clusts){
  kc.org <- kmeans(troughs,clusts)
  kc.org.cluster <- kc.org$cluster
  print(table(kc.org.cluster))
  #numclusters <- unique(kc.org.cluster)
  clusters <- list()
  for (i in 1:clusts){
    clusters[[i]] <- names(kc.org.cluster[kc.org.cluster==i])
  }
  return(clusters)
}

# This function helps to create a vector that can be used to discard troughs that are unnecessary
clusters.to.discard <- function(troughnamesToDiscard){
  to.discard.num <- rep(0,length(troughnamesToDiscard))
  for (i in 1:length(to.discard)){
    to.discard.num[i] <- as.numeric(unlist(regmatches(troughnamesToDiscard[i], gregexpr('\\(?[0-9,.]+', troughnamesToDiscard[i]))))
  }
  return(to.discard.num)
}

# This function plots individual plot for a feature and computes the p-value based on a wilcoxon rank sum test
boxplot.individual.feat <- function(feature,troughs,cluster1,cluster2,showpval,alt=NULL,sideval=NULL,lineval=NULL){
  x = abs(troughs[cluster1,feature])
  y = abs(troughs[cluster2,feature])
  # Plot boxplot
  boxplot(x,y,names=c("Dissimilar","Similar"),ylab=feature)
  if (showpval){
    # Perform statistical test
    w <- wilcox.test(x,y,alternative=alt,correct=FALSE)
    pval <- format(w$p.value,digits=2)
    mtext(paste("p-val=",pval,sep=""),side=sideval,line=lineval)
  }
}

# This function plots boxplots for a feature for all the clusters in variable clusters
showboxplotsForFeat <- function(troughs,clusters,feature,absval=FALSE){
  numclusts <- length(clusters)
  valstoplot <- list()
  for (i in 1:numclusts){
    if (absval){
      df <- data.frame(abs(troughs[clusters[[i]],feature]))
    }
    else{
      df <- data.frame(troughs[clusters[[i]],feature])
    }
    
    colnames(df) <- paste("Cluster",toString(i),sep=" ")
    valstoplot[[i]] <- df
  }
  alldf <- melt(valstoplot)
  # Plot boxplot
  qplot(factor(variable), value, data = alldf, geom = "boxplot")
}

gettRNAs <- function(strain.trnas,trnas,strand,posorneg,window){
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