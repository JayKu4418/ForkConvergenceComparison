troughComparisons <- read.delim("~/Documents/SeqData/DSBioinformaticsSpyder/CompareOkazakiFragDist/ForkConvergenceComparison/DataToTest/troughComparisonsMoreInfoLatest2.txt")
View(troughComparisons)
tn <- rep(NA,nrow(troughComparisons))
for (i in 1:nrow(troughComparisons)){
tn[i] <- paste('Trough',i,sep="")
}
row.names(troughComparisons) <- tn
tc <- troughComparisons[,c(6:ncol(troughComparisons))]
WTMUTdiff <- tc$WT.Term.Ratio - tc$MUT.Term.Ratio
tc <- cbind(WTMUTdiff,tc[,3:ncol(tc)])
row.names(tc) <- tn
View(tc)
tc.t <- t(tc)
colnames(tc.t) <- tn
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
pdf("sigsilplot.pdf")
plot(avg.sil.values)
dev.off()
siggroups = cutree(tc.clust, k=3)
table(siggroups)

kc.org <- kmeans(tc,2)
kc.org.cluster <- kc.org$cluster
table(kc.org.cluster)
to.discard <- names(kc.org.cluster[kc.org.cluster==2])
to.discard.num <- rep(0,length(to.discard))
for (i in 1:length(to.discard)){
  to.discard.num[i] <- as.numeric(unlist(regmatches(to.discard[i], gregexpr('\\(?[0-9,.]+', to.discard[i]))))
}
tc.new <- tc[-to.discard.num,]
kc.new <- kmeans(tc.new,2)
kclust.new <- kc.new$cluster
clust1 <- names(kclust.new[kclust.new==1])
clust2 <- names(kclust.new[kclust.new==2])
                
tc.clust1 <- tc[clust1,]
tc.clust2 <- tc[clust2,]

boxplot(tc.clust1$Corr.Coeff,tc.clust2$Corr.Coeff)
boxplot(tc.clust1$Euclidean.Dist,tc.clust2$Euclidean.Dist)
boxplot(tc.clust1$FirstDer.Corr,tc.clust2$FirstDer.Corr)
boxplot(tc.clust1$FirstDer.Dist,tc.clust2$FirstDer.Dist)
boxplot(tc.clust1$SecDer.Corr,tc.clust2$SecDer.Corr)
boxplot(tc.clust1$SecDer.Dist,tc.clust2$SecDer.Dist)

write.table(troughComparisons[clust1,],file = "TroughsInClust1Sim.txt",sep = "\t",col.names = TRUE,row.names=FALSE)
write.table(troughComparisons[clust2,],file = "TroughsInClust2Dissim.txt",sep = "\t",col.names = TRUE,row.names=FALSE)
