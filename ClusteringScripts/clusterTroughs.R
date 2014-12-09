troughComparisons <- read.delim("~/Documents/SeqData/DSBioinformaticsSpyder/CompareOkazakiFragDist/ForkConvergenceComparison/DataCreated/troughComparisonsWOtransposons.txt")
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
clust1 <- names(kc.org.cluster[kc.org.cluster==2])
clust2 <- names(kc.org.cluster[kc.org.cluster==1])
#to.discard <- names(kc.org.cluster[kc.org.cluster==2])
#to.discard.num <- rep(0,length(to.discard))
#for (i in 1:length(to.discard)){
  #to.discard.num[i] <- as.numeric(unlist(regmatches(to.discard[i], gregexpr('\\(?[0-9,.]+', to.discard[i]))))
#}
#tc.new <- tc[-to.discard.num,]
#kc.new <- kmeans(tc.new,2)
#kclust.new <- kc.new$cluster
#clust1 <- names(kclust.new[kclust.new==1])
#clust2 <- names(kclust.new[kclust.new==2])
                
tc.clust1.org <- tc[clust1,]
#tc.clust1 <- tc[clust1,]
tc.clust2 <- tc[clust2,]
# I want to grab the same number of samples from tc.clust1 randomly as there are in tc.clust2
rs <- sample(x= 1:nrow(tc.clust1.org),size = length(clust2), replace=FALSE)
tc.clust1 <- tc.clust1.org[rs,]

# Grab different random samples from tc.clust1
rs2 <- sample(x= 1:nrow(tc.clust1.org),size = length(clust2), replace=FALSE)
tc.clust1 <- tc.clust1.org[rs2,]
boxplot(abs(tc.clust1$WTMUTdiff),abs(tc.clust2$WTMUTdiff))
boxplot(tc.clust1$Corr.Coeff,tc.clust2$Corr.Coeff)
boxplot(tc.clust1$Euclidean.Dist,tc.clust2$Euclidean.Dist)
boxplot(tc.clust1$FirstDer.Corr,tc.clust2$FirstDer.Corr)
boxplot(tc.clust1$FirstDer.Dist,tc.clust2$FirstDer.Dist)
boxplot(tc.clust1$SecDer.Corr,tc.clust2$SecDer.Corr)
boxplot(tc.clust1$SecDer.Dist,tc.clust2$SecDer.Dist)
tc.WT.diff.clust1 <- troughComparisons[clust1[rs],]$WT.O2 - troughComparisons[clust1[rs],]$WT.O1
tc.WT.diff.clust2 <- troughComparisons[clust2,]$WT.O2 - troughComparisons[clust2,]$WT.O1
boxplot(tc.WT.diff.clust1,tc.WT.diff.clust2)

write.table(troughComparisons[clust1,],file = "TroughsInClust1Sim.txt",sep = "\t",col.names = TRUE,row.names=FALSE)
write.table(troughComparisons[clust2,],file = "TroughsInClust2Dissim.txt",sep = "\t",col.names = TRUE,row.names=FALSE)
write.table(troughComparisons[clust1[rs],],file = "TroughsInClust1Limit70Sim.txt",sep = "\t",col.names = TRUE,row.names=FALSE)

write.table(troughComparisons[clust1,],file = "TroughsInClust1WOTransposonsSim.txt",sep = "\t",col.names = TRUE,row.names=FALSE)
write.table(troughComparisons[clust2,],file = "TroughsInClust2WOTransposonsDissim.txt",sep = "\t",col.names = TRUE,row.names=FALSE)
write.table(troughComparisons[clust1[rs],],file = "TroughsInClust1LimitWOTransposons46Sim.txt",sep = "\t",col.names = TRUE,row.names=FALSE)
write.table(troughComparisons[clust1[rs2],],file = "TroughsInClust1LimitWOTransposonsDiff46Sim.txt",sep = "\t",col.names = TRUE,row.names=FALSE)

tc.1 <- cbind(tc.clust1,rep(x = "Similar",times=nrow(tc.clust1)))
colnames(tc.1) <- c(colnames(tc.1)[1:7],"Trough.Type")

tc.2 <- cbind(tc.clust2,rep(x = "Dissimilar",times=nrow(tc.clust2)))
colnames(tc.2) <- c(colnames(tc.2)[1:7],"Trough.Type")

tc.fin <- rbind(tc.1,tc.2)

# tau is the true shift in location of Similar Troughs vs Dissimilar troughs 
# perform a upper-tailed test, corr coeff of similar troughs greater than dissimilar
# H0: tau <= 0
# HA: tau > 0
r <- with(tc.fin, rank(Corr.Coeff))
T.1 <- print(sum(r[tc.fin$Trough.Type=="Similar"]))
T.2 <- print(sum(r[tc.fin$Trough.Type=="Dissimilar"]))
n.1 <- sum(tc.fin$Trough.Type=="Similar")
n.2 <- sum(tc.fin$Trough.Type=="Dissimilar")
W.1 <- T.1 - n.1*(n.1+1)/2
W.2 <- T.2 - n.2*(n.2+1)/2
W.star <- W.1
pwilcox(W.star,n.1,n.2,lower.tail = FALSE)
# p-value is small

wilcox.test(abs(WTMUTdiff) ~ Trough.Type, alternative = "less", data = tc.fin, correct = FALSE)
wilcox.test(Corr.Coeff ~ Trough.Type, alternative = "greater", data = tc.fin, correct = FALSE)
wilcox.test(Euclidean.Dist ~ Trough.Type, alternative = "less", data = tc.fin, correct = FALSE)
wilcox.test(FirstDer.Corr ~ Trough.Type, alternative = "greater", data = tc.fin, correct = FALSE)
wilcox.test(FirstDer.Dist ~ Trough.Type, alternative = "less", data = tc.fin, correct = FALSE)
wilcox.test(SecDer.Corr ~ Trough.Type, alternative = "greater", data = tc.fin, correct = FALSE)
wilcox.test(SecDer.Dist ~ Trough.Type, alternative = "less", data = tc.fin, correct = FALSE)


# Copied over from python
# Number of tRNA genes found per base pair With Transposons
zd <- c(1.2965123816932451e-05, 1.4810755570695439e-05, 7.108499395777551e-05, 2.7941546285171422e-05, 6.117031037815486e-05, 0.0, 0.0, 2.446213872478871e-05, 8.647526807333103e-05, 0.0, 2.44143605268619e-05, 3.9197240514267794e-05, 9.311420457190745e-06, 3.801341873681409e-05, 2.6090586516384887e-05, 2.8372282171803625e-05, 5.755395683453237e-05, 7.755094127454972e-05, 2.5959866047091195e-05, 9.355995633868704e-05, 1.994614540740002e-05, 2.4103645676408555e-05, 2.910840941948129e-05, 3.1282749127993365e-05, 6.181615874389565e-05, 7.162682424568e-05, 6.297923259805079e-05, 7.484609771158056e-05, 1.1052899175453721e-05, 8.659758104090293e-05, 1.593447742881272e-05, 3.666630000366663e-05, 9.861689800547324e-05, 3.112501361719346e-05, 1.5954815961197887e-05, 4.6576618537494175e-05, 1.3909947003101918e-05, 6.60414740457007e-05, 0.00016492707867021652, 7.642922653622745e-05, 8.667013347200555e-05, 2.7575176825821397e-05, 1.4094234048850615e-05, 0.00011185682326621924, 5.6902241948332766e-05, 2.104244260673779e-05, 2.0158851751804218e-05, 2.637409009389176e-05, 0.0, 1.1640494953845437e-05, 2.170633499386796e-05, 5.229487335924835e-05, 4.075228722212034e-05, 0.0, 3.298914657077821e-05, 4.8214679303504516e-05, 3.5575161422294955e-05, 3.732875433946769e-05, 1.0967557963543837e-05, 6.32657795398469e-05, 2.3102424599461713e-05, 0.0, 3.3315010077790546e-05, 5.947778504728484e-05, 2.003606491685033e-05, 2.7578979301976035e-05, 1.043231513937573e-05, 2.6736538153039946e-05, 3.2712879060486114e-05, 4.3505688368754214e-05)
zs <- c(2.088075003654131e-05, 0.00010990822663076331, 3.540031860286742e-05, 0.0, 0.0, 0.0, 1.8742034635280007e-05, 4.591930447560154e-05, 0.0, 3.909533397189046e-05, 0.0, 0.0, 9.63020030816641e-05, 0.0, 1.1402248523408816e-05, 5.203590477429426e-05, 1.7203413157170383e-05, 2.491125365884038e-05, 0.0, 2.788155913678693e-05, 2.5801124929046907e-05, 3.6976778583049845e-05, 9.382916836080442e-05, 0.0, 2.6404731727925644e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 1.9893371528606668e-05, 4.1531257462647826e-05, 0.0, 2.483608185972581e-05, 0.0, 4.169620147604553e-05, 8.614748449345279e-05, 1.2663838409421896e-05, 0.0, 2.11891342119761e-05, 3.825042553598409e-05, 6.0201071579074104e-05, 2.5519969376036748e-05, 1.6601367952719304e-05, 0.0, 7.880841673890771e-05, 2.2661862351848074e-05, 2.2012370952475292e-05, 6.948546016745995e-05, 0.0, 0.0, 1.7152364453439908e-05, 2.2404947012300316e-05, 1.8270179413161837e-05, 2.3379234563860378e-05, 0.0, 2.078353943676608e-05, 2.4945120734384356e-05, 0.0, 2.481143310837634e-05, 0.0, 9.358914365933551e-05, 3.0353619669145544e-05, 3.2196786760681284e-05, 7.224650507531698e-05, 0.0, 1.2662074554294975e-05, 3.895749737036893e-05, 0.0, 0.0)
boxplot(zs,zd)

# Number of tRNAgenes found per base pair without Transposons
zs <- c(9.358914365933551e-05, 0.0, 0.0, 1.2662074554294975e-05, 0.0, 3.112501361719346e-05, 4.169620147604553e-05, 0.0, 1.043231513937573e-05, 1.657302905251993e-05, 2.318195516609871e-05, 6.297923259805079e-05, 3.540031860286742e-05, 0.0, 0.0, 0.0, 9.382916836080442e-05, 3.0353619669145544e-05, 0.0, 2.483608185972581e-05, 0.0, 0.0, 2.481143310837634e-05, 0.0, 0.0, 5.318583129454313e-05, 7.224650507531698e-05, 2.2012370952475292e-05, 2.2661862351848074e-05, 8.614748449345279e-05, 1.2797870434359723e-05, 0.0, 0.0, 0.0, 1.3737395939225761e-05, 0.0, 0.0, 2.6736538153039946e-05, 2.689618074233459e-05, 4.476275738585497e-05, 1.0967557963543837e-05, 2.6090586516384887e-05, 0.0, 2.078353943676608e-05, 1.8742034635280007e-05, 3.895749737036893e-05)
# Diff randoms from trough1
zs <- c(0.0, 1.7203413157170383e-05, 3.895749737036893e-05, 0.0, 2.689618074233459e-05, 0.0, 0.0, 2.318195516609871e-05, 7.452675510508272e-05, 0.0, 6.948546016745995e-05, 3.112501361719346e-05, 2.2661862351848074e-05, 0.0, 0.00010990822663076331, 2.078353943676608e-05, 5.52008390527536e-05, 1.2797870434359723e-05, 0.0, 2.3379234563860378e-05, 2.5801124929046907e-05, 0.0, 2.6404731727925644e-05, 0.0, 2.5519969376036748e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 2.6736538153039946e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 1.043231513937573e-05, 5.460154522372983e-05, 0.0, 4.476275738585497e-05, 7.880841673890771e-05, 7.822889775483064e-05, 0.0, 5.9591204338239674e-05, 0.0, 9.313299391531106e-05)
zd <- c(1.2965123816932451e-05, 1.4810755570695439e-05, 6.117031037815486e-05, 0.0, 0.0, 2.446213872478871e-05, 8.647526807333103e-05, 0.0, 2.44143605268619e-05, 3.9197240514267794e-05, 5.755395683453237e-05, 7.755094127454972e-05, 9.355995633868704e-05, 1.994614540740002e-05, 2.910840941948129e-05, 3.1282749127993365e-05, 7.162682424568e-05, 7.484609771158056e-05, 1.1052899175453721e-05, 5.683432793407218e-05, 8.659758104090293e-05, 3.666630000366663e-05, 9.861689800547324e-05, 1.5954815961197887e-05, 4.6576618537494175e-05, 1.3909947003101918e-05, 6.60414740457007e-05, 0.00016492707867021652, 7.642922653622745e-05, 8.667013347200555e-05, 2.7575176825821397e-05, 0.00011185682326621924, 5.6902241948332766e-05, 2.104244260673779e-05, 2.0158851751804218e-05, 2.637409009389176e-05, 0.0, 1.1640494953845437e-05, 4.075228722212034e-05, 0.0, 3.5575161422294955e-05, 3.3315010077790546e-05, 5.947778504728484e-05, 2.003606491685033e-05, 3.2712879060486114e-05, 4.3505688368754214e-05)
boxplot(zs,zd)
wilcox.test(zs,zd,alternative='less',correct=FALSE)

# Number of G4s found per base pair without Transposons
zgs <- c(6.239276243955701e-05, 3.4298257648511456e-05, 0.00010487402008337485, 2.532414910858995e-05, 4.7885840157065554e-05, 3.112501361719346e-05, 4.169620147604553e-05, 5.8343057176196036e-05, 7.302620597563011e-05, 8.286514526259965e-05, 4.636391033219742e-05, 3.148961629902539e-05, 3.540031860286742e-05, 8.1113965120995e-05, 6.570302233902759e-05, 5.7520851308599364e-05, 0.0002814875050824133, 9.106085900743663e-05, 4.9331557397267034e-05, 2.483608185972581e-05, 0.0, 9.941346058256288e-05, 7.443429932512901e-05, 0.0, 0.0, 5.318583129454313e-05, 3.612325253765849e-05, 6.603711285742587e-05, 2.2661862351848074e-05, 8.614748449345279e-05, 2.5595740868719446e-05, 6.845564074479737e-05, 3.091094556582486e-05, 0.00010121457489878542, 2.7474791878451522e-05, 0.0001336719689881032, 3.9207998431680065e-05, 2.6736538153039946e-05, 5.379236148466918e-05, 8.952551477170994e-05, 4.387023185417535e-05, 0.0, 8.899172376968942e-05, 2.078353943676608e-05, 9.371017317640003e-05, 7.791499474073786e-05)
zgd <- c(3.889537145079735e-05, 7.405377785347719e-05, 7.340437245378583e-05, 1.2612249016244577e-05, 7.308979812597757e-05, 4.892427744957742e-05, 5.765017871555402e-05, 3.059008269519022e-05, 3.662154079029285e-05, 1.3065746838089265e-05, 0.00011510791366906474, 0.00011632641191182458, 0.0, 3.989229081480004e-05, 3.881121255930838e-05, 3.1282749127993365e-05, 5.372011818426001e-05, 5.6134573283685424e-05, 4.4211596701814886e-05, 3.788955195604812e-05, 2.886586034696764e-05, 1.8333150001833314e-05, 0.00014792534700820985, 7.977407980598944e-05, 4.6576618537494175e-05, 4.1729841009305755e-05, 1.320829480914014e-05, 7.068303371580708e-05, 0.0, 1.7334026694401108e-05, 5.515035365164279e-05, 0.00011185682326621924, 1.4225560487083191e-05, 9.469099173032006e-05, 0.00016127081401443375, 6.59352252347294e-05, 8.662508662508663e-05, 4.656197981538175e-05, 4.075228722212034e-05, 6.410667350471184e-05, 8.893790355573738e-05, 6.663002015558109e-05, 5.947778504728484e-05, 5.009016229212583e-05, 4.089109882560764e-05, 2.1752844184377107e-05)
boxplot(zgs,zgd)
wilcox.test(zgs,zgd,alternative='less',correct=FALSE)


sim.t.trna <- data.frame(trna.ratio=zs,trough.type=rep("Similar",times=length(zs)))
dissim.t.trna <- data.frame(trna.ratio=zd,trough.type=rep("Dissimilar",times=length(zd)))

all.trna <- rbind(sim.t.trna,dissim.t.trna)

# tau is the true shift in location of tRNA Ratios for Similar Troughs vs Dissimilar troughs 
# perform a lower-tailed test, tRNA Ratios of similar troughs lower than dissimilar
# H0: tau >= 0
# HA: tau < 0
r <- with(all.trna, rank(trna.ratio))
T.1 <- print(sum(r[tc.fin$Trough.Type=="Similar"]))
T.2 <- print(sum(r[tc.fin$Trough.Type=="Dissimilar"]))
n.1 <- sum(tc.fin$Trough.Type=="Similar")
n.2 <- sum(tc.fin$Trough.Type=="Dissimilar")
W.1 <- T.1 - n.1*(n.1+1)/2
W.2 <- T.2 - n.2*(n.2+1)/2
W.star <- W.1
pwilcox(W.star,n.1,n.2,lower.tail = FALSE)
mu.W <- n.1*n.2/2
sig2.W <- (n.1*n.2/12)*( (n.1+n.2+1) - (14*(14^2 - 1)/((n.1+n.2)*(n.1+n.2-1))))
Z.star <- (W.star - mu.W)/sqrt(sig2.W)
pnorm(Z.star,lower.tail = FALSE)



wilcox.test(trna.ratio ~ trough.type, alternative = "greater", data = all.trna, correct = FALSE)

wilcox.test(zd,zs,alternative = "greater",correct=FALSE)
