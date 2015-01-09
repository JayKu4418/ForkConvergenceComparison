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
  pdf("sigsilplot.pdf")
  plot(avg.sil.values)
  dev.off()
  return(tc.clust)
}

siggroups = cutree(tc.clust, k=3)
table(siggroups)

# This function uses kmeans clustering to cluster troughs. Need to provide the number of clusters you want.
kmeans.cluster.grps <- function(troughs,clusts){
  kc.org <- kmeans(troughs,clusts)
  kc.org.cluster <- kc.org$cluster
  print(table(kc.org.cluster))
  numclusters <- unique(kc.org.cluster)
  clusters <- list()
  for (i in 1:length(numclusters)){
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
  boxplot(x,y,names=c("Similar","Dissimilar"),ylab=feature)
  if (showpval){
    # Perform statistical test
    w <- wilcox.test(x,y,alternative=alt,correct=FALSE)
    pval <- format(w$p.value,digits=2)
    mtext(paste("p-val=",pval,sep=""),side=sideval,line=lineval)
  }
}
#################TOEXECUTE##############################
setwd("~/Documents/SeqData/DSBioinformaticsSpyder/CompareOkazakiFragDist/ForkConvergenceComparison/DataCreated")
alltroughs <- readTroughComparisonsAndAddTroughNames("troughComparisonsWOtransposons.txt",F,"")
troughsJustFeats <- minimizeColsInTroughs(alltroughs,F)
clusts <- kmeans.cluster.grps(troughsJustFeats,2)
# Take clusters to be randomly picked 80 and then another 80
#clust1 <- sample(x=c(1:160),size = 80,replace = FALSE)
#clust2 <- setdiff(c(1:160),clust1)
#names.clust1 <- rownames(alltroughs[clust1,])
#names.clust2 <- rownames(alltroughs[clust2,])
#clusts = list()
#clusts[[1]] <- names.clust1
#clusts[[2]] <- names.clust2
boxplot.individual.feat("WTMUTdiff",troughsJustFeats,clusts[[1]],clusts[[2]],TRUE,alt="less",sideval=3,lineval=-5)
boxplot.individual.feat("Corr.Coeff",troughsJustFeats,clusts[[1]],clusts[[2]],TRUE,alt="greater",sideval=1,lineval=-12)
boxplot.individual.feat("Euclidean.Dist",troughsJustFeats,clusts[[1]],clusts[[2]],TRUE,alt="less",sideval=3,lineval=-12)
boxplot.individual.feat("FirstDer.Corr",troughsJustFeats,clusts[[1]],clusts[[2]],TRUE,alt="greater",sideval=1,lineval=-12)
boxplot.individual.feat("FirstDer.Dist",troughsJustFeats,clusts[[1]],clusts[[2]],TRUE,alt="less",sideval=3,lineval=-12)
boxplot.individual.feat("SecDer.Corr",troughsJustFeats,clusts[[1]],clusts[[2]],TRUE,alt="greater",sideval=1,lineval=-12)
boxplot.individual.feat("SecDer.Dist",troughsJustFeats,clusts[[1]],clusts[[2]],TRUE,alt="less",sideval=3,lineval=-12)


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
boxplot(zs,zd,names=c("Similar","Dissimilar"),ylab="tRNAGenesPerBase")
pval <- wilcox.test(zs,zd,alternative='less',correct=FALSE)$p.value
mtext(paste("p-val=",pval,sep=""),side=3,line=-5)

# Number of G4s found per base pair without Transposons
zgs <- c(6.239276243955701e-05, 3.4298257648511456e-05, 0.00010487402008337485, 2.532414910858995e-05, 4.7885840157065554e-05, 3.112501361719346e-05, 4.169620147604553e-05, 5.8343057176196036e-05, 7.302620597563011e-05, 8.286514526259965e-05, 4.636391033219742e-05, 3.148961629902539e-05, 3.540031860286742e-05, 8.1113965120995e-05, 6.570302233902759e-05, 5.7520851308599364e-05, 0.0002814875050824133, 9.106085900743663e-05, 4.9331557397267034e-05, 2.483608185972581e-05, 0.0, 9.941346058256288e-05, 7.443429932512901e-05, 0.0, 0.0, 5.318583129454313e-05, 3.612325253765849e-05, 6.603711285742587e-05, 2.2661862351848074e-05, 8.614748449345279e-05, 2.5595740868719446e-05, 6.845564074479737e-05, 3.091094556582486e-05, 0.00010121457489878542, 2.7474791878451522e-05, 0.0001336719689881032, 3.9207998431680065e-05, 2.6736538153039946e-05, 5.379236148466918e-05, 8.952551477170994e-05, 4.387023185417535e-05, 0.0, 8.899172376968942e-05, 2.078353943676608e-05, 9.371017317640003e-05, 7.791499474073786e-05)
zgd <- c(3.889537145079735e-05, 7.405377785347719e-05, 7.340437245378583e-05, 1.2612249016244577e-05, 7.308979812597757e-05, 4.892427744957742e-05, 5.765017871555402e-05, 3.059008269519022e-05, 3.662154079029285e-05, 1.3065746838089265e-05, 0.00011510791366906474, 0.00011632641191182458, 0.0, 3.989229081480004e-05, 3.881121255930838e-05, 3.1282749127993365e-05, 5.372011818426001e-05, 5.6134573283685424e-05, 4.4211596701814886e-05, 3.788955195604812e-05, 2.886586034696764e-05, 1.8333150001833314e-05, 0.00014792534700820985, 7.977407980598944e-05, 4.6576618537494175e-05, 4.1729841009305755e-05, 1.320829480914014e-05, 7.068303371580708e-05, 0.0, 1.7334026694401108e-05, 5.515035365164279e-05, 0.00011185682326621924, 1.4225560487083191e-05, 9.469099173032006e-05, 0.00016127081401443375, 6.59352252347294e-05, 8.662508662508663e-05, 4.656197981538175e-05, 4.075228722212034e-05, 6.410667350471184e-05, 8.893790355573738e-05, 6.663002015558109e-05, 5.947778504728484e-05, 5.009016229212583e-05, 4.089109882560764e-05, 2.1752844184377107e-05)
boxplot(zgs,zgd,names=c("Similar","Dissimilar"),ylab="G4QuadsPerBase")
pval <- wilcox.test(zgs,zgd,alternative='less',correct=FALSE)$p.value
mtext(paste("p-val=",pval,sep=""),side=3,line=-8)


# Number of C[AT]G repeats found per base pair without Transposons
zcs <- c(0.0, 3.3137820194187624e-05, 2.7037988373664998e-05, 0.0, 0.0, 0.0, 2.6090586516384887e-05, 0.0, 3.7484069270560014e-05, 0.0, 0.0, 3.119638121977851e-05, 4.7885840157065554e-05, 1.5744808149512697e-05, 0.0, 4.311645755184754e-05, 0.0, 6.070723933829109e-05, 0.0, 7.399733609590054e-05, 0.0, 0.0, 0.0, 0.0, 3.4298257648511456e-05, 3.314605810503986e-05, 2.2012370952475292e-05, 0.0, 3.3738191632928474e-05, 1.180010620095581e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.318583129454313e-05, 0.0, 0.0, 0.0, 1.2662074554294975e-05, 0.0, 1.9603999215840032e-05, 0.0)
# Diff set of similar troughs
zcs2 <- c(0.0, 3.4406826314340766e-05, 0.0, 1.0959744857139726e-05, 7.452675510508272e-05, 2.3379234563860378e-05, 0.0, 4.243131430996075e-05, 3.4742730083729977e-05, 0.0, 0.0, 0.0, 4.311645755184754e-05, 0.0, 1.8400279684251202e-05, 4.7885840157065554e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.3202365863962822e-05, 0.0, 0.0, 3.797804868785842e-05, 0.0, 3.8669760247486464e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.2082123416987777e-05, 0.0, 2.962962962962963e-05)
zcd <- c(1.2965123816932451e-05, 7.405377785347719e-06, 3.670218622689292e-05, 1.2612249016244577e-05, 1.4617959625195516e-05, 3.669320808718306e-05, 2.882508935777701e-05, 0.0, 2.44143605268619e-05, 3.9197240514267794e-05, 8.633093525179857e-05, 0.0, 3.118665211289568e-05, 3.989229081480004e-05, 9.702803139827096e-06, 4.6924123691990054e-05, 0.0, 0.0, 1.1052899175453721e-05, 3.788955195604812e-05, 2.886586034696764e-05, 1.8333150001833314e-05, 0.0, 1.5954815961197887e-05, 0.0, 0.0, 0.0, 2.356101123860236e-05, 0.0, 0.0, 0.0, 0.0, 1.4225560487083191e-05, 4.208488521347558e-05, 0.0, 0.0, 5.197505197505198e-05, 0.0, 0.0, 1.602666837617796e-05, 0.0, 1.6657505038895273e-05, 1.9825928349094945e-05, 3.0054097375275497e-05, 8.178219765121528e-06, 0.0)
# With 3 repeats AND exactly within trough, so window=0
zcd0 <- c(9.075586671852716e-05, 3.7026888926738595e-05, 4.893624830252389e-05, 0.00012612249016244576, 4.385387887558655e-05, 6.115534681197177e-05, 8.647526807333103e-05, 9.177024808557065e-05, 0.0001220718026343095, 0.00013065746838089264, 0.00014388489208633093, 9.693867659318716e-05, 6.237330422579136e-05, 0.00011967687244440011, 6.791962197878968e-05, 0.0002033378693319569, 3.581341212284e-05, 5.6134573283685424e-05, 6.631739505272233e-05, 0.00015155820782419247, 8.659758104090293e-05, 0.00010999890001099989, 4.930844900273662e-05, 0.0001276385276895831, 9.315323707498835e-05, 9.736962902171343e-05, 0.00011887465328226125, 4.712202247720472e-05, 0.0001528584530724549, 0.0, 6.893794206455349e-05, 7.457121551081282e-05, 9.957892340958233e-05, 7.364854912358227e-05, 0.0, 2.637409009389176e-05, 8.662508662508663e-05, 4.656197981538175e-05, 6.112843083318051e-05, 0.00017629335213795757, 3.5575161422294955e-05, 8.328752519447637e-05, 0.00019825928349094947, 0.00011019835704267681, 4.089109882560764e-05, 8.701137673750843e-05)
zcs0 <- c(0.0, 3.3137820194187624e-05, 8.1113965120995e-05, 8.628127696289905e-05, 4.3073742246726395e-05, 0.0001240571655418817, 7.827175954915466e-05, 7.791499474073786e-05, 0.00011245220781168004, 5.2437010041687425e-05, 9.382916836080442e-05, 0.00012478552487911403, 0.0002873150409423933, 0.00011021365704658889, 4.967216371945162e-05, 0.0001293493726555426, 6.235061831029825e-05, 0.00018212171801487327, 0.0001089375476601771, 0.0001479946721918011, 0.00017113910186199342, 0.00012364378226329943, 8.339240295209106e-05, 0.00015549681231534753, 3.4298257648511456e-05, 0.00019887634863023914, 8.804948380990117e-05, 4.532372470369615e-05, 6.747638326585695e-05, 5.9000531004779045e-05, 2.318195516609871e-05, 4.121218781767728e-05, 0.00013140604467805518, 4.45573229960344e-05, 0.0, 5.4837789817719185e-05, 0.00014238675803150306, 5.8343057176196036e-05, 0.00021274332517817252, 3.612325253765849e-05, 5.119148173743889e-05, 0.00013428827215756492, 1.2662074554294975e-05, 5.216157569687865e-05, 0.00011762399529504019, 2.6736538153039946e-05)
boxplot(zcs0,zcd0,names=c("Similar","Dissimilar"),ylab="C[AT]GRepeatsPerBase")
pval <- wilcox.test(zcs0,zcd0,alternative='less',correct=FALSE)$p.value
# With 3 repeats and more AND exactly within trough, so window = 0, for C[CG]G Repeats
zcdCG0 <- c(2.5930247633864902e-05, 2.221613335604316e-05, 1.2234062075630972e-05, 0.0, 2.1926939437793273e-05, 2.446213872478871e-05, 0.0, 1.0196694231730073e-05, 1.220718026343095e-05, 2.613149367617853e-05, 0.0, 1.938773531863743e-05, 3.118665211289568e-05, 5.983843622220006e-05, 0.0, 1.5641374563996683e-05, 0.0, 0.0, 1.1052899175453721e-05, 0.0, 0.0, 0.0, 4.930844900273662e-05, 6.381926384479155e-05, 0.0, 0.0, 1.320829480914014e-05, 2.356101123860236e-05, 3.8214613268113725e-05, 1.7334026694401108e-05, 4.1362765238732097e-05, 0.0, 2.8451120974166383e-05, 3.156366391010668e-05, 2.0158851751804218e-05, 0.0, 1.7325017325017324e-05, 0.0, 0.0, 0.0, 1.7787580711147478e-05, 3.3315010077790546e-05, 0.0, 1.0018032458425165e-05, 8.178219765121528e-06, 3.262926627656566e-05)
zcsCG0 <- c(0.0, 3.3137820194187624e-05, 0.0, 5.7520851308599364e-05, 0.0, 0.0, 0.0, 0.0, 1.8742034635280007e-05, 0.0, 0.0, 0.0, 4.7885840157065554e-05, 1.5744808149512697e-05, 2.483608185972581e-05, 0.0, 0.0, 0.0, 1.556250680859673e-05, 2.4665778698633517e-05, 0.0, 3.091094556582486e-05, 0.0, 0.0, 0.0, 3.314605810503986e-05, 0.0, 0.0, 0.0, 0.0, 4.636391033219742e-05, 0.0, 0.0, 4.45573229960344e-05, 0.0, 1.0967557963543837e-05, 1.7798344753937883e-05, 2.9171528588098018e-05, 0.0, 3.612325253765849e-05, 1.2797870434359723e-05, 4.476275738585497e-05, 2.532414910858995e-05, 1.043231513937573e-05, 3.9207998431680065e-05, 2.6736538153039946e-05)
boxplot(zcsCG0,zcdCG0,names=c("Similar","Dissimilar"),ylab="C[CG]GRepeatsPerBase")
pval <- wilcox.test(zcsCG0,zcdCG0,alternative='less',correct=FALSE)$p.value


troughsInWhichtRNAsFound = unique(c('Trough3', 'Trough80', 'Trough81', 'Trough81', 'Trough82', 'Trough83', 'Trough84', 'Trough84', 'Trough84', 'Trough84', 'Trough84', 'Trough85', 'Trough85', 'Trough85', 'Trough85', 'Trough85', 'Trough86', 'Trough86', 'Trough88', 'Trough88', 'Trough88', 'Trough88', 'Trough89', 'Trough89', 'Trough90', 'Trough91', 'Trough91', 'Trough92', 'Trough93', 'Trough94', 'Trough95', 'Trough97', 'Trough97', 'Trough97', 'Trough98', 'Trough101', 'Trough101', 'Trough101', 'Trough102', 'Trough102', 'Trough103', 'Trough105', 'Trough105', 'Trough108', 'Trough109', 'Trough109', 'Trough109', 'Trough110', 'Trough111', 'Trough111', 'Trough111', 'Trough111', 'Trough113', 'Trough114', 'Trough116', 'Trough117', 'Trough117', 'Trough119', 'Trough119', 'Trough122', 'Trough123', 'Trough123', 'Trough129', 'Trough130', 'Trough132', 'Trough132', 'Trough135', 'Trough137', 'Trough138', 'Trough138', 'Trough139', 'Trough139', 'Trough139', 'Trough141', 'Trough141', 'Trough142', 'Trough143', 'Trough143', 'Trough146', 'Trough147', 'Trough150', 'Trough152', 'Trough154', 'Trough156', 'Trough158', 'Trough158', 'Trough158', 'Trough158', 'Trough160', 'Trough160', 'Trough160', 'Trough160', 'Trough5', 'Trough5', 'Trough7', 'Trough7', 'Trough7', 'Trough7', 'Trough7', 'Trough10', 'Trough10', 'Trough13', 'Trough13', 'Trough14', 'Trough14', 'Trough14', 'Trough15', 'Trough16', 'Trough17', 'Trough21', 'Trough21', 'Trough22', 'Trough22', 'Trough23', 'Trough24', 'Trough25', 'Trough25', 'Trough25', 'Trough27', 'Trough29', 'Trough30', 'Trough30', 'Trough33', 'Trough34', 'Trough34', 'Trough35', 'Trough35', 'Trough35', 'Trough35', 'Trough37', 'Trough37', 'Trough38', 'Trough39', 'Trough39', 'Trough41', 'Trough42', 'Trough42', 'Trough42', 'Trough44', 'Trough46', 'Trough46', 'Trough47', 'Trough47', 'Trough47', 'Trough48', 'Trough49', 'Trough50', 'Trough51', 'Trough51', 'Trough51', 'Trough52', 'Trough52', 'Trough52', 'Trough53', 'Trough53', 'Trough55', 'Trough55', 'Trough55', 'Trough55', 'Trough56', 'Trough56', 'Trough56', 'Trough56', 'Trough57', 'Trough57', 'Trough57', 'Trough57', 'Trough58', 'Trough60', 'Trough61', 'Trough61', 'Trough62', 'Trough62', 'Trough63', 'Trough65', 'Trough66', 'Trough68', 'Trough68', 'Trough73', 'Trough74', 'Trough74', 'Trough74', 'Trough74', 'Trough76', 'Trough76', 'Trough77', 'Trough77', 'Trough77'))
to.discard.num <- rep(0,length(troughsInWhichtRNAsFound))
for (i in 1:length(troughsInWhichtRNAsFound)){
to.discard.num[i] <- as.numeric(unlist(regmatches(troughsInWhichtRNAsFound[i], gregexpr('\\(?[0-9,.]+', troughsInWhichtRNAsFound[i]))))
}
tc.withtRNAs = troughComparisons[troughsInWhichtRNAsFound,]
tc.WithouttRNAs = troughComparisons[-to.discard.num,]
boxplot(tc.WithouttRNAs$Corr.Coeff,tc.withtRNAs$Corr.Coeff)
wilcox.test(tc.WithouttRNAs$Corr.Coeff,tc.withtRNAs$Corr.Coeff,alternative = "greater",correct = FALSE)
boxplot(tc.WithouttRNAs$Euclidean.Dist,tc.withtRNAs$Euclidean.Dist)
wilcox.test(tc.WithouttRNAs$Euclidean.Dist,tc.withtRNAs$Euclidean.Dist,alternative = "less",correct = FALSE)

# This function will calculate the probability of seeing x number of CNG repeats in a certain length of troughs 
calcpvalForNumCNGreps <- function(totalnumbases,samplenumbases,successesinpop,successesinsamp,numrepeat){
  popsize = totalnumbases/(numrepeat*3)
  samplesize = samplenumbases/(numrepeat*3)
  return(phyper(successesinsamp-1,successesinpop,popsize-successesinpop,samplesize,lower.tail = FALSE))
}