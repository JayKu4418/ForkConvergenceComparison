# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 17:14:02 2014

@author: jayashreekumar
"""
#######################################################################################################################
# IMPORT
#import numpy as np
#import matplotlib.pyplot as plt
import watsoncrickRatio as wcr
import random
import ManipulateSeqData.findrepeats as fr
yeastsize = {'1':230218,'2':813184,'3':316620,'4':1531933,'5':576874,'6':270161,'7':1090940,'8':562643,'9':439888,'10':745751,'11':666816,'12':1078177,'13':924431,'14':784333,'15':1091291,'16':948066,'Total':12157105}
trnatrans = {'tA':'Alanine','tF':'Phenylalanine','tL':'Leucine','tI':'Isoleucine','tM':'Methionine','tV':'Valine','tS':'Serine','tP':'Proline','tT':'Threonine','tY':'Tyrosine','tH':'Histidine','tQ':'Glutamine','tN':'Asparagine','tK':'Lysine','tD':'Aspartic acid','tE':'Glutamic acid','tC':'Cysteine','tW':'Tryptophan','tR':'Arginine','tG':'Glycine','tX':'Undetermined'}
#import matplotlib.pyplot as plt
#######################################################################################################################
# FUNCTIONS
# This function grabs wcratio vals from the WT and MUT files for a particular region that has a chromosome and start
# and end coordinate and calculates vals for the feat indicated in each region
def featureOfRegionsOfInterestForChromosome(interestRegions,feat,valsforWT,valsforMUT,window):

    allregions = []
    
    for i in interestRegions:
        startend = [i[1]-window,i[2]+window]        
        #valsforWT = [float(m[2]) for m in valsWT if m[0]==i[0]]
        #valsforMUT = [float(m[2]) for m in valsMUT if m[0]==i[0]]
        if feat=="Corr":
            featval = wcr.getCorrCoeffBetweenPairsOfCoordsForChromosome(valsforWT,startend,valsforMUT,startend)[0]
        elif feat=="EucDist":
            featval = wcr.getEuclideanDistanceBetweenPairsOfCoordsForChromosome(valsforWT,startend,valsforMUT,startend)
        elif feat=="FirstDerCorr":
            firstder = wcr.getCorrCoeffAndEucDistBetweenFirstDerOfPairsOfCoordsForChromosome(valsforWT,startend,valsforMUT,startend)
            featval = firstder[0][0]
        elif feat=="FirstDerEucDist":
            firstder = wcr.getCorrCoeffAndEucDistBetweenFirstDerOfPairsOfCoordsForChromosome(valsforWT,startend,valsforMUT,startend)
            featval = firstder[1]
        elif feat=="SecDerCorr":
            secder = wcr.getCorrCoeffAndEucDistBetweenSecDerOfPairsOfCoordsForChromosome(valsforWT,startend,valsforMUT,startend)
            featval = secder[0][0]
        elif feat=="SecDerEucDist":
            secder = wcr.getCorrCoeffAndEucDistBetweenSecDerOfPairsOfCoordsForChromosome(valsforWT,startend,valsforMUT,startend)
            featval = secder[1]
        
        allregions.append([[i[0],i[1]-window,i[2]+window],featval])

    return allregions

################################TRNA GENES##################################
# This function takes a tRNAfile which has order tRNA name, chromosome, 
# start, end, strand. Note that for crick, the start is larger than the end
def extractStartEndChromStrandFortRNAs(tRNAfile):
    with open(tRNAfile) as f:
        trnas = [line.strip().split('\t') for line in f]
        
    tRNAsToRet = [[i[1],int(i[2]),int(i[3]),i[4]] if i[4]=='W' else [i[1],int(i[3]),int(i[2]),i[4]] for i in trnas]
    
    return tRNAsToRet
    
# This function returns troughs within which tRNA gene are found
def troughsInWhichtRNAGeneFoundForChrom(troughsfile,tRNAGene):
    
    start = tRNAGene[1]
    end = tRNAGene[2]
    chromosome = tRNAGene[0]
    
    with open(troughsfile) as f:
        troughs = [line.strip().split('\t') for line in f]
        
    justchrmstartend = [[i[0],i[1],min(int(i[2]),int(i[4])),max(int(i[3]),int(i[5]))] for i in troughs]
    
    troughsFound = [j[0].split('"')[1] for j in justchrmstartend if chromosome == j[1] and start >= j[2] and end <= j[3]]
    
    return troughsFound
    
# This function returns all troughs in which all tRNA genes are found
def troughsWhichAllContaintRNAGenes(troughsfile,tRNAfile):
    
    trnas = extractStartEndChromStrandFortRNAs(tRNAfile)
   
    alltroughs = []
    
    for t in trnas:
        alltroughs.extend(troughsInWhichtRNAGeneFoundForChrom(troughsfile,t))
        
    return alltroughs
    
# This function returns features around a certain window around tRNA genes
def maintRNAgenefeats(tRNAfile,feat,wcratiofileWT,wcratiofileMUT,window):
    trnas = extractStartEndChromStrandFortRNAs(tRNAfile)
    with open(wcratiofileWT) as f:
        valsforWT = [line.strip().split('\t') for line in f]
    with open(wcratiofileMUT) as f:
        valsforMUT = [line.strip().split('\t') for line in f]
    feats = featureOfRegionsOfInterestForChromosome(trnas,feat,valsforWT,valsforMUT,window)
    return feats
    
def maintRNAgenefeatOverWindowSizes(chromosome,tRNAfile,inc,maxwindow,feat,wcratiofileWT,wcratiofileMUT):
    trnasforchrom = [i for i in extractStartEndChromStrandFortRNAs(tRNAfile) if i[0]==chromosome]
    if len(trnasforchrom)!=0:
        print(len(trnasforchrom))
        with open(wcratiofileWT) as f:
            valsforWT = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
        with open(wcratiofileMUT) as f:
            valsforMUT = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
        featsThroughAllWindow = {}   
        for r in range(inc,maxwindow+inc,inc):
            feats = featureOfRegionsOfInterestForChromosome(trnasforchrom,feat,valsforWT,valsforMUT,r)
            featsThroughAllWindow[r] = [i[1] for i in feats]
        return featsThroughAllWindow
    else:
        return {}

def maintRNAgenefeatOverWindowSizesAllChrom(tRNAfile,inc,maxwindow,feat,wcratiofileWT,wcratiofileMUT):
    featforallchrom = {}
    for c in range(1,17):
        chromosome = str(c)
        print(chromosome)
        featforallchrom[chromosome] = maintRNAgenefeatOverWindowSizes(chromosome,tRNAfile,inc,maxwindow,feat,wcratiofileWT,wcratiofileMUT)
    return featforallchrom
##############################################################################

################################CAG REPEATS##################################
def mainCNGrepeatsfeatForChromosome(fastafile,repeat,numrepeats,chromosome,window,feat,wcratiofileWT,wcratiofileMUT):
    locsOfcngrepeats = fr.getRepeatsForChromsome(fastafile,chromosome,repeat,numrepeats)
    with open(wcratiofileWT) as f:
        valsforWT = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
    with open(wcratiofileMUT) as f:
        valsforMUT = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
    properlocscngrepeats = [[chromosome,i[0],i[1]] for i in locsOfcngrepeats]
    feats = featureOfRegionsOfInterestForChromosome(properlocscngrepeats,feat,valsforWT,valsforMUT,window)
    return feats
            
def mainCNGrepeatsfeatForAll(fastafile,repeat,numrepeats,window,feat,wcratiofileWT,wcratiofileMUT):
    featforallchrom = []
    for c in range(1,17):
        chromosome = str(c)
        print(chromosome)
        featforallchrom.extend(mainCNGrepeatsfeatForChromosome(fastafile,repeat,numrepeats,chromosome,window,feat,wcratiofileWT,wcratiofileMUT))
    return featforallchrom

def mainCNGrepeatsfeatOverWindowSizes(chromosome,fastafile,repeat,numrepeats,inc,maxwindow,feat,wcratiofileWT,wcratiofileMUT):
    locsOfcngrepeats = fr.getRepeatsForChromsome(fastafile,chromosome,repeat,numrepeats)
    if len(locsOfcngrepeats)!=0:
        print(len(locsOfcngrepeats))
        with open(wcratiofileWT) as f:
            valsforWT = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
        with open(wcratiofileMUT) as f:
            valsforMUT = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
        properlocscngrepeats = [[chromosome,i[0],i[1]] for i in locsOfcngrepeats]
        featsThroughAllWindow = {}   
        for r in range(inc,maxwindow+inc,inc):
            feats = featureOfRegionsOfInterestForChromosome(properlocscngrepeats,feat,valsforWT,valsforMUT,r)
            featsThroughAllWindow[r] = [i[1] for i in feats]
        return featsThroughAllWindow
    else:
        return {}
def mainCNGrepeatsfeatOverWindowSizesAllChrom(fastafile,repeat,numrepeats,inc,maxwindow,feat,wcratiofileWT,wcratiofileMUT):
    featforallchrom = {}
    for c in range(1,17):
        chromosome = str(c)
        print(chromosome)
        featforallchrom[chromosome] = mainCNGrepeatsfeatOverWindowSizes(chromosome,fastafile,repeat,numrepeats,inc,maxwindow,feat,wcratiofileWT,wcratiofileMUT)
    return featforallchrom
    
def grabJustDataForFeatOverWindowSizesAllChrom(feats,inc,maxwindow):
    featsforallr = []
    for r in range(inc,maxwindow+inc,inc):
        print r
        featsforr = []
        for c in range(1,17):
            chrom = str(c)
            if len(feats[chrom])!=0:
                featsforr.extend(feats[chrom][r])
        featsforallr.append(featsforr)
    
    return featsforallr
    #plt.boxplot(featsforallr)
    #plt.show()
def grabJustDataFor2FeatOverWindowSizesAllChrom(feats1,feats2,inc,maxwindow):
    featsforallr = []
    for r in range(inc,maxwindow+inc,inc):
        print r
        featsforr = []
        for c in range(1,17):
            chrom = str(c)
            if len(feats1[chrom])!=0:
                featsforr.extend(feats1[chrom][r])
            if len(feats2[chrom])!=0:
                featsforr.extend(feats2[chrom][r]) 
        featsforallr.append(featsforr)
    
    return featsforallr

# This function combines the two functions above to grab feats for both repeats
# and combines them into one dataset
def grabfeatForCNGrepeatsOverWindows(fastafile,repeat1,repeat2,numrepeats,inc,maxwindow,feat,wcratiofileWT,wcratiofileMUT):
    feats1 =  mainCNGrepeatsfeatOverWindowSizesAllChrom(fastafile,repeat1,numrepeats,inc,maxwindow,feat,wcratiofileWT,wcratiofileMUT)  
    feats2 =  mainCNGrepeatsfeatOverWindowSizesAllChrom(fastafile,repeat2,numrepeats,inc,maxwindow,feat,wcratiofileWT,wcratiofileMUT)
    return grabJustDataFor2FeatOverWindowSizesAllChrom(feats1,feats2,inc,maxwindow)
##############################################################################
################################RANDOM########################################
# Collect random coordinates for chromosome
def randomCoordinates(numpoints,size):
    randCoordStartAndEnd = []
    for i in range(numpoints):
        chromosome = str(random.randint(1,16))
        randCoord = random.randint(1,yeastsize[chromosome]+1)
    
        randCoordStartAndEnd.append([chromosome, randCoord,randCoord+size])
    
    return randCoordStartAndEnd
 
#def mainrandomfeatforchromosome(chromosome,numpoints,size,feat,wcratiofileWT,wcratiofileMUT,window):
    #randcoords = randomCoordinatesforchromosome(chromosome,numpoints,size)
    #with open(wcratiofileWT) as f:
        #valsforWT = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
    #with open(wcratiofileMUT) as f:
        #valsforMUT = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
    #feats = featureOfRegionsOfInterestForChromosome(randcoords,feat,valsforWT,valsforMUT,window)
    #return feats

def mainrandomfeatOverWindowSizes(numpoints,size,inc,maxwindow,feat,wcratiofileWT,wcratiofileMUT):
    featforallchrom = {}
    randcoords = randomCoordinates(numpoints,size)
    for c in range(1,17):
        chromosome = str(c)
        print chromosome
        coords = [i for i in randcoords if i[0]==chromosome]
        if len(coords)!=0:
            print "Here"
            with open(wcratiofileWT) as f:
                valsforWT = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
            with open(wcratiofileMUT) as f:
                valsforMUT = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
            featsThroughAllWindow = {}   
            for r in range(inc,maxwindow+inc,inc):
                feats = featureOfRegionsOfInterestForChromosome(coords,feat,valsforWT,valsforMUT,r)
                featsThroughAllWindow[r] = [i[1] for i in feats]
            
            featforallchrom[chromosome] = featsThroughAllWindow
        else:
            featforallchrom[chromosome] = []
    return featforallchrom