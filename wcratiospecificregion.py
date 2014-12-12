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
yeastsize = {'1':230218,'2':813184,'3':316620,'4':1531933,'5':576874,'6':270161,'7':1090940,'8':562643,'9':439888,'10':745751,'11':666816,'12':1078177,'13':924431,'14':784333,'15':1091291,'16':948066,'Total':12157105}

#######################################################################################################################
# FUNCTIONS
# This function grabs wcratio vals from the WT and MUT files for a particular region that has a chromosome and start
# and end coordinate and calculates vals for the feat indicated in each region
def featuresOfRegionsOfInterest(interestRegions,feat,valsWT,valsMUT,window):

    allregions = []
    
    for i in interestRegions:
        print i
        startend = [i[1]-window,i[2]+window]        
        valsforWT = [float(m[2]) for m in valsWT if m[0]==i[0]]
        valsforMUT = [float(m[2]) for m in valsMUT if m[0]==i[0]]
        print("Here")
        if feat=="Corr":
            featval = wcr.getCorrCoeffBetweenPairsOfCoordsForChromosome(valsforWT,startend,valsforMUT,startend)
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
        
        allregions.append([i,featval])

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
    feats = featuresOfRegionsOfInterest(trnas,feat,valsforWT,valsforMUT,window)
    return feats

##############################################################################

################################RANDOM########################################
# Collect random coordinates for chromosome
def randomCoordinatesforchromosome(chromosome,numpoints,window):
    coordRangeToSelect = range(1,yeastsize[chromosome]+1)
    
    randCoord = random.sample(coordRangeToSelect,numpoints)
    
    randCoordStartAndEnd = [[chromosome, i,i+window] for i in randCoord]
    
    return randCoordStartAndEnd
    