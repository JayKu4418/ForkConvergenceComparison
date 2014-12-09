# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 17:14:02 2014

@author: jayashreekumar
"""

#import numpy as np
#import matplotlib.pyplot as plt
import watsoncrickRatio as wcr
import random
yeastsize = {'1':230218,'2':813184,'3':316620,'4':1531933,'5':576874,'6':270161,'7':1090940,'8':562643,'9':439888,'10':745751,'11':666816,'12':1078177,'13':924431,'14':784333,'15':1091291,'16':948066,'Total':12157105}
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
    
def featuresOftRNAGenesForChromosome(tRNAs,valsforWT,valsforMUT):

    allfeatures = []
    
    for i in tRNAs:
        
        startend = [i[1]-1000,i[2]+1000]        
        
        featureperpair = {}
        corr = wcr.getCorrCoeffBetweenPairsOfCoordsForChromosome(valsforWT,startend,valsforMUT,startend)
        dist = wcr.getEuclideanDistanceBetweenPairsOfCoordsForChromosome(valsforWT,startend,valsforMUT,startend)
        
        firstder = wcr.getCorrCoeffAndEucDistBetweenFirstDerOfPairsOfCoordsForChromosome(valsforWT,startend,valsforMUT,startend)
        
    
        secder = wcr.getCorrCoeffAndEucDistBetweenSecDerOfPairsOfCoordsForChromosome(valsforWT,startend,valsforMUT,startend)
        
        featureperpair['Corr'] = corr
        featureperpair['Dist'] = dist
        featureperpair['FirstDer'] = [firstder[0][0],firstder[1]]
        featureperpair['SecDer'] = [secder[0][0],secder[1]]

        allfeatures.append([i,featureperpair])

    return allfeatures
    
##############################################################################

################################RANDOM########################################
# Collect random coordinates for chromosome
def randomCoordinatesforchromosome(chromosome,numpoints,window):
    coordRangeToSelect = range(1,yeastsize[chromosome]+1)
    
    randCoord = random.sample(coordRangeToSelect,numpoints)
    
    randCoordStartAndEnd = [[chromosome, i,i+window] for i in randCoord]
    
    return randCoordStartAndEnd
    