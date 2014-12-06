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
    
# This function returns wcratio dictionary near a single tRNA gene for both WT and MUT. The tRNA gene contains
# chromosome, start, end, and strand
def wcratioNeartRNAGeneForChrom(wcratiofileWT,wcratiofileMUT,tRNAGene):
    
    start = tRNAGene[1] - 1000
    end = tRNAGene[2] + 1000
    chromosome = tRNAGene[0]
    
    with open(wcratiofileWT) as f:
        valsforchromWT = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome][start:end]
        
    with open(wcratiofileMUT) as f:
        valsforchromMUT = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome][start:end]
        
    wcratio = {'WT':valsforchromWT,'MUT':valsforchromMUT}
    
    return wcratio
    
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
    