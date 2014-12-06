# -*- coding: utf-8 -*-
"""
Created on Sun Aug 17 22:06:56 2014

@author: jayashreekumar
"""

#######################################################################################################################
# IMPORT

#import smoothing as s
import replicationforktermination as rft
import smoothing as s
from scipy.stats import pearsonr
from operator import itemgetter

# dictionary that converts values from the form 'chrmIII' to '3'
chrmconvert = {'chrI':'1','chrII':'2','chrIII':'3','chrIV':'4','chrV':'5','chrVI':'6','chrVII':'7','chrVIII':'8','chrIX':'9','chrX':'10','chrXI':'11','chrXII':'12','chrXIII':'13','chrXIV':'14','chrXV':'15','chrXVI':'16'}

#######################################################################################################################
# FUNCTIONS

# This function grabs all the similar troughs for a chromosome between two oemfiles
def troughsForChromosomeForTwoOEMfiles(oemfile1,oemfile2,chromosome,window,withoem):
    
    # Grab origins and terminatons for oemfile1 for chromosome
    oandt1 = [i for i in rft.getListOfOriginsAndTerminationPoints(oemfile1,chromosome) if (i[2]=='O' and i[1] > 0.1) or (i[2]=='T' and i[1] < -0.1)]
    #print oandt1
    
    # Grab origins and terminatons for oemfile2 for chromosome
    oandt2 = [i for i in rft.getListOfOriginsAndTerminationPoints(oemfile2,chromosome) if (i[2]=='O' and i[1] > 0.1) or (i[2]=='T' and i[1] < -0.1)]
    #print oandt2
    
    # Grab origins that are close to each other in the two oemfiles    
    closeorigins = rft.graboriginsThatAreClose(oemfile1,oemfile2,chromosome,window,0.1)
    #print(closeorigins)
    
    # Number of common origins
    numcomorgs = len(closeorigins[0])
    
    # Empty list to store common start and end points of troughs for each oemfile
    startendforboth = []
    
    # Grab start and end points of troughs to compare for each oemfile   
    for i in range(0,numcomorgs):
        org1 = closeorigins[0][i]
        org2 = closeorigins[1][i]
        #print org1
        #print org2
        # Look through origins and terminations for oemfile1 and grab appropriate start and end points
        if withoem:
            startendoem1 = [[oandt1[i][0],oandt1[i][1],oandt1[i+1][0],oandt1[i+1][1],oandt1[i+2][0],oandt1[i+2][1]] for i in range(len(oandt1)) if i<=len(oandt1)-3 and oandt1[i][0]==org1 and oandt1[i][2]=='O' and oandt1[i+1][2]=='T' and oandt1[i+2][2]=='O']
            startendoem2 = [[oandt2[i][0],oandt2[i][1],oandt2[i+1][0],oandt2[i+1][1],oandt2[i+2][0],oandt2[i+2][1]] for i in range(len(oandt2)) if i<=len(oandt2)-3  and oandt2[i][0]==org2 and oandt2[i][2]=='O' and oandt2[i+1][2]=='T' and oandt2[i+2][2]=='O']
            if len(startendoem1)!=0 and len(startendoem2)!=0:        
                if (startendoem1[0][4] >= startendoem2[0][4]-window and startendoem1[0][4] <= startendoem2[0][4]+window) or (startendoem2[0][4] >= startendoem1[0][4]-window and startendoem2[0][4] <= startendoem1[0][4]+window):
                    startendforboth.append([startendoem1[0],startendoem2[0]])
        else:
            startendoem1 = [[oandt1[i][0],oandt1[i+1][0],oandt1[i+2][0]] for i in range(len(oandt1)) if i<=len(oandt1)-3 and oandt1[i][0]==org1 and oandt1[i][2]=='O' and oandt1[i+1][2]=='T' and oandt1[i+2][2]=='O']
            #print startendoem1
            # Look through origins and terminations for oemfile2 and grab appropriate start and end points
            startendoem2 = [[oandt2[i][0],oandt2[i+1][0],oandt2[i+2][0]] for i in range(len(oandt2)) if i<=len(oandt2)-3  and oandt2[i][0]==org2 and oandt2[i][2]=='O' and oandt2[i+1][2]=='T' and oandt2[i+2][2]=='O']
            #print startendoem2
            if len(startendoem1)!=0 and len(startendoem2)!=0:        
                if (startendoem1[0][2] >= startendoem2[0][2]-window and startendoem1[0][2] <= startendoem2[0][2]+window) or (startendoem2[0][2] >= startendoem1[0][2]-window and startendoem2[0][2] <= startendoem1[0][2]+window):
                    startendforboth.append([startendoem1[0],startendoem2[0]])
        
    return startendforboth

# This function calculates the gradient for each trough for chromosome
def gradientCorrelationOfTroughsForChromosome(troughsInput,oemfile1,oemfile2,chromosome,window,freq):
    
    # Pearson r gradient for trough 
    troughPlusGradCorr = []
    # Troughs for chromosome    
    if troughsInput == []:    
        troughs = troughsForChromosomeForTwoOEMfiles(oemfile1,oemfile2,chromosome,window,False)
    else:
        troughs = troughsInput
    for i in troughs:
        grad1 = s.gradientOfStrain(oemfile1,i[0][0],i[0][2],chromosome,freq)
        grad2 = s.gradientOfStrain(oemfile2,i[1][0],i[1][2],chromosome,freq)
        lenofdatatograb = min(len(grad1),len(grad2))
        p = pearsonr(grad1[0:lenofdatatograb],grad2[0:lenofdatatograb])
        troughPlusGradCorr.append([i[0],i[1],p])
        
    return troughPlusGradCorr
    
# This function calculates the Midpoint difference between each pair of troughs for chromosome
def midpointDiffOfTroughsForChromosome(troughsInput,oemfile1,oemfile2,chromosome,window,freq):
    
    # Midpt diff for troughs
    troughPlusMidPtDiff = []
    # Troughs for chromosome    
    if troughsInput == []:    
        troughs = troughsForChromosomeForTwoOEMfiles(oemfile1,oemfile2,chromosome,window,False)
    else:
        troughs = troughsInput
    for i in troughs:
        midpt1 = s.midpointOfStrain(oemfile1,i[0][0],i[0][2],chromosome,freq)
        midpt2 = s.midpointOfStrain(oemfile2,i[1][0],i[1][2],chromosome,freq)
        midptdiff = abs(midpt1-midpt2)
        troughPlusMidPtDiff.append([i[0],i[1],[midpt1,midpt2,midptdiff]])
        
    return troughPlusMidPtDiff

# This function calculates the difference in OEM of peaks between each pair of troughs for chromosome
def oemPeakDiffOfTroughsForChromosome(troughsInput,oemfile1,oemfile2,chromosome,window,freq):
    
    # oem Peak diff for troughs
    troughPlusOEMPeakDiff = []
    # Troughs for chromosome    
    if troughsInput == []:    
        troughs = troughsForChromosomeForTwoOEMfiles(oemfile1,oemfile2,chromosome,window,True)
    else:
        troughs = troughsInput
    for i in troughs:
        peak1oemdiff = abs(i[0][1]-i[1][1])
        troughoemdiff = abs(i[0][3]-i[1][3])
        peak2oemdiff = abs(i[0][5]-i[1][5])
        troughPlusOEMPeakDiff.append([[i[0][0],i[0][2],i[0][4]],[i[1][0],i[1][2],i[1][4]],[peak1oemdiff,troughoemdiff,peak2oemdiff]])
        
    return troughPlusOEMPeakDiff
   
# This function assign labels for pairs of troughs for each chromosome: S for similar, D for different, NA for Not Applicable
def assignLabelsTroughsForChromosome(oemfile1,oemfile2,chromosome,window,freq,additionaldet):
    
    # List to assign troughs if they are similar, different or not applicable
    finaltroughChromosome = []
    # Get troughs for chromosome without oems
    troughsWithoutOEM = troughsForChromosomeForTwoOEMfiles(oemfile1,oemfile2,chromosome,window,False)
    # Get troughs for chromosome with oem
    troughsWithOEM = troughsForChromosomeForTwoOEMfiles(oemfile1,oemfile2,chromosome,window,True)
    
    # Gradient correlation of troughs
    grad = gradientCorrelationOfTroughsForChromosome(troughsWithoutOEM,oemfile1,oemfile2,chromosome,window,freq)
    
    # MidPoint Diff
    midptdiff = midpointDiffOfTroughsForChromosome(troughsWithoutOEM,oemfile1,oemfile2,chromosome,window,freq)
    
    # OEM Peak Diffs 
    oemdiff = oemPeakDiffOfTroughsForChromosome(troughsWithOEM,oemfile1,oemfile2,chromosome,window,freq)
    
    # Number of troughs 
    numtroughs = len(troughsWithOEM)
    
    for i in range(numtroughs):
        if oemdiff[i][2][0] < 0.15 and oemdiff[i][2][2] < 0.15 and (oemdiff[i][2][0] < 0.1 or oemdiff[i][2][2] < 0.1):
            if grad[i][2][0] > 0.95 and midptdiff[i][2][2] < 0.1 and oemdiff[i][2][1] < 0.1 :
                if additionaldet:
                    score = (1- grad[i][2][0]) + midptdiff[i][2][2] + oemdiff[i][2][1]
                    finaltroughChromosome.append([chromosome,troughsWithoutOEM[i],grad[i][2][0],midptdiff[i][2][2],oemdiff[i][2][1],score,'S'])
                else:
                    finaltroughChromosome.append([troughsWithoutOEM[i],'S'])
            else:
                if additionaldet:
                    score = (1- grad[i][2][0]) + midptdiff[i][2][2] + oemdiff[i][2][1]
                    finaltroughChromosome.append([chromosome,troughsWithoutOEM[i],grad[i][2][0],midptdiff[i][2][2],oemdiff[i][2][1],score,'D'])
                else:
                    finaltroughChromosome.append([troughsWithoutOEM[i],'D'])
        else:
            if additionaldet:
                score = (1- grad[i][2][0]) + midptdiff[i][2][2] + oemdiff[i][2][1]
                finaltroughChromosome.append([chromosome,troughsWithoutOEM[i],grad[i][2][0],midptdiff[i][2][2],oemdiff[i][2][1],score,'NA'])
            else:
                finaltroughChromosome.append([troughsWithoutOEM[i],'NA'])

    return finaltroughChromosome
                
# Rank troughs based on label
def rankTroughsBasedOnLabel(label,oemfile1,oemfile2,window,freq,writefile=''):
    # assign labels to troughs for oemfile1 and oemfile2
    labeledtroughs = []
    for c in range(1,17):
        chromosome = str(c)
        print chromosome
        labeledtroughs.extend(assignLabelsTroughsForChromosome(oemfile1,oemfile2,chromosome,window,freq,True))
        
    troughsToRank = [i for i in labeledtroughs if i[6]==label]
    #print troughsToRank
    if label == 'S':    
        rankedTroughs = sorted(troughsToRank,key=itemgetter(5))
    else:
        rankedTroughs = sorted(troughsToRank,key=itemgetter(5),reverse=True)
    
    if writefile:
        with open(writefile,'w') as fw:
            for i in rankedTroughs:
                fw.write(i[0] + '\t' + str(i[1][0]) + '\t' + str(i[1][1]) + '\t' + str(i[2]) + '\t' + str(i[3]) + '\t' + str(i[4]) + '\t' + str(i[5]) + '\t' + i[6] + '\n')
    else:
        return rankedTroughs
    
# For Troughs that are different, find difference and cluster the troughs
"""
def clusterDissimilarTroughsBasedOnDifferences(oemfile1,oemfile2,window,freq):
    dissimilartroughs = []
    for c in range(1,17):
        chromosome = str(c)
        print chromosome
        alltroughs = assignLabelsTroughsForChromosome(oemfile1,oemfile2,chromosome,window,freq,False)
        dissimilartroughs = [[chromosome,i[0]] for i in alltroughs if i[1]=='D']
        
    
"""
# tRNA genes found within a single trough
def tRNAGeneFoundWithinTrough(trnagenes,start,end,chromosome,window):
    with open(trnagenes) as f:
        trnas = [line.strip().split('\t') for line in f if line.strip().split('\t')[1]==chromosome]
        
    trnasFoundWithinTrough = [i for i in trnas if (int(i[2]) >= start-window and int(i[2]) <= end+window) or (int(i[3]) >= start-window and int(i[3]) <= end+window)]
    
    return trnasFoundWithinTrough
    
# tRNA genes found within troughs of a certain label
def tRNAGenesFoundwithinLabelledTroughs(trnagenes,labelledtroughsfile,window):
    
    # Variable to store tRNAs found in all troughs
    trnasAllTroughs = []
    
    # Grab the first 50 troughs and store in variable troughs     
    with open(labelledtroughsfile) as f:
        troughs = [line.strip().split('\t') for line in f][0:48]
        
    for i in troughs:
        # trna genes Found in first Trough
        firsttroughStart = int(i[1].split('[')[1].split(']')[0].split(',')[0])
        secondtroughStart = int(i[2].split('[')[1].split(']')[0].split(',')[0])
        firsttroughEnd = int(i[1].split('[')[1].split(']')[0].split(',')[2])
        secondtroughEnd = int(i[2].split('[')[1].split(']')[0].split(',')[2])
        trnasFoundInTrough = tRNAGeneFoundWithinTrough(trnagenes,min(firsttroughStart,secondtroughStart),max(firsttroughEnd,secondtroughEnd),i[0],window)
        trnasAllTroughs.extend(trnasFoundInTrough)
        
    return trnasAllTroughs
        
# tRNA Genes found within troughs of a certain label over a number of windows
def tRNAGenesFoundWithinLabelledTroughsOverRangeWindow(trnagenes,labelledtroughsfile,interval,maxwindow):
    totalnumOftRNAGenesFoundPerInterval = []    
    for w in range(0,maxwindow+interval,interval):
        #print w
        totalnumOftRNAGenesFoundPerInterval.append(len(tRNAGenesFoundwithinLabelledTroughs(trnagenes,labelledtroughsfile,w)))
        
    return totalnumOftRNAGenesFoundPerInterval
    
def comparelengthOfTroughInterval(labelledtroughsfile):
     # Grab the first 50 troughs and store in variable troughs
    firsttroughs = []
    secondtroughs = []     
    with open(labelledtroughsfile) as f:
        troughs = [line.strip().split('\t') for line in f][0:50]
        
    for i in troughs:
        firsttroughStart = int(i[1].split('[')[1].split(']')[0].split(',')[0])
        secondtroughStart = int(i[2].split('[')[1].split(']')[0].split(',')[0])
        firsttroughEnd = int(i[1].split('[')[1].split(']')[0].split(',')[2])
        secondtroughEnd = int(i[2].split('[')[1].split(']')[0].split(',')[2])
        firsttroughs.append(firsttroughEnd-firsttroughStart)
        secondtroughs.append(secondtroughEnd-secondtroughStart)
        
    return [firsttroughs,secondtroughs]
    
# Remove troughs found within Ty1 elements
def removeTroughsFoundWithinTransposableElements(labelledtroughsfile,transfile,writefile):
    with open(labelledtroughsfile) as f:
        troughs = [line.strip().split('\t') for line in f]
    with open(writefile,'w') as fw:
        for i in troughs:
            chrm = i[0]
            firsttroughStart = int(i[1].split('[')[1].split(']')[0].split(',')[0])
            secondtroughStart = int(i[2].split('[')[1].split(']')[0].split(',')[0])
            firsttroughEnd = int(i[1].split('[')[1].split(']')[0].split(',')[2])
            secondtroughEnd = int(i[2].split('[')[1].split(']')[0].split(',')[2])
            with open(transfile) as f:
                transValid = [line.strip().split(' ') for line in f if chrmconvert[line.strip().split(' ')[1]]==chrm]
            #transFoundInTrough = [j for j in transValid if (j[4]=='W' and ((int(j[2]) >= firsttroughStart and int(j[2]) <= firsttroughEnd) or (int(j[2]) >= secondtroughStart and int(j[2]) <= secondtroughEnd) )) or (j[4]=='C' and ((int(j[3]) >= firsttroughStart and int(j[3]) <= firsttroughEnd) or (int(j[3]) >= secondtroughStart and int(j[3]) <= secondtroughEnd)))]
            transFoundInTrough = [j for j in transValid if ((int(j[2]) >= firsttroughStart and int(j[2]) <= firsttroughEnd) or (int(j[2]) >= secondtroughStart and int(j[2]) <= secondtroughEnd) )]
            if len(transFoundInTrough)==0:
                fw.write(i[0] + '\t' + i[1] + '\t' + i[2] + '\t' + i[3] + '\t' + i[4] + '\t' + i[5] + '\t' + i[6] + '\t' + i[7] + '\n')