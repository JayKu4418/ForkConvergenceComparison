# -*- coding: utf-8 -*-
"""
Created on Mon Jul 14 14:24:23 2014

@author: jayashreekumar
"""

#######################################################################################################################
# IMPORTS
#import matplotlib.pyplot as plt
cerevisiaechrmconvert = {'1':'chrI','2':'chrII','3':'chrIII','4':'chrIV','5':'chrV','6':'chrVI','7':'chrVII','8':'chrVIII','9':'chrIX','10':'chrX','11':'chrXI','12':'chrXII','13':'chrXIII','14':'chrXIV','15':'chrXV','16':'chrXVI'}
#import numpy as np
from operator import itemgetter
#import bokeh.plotting as bk
import matplotlib.pyplot as plt
import ManipulateSeqData.originscomparison as oc
import ManipulateSeqData.replicationtimes as rt
#######################################################################################################################
# FUNCTIONS

# This is a dictionary of yeast chromosome sizes
yeastsize = {'1':230218,'2':813184,'3':316620,'4':1531933,'5':576874,'6':270161,'7':1090940,'8':562643,'9':439888,'10':745751,'11':666816,'12':1078177,'13':924431,'14':784333,'15':1091291,'16':948066,'Total':12157105}

# This function returns the key to the maximum value in the dictionary
def findminInDict(min2finddict):
    return min(min2finddict.iteritems(), key=itemgetter(1))

# Get a list of bases which have min OEM in this dictionary as long as OEM is negative
def basesWminoem(splitlist):
    # Empty list to store all the bases with max positive OEM values
    minoembases = []
    # Go through every dictionary in the split list
    for dictoem in splitlist:
        # If first value in dictionary is positive, means all the values in that dictionary are positive since they
        # should all have the same sign
        if dictoem.values()[0] < 0:
            # Multiply each the base with the max OEM value by perbasebox since OEMs were calculated for every perbasebox value
            minoembases.append(findminInDict(dictoem))

    return minoembases
    
    
# This function returns the predicted termination points for the oemfile
def getPredictedTerminationPointsForOEMList(oemfile,chromosome):
    
    #Grab the OEM for sgrcomblist
    with open(oemfile) as f:
        oemc = [[int(line.strip().split('\t')[1]),float(line.strip().split('\t')[2])] for line in f.readlines() if line.strip().split('\t')[0]==chromosome]
    
    # Split up the OEM values based on sign, so consecutive positive numbers and negative numbers get grouped
    s = oc.splitListbasedOnSign(oemc)
    # Grab bases which have the max positive OEM value per list
    maxBases = basesWminoem(s)

    return maxBases
    
def getListOfOriginsAndTerminationPoints(oemfile,chromosome):
    predictedorigins = oc.getPredictedOriginsForOEMList(oemfile,chromosome)
    predictedterms = getPredictedTerminationPointsForOEMList(oemfile,chromosome)
    
    po = [[i[0],i[1],'O'] for i in predictedorigins]
    pt = [[i[0],i[1],'T'] for i in predictedterms]
    
    allp = po + pt
    
    allp.sort(key=lambda x: x[0])    
    
    return allp
    
def calculateTermBaseDistFromNeighbouringOrigins(oemfile,chromosome,oemlim):
    x = getListOfOriginsAndTerminationPoints(oemfile,chromosome)
    
    distfromorg = []
    for i in range(1,len(x)-1):
        if x[i][2]=='T' and x[i-1][2]=='O' and x[i+1][2]=='O' and x[i][1]<=oemlim:
            rightfork = x[i][0]-x[i-1][0]
            leftfork = x[i+1][0]-x[i][0]
            avg = (rightfork+leftfork)/2
            distfromorg.append([x[i],rightfork,leftfork,avg])
            
    return distfromorg
            
# Plot just predicted origins and termination points above and below an oemlim
def plotOEMForPredictedOriginsTerminations(oemfile,chromosome,oemlim,strain,col):
    
    orgsandterm = getListOfOriginsAndTerminationPoints(oemfile,chromosome)
    orgsandtermLim = [i for i in orgsandterm if (i[2]=='O' and i[1]>=oemlim) or (i[2]=='T' and i[1]<=-oemlim)]
    cood = [i[0] for i in orgsandtermLim]
    oems = [i[1] for i in orgsandtermLim]
    plt.plot(cood,oems,'.-'+col,label=strain)
    plt.axhline(xmax=yeastsize[chromosome],color='black')
    plt.show()
    
def ratioOrDiffOfOriginDistancefromTerm(oemfile,reptimefile,giveratio,option,chromosome):
    
    allpoints = getListOfOriginsAndTerminationPoints(oemfile,chromosome)
    
    splitTerm = [[allpoints[i-1],allpoints[i],allpoints[i+1]] for i in range(1,len(allpoints)-1) if allpoints[i][2]=='T']
    
    justorg = [[i[0],i[1]] for i in allpoints if i[2]=='O']
    
    orgwithreptime = rt.assignReplicationTimesForPredictedOrigins(reptimefile,justorg,chromosome)
        
    ratiosordiffs = []
    if option=='early-early':
        earlyorg = [[i[0],i[1],'O'] for i in orgwithreptime if i[2]<30]        
        for i in splitTerm:
            if i[0] in earlyorg and i[2] in earlyorg:
                rightfork = i[1][0]-i[0][0]
                leftfork = i[2][0]-i[1][0]
                if giveratio:
                    ratiosordiffs.append(float(rightfork)/leftfork)
                else:
                    ratiosordiffs.append(abs(rightfork-leftfork))
    elif option=='late-late':
        lateorg = [[i[0],i[1],'O'] for i in orgwithreptime if i[2]>=30]
        for i in splitTerm:
            if i[0] in lateorg and i[2] in lateorg:
                rightfork = i[1][0]-i[0][0]
                leftfork = i[2][0]-i[1][0]
                if giveratio:
                    ratiosordiffs.append(float(rightfork)/leftfork)
                else:
                    ratiosordiffs.append(abs(rightfork-leftfork))
    elif option=='early-late':
        earlyorg = [[i[0],i[1],'O'] for i in orgwithreptime if i[2]<30]
        lateorg = [[i[0],i[1],'O'] for i in orgwithreptime if i[2]>=30]
        for i in splitTerm:
            if (i[0] in earlyorg and i[2] in lateorg) or (i[2] in earlyorg and i[0] in lateorg) :
                rightfork = i[1][0]-i[0][0]
                leftfork = i[2][0]-i[1][0]
                if giveratio:
                    ratiosordiffs.append(float(rightfork)/leftfork)
                else:
                    ratiosordiffs.append(abs(rightfork-leftfork))
    elif option=='same time':
        for i in splitTerm:
            time1 = [j[2] for j in orgwithreptime if i[0][0]==j[0]]
            time2 = [j[2] for j in orgwithreptime if i[2][0]==j[0]]
            if time1 == time2:
                rightfork = i[1][0]-i[0][0]
                leftfork = i[2][0]-i[1][0]
                if giveratio:
                    ratiosordiffs.append(float(rightfork)/leftfork)
                else:
                    ratiosordiffs.append(abs(rightfork-leftfork))
    elif option=="all":
        for i in splitTerm:
            rightfork = i[1][0]-i[0][0]
            leftfork = i[2][0]-i[1][0]
            if giveratio:
                ratiosordiffs.append(float(rightfork)/leftfork)
            else:
                ratiosordiffs.append(abs(rightfork-leftfork))        
    return ratiosordiffs
    
def ratioOrDiffOfOriginDistFromTermForAllChromosomes(oemfile,reptimefile,getratio,option):
    allratios = []    
    for c in range(1,17):
        chromosome= str(c)
        print chromosome
        ratioperchrom = ratioOrDiffOfOriginDistancefromTerm(oemfile,reptimefile,getratio,option,chromosome)
        allratios.extend(ratioperchrom)
        
    return allratios
    
def graboriginsThatAreClose(oemfile1,oemfile2,chromosome,window,oemlimit):
    po1 = oc.getPredictedOriginsForOEMList(oemfile1,chromosome)
    po2 = oc.getPredictedOriginsForOEMList(oemfile2,chromosome)
    morg = oc.matchingOriginsBasedOnWindow([i[0] for i in po1 if i[1]>=oemlimit],[i[0] for i in po2 if i[1]>=oemlimit],window)
    o1 = [i[0] for i in morg]
    o2 = [i[1] for i in morg]
    return [o1,o2]

def grabWTOriginsCommonForAll(oemfileWT,oemfile1,oemfile2,oemfile3,chromosome,window,oemlimit):
    bothorigins1 = graboriginsThatAreClose(oemfileWT,oemfile1,chromosome,window,oemlimit)
    bothorigins2 = graboriginsThatAreClose(oemfileWT,oemfile2,chromosome,window,oemlimit)
    bothorigins3 = graboriginsThatAreClose(oemfileWT,oemfile3,chromosome,window,oemlimit)
    commonorigins = (set(bothorigins3[0])).intersection((set(bothorigins2[0])).intersection(set(bothorigins1[0])))
    
    return commonorigins
    
def grabWTOriginsCommonForAllForAllChromosomes(oemfileWT,oemfile1,oemfile2,oemfile3,window,oemlimit):
    commorgsallchrom = {}    
    for c in range(1,17):
        chromosome = str(c)
        commorgsallchrom[chromosome]=(list(grabWTOriginsCommonForAll(oemfileWT,oemfile1,oemfile2,oemfile3,chromosome,window,oemlimit)))
    return commorgsallchrom
    
def ratioOrDiffOfOriginDistancefromTermForCommonOrigins(oemfile,reptimefile,giveratio,option,chromosome,commonorigins,window,oemlimit):
    
    allpoints = [i for i in getListOfOriginsAndTerminationPoints(oemfile,chromosome) if (i[2]=='O' and i[1]>=oemlimit) or (i[2]=='T' and i[1]<=-oemlimit)]
    predictedorgs = [i[0] for i in allpoints if i[2]=='O']
    predictedorgNearCommonOrigins = [i[1] for i in oc.matchingOriginsBasedOnWindow(commonorigins,predictedorgs,window)]
    
    splitTerm = [[allpoints[i-1],allpoints[i],allpoints[i+1]] for i in range(1,len(allpoints)-1) if allpoints[i][2]=='T' and allpoints[i-1][0] in predictedorgNearCommonOrigins and allpoints[i+1][0] in predictedorgNearCommonOrigins]
    
    
    justorg = [[i[0],i[1]] for i in allpoints if i[2]=='O']
    
    orgwithreptime = rt.assignReplicationTimesForPredictedOrigins(reptimefile,justorg,chromosome)
        
    ratiosordiffs = []
    if option=='early-early':
        earlyorg = [[i[0],i[1],'O'] for i in orgwithreptime if i[2]<30]        
        for i in splitTerm:
            if i[0] in earlyorg and i[2] in earlyorg:
                if giveratio:
                    interorigindist = i[2][0]-i[0][0]
                    distoftermfromfirstorg = i[1][0]-i[0][0]
                    ratiosordiffs.append(float(distoftermfromfirstorg)/interorigindist)
                else:
                    rightfork = i[1][0]-i[0][0]
                    leftfork = i[2][0]-i[1][0]
                    ratiosordiffs.append(abs(rightfork-leftfork))
    elif option=='late-late':
        lateorg = [[i[0],i[1],'O'] for i in orgwithreptime if i[2]>=30]
        for i in splitTerm:
            if i[0] in lateorg and i[2] in lateorg:
                if giveratio:
                    interorigindist = i[2][0]-i[0][0]
                    distoftermfromfirstorg = i[1][0]-i[0][0]
                    ratiosordiffs.append(float(distoftermfromfirstorg)/interorigindist)
                else:
                    rightfork = i[1][0]-i[0][0]
                    leftfork = i[2][0]-i[1][0]
                    ratiosordiffs.append(abs(rightfork-leftfork))
    elif option=='early-late':
        earlyorg = [[i[0],i[1],'O'] for i in orgwithreptime if i[2]<30]
        lateorg = [[i[0],i[1],'O'] for i in orgwithreptime if i[2]>=30]
        for i in splitTerm:
            if (i[0] in earlyorg and i[2] in lateorg) or (i[2] in earlyorg and i[0] in lateorg) :
                if giveratio:
                    interorigindist = i[2][0]-i[0][0]
                    distoftermfromfirstorg = i[1][0]-i[0][0]
                    ratiosordiffs.append(float(distoftermfromfirstorg)/interorigindist)
                else:
                    rightfork = i[1][0]-i[0][0]
                    leftfork = i[2][0]-i[1][0]
                    ratiosordiffs.append(abs(rightfork-leftfork))
    elif option=='same time':
        for i in splitTerm:
            time1 = [j[2] for j in orgwithreptime if i[0][0]==j[0]]
            time2 = [j[2] for j in orgwithreptime if i[2][0]==j[0]]
            if time1 == time2:
                if giveratio:
                    interorigindist = i[2][0]-i[0][0]
                    distoftermfromfirstorg = i[1][0]-i[0][0]
                    ratiosordiffs.append(float(distoftermfromfirstorg)/interorigindist)
                else:
                    rightfork = i[1][0]-i[0][0]
                    leftfork = i[2][0]-i[1][0]
                    ratiosordiffs.append(abs(rightfork-leftfork))
    elif option=="all":
        for i in splitTerm:
            if giveratio:
                interorigindist = i[2][0]-i[0][0]
                distoftermfromfirstorg = i[1][0]-i[0][0]
                ratiosordiffs.append(float(distoftermfromfirstorg)/interorigindist)
            else:
                rightfork = i[1][0]-i[0][0]
                leftfork = i[2][0]-i[1][0]
                ratiosordiffs.append(abs(rightfork-leftfork))
                 
    return ratiosordiffs

def ratioOrDiffOfOriginDistFromTermForAllChromosomesForCommonOrigins(oemfile,reptimefile,getratio,option,window,commonorgs,oemlimit):
    allratios = []    
    for c in range(1,17):
        chromosome= str(c)
        ratioperchrom = ratioOrDiffOfOriginDistancefromTermForCommonOrigins(oemfile,reptimefile,getratio,option,chromosome,commonorgs[chromosome],window,oemlimit)
        if len(ratioperchrom)!=0:        
            allratios.extend(ratioperchrom)
        
    return allratios
    
def terminationsFoundNeartRNAGenes(tRNAgenesfile,oemfile,window):
    relevantTermpoints = []    
    for c in range(1,17):
        chromosome = str(c)
        print chromosome        
        terms = getPredictedTerminationPointsForOEMList(oemfile,chromosome)
        with open(tRNAgenesfile) as f:
            genesforchrom = [line.strip().split('\t') for line in f if line.strip().split('\t')[1]==chromosome]
            
        for i in terms:
            for j in genesforchrom:
                if (j[4]=='W' and i[0] >= int(j[2])-window and i[0]<=int(j[2])+window)  or (j[4]=='C' and i[0] >= int(j[3])-window and i[0]<=int(j[3])+window):
                    relevantTermpoints.append([chromosome,i[0],i[1]])
    return relevantTermpoints
    
def distanceTermPointsFromtRNAgenes(tRNAgenesfile,oemfile,chromosome):
    alldistforallterm = []
    terms = getPredictedTerminationPointsForOEMList(oemfile,chromosome)
    with open(tRNAgenesfile) as f:
            genesforchrom = [line.strip().split('\t') for line in f if line.strip().split('\t')[1]==chromosome]
    for i in terms:
        for j in range(len(genesforchrom    )):
            if j==0 and genesforchrom[j][4]=='W':            
                distforeachi = abs(int(genesforchrom[j][2])-i[0])
            elif j==0 and genesforchrom[j][4]=='C':
                distforeachi = abs(int(genesforchrom[j][3])-i[0])
            else:
                if genesforchrom[j][4]=='W' and abs(int(genesforchrom[j][2])-i[0]) <= distforeachi:
                    distforeachi = abs(int(genesforchrom[j][2])-i[0])
                elif genesforchrom[j][4]=='C' and abs(int(genesforchrom[j][3])-i[0]) <= distforeachi:
                    distforeachi = abs(int(genesforchrom[j][2])-i[0])
            print distforeachi
        alldistforallterm.append([i[0],i[1],distforeachi])
        
    return alldistforallterm
            
def main(): 
    
    commonorigins = grabWTOriginsCommonForAllForAllChromosomes('../cdc9deg/Raw/cdc9degRawOEM.txt','../rrm3d_pif1m2/Raw/rrm3d_pif1m2RawOEM.txt','../rrm3d/Raw/rrm3dRawOEM.txt','../pif1m2/Raw/pif1m2RawOEM.txt',1000,0.1)
    print('Plotting BoxPlots for two early replication origins')    
    ee_WT = ratioOrDiffOfOriginDistFromTermForAllChromosomesForCommonOrigins('../cdc9deg/Raw/cdc9degRawOEM.txt','DataFiles/ragutiming.txt',True,'early-early',1000,commonorigins,0.1)
    ee_rrm3d = ratioOrDiffOfOriginDistFromTermForAllChromosomesForCommonOrigins('../rrm3d/Raw/rrm3dRawOEM.txt','DataFiles/ragutiming.txt',True,'early-early',1000,commonorigins,0.1)
    ee_rrm3dpif1m2 = ratioOrDiffOfOriginDistFromTermForAllChromosomesForCommonOrigins('../rrm3d_pif1m2/Raw/rrm3d_pif1m2RawOEM.txt','DataFiles/ragutiming.txt',True,'early-early',1000,commonorigins,0.1)
    ee_pif1m2 = ratioOrDiffOfOriginDistFromTermForAllChromosomesForCommonOrigins('../pif1m2/Raw/pif1m2RawOEM.txt','DataFiles/ragutiming.txt',True,'early-early',1000,commonorigins,0.1)
    plt.figure(1)    
    plt.clf()
    bp = plt.boxplot([ee_WT,ee_rrm3d,ee_pif1m2,ee_rrm3dpif1m2])
    for line in bp['medians']:
        x1, y = line.get_xydata()[0]
        x2 = line.get_xydata()[1][0]
        x = (x1+x2)/2
        plt.text(x, y, '%.3f' % y,horizontalalignment='center')
    plt.xticks(range(1,5),('WT','rrm3d','pif1m2','rrm3d-pif1m2'))
    plt.title('Difference in Distance between origins and fork termination point - both early replication origins')

    print('Plotting BoxPlots for one early and one late replication origin') 
    el_WT = ratioOrDiffOfOriginDistFromTermForAllChromosomesForCommonOrigins('../cdc9deg/Raw/cdc9degRawOEM.txt','DataFiles/ragutiming.txt',True,'early-late',1000,commonorigins,0.1)
    el_rrm3d = ratioOrDiffOfOriginDistFromTermForAllChromosomesForCommonOrigins('../rrm3d/Raw/rrm3dRawOEM.txt','DataFiles/ragutiming.txt',True,'early-late',1000,commonorigins,0.1)
    el_rrm3dpif1m2 = ratioOrDiffOfOriginDistFromTermForAllChromosomesForCommonOrigins('../rrm3d_pif1m2/Raw/rrm3d_pif1m2RawOEM.txt','DataFiles/ragutiming.txt',True,'early-late',1000,commonorigins,0.1)
    el_pif1m2 = ratioOrDiffOfOriginDistFromTermForAllChromosomesForCommonOrigins('../pif1m2/Raw/pif1m2RawOEM.txt','DataFiles/ragutiming.txt',True,'early-late',1000,commonorigins,0.1)
    plt.figure(2)    
    plt.clf()
    bp = plt.boxplot([el_WT,el_rrm3d,el_pif1m2,el_rrm3dpif1m2])
    for line in bp['medians']:
        x1, y = line.get_xydata()[0]
        x2 = line.get_xydata()[1][0]
        x = (x1+x2)/2
        plt.text(x, y, '%.3f' % y,horizontalalignment='center')
    plt.xticks(range(1,5),('WT','rrm3d','pif1m2','rrm3d-pif1m2'))
    plt.title('Difference in Distance between origins and fork termination point - one early and one late replication origin')
    
    print('Plotting BoxPlots for two late replication origins')
    ll_WT = ratioOrDiffOfOriginDistFromTermForAllChromosomesForCommonOrigins('../cdc9deg/Raw/cdc9degRawOEM.txt','DataFiles/ragutiming.txt',True,'late-late',1000,commonorigins,0.1)
    ll_rrm3d = ratioOrDiffOfOriginDistFromTermForAllChromosomesForCommonOrigins('../rrm3d/Raw/rrm3dRawOEM.txt','DataFiles/ragutiming.txt',True,'late-late',1000,commonorigins,0.1)
    ll_rrm3dpif1m2 = ratioOrDiffOfOriginDistFromTermForAllChromosomesForCommonOrigins('../rrm3d_pif1m2/Raw/rrm3d_pif1m2RawOEM.txt','DataFiles/ragutiming.txt',True,'late-late',1000,commonorigins,0.1)
    ll_pif1m2 = ratioOrDiffOfOriginDistFromTermForAllChromosomesForCommonOrigins('../pif1m2/Raw/pif1m2RawOEM.txt','DataFiles/ragutiming.txt',True,'late-late',1000,commonorigins,0.1)
    plt.figure(3)    
    plt.clf()
    bp = plt.boxplot([ll_WT,ll_rrm3d,ll_pif1m2,ll_rrm3dpif1m2])
    for line in bp['medians']:
        x1, y = line.get_xydata()[0]
        x2 = line.get_xydata()[1][0]
        x = (x1+x2)/2
        plt.text(x, y, '%.3f' % y,horizontalalignment='center')
    plt.xticks(range(1,5),('WT','rrm3d','pif1m2','rrm3d-pif1m2'))
    plt.title('Difference in Distance between origins and fork termination point - both late replication origins')
    
    print('Plotting BoxPlots for replication origins firing at the same time')
    same_WT = ratioOrDiffOfOriginDistFromTermForAllChromosomesForCommonOrigins('../cdc9deg/Raw/cdc9degRawOEM.txt','DataFiles/ragutiming.txt',True,'same time',1000,commonorigins,0.1)
    same_rrm3d = ratioOrDiffOfOriginDistFromTermForAllChromosomesForCommonOrigins('../rrm3d/Raw/rrm3dRawOEM.txt','DataFiles/ragutiming.txt',True,'same time',1000,commonorigins,0.1)
    same_rrm3dpif1m2 = ratioOrDiffOfOriginDistFromTermForAllChromosomesForCommonOrigins('../rrm3d_pif1m2/Raw/rrm3d_pif1m2RawOEM.txt','DataFiles/ragutiming.txt',True,'same time',1000,commonorigins,0.1)
    same_pif1m2 = ratioOrDiffOfOriginDistFromTermForAllChromosomesForCommonOrigins('../pif1m2/Raw/pif1m2RawOEM.txt','DataFiles/ragutiming.txt',True,'same time',1000,commonorigins,0.1)
    plt.figure(4)    
    plt.clf()
    bp = plt.boxplot([same_WT,same_rrm3d,same_pif1m2,same_rrm3dpif1m2])
    for line in bp['medians']:
        x1, y = line.get_xydata()[0]
        x2 = line.get_xydata()[1][0]
        x = (x1+x2)/2
        plt.text(x, y, '%.3f' % y,horizontalalignment='center')
    plt.xticks(range(1,5),('WT','rrm3d','pif1m2','rrm3d-pif1m2'))
    plt.title('Difference in Distance between origins and fork termination point - replication origins firing at same time')
    
    print('Plotting BoxPlots for all replication origins')
    all_WT = ratioOrDiffOfOriginDistFromTermForAllChromosomesForCommonOrigins('../cdc9deg/Raw/cdc9degRawOEM.txt','DataFiles/ragutiming.txt',True,'all',1000,commonorigins,0.1)
    all_rrm3d = ratioOrDiffOfOriginDistFromTermForAllChromosomesForCommonOrigins('../rrm3d/Raw/rrm3dRawOEM.txt','DataFiles/ragutiming.txt',True,'all',1000,commonorigins,0.1)
    all_rrm3dpif1m2 = ratioOrDiffOfOriginDistFromTermForAllChromosomesForCommonOrigins('../rrm3d_pif1m2/Raw/rrm3d_pif1m2RawOEM.txt','DataFiles/ragutiming.txt',True,'all',1000,commonorigins,0.1)
    all_pif1m2 = ratioOrDiffOfOriginDistFromTermForAllChromosomesForCommonOrigins('../pif1m2/Raw/pif1m2RawOEM.txt','DataFiles/ragutiming.txt',True,'all',1000,commonorigins,0.1)
    plt.figure(5)    
    plt.clf()
    bp = plt.boxplot([all_WT,all_rrm3d,all_pif1m2,all_rrm3dpif1m2])
    for line in bp['medians']:
        x1, y = line.get_xydata()[0]
        x2 = line.get_xydata()[1][0]
        x = (x1+x2)/2
        plt.text(x, y, '%.3f' % y,horizontalalignment='center')
    plt.xticks(range(1,5),('WT','rrm3d','pif1m2','rrm3d-pif1m2'))
    plt.title('Difference in Distance between origins and fork termination point - all replication origins')