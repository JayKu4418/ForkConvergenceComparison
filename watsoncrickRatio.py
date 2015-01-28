# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 15:14:50 2014

@author: jayashreekumar
"""
import numpy as np
import statsmodels.api as sm
import pyfftw
import matplotlib.pyplot as plt
import ManipulateSeqData.originscomparison as oc
import ManipulateSeqData.oem as oem
#from operator import itemgetter
from scipy.stats import pearsonr
from scipy.spatial import distance
import ManipulateSeqData.findrepeats as fr
import wcratiospecificregion as wcs
import random

#import matplotlib.pyplot as plt
#import Bio.Statistics as bs
#from datetime import datetime

# This is a dictionary of yeast chromosome sizes
yeastsize = {'1':230218,'2':813184,'3':316620,'4':1531933,'5':576874,'6':270161,'7':1090940,'8':562643,'9':439888,'10':745751,'11':666816,'12':1078177,'13':924431,'14':784333,'15':1091291,'16':948066,'Total':12157105}

# dictionary that converts values from the form 'chrmIII' to '3'
chrmconvert = {'chrI':'1','chrII':'2','chrIII':'3','chrIV':'4','chrV':'5','chrVI':'6','chrVII':'7','chrVIII':'8','chrIX':'9','chrX':'10','chrXI':'11','chrXII':'12','chrXIII':'13','chrXIV':'14','chrXV':'15','chrXVI':'16'}

# This function calcuates the 2log ratio of watson to crick hits
def watsoncricklogratio(sgrfileW,sgrfileC,chrom):
    with open(sgrfileW) as f:
        w = [int(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chrom]
        
    with open(sgrfileC) as f:
        c = [int(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chrom]
    smalle = 10**(-5)
    # Add 1 to  is to get rid of zero values an
    x = [(float(w[i]+smalle)/(c[i]+smalle)) for i in range(0,len(w))]
       
    y =  [np.log2(i) for i in x]
#def smoothWatsonCricklogratio()
    #print(str(datetime.now()))
    #yl = sm.nonparametric.lowess(y,range(1,(len(w)+1)))
    #print(str(datetime.now()))
    return y

# This function writes the watson and crick log2 ratio unsmoothed as file for for all chromosomes
def watsoncricklogratioWriteFile(sgrfileW,sgrfileC,writefile):
    with open(writefile,'w') as fw:
        for c in range(1,17):
            chrom = str(c)
            print(chrom)
            vals = watsoncricklogratio(sgrfileW,sgrfileC,chrom)
            for i in range(1,len(vals)+1):
                fw.write(chrom + '\t' + str(i) + '\t' + str(vals[i-1]) + '\n')
###############################SMOOTHING FUNCTIONS##################################
# this function calculates a moving average for vals using a window specified
def movingavg(vals,window):
    movavgforvals = []
    for i in range(len(vals)):
        if (i < window) or (i > len(vals)-window):
            movavgforvals.append(vals[i])
        else:
            movavgforvals.append(np.mean(vals[i-window:i+window]))
    
    return movavgforvals
    
# This function runs lowess function for vals using frac f  and delta d specified, used f = 0.01 and d= 0.001
def lowessforWC(vals,f,d):
    valslow = sm.nonparametric.lowess(vals,range(1,len(vals)+1),frac=f,delta=d*len(vals))
    return [i[1] for i in valslow]
    
# This function writes the watson and crick log2 ratio unsmoothed as file for for all chromosomes
def watsoncricklogratioLowessSmoothedWriteFile(sgrfileW,sgrfileC,f,d,writefile):
    with open(writefile,'w') as fw:
        for c in range(1,17):
            chrom = str(c)
            print(chrom)
            vals = watsoncricklogratio(sgrfileW,sgrfileC,chrom)
            vals_smooth = lowessforWC(vals,f,d)
            for i in range(1,len(vals_smooth)+1):
                fw.write(chrom + '\t' + str(i) + '\t' + str(vals_smooth[i-1]) + '\n')

# Use Fourier transform to smooth vals
def ftSmoothforWC(vals,nfreq):
    vals = np.array(vals)
    f_vals = pyfftw.interfaces.numpy_fft.rfft(vals)
    f_vals[nfreq:] = 0
    # Inverse fourier transform to get back smoothed signal
    fvals_smooth = pyfftw.interfaces.numpy_fft.irfft(f_vals)
    
    return fvals_smooth

###################################################################################

##########################EFFICIENCY##########################################

"""
def calculateDer(vals,window):
    gradientOfvals = []
    for i in range(len(vals)):
        if i > window/2 and i < (len(vals)-window/2):
            grad = -1000*((vals[i+(window/2)]-vals[i-(window/2)])/window)
            gradientOfvals.append(grad)
        else:
            gradientOfvals.append(0)
    return gradientOfvals
"""
    
from scipy import interpolate

from pylab import plot
from numpy import arange

def draw_tangent(x,y,a):
 # interpolate the data with a spline
    spl = interpolate.splrep(x,y)
    small_t = arange(a-1000,a+1000)
    fa = interpolate.splev(a,spl,der=0)     # f(a)
    fprime = interpolate.splev(a,spl,der=1) # f'(a)
    tan = fa+fprime*(small_t-a) # tangent
    #plot(x,y,alpha=0.5)
    plot(a,fa,'om',small_t,tan,'--r')
"""
def calculateDer(vals,d=1):
    spl = interpolate.splrep(range(len(vals)),vals)
    fprime = [1000*(interpolate.splev(i,spl,der=d)) for i in range(len(vals))]
    return fprime
"""

def calculateDer(vals,d=1):
    if d==1:
        return np.gradient(vals,1)
    else:
        firstder = np.gradient(vals,1)
        return np.gradient(firstder,1)
        
def writeDerFile(wcratiofile,d,writefile):
    with open(writefile,'w') as fw:
        for c in range(1,17):
            chromosome = str(c)
            print chromosome
            with open(wcratiofile) as f:
                vals = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
            
            dercalculated = calculateDer(vals,d)
            for i in range(1,len(vals)+1):
                fw.write(chromosome + '\t' + str(i) + '\t' + str(dercalculated[i-1]) + '\n')
    
###################################################################################  
    
#############################ORIGINS-PREDICTION#########################################
# Function to identify origins using log2 ratio of Watson Crick. Identify
# where 0's are and then check first derivative at 0 to see if its negatve
def identifyOrigins(svals):
    origins = []
    for i in range(30000,len(svals)-30000):
        if svals[i]>0 and svals[i+1]<0:
            fd = (svals[i+50] - svals[i-50])/100
            if fd < 0:             
                origins.append([i,fd])
    return origins
# Function that compares predicted Origins with OriDB origins
def comparePredictedOrigins(svals,chromosome,origintype,comparewindow):
    comparableOrigins = [[i[0],i[1]] for i in oem.readOriginsfile('S-cerevisiae_origins.txt',chromosome,origintype)]
    predictedOrigins = identifyOrigins(svals)
    
    return oc.predictedComparedToOriginsList(predictedOrigins,comparableOrigins,comparewindow)
# This function gets number of predicted origins in an oem file that are within originwindow bases of confirmed,likely or
# dubious origins. Origins not fitting any of these categories are also recorded
def getnumComparedOrigins(wcratiofile,originwindow):
    # Go through every chromosome in yeast

    # Intialize variables to store total num of predicted matching origins to be zero
    totalConfirmed = 0
    totalLikely = 0
    totalDubious = 0
    totalNoClass = 0

    for i in range(1,17):
        # Need to convert number into string to access it in all lists
        chromosome = str(i)
        print(chromosome)
        with open(wcratiofile) as f:
            wcratiovals = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
        # Return predicted origins
        po = identifyOrigins(wcratiovals)

        # This is a list of confirmed origins from OriDB
        confirmedOrigins = [[i[0],i[1]] for i in oem.readOriginsfile('S-cerevisiae_origins.txt',chromosome,'Confirmed')]
        # This is a list of likely origins from OriDB
        likelyOrigins = [[i[0],i[1]] for i in oem.readOriginsfile('S-cerevisiae_origins.txt',chromosome,'Likely')]
        # This is a list of dubious origins from OriDB
        dubiousOrigins = [[i[0],i[1]] for i in oem.readOriginsfile('S-cerevisiae_origins.txt',chromosome,'Dubious')]

        # Number of origins that are confirmed,likely,dubious and None of those among predicted origins for chromosome
        originsConfirmedInAndOut = oc.predictedComparedToOriginsList(po,confirmedOrigins,originwindow)
        originsLikelyInAndOut = oc.predictedComparedToOriginsList(originsConfirmedInAndOut['Not'],likelyOrigins,originwindow)
        originsDubiousInAndOut = oc.predictedComparedToOriginsList(originsLikelyInAndOut['Not'],dubiousOrigins,originwindow)
        originsNone = originsDubiousInAndOut['Not']

        # Add it to the total number of confirmed, likely, dubious or none of those among predicted origins
        totalConfirmed = totalConfirmed + len(originsConfirmedInAndOut['Within'])
        totalLikely = totalLikely + len(originsLikelyInAndOut['Within'])
        totalDubious = totalDubious + len(originsDubiousInAndOut['Within'])
        totalNoClass = totalNoClass + len(originsNone)



    return {'C':totalConfirmed,'L':totalLikely,'D':totalDubious,'N':totalNoClass}

# This function grabs pairs of origins in WT and MUT that are close to each other within a certain window
def graboriginsThatAreClose(valsWT,valsMUT,window):
    poWT = identifyOrigins(valsWT)
    poMUT = identifyOrigins(valsMUT)
    morg = oc.matchingOriginsBasedOnWindow([i[0] for i in poWT],[i[0] for i in poMUT],window)
    o1 = [i[0] for i in morg]
    o2 = [i[1] for i in morg]
    return [o1,o2]

# This function groups Origins together so that matching pairs of WT and MUT origins that are close to each other are grouped tgether
def groupOriginsTogetherForWTandMUT(poWT,poMUT):
    
    if (len(poWT)==len(poMUT)):
        grpOrgs = [[[poWT[i],poWT[i+1]],[poMUT[i],poMUT[i+1]]] for i in range(len(poWT)-1)]
    else:
        grpOrgs = []
        
    return grpOrgs

# This function looks at the watson crick ratio at each known origin and writes it
# into a file about 10000 bases away
def originsWCRatioWriteFile(originsfile,wcratiofile,window,writefile):
    with open(originsfile) as f:
        origins = [line.strip().split('\t') for line in f]
        
    originsToRet = [[i[0],int(i[1]),int(i[2])] for i in origins if i[5]=='Confirmed' and int(i[1])-window > 0 and int(i[2])+window < yeastsize[i[0]]]
    
    with open(writefile,'w') as fw:
        for c in range(1,17):
            chromosome = str(c)
            with open(wcratiofile) as f:
                wcratiosforchrom = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
            print(len(wcratiosforchrom))
            originsforchrom = [i for i in originsToRet if i[0]==chromosome]
            for t in originsforchrom:
                start = t[1]-window
                end = t[1]+window
                fw.write(t[0]+'\t'+str(t[1]) + '\t' + str(t[2]) + '\t')
                listToWrite = [str(i) for i in wcratiosforchrom[start:end]]
                fw.write('\t'.join(listToWrite))
                fw.write('\n')

# This function calculates the average watson crick ratio for every base pair 
# away +window and -window away from the origin 
def originsWCRatioMean(originsWCRatiofile,originsToExclude,window):
    with open(originsWCRatiofile) as f:
        wcratios = [line.strip().split('\t') for line in f]
    print(len(wcratios))
    wcratiosToInclude = [wcratios[i] for i in range(len(wcratios)) if i+1 not in originsToExclude]
    print(len(wcratiosToInclude))
    avgratios = []
    for r in range(3,2*window+3):
        ratiosforr = [float(i[r]) for i in wcratiosToInclude]
        avgratios.append(np.mean(ratiosforr))
    return avgratios

###################################################################################

##########################TERMINATIONS##########################################
# This function gets positions of fork mergers between origins using watson crick vals and predicted origins
# This is for a single origin pair
def forkmergersBtwnOriginPairsForChromosome(vals,originpair,withfd,getnum):
    valsplusInd = [[i+1,vals[i]] for i in range(len(vals))]
    valsbtwnOrgs = valsplusInd[originpair[0]:originpair[1]]
    terminationsPerOriginPair = []
    for j in range(len(valsbtwnOrgs)-1):
        if valsbtwnOrgs[j][1]<0 and valsbtwnOrgs[j+1][1]>0:
            fd = (valsbtwnOrgs[j+50][1] - valsbtwnOrgs[j-50][1])/100
            if fd > 0:
                if withfd:
                    terminationsPerOriginPair.append([valsbtwnOrgs[j][0],fd])
                else:
                    terminationsPerOriginPair.append(valsbtwnOrgs[j][0])
    if getnum:
        return len(terminationsPerOriginPair)
    else:
        return terminationsPerOriginPair

# This function gets positions of fork mergers between origins using watson crick vals and predicted origins
# This is for a single chromosome for all origin pairs
def forkmergersbetweenAllOriginsForChromosome(vals,po,withfd,getnum):    
    allterminations = []
    predictedorigins= sorted(po)
    for i in range(0,len(predictedorigins)-1):
        terminationsPerOriginPair = forkmergersBtwnOriginPairsForChromosome(vals,[predictedorigins[i],predictedorigins[i+1]],withfd,getnum)
        allterminations.append([predictedorigins[i],terminationsPerOriginPair,predictedorigins[i+1]])

    return allterminations
    

###############################################################################
    
##########################FEATURES#############################################
# This function returns a ratio of distance of termination from the first origin and the total distance between the origins
def ratioForkMergerForPairsOfOrigins(originsplusterms):
    if len(originsplusterms[1])==1:
        fract = (originsplusterms[1][0] - originsplusterms[0])/float(originsplusterms[2]-originsplusterms[0])
    else:
        termavg = np.mean(originsplusterms[1])
        fract = (termavg - originsplusterms[0])/float(originsplusterms[2]-originsplusterms[0])
        
    return [fract,originsplusterms[0],originsplusterms[1],originsplusterms[2]]
    
# This function returns ratios of distance of termination from the first origin and the total distance between the origins, for all
# origins in a chromosome
def ratioForkMergerForChromosome(originsplusterms):
    fracsPlusOrgsAndterm = []
    for i in originsplusterms:
        fracsPlusOrgsAndterm.append(ratioForkMergerForPairsOfOrigins(i))
        
    return fracsPlusOrgsAndterm

# This function returns ratios of distance of termination from the first origin and the total distance between the origins, for all
# origins in all chromosome       
def fractionOfForkMergerAllChromosomes(wcratiofile,window,justfracs):
    if justfracs:
        allfracs = []
    else:
        allfracs = {}
    
    for c in range(1,17):
        chromosome = str(c)
        print chromosome
        with open(wcratiofile) as f:
            valsforchrom = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
        
        comparepo = comparePredictedOrigins(valsforchrom,chromosome,'Confirmed',window)
        
        poWithin = [i[0] for i in comparepo['Within']]
        
        poplusterms = forkmergersbetweenAllOriginsForChromosome(valsforchrom,poWithin,False)
        
        fracsForMerger = ratioForkMergerForChromosome(poplusterms)
        if justfracs:
            allfracs.extend([i[0] for i in fracsForMerger])
        else:
            allfracs[chromosome] = fracsForMerger

    return allfracs


# This function calculates the correlation coefficient between the troughs of WT and MUT
def getCorrCoeffBetweenPairsOfCoordsForChromosome(valsforWT,coordWT,valsforMUT,coordMUT):
     valsbetweenCoordWT = valsforWT[coordWT[0]:coordWT[1]]
     valsbetweenCoordMUT = valsforMUT[coordMUT[0]:coordMUT[1]]
     minnumvals = min(len(valsbetweenCoordWT),len(valsbetweenCoordMUT))
     p = pearsonr(valsbetweenCoordWT[0:minnumvals],valsbetweenCoordMUT[0:minnumvals])
     
     return p

# This function calculates the euclidean distance between the troughs of WT and MUT
def getEuclideanDistanceBetweenPairsOfCoordsForChromosome(valsforWT,coordWT,valsforMUT,coordMUT):
    maxStartCoord = max(coordWT[0],coordMUT[0])
    minEndCoord = min(coordWT[1],coordMUT[1])
    valsbetweenCoordWT = valsforWT[maxStartCoord:minEndCoord]
    valsbetweenCoordMUT = valsforMUT[maxStartCoord:minEndCoord]
    
    dst = distance.euclidean(valsbetweenCoordWT,valsbetweenCoordMUT)
    
    return dst
    
# This function calculates the first derivative beween the troughs of WT and MUT and spits out the corr coeff and euclidean dist of vector
def getCorrCoeffAndEucDistBetweenFirstDerOfPairsOfCoordsForChromosome(valsforWT,coordWT,valsforMUT,coordMUT):
    maxStartCoord = max(coordWT[0],coordMUT[0])
    minEndCoord = min(coordWT[1],coordMUT[1])
    valsbetweenCoordWT = valsforWT[maxStartCoord:minEndCoord]
    valsbetweenCoordMUT = valsforMUT[maxStartCoord:minEndCoord]
    firstderWT = calculateDer(valsbetweenCoordWT,1)
    firstderMUT = calculateDer(valsbetweenCoordMUT,1)
    
    p = pearsonr(firstderWT,firstderMUT)
    
    dist = distance.euclidean(firstderWT,firstderMUT)

    return [p,dist]

# This function calculates the second derivative beween the troughs of WT and MUT and spits out the corr coeff and euclidean dist of vector
def getCorrCoeffAndEucDistBetweenSecDerOfPairsOfCoordsForChromosome(valsforWT,coordWT,valsforMUT,coordMUT):
    maxStartCoord = max(coordWT[0],coordMUT[0])
    minEndCoord = min(coordWT[1],coordMUT[1])
    valsbetweenCoordWT = valsforWT[maxStartCoord:minEndCoord]
    valsbetweenCoordMUT = valsforMUT[maxStartCoord:minEndCoord]
    #firstderWT = calculateDer(valsbetweenOrgWT,window)
    #firstderMUT = calculateDer(valsbetweenOrgMUT,window)
    secderWT = calculateDer(valsbetweenCoordWT,2)
    secderMUT = calculateDer(valsbetweenCoordMUT,2)
    
    p = pearsonr(secderWT,secderMUT)
    
    dist = distance.euclidean(secderWT,secderMUT)

    return [p,dist]    
    
def featuresOfOriginPairsForChromosome(valsforWT,valsforMUT,window):
    closeorigins = graboriginsThatAreClose(valsforWT,valsforMUT,window)
    groupedorigins = groupOriginsTogetherForWTandMUT(closeorigins[0],closeorigins[1])
    
    allfeatures = []
    
    for i in groupedorigins:
        
        featureperpair = {}
        corr = getCorrCoeffBetweenPairsOfCoordsForChromosome(valsforWT,i[0],valsforMUT,i[1])
        dist = getEuclideanDistanceBetweenPairsOfCoordsForChromosome(valsforWT,i[0],valsforMUT,i[1])
        termsforWT = forkmergersBtwnOriginPairsForChromosome(valsforWT,i[0],False,False)
        numtersforWT = len(termsforWT)
        fractermWT = ratioForkMergerForPairsOfOrigins([i[0][0],termsforWT,i[0][1]])[0]
        
        termsforMUT = forkmergersBtwnOriginPairsForChromosome(valsforMUT,i[1],False,False)
        numtersforMUT = len(termsforMUT)
        fractermMUT = ratioForkMergerForPairsOfOrigins([i[1][0],termsforMUT,i[1][1]])[0]
        
        firstder = getCorrCoeffAndEucDistBetweenFirstDerOfPairsOfCoordsForChromosome(valsforWT,i[0],valsforMUT,i[1])
        
    
        secder = getCorrCoeffAndEucDistBetweenSecDerOfPairsOfCoordsForChromosome(valsforWT,i[0],valsforMUT,i[1])
        
        featureperpair['Corr'] = corr
        featureperpair['Dist'] = dist
        featureperpair['NumberOfCrosses'] = [numtersforWT,numtersforMUT]
        featureperpair['TermOccurs'] = [fractermWT,fractermMUT]
        featureperpair['FirstDer'] = [firstder[0][0],firstder[1]]
        featureperpair['SecDer'] = [secder[0][0],secder[1]]

        allfeatures.append([i,featureperpair])

    return allfeatures
    
def writefeaturesOfOriginPairsForAllChromosomes(wcratiofileWT,wcratiofileMUT,writefile,window=2500):
    # Go through all chromosomes
    with open(writefile,'w') as fw:
        fw.write('Chr' + '\t' + 'WT O1' + '\t' + 'WT O2' + '\t' + 'MUT O1' + '\t' + 'MUT O2' + '\t' + 'WT Term Ratio' + '\t' + 'MUT Term Ratio' + '\t' + 'Corr Coeff' + '\t' + 'Euclidean Dist' + '\t' + 'FirstDer Corr' + '\t' + 'FirstDer Dist' + '\t' + 'SecDer Corr' + '\t' + 'SecDer Dist' + '\n')
        for c in range(1,17):
            chromosome = str(c)
            print chromosome
            with open(wcratiofileWT) as f:
                valsforWT = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
            with open(wcratiofileMUT) as f:
                valsforMUT = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
            
            feats = featuresOfOriginPairsForChromosome(valsforWT,valsforMUT,window)
            
            for i in feats:
                wto1 = str(i[0][0][0])
                wto2 = str(i[0][0][1])
                muto1 = str(i[0][1][0])
                muto2 = str(i[0][1][1])
                wttermratio = str(i[1]['TermOccurs'][0])
                muttermratio = str(i[1]['TermOccurs'][1])
                corr = str(i[1]['Corr'][0])
                dist = str(i[1]['Dist'])
                firstdercorr = str(i[1]['FirstDer'][0])
                firstderdist = str(i[1]['FirstDer'][1])
                secdercorr = str(i[1]['SecDer'][0])
                secderdist = str(i[1]['SecDer'][1])
                fw.write(chromosome + '\t' + wto1 + '\t' + wto2 + '\t' + muto1 + '\t' + muto2 + '\t' + wttermratio + '\t' + muttermratio + '\t' + corr + '\t' + dist + '\t' + firstdercorr + '\t' + firstderdist + '\t' + secdercorr + '\t' + secderdist + '\n')
        
# thos function gets rid of certain regions that overlap with the transposon regions where you have sequence gaps
# Remove troughs found within Ty1 elements
def removeRegionsFoundWithinTransposableElements(labelledtroughsfile,transfile,writefile):
    with open(labelledtroughsfile) as f:
        troughs = [line.strip().split('\t') for line in f][1:]
    with open(writefile,'w') as fw:
        fw.write('Chr' + '\t' + 'WT O1' + '\t' + 'WT O2' + '\t' + 'MUT O1' + '\t' + 'MUT O2' + '\t' + 'WT Term Ratio' + '\t' + 'MUT Term Ratio' + '\t' + 'Corr Coeff' + '\t' + 'Euclidean Dist' + '\t' + 'FirstDer Corr' + '\t' + 'FirstDer Dist' + '\t' + 'SecDer Corr' + '\t' + 'SecDer Dist' + '\n')
        for i in troughs:
            chrm = i[0]
            firsttroughStart = int(i[1])
            secondtroughStart = int(i[3])
            firsttroughEnd = int(i[2])
            secondtroughEnd = int(i[4])
            with open(transfile) as f:
                transValid = [line.strip().split(' ') for line in f if chrmconvert[line.strip().split(' ')[1]]==chrm]
            #transFoundInTrough = [j for j in transValid if (j[4]=='W' and ((int(j[2]) >= firsttroughStart and int(j[2]) <= firsttroughEnd) or (int(j[2]) >= secondtroughStart and int(j[2]) <= secondtroughEnd) )) or (j[4]=='C' and ((int(j[3]) >= firsttroughStart and int(j[3]) <= firsttroughEnd) or (int(j[3]) >= secondtroughStart and int(j[3]) <= secondtroughEnd)))]
            transFoundInTrough = [j for j in transValid if ((int(j[2]) >= firsttroughStart and int(j[2]) <= firsttroughEnd) or (int(j[2]) >= secondtroughStart and int(j[2]) <= secondtroughEnd) )]
            if len(transFoundInTrough)==0:
                fw.write(i[0] + '\t' + i[1] + '\t' + i[2] + '\t' + i[3] + '\t' + i[4] + '\t' + i[5] + '\t' + i[6] + '\t' + i[7] + '\t' + i[8] + '\t' + i[9] + '\t' + i[10] + '\t' + i[11] + '\t' + i[12] + '\n')
               
###############################################################################
##########################NORMALIZE WATSON CRICK VALS##########################
"""
# This function grabs the values between the originpair and normalizes them by
# subtracting the mean of the vector and dividing by the standard deviation
def normalizeValsBtwnOriginPairs(vals,originpair):
    # Values to Normalize grab btwn the origins    
    valsToNormalize = vals[originpair[0]:originpair[1]]
    muVals = np.mean(valsToNormalize)
    stdVals = np.std(valsToNormalize)
    normalizedVals = (((np.array(valsToNormalize))-muVals)/stdVals)
    
    return list(normalizedVals)
"""    
###############################################################################
##########################IDENTIFY DIFF TROUGHS################################
# This function identifies troughs in which the number of crossovers between predicted origins are different
def differentiateTroughsUsingCrossoverNums(valsforWT,valsforMUT,window):
    feats = featuresOfOriginPairsForChromosome(valsforWT,valsforMUT,window)
    
    numorigins = [i for i in feats if i[1]['NumberOfCrosses'][0]!=i[1]['NumberOfCrosses'][1]]


    return numorigins

# This function identifies troughs in which corr coeff is less than a certain lim
def differentiateTroughsUsingCorrCoeff(valsforWT,valsforMUT,window,corrcoefflim):
    feats = featuresOfOriginPairsForChromosome(valsforWT,valsforMUT,window)
    
    numorigins = [i for i in feats if i[1]['Corr'][0]<corrcoefflim]


    return numorigins
    
# This function identifies troughs in which midpt difference is greater than a certain limit
def differentiateTroughsUsingMidPts(valsforWT,valsforMUT,window,midptlim):
    feats = featuresOfOriginPairsForChromosome(valsforWT,valsforMUT,window)
    
    #numorigins = [i for i in feats if abs(i[1]['TermOccurs'][0]-i[1]['TermOccurs'][1])>midptlim]
    orgs = [i for i in feats if abs(i[1]['TermOccurs'][0]-i[1]['TermOccurs'][1])>midptlim]

    #return numorigins
    return orgs


"""
# This function grabs the overlaps between the differentiated troughs using features
def overlappingDiffTroughs(valsforWT,valsforMUT,window,corrcoefflim,midptlim):
    featsCrossOver = differentiateTroughsUsingCrossoverNums(valsforWT,valsforMUT,window)
    featsCorrCoeff = differentiateTroughsUsingCorrCoeff(valsforWT,valsforMUT,window,corrcoefflim)
    featsMidpoint = differentiateTroughsUsingMidPts(valsforWT,valsforMUT,window,midptlim)
    
    orgsWTCrossOver = [i[0][0][0] for i in featsCrossOver]
    orgsWTCorrCoeff = [i[0][0][0] for i in featsCorrCoeff]
    orgsWTMidpoint = [i[0][0][0] for i in featsMidpoint]
    
    orgsMUTCrossOver = [i[0][1][0] for i in featsCrossOver]
    orgsMUTCorrCoeff = [i[0][1][0] for i in featsCorrCoeff]
    orgsMUTMidpoint = [i[0][1][0] for i in featsMidpoint]
    
    orgsInCommonWT = [i for i in orgsWTCrossOver if i in [i for i in orgsWTCorrCoeff if i in orgsWTMidpoint]]
    
    orgsInCommonMUT = [i for i in orgsMUTCrossOver if i in [i for i in orgsMUTCorrCoeff if i in orgsMUTMidpoint]]
"""    
#############################GENOMIC ELEMENTS OF INTEREST #########################################
# This function grabs the tRNA genes found within a window around a region. 
# The region contains the chromosome and the start and end coordinate 
def tRNAGenesFoundInRegion(trnagenesfile,region,window):
    trnas = wcs.extractStartEndChromStrandFortRNAs(trnagenesfile)
    
    chromosome = region[0]
    start = region[1]
    end = region[2]
    
    validtrnasForRegion = [i for i in trnas if i[0]==chromosome and ( (int(i[1]) <= end+window and int(i[1]) >= start-window) or (int(i[2]) <= end+window and int(i[2]) >= start-window) ) ]
    
    return validtrnasForRegion
    
# This function counts number of tRNA genes in a bunch of regions
def tRNAGenesFoundInMultipleRegions(trnagenesfile,regionfile,window):
    with open(regionfile) as f:
        regsBothWTandMUT = [line.strip().split('\t') for line in f][1:]
        
    tRNAsFound = []
    
    for i in regsBothWTandMUT:
        reg = [i[0],min(int(i[1]),int(i[3])),max(int(i[2]),int(i[4]))]
        tRNASInReg = tRNAGenesFoundInRegion(trnagenesfile,reg,window)        
        tRNAsFound.append([reg,tRNASInReg])
        
    return tRNAsFound

# This function looks at the watson crick ratio at each tRNA gene and writes it
# into a file about 10000 bases away
def tRNAGenesWCRatioWriteFile(trnagenesfile,wcratiofile,window,writefile):
    trnas = wcs.extractStartEndChromStrandFortRNAs(trnagenesfile)
    with open(writefile,'w') as fw:
        for c in range(1,17):
            chromosome = str(c)
            with open(wcratiofile) as f:
                wcratiosforchrom = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
            print(len(wcratiosforchrom))
            trnasforchrom = [i for i in trnas if i[0]==chromosome]
            for t in trnasforchrom:
                start = t[1]-window
                end = t[1]+window
                fw.write(t[0]+'\t'+str(t[1]) + '\t' + str(t[2]) + '\t' + t[3] + '\t')
                listToWrite = [str(i) for i in wcratiosforchrom[start:end]]
                fw.write('\t'.join(listToWrite))
                fw.write('\n')
                
# This function calculates the average watson crick ratio for every base pair 
# away +500 and =500 away from the tRNA gene 
def tRNAGenesWCRatioMeanForStrand(trnageneWCRatiofile,strand,tRNAsToExclude,window):
    with open(trnageneWCRatiofile) as f:
        wcratios = [line.strip().split('\t') for line in f if line.strip().split('\t')[3]==strand]
    print(len(wcratios))
    wcratiosToInclude = [wcratios[i] for i in range(len(wcratios)) if i+1 not in tRNAsToExclude]
    print(len(wcratiosToInclude))
    avgratios = []
    for r in range(4,2*window+4):
        ratiosforr = [float(i[r]) for i in wcratiosToInclude]
        avgratios.append(np.mean(ratiosforr))
    return avgratios
    
# This function calculates the average watson crick ratio for every base pair 
# away +500 and =500 away from the tRNA gene 
def tRNAGenesDiffTwoStrainsWCRatioMeanForStrand(trnageneWCRatiofile1,trnageneWCRatiofile2,strand,tRNAsToExclude,window):
    with open(trnageneWCRatiofile1) as f:
        wcratios1 = [line.strip().split('\t') for line in f if line.strip().split('\t')[3]==strand]
    with open(trnageneWCRatiofile2) as f:
        wcratios2 = [line.strip().split('\t') for line in f if line.strip().split('\t')[3]==strand]
    print(len(wcratios1))
    print(len(wcratios2))
    wcratiosToInclude1 = [wcratios1[i] for i in range(len(wcratios1)) if i+1 not in tRNAsToExclude]
    wcratiosToInclude2 = [wcratios2[i] for i in range(len(wcratios2)) if i+1 not in tRNAsToExclude]
    print(len(wcratiosToInclude1))
    print(len(wcratiosToInclude2))
    oneminustwo = []
    for r in range(len(wcratiosToInclude1)):
        w = wcratiosToInclude1[r][4:]
        m = wcratiosToInclude2[r][4:]
        oneminustwo.append([float(i)-float(j) for i,j in zip(w,m)])
    avgratios = []
    for r in range(2*window):
        ratiosforr = [float(i[r]) for i in oneminustwo]
        avgratios.append(np.mean(ratiosforr))
    return avgratios

# This function spits out a plot of wc ratio of trnas found in Similar troughs
# and Dissimilar troughs.
def tRNAGenesWCRatioWithinTroughs(trnagenesfile,trnageneWCRatiofile1,trnageneWCRatiofile2,strand,tRNAsToExclude,troughpairfile,window1,window2):
    trnasFound = tRNAGenesFoundInMultipleRegions(trnagenesfile,troughpairfile,window1)
    trnas = []
    for i in trnasFound:
        trnas.extend([[str(j[0]),str(j[1]),str(j[2]),str(j[3])] for j in i[1]])
    print(len(trnas))  
    with open(trnageneWCRatiofile1) as f:
        wcratios1 = [line.strip().split('\t') for line in f if line.strip().split('\t')[3]==strand]
    with open(trnageneWCRatiofile2) as f:
        wcratios2 = [line.strip().split('\t') for line in f if line.strip().split('\t')[3]==strand]
    wcratiosToInclude1 = [wcratios1[i] for i in range(len(wcratios1)) if i+1 not in tRNAsToExclude]
    wcratiosToInclude2 = [wcratios2[i] for i in range(len(wcratios2)) if i+1 not in tRNAsToExclude]
    print(len(wcratiosToInclude1))
    print(len(wcratiosToInclude2))
    wcratiosFortrnas1 = [i for i in wcratiosToInclude1 if i[0:4] in trnas]
    wcratiosFortrnas2 = [i for i in wcratiosToInclude2 if i[0:4] in trnas]
    print(len(wcratiosFortrnas1))
    print(len(wcratiosFortrnas2))
    avgratios1 = []
    avgratios2 = []
    for r in range(4,2*window2+4):
        ratiosforr1 = [float(i[r]) for i in wcratiosFortrnas1]
        avgratios1.append(np.mean(ratiosforr1))
        ratiosforr2 = [float(i[r]) for i in wcratiosFortrnas2]
        avgratios2.append(np.mean(ratiosforr2))
    
    return [avgratios1,avgratios2]
    #return np.array(avgratios2) - np.array(avgratios1)
    
# This function calculate the number of tRNAs found per base
def calculateNumtRNAsFoundPerBase(tRNAsgenesfile,troughpairfile,window,NumPerBase):
    trnasFound = tRNAGenesFoundInMultipleRegions(tRNAsgenesfile,troughpairfile,window)
    if NumPerBase:    
        return [len(i[1])/float(i[0][2]-i[0][1]) for i in trnasFound]
    else:
        trnas = []
        for i in trnasFound:
            trnas.extend([[j[1],j[2],j[3],j[4]] for j in i[1]])
        return trnas

# This function looks at the watson crick ratio at each tRNA gene and writes it
# into a file about 10000 bases away
def tRNAGenesStrandWriteFile(trnagenesfile,sgrfile,window,writefile):
    trnas = wcs.extractStartEndChromStrandFortRNAs(trnagenesfile)
    with open(writefile,'w') as fw:
        for c in range(1,17):
            chromosome = str(c)
            with open(sgrfile) as f:
                sgrforchrom = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
            print(len(sgrforchrom))
            trnasforchrom = [i for i in trnas if i[0]==chromosome]
            for t in trnasforchrom:
                start = t[1]-window
                end = t[1]+window
                fw.write(t[0]+'\t'+str(t[1]) + '\t' + str(t[2]) + '\t' + t[3] + '\t')
                listToWrite = [str(i) for i in sgrforchrom[start:end]]
                fw.write('\t'.join(listToWrite))
                fw.write('\n')
# This function calculates the average watson crick ratio for every base pair 
# away +500 and =500 away from the tRNA gene 
def tRNAGenesStrandDensityMeanForStrand(trnagenesgrfile,tRNAsToExclude,window):
    with open(trnagenesgrfile) as f:
        sgrdata = [line.strip().split('\t') for line in f]
    print(len(sgrdata))
    sgrdataToInclude = [sgrdata[i] for i in range(len(sgrdata)) if i+1 not in tRNAsToExclude]
    print(len(sgrdataToInclude))
    avgsgr = []
    for r in range(4,2*window+4):
        ratiosforr = [float(i[r]) for i in sgrdataToInclude]
        avgsgr.append(np.mean(ratiosforr))
    return avgsgr
# This function spits out a plot of tRNAs found in similar and dissimilar troughs
# of two strains
def plotwcratiotRNAsWithinTroughs(strain1,strain2,strand,window1,window2,plotall):
    trnageneWCRatiofile1 = '../TestData/' + strain1 + 'tRNAGenesRawWCRatio-' + str(window2) + '.txt'
    trnageneWCRatiofile2 = '../TestData/' + strain2 + 'tRNAGenesRawWCRatio-' + str(window2) + '.txt'
    if strand=='W':
        tRNAsToExclude = [1,2,4,13,15,17,66,70,114,115,119,120]
    else:
        tRNAsToExclude = [5,7,8,31,33,35,39,40,57,59,64,111,113,114,124]
    troughpairfileSim = '../TestData/SimilarTroughs' + strain1 + strain2 + 'WOTransposons.txt'
    troughpairfileDissim = '../TestData/DissimilarTroughs' + strain1 + strain2 + 'WOTransposons.txt'
    wcratiosSim = tRNAGenesWCRatioWithinTroughs('DataToUse/tRNAGenes.txt',trnageneWCRatiofile1,trnageneWCRatiofile2,strand,tRNAsToExclude,troughpairfileSim,window1,window2)
    wcratiosDissim = tRNAGenesWCRatioWithinTroughs('DataToUse/tRNAGenes.txt',trnageneWCRatiofile1,trnageneWCRatiofile2,strand,tRNAsToExclude,troughpairfileDissim,window1,window2)
    wc_strain1Sim = wcratiosSim[0]
    wc_strain2Sim = wcratiosSim[1]
    wc_strain1Dissim = wcratiosDissim[0]
    wc_strain2Dissim = wcratiosDissim[1]
    if plotall:
        plt.plot(range(-window2,window2),wc_strain1Sim,label=strain1+'-Similar')
        plt.plot(range(-window2,window2),wc_strain2Sim,label=strain2+'-Similar')
        plt.plot(range(-window2,window2),np.array(wc_strain2Sim)-np.array(wc_strain1Sim),label=strain1 + strain2 + '-Diff-Similar')
        plt.plot(range(-window2,window2),wc_strain1Dissim,label=strain1+'-Dissimilar')
        plt.plot(range(-window2,window2),wc_strain2Dissim,label=strain2+'-Dissimilar')
        plt.plot(range(-window2,window2),np.array(wc_strain2Dissim)-np.array(wc_strain1Dissim),label=strain1 + strain2 + '-Diff-Dissimilar')
    else:
        plt.plot(range(-window2,window2),np.array(wc_strain2Sim)-np.array(wc_strain1Sim),label=strain1 + strain2 + '-Diff-Similar')
        plt.plot(range(-window2,window2),np.array(wc_strain2Dissim)-np.array(wc_strain1Dissim),label=strain1 + strain2 + '-Diff-Dissimilar')
    plt.xlim(-window2,window2)
    plt.xlabel('Bases')
    plt.ylabel('Watson Crick Log 2 Ratio')
    plt.axhline(xmax=window2,color='black')
    plt.axvline(x=0,color='black')
# This function counts number of CNG repeats found within a window around a 
# region. The region should contain the chromosome and start and end coordinate
def cngrepeatsFoundInRegion(cngreps,region,window):
    
    chromosome = region[0]
    start = region[1]
    end = region[2]
    
    validcngsForRegion = [i for i in cngreps if i[0]==chromosome and ( (int(i[1]) <= end+window and int(i[1]) >= start-window) or (int(i[2]) <= end+window and int(i[2]) >= start-window) ) ]
    return validcngsForRegion
    
# This function counts number of cng repeats in a bunch of regions
def cngRepsFoundInMultipleRegions(fastafile,repeat1,repeat2,regionfile,window):
    with open(regionfile) as f:
        regsBothWTandMUT = [line.strip().split('\t') for line in f][1:]
        
    cngRepsFound = []
    for c in range(1,17):
        chromosome = str(c)
        print chromosome
        regsforchrom = [i for i in regsBothWTandMUT if i[0]==chromosome]
        cngrepsforchrom = []   
        for r in range(3,14):
            cngrepsforchrom.extend(fr.getRepeatsForChromsome(fastafile,chromosome,repeat1,r))
            cngrepsforchrom.extend(fr.getRepeatsForChromsome(fastafile,chromosome,repeat2,r))
        #print cngrepsforchrom
        for i in regsforchrom:
            reg = [i[0],min(int(i[1]),int(i[3])),max(int(i[2]),int(i[4]))]
            cngRepsInReg = cngrepeatsFoundInRegion(cngrepsforchrom,reg,window)        
            cngRepsFound.append([reg,cngRepsInReg])
        
    return cngRepsFound
    
# This function counts the number of x repeats of CNG found in the regions of interest
def numberOfcngrepeatsFoundInRegion(allcngreps,numrepeats):
    xreps = [i[1] for i in allcngreps if len(i[1])>0]
    xrepslen = []
    for i in xreps:
        xrepslen.extend([j[2]-j[1] for j in i])
        
    xdrepslenOfinterest = [i for i in xrepslen if i==numrepeats*3]
    
    return len(xdrepslenOfinterest)
# This function grabs the G4 quads found within a window around a region. 
# The region contains the chromosome and the start and end coordinate    
def g4QuadsFoundInRegion(g4quadsfile,region,window):
    with open(g4quadsfile) as f:
        g4s = [line.strip().split('\t') for line in f]
    
    chromosome = region[0]
    start = region[1]
    end = region[2]
    
    validg4sForRegion = [i for i in g4s if i[0]==chromosome and ( (int(i[1]) <= end+window and int(i[1]) >= start-window) or (int(i[2]) <= end+window and int(i[2]) >= start-window) ) ]
    
    return validg4sForRegion
    
# This function counts number of tRNA genes in a bunch of regions
def g4QuadsFoundInMultipleRegions(g4quadsfile,regionfile,window):
    with open(regionfile) as f:
        regsBothWTandMUT = [line.strip().split('\t') for line in f][1:]
        
    g4sFound = []
    
    for i in regsBothWTandMUT:
        reg = [i[0],min(int(i[1]),int(i[3])),max(int(i[2]),int(i[4]))]
        g4sInReg = g4QuadsFoundInRegion(g4quadsfile,reg,window)        
        g4sFound.append([reg,g4sInReg])
        
    return g4sFound
    
# This function grabs the tRNA genes found within a window around a region. 
# The region contains the chromosome and the start and end coordinate 
def RNAPolIIGenesFoundInRegion(rnapol2genesfile,region,window):
    with open(rnapol2genesfile) as f:
        rnapol2genes = [line.strip().split('\t')[2:5] for line in f][1:]

    chromosome = region[0]
    start = region[1]
    end = region[2]
    
    validrnapol2genesForRegion = [i for i in rnapol2genes if i[0]==chromosome and ( (int(i[1]) <= end+window and int(i[1]) >= start-window) or (int(i[2]) <= end+window and int(i[2]) >= start-window) ) ]
    
    return validrnapol2genesForRegion
    
# This function counts number of tRNA genes in a bunch of regions
def RNAPolIIGenesFoundInMultipleRegions(rnapol2genesfile,regionfile,window):
    with open(regionfile) as f:
        regsBothWTandMUT = [line.strip().split('\t') for line in f][1:]
        
    RNAPolIIGenesFound = []
    
    for i in regsBothWTandMUT:
        reg = [i[0],min(int(i[1]),int(i[3])),max(int(i[2]),int(i[4]))]
        RNAPolIIGenesInReg = RNAPolIIGenesFoundInRegion(rnapol2genesfile,reg,window)        
        RNAPolIIGenesFound.append([reg,RNAPolIIGenesInReg])
        
    return RNAPolIIGenesFound

# This function looks at the watson crick ratio at each tRNA gene and writes it
# into a file about 10000 bases away
def RNAPolIIGenesWCRatioWriteFile(rnapol2genesfile,wcratiofile,window,writefile):
    with open(rnapol2genesfile) as f:
        rnapol2genes = [line.strip().split('\t')[2:6] for line in f][1:]
    with open(writefile,'w') as fw:
        for c in range(1,17):
            chromosome = str(c)
            with open(wcratiofile) as f:
                wcratiosforchrom = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
            print(len(wcratiosforchrom))
            rnapol2genessforchrom = [i for i in rnapol2genes if i[0]==chromosome]
            for t in rnapol2genessforchrom:
                start = int(t[1])-window
                end = int(t[1])+window
                fw.write(t[0]+'\t'+str(t[1]) + '\t' + str(t[2]) + '\t' + t[3] + '\t')
                listToWrite = [str(i) for i in wcratiosforchrom[start:end]]
                fw.write('\t'.join(listToWrite))
                fw.write('\n')
                
# This function calculates the average watson crick ratio for every base pair 
# away +500 and =500 away from the highly transcribed RNA Pol II gene 
def RNAPolIIGenesWCRatioMeanForStrand(rnapol2geneWCRatiofile,strand,RNAPol2GenesToExclude,window):
    with open(rnapol2geneWCRatiofile) as f:
        wcratios = [line.strip().split('\t') for line in f if line.strip().split('\t')[3]==strand]
    print(len(wcratios))
    wcratiosToInclude = [wcratios[i] for i in range(len(wcratios)) if i+1 not in RNAPol2GenesToExclude]
    print(len(wcratiosToInclude))
    avgratios = []
    for r in range(4,2*window+4):
        ratiosforr = [float(i[r]) for i in wcratiosToInclude]
        avgratios.append(np.mean(ratiosforr))
    return avgratios
    
# This function calculates the average watson crick ratio for every base pair 
# away +500 and =500 away from the tRNA gene 
def RNAPolIIGenesDiffTwoStrainsWCRatioMeanForStrand(rnapol2geneWCRatiofile1,rnapol2geneWCRatiofile2,strand,RNAPol2GenesToExclude,window):
    with open(rnapol2geneWCRatiofile1) as f:
        wcratios1 = [line.strip().split('\t') for line in f if line.strip().split('\t')[3]==strand]
    with open(rnapol2geneWCRatiofile2) as f:
        wcratios2 = [line.strip().split('\t') for line in f if line.strip().split('\t')[3]==strand]
    print(len(wcratios1))
    print(len(wcratios2))
    wcratiosToInclude1 = [wcratios1[i] for i in range(len(wcratios1)) if i+1 not in RNAPol2GenesToExclude]
    wcratiosToInclude2 = [wcratios2[i] for i in range(len(wcratios2)) if i+1 not in RNAPol2GenesToExclude]
    print(len(wcratiosToInclude1))
    print(len(wcratiosToInclude2))
    oneminustwo = []
    for r in range(len(wcratiosToInclude1)):
        w = wcratiosToInclude1[r][4:]
        m = wcratiosToInclude2[r][4:]
        oneminustwo.append([float(i)-float(j) for i,j in zip(w,m)])
    avgratios = []
    for r in range(2*window):
        ratiosforr = [float(i[r]) for i in oneminustwo]
        avgratios.append(np.mean(ratiosforr))
    return avgratios
 
# This function looks at the watson crick ratio at each tRNA gene and writes it
# into a file about 10000 bases away
def RandomSpotsWCRatioWriteFile(numrandomspots,size,wcratiofile,window,writefile):
    rs = []
    for i in range(numrandomspots):
        chromosome = str(random.randint(1,16))
        randCoord = random.randint(window,yeastsize[chromosome]+1-window)
    
        rs.append([chromosome, randCoord,randCoord+size])
        
    with open(writefile,'w') as fw:
        for c in range(1,17):
            chromosome = str(c)
            with open(wcratiofile) as f:
                wcratiosforchrom = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
            print(len(wcratiosforchrom))
            rsforchrom = [i for i in rs if i[0]==chromosome]
            for t in rsforchrom:
                start = int(t[1])-window
                end = int(t[1])+window
                fw.write(t[0]+'\t'+str(t[1]) + '\t' + str(t[2]) + '\t')
                listToWrite = [str(i) for i in wcratiosforchrom[start:end]]
                fw.write('\t'.join(listToWrite))
                fw.write('\n')
# This function calculates the average watson crick ratio for every base pair 
# away +500 and =500 away from the highly transcribed RNA Pol II gene 
def RandomSpotsWCRatioMeanForStrand(randomspotsWCRatiofile,window):
    with open(randomspotsWCRatiofile) as f:
        wcratios = [line.strip().split('\t') for line in f]
    avgratios = []
    for r in range(3,2*window+3):
        ratiosforr = [float(i[r]) for i in wcratios]
        avgratios.append(np.mean(ratiosforr))
    return avgratios
#############################PLOT#########################################
# Function to plot watson-crick ratios


def plotWatsonCrickRatioForOriginPair(vals,originpair,window,strain):
    valsbetweenOrg = vals[originpair[0]:originpair[1]]
    print("Calculate first derivative")
    firstder = calculateDer(valsbetweenOrg,1)
    print("Calculate second derivative")
    secder = calculateDer(valsbetweenOrg,2)
    plt.plot(valsbetweenOrg,label=strain + '-Original')
    plt.plot(firstder,label=strain + '-First Der')
    plt.plot(secder,label=strain+'-Second Der')
    plt.xlim(xmax=len(valsbetweenOrg))
    # Set a clear x axis
    plt.axhline(xmax=len(vals),color='black')
    plt.ylabel('Log2 Watson Crick Ratio')
    plt.xlabel('Base')
    plt.legend(loc='best')
def plotWatsonCrickRatioFromFile(wcratiofile,chromosome,sbnum,color,strain):
    with open(wcratiofile) as f:
        vals = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
    # Start up a figure
    fig = plt.figure(1)
    # Add a subplot
    ax = fig.add_subplot(sbnum)
    # Plot reference base vs OEM
    ax.plot(vals,color=color,label=strain)
    # Set the x limit to be the number of bases for that chromosome
    plt.xlim(xmax=len(vals))
    # Set y limits to be between -1 and 1 since that's max and min of OEM
    #plt.ylim(ymin=-1,ymax=1)
    # Set a clear x axis
    plt.axhline(xmax=len(vals),color='black')
    # Fill in the curves so we get filled in graph for WT in black
    #plt.fill_between(ind,just_oem,color=color)
    plt.ylabel('Log2 Watson Crick Ratio')
    plt.xlabel('Base')
    plt.show()

def plotTwoWatsonCrickRatioFromFile(wcratiofile1,wcratiofile2,chromosome,strain1,strain2):
    with open(wcratiofile1) as f:
        vals1 = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
    with open(wcratiofile2) as f:
        vals2 = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
    
    # Plot first file
    plt.plot(vals1,'b',label=strain1)
    # Plot second file
    plt.plot(vals2,'g',label=strain2)
    # Set the x limit to be the number of bases for that chromosome
    plt.xlim(xmax=len(vals1))
    # Set a clear x axis
    plt.axhline(xmax=len(vals1),color='black')
    
    plt.legend(loc='best')
    plt.ylabel('Log2 Watson Crick Ratio')
    plt.xlabel('Base')
    
    plt.show()
    
def plotTwoGradientWatsonCrickRatioFromFile(wcratiofile1,wcratiofile2,chromosome,window,strain1,strain2):
    with open(wcratiofile1) as f:
        vals1 = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
    with open(wcratiofile2) as f:
        vals2 = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
    
    # Calculate Gradient of vals1 and vals2
    gvals1 = calculateDer(vals1,window)
    gvals2 = calculateDer(vals2,window)
    
    # Plot first file
    plt.plot(gvals1,label=strain1)
    # Plot second file
    plt.plot(gvals2,label=strain2)
    # Set the x limit to be the number of bases for that chromosome
    plt.xlim(xmax=len(gvals1))
    # Set a clear x axis
    plt.axhline(xmax=len(gvals1),color='black')
    
    plt.legend(loc='best')
    plt.ylabel('First Derivative Log2 Watson Crick Ratio')
    plt.xlabel('Base')
    
    plt.show()
    
def plotTwoSecDerWatsonCrickRatioFromFile(wcratiofile1,wcratiofile2,chromosome,window,strain1,strain2):
    with open(wcratiofile1) as f:
        vals1 = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
    with open(wcratiofile2) as f:
        vals2 = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
    
    # Calculate Gradient of vals1 and vals2
    g1vals1 = calculateDer(vals1,window)
    g1vals2 = calculateDer(vals2,window)
    gvals1 = calculateDer(g1vals1,window)
    gvals2 = calculateDer(g1vals2,window)
    # Plot first file
    plt.plot(np.array(gvals1)-np.array(gvals2),label=strain1)
    # Plot second file
    #plt.plot(gvals2,label=strain2)
    # Set the x limit to be the number of bases for that chromosome
    plt.xlim(xmax=len(gvals1))
    # Set a clear x axis
    plt.axhline(xmax=len(gvals1),color='black')
    
    plt.legend(loc='best')
    plt.ylabel('Diff Second Derivative Log2 Watson Crick Ratio')
    plt.xlabel('Base')
    
    plt.show()

def plotWCRatiotRNAGenesSpecificallyForAminoAcid(aminoacid):
    aminoacid_cdc9deg_watson = tRNAGenesWCRatioMeanForStrand('../TestData/cdc9degtRNAGenesRawWCRatio-' + aminoacid + '-5000.txt','W',[],5000)
    aminoacid_cdc9deg_crick = tRNAGenesWCRatioMeanForStrand('../TestData/cdc9degtRNAGenesRawWCRatio-' + aminoacid + '-5000.txt','C',[],5000)
    aminoacid_rrm3d_watson = tRNAGenesWCRatioMeanForStrand('../TestData/rrm3dtRNAGenesRawWCRatio-' + aminoacid + '-5000.txt','W',[],5000)
    aminoacid_rrm3d_crick = tRNAGenesWCRatioMeanForStrand('../TestData/rrm3dtRNAGenesRawWCRatio-' + aminoacid + '-5000.txt','C',[],5000)
    aminoacid_pif1m2_watson = tRNAGenesWCRatioMeanForStrand('../TestData/pif1m2tRNAGenesRawWCRatio-' + aminoacid + '-5000.txt','W',[],5000)
    aminoacid_pif1m2_crick = tRNAGenesWCRatioMeanForStrand('../TestData/pif1m2tRNAGenesRawWCRatio-' + aminoacid + '-5000.txt','C',[],5000)
    aminoacid_rrm3dpif1m2_watson = tRNAGenesWCRatioMeanForStrand('../TestData/rrm3d_pif1m2tRNAGenesRawWCRatio-' + aminoacid + '-5000.txt','W',[],5000)
    aminoacid_rrm3dpif1m2_crick = tRNAGenesWCRatioMeanForStrand('../TestData/rrm3d_pif1m2tRNAGenesRawWCRatio-' + aminoacid + '-5000.txt','C',[],5000)
    
    fig = plt.figure(1)
    ax = fig.add_subplot(121)
    ax.plot(range(-5000,5000),np.array(aminoacid_rrm3d_watson)-np.array(aminoacid_cdc9deg_watson),label='rrm3d')
    ax.plot(range(-5000,5000),np.array(aminoacid_pif1m2_watson)-np.array(aminoacid_cdc9deg_watson),label='pif1m2')
    ax.plot(range(-5000,5000),np.array(aminoacid_rrm3dpif1m2_watson)-np.array(aminoacid_cdc9deg_watson),label='rrm3dpif1m2')
    plt.legend(loc='best')
    plt.xlabel('Bases')
    plt.xlim(-5000,5000)
    plt.axhline(xmax=5000,color='black')
    plt.axvline(x=0,color='black') 
    plt.ylabel('Log2 Watson-Crick Ratio')
    plt.title('Watson-Crick Ratio Around ' + aminoacid + ' tRNA Watson Genes')

    ax = fig.add_subplot(122)
    ax.plot(range(-5000,5000),np.array(aminoacid_rrm3d_crick)-np.array(aminoacid_cdc9deg_crick),label='rrm3d')
    ax.plot(range(-5000,5000),np.array(aminoacid_pif1m2_crick)-np.array(aminoacid_cdc9deg_crick),label='pif1m2')
    ax.plot(range(-5000,5000),np.array(aminoacid_rrm3dpif1m2_crick)-np.array(aminoacid_cdc9deg_crick),label='rrm3dpif1m2')
    plt.legend(loc='best')
    plt.xlabel('Bases')
    plt.xlim(-5000,5000)
    plt.axhline(xmax=5000,color='black')
    plt.axvline(x=0,color='black') 
    plt.ylabel('Log2 Watson-Crick Ratio')
    plt.title('Watson-Crick Ratio Around ' + aminoacid + ' tRNA Crick Genes')
    
    plt.show()

def plotWCRatiotRNAGenesSpecificallyForAminoAcidGroups(aminoacidgroup):
    if aminoacidgroup=='non-polar':
        aminoacids = ['Alanine','Glycine','Valine','Leucine','Methionine','Isoleucine']
    elif aminoacidgroup == 'polar':
        aminoacids = ['Serine','Threonine','Cysteine','Proline','Asparagine','Glutamine']
    elif aminoacidgroup == 'aromatic':
        aminoacids = ['Phenylalanine','Tyrosine','Tryptophan']
    elif aminoacidgroup == 'positive-R-group':
        aminoacids = ['Lysine','Arginine','Histidine']
    elif aminoacidgroup == 'negative-R-group':
        aminoacids = ['Aspartic_acid','Glutamic_acid']
    
    aminoacid_cdc9deg_watson = []
    aminoacid_cdc9deg_crick = []
    aminoacid_rrm3d_watson = []
    aminoacid_rrm3d_crick = []
    aminoacid_pif1m2_watson = []
    aminoacid_pif1m2_crick = []
    aminoacid_rrm3dpif1m2_watson = []
    aminoacid_rrm3dpif1m2_crick = []
    
    for i in aminoacids:    
        aminoacid_cdc9deg_watson.append(tRNAGenesWCRatioMeanForStrand('../TestData/cdc9degtRNAGenesRawWCRatio-' + i + '-5000.txt','W',[],5000))
        aminoacid_cdc9deg_crick.append(tRNAGenesWCRatioMeanForStrand('../TestData/cdc9degtRNAGenesRawWCRatio-' + i + '-5000.txt','C',[],5000))
        aminoacid_rrm3d_watson.append(tRNAGenesWCRatioMeanForStrand('../TestData/rrm3dtRNAGenesRawWCRatio-' + i + '-5000.txt','W',[],5000))
        aminoacid_rrm3d_crick.append(tRNAGenesWCRatioMeanForStrand('../TestData/rrm3dtRNAGenesRawWCRatio-' + i + '-5000.txt','C',[],5000))
        aminoacid_pif1m2_watson.append(tRNAGenesWCRatioMeanForStrand('../TestData/pif1m2tRNAGenesRawWCRatio-' + i + '-5000.txt','W',[],5000))
        aminoacid_pif1m2_crick.append(tRNAGenesWCRatioMeanForStrand('../TestData/pif1m2tRNAGenesRawWCRatio-' + i + '-5000.txt','C',[],5000))
        aminoacid_rrm3dpif1m2_watson.append(tRNAGenesWCRatioMeanForStrand('../TestData/rrm3d_pif1m2tRNAGenesRawWCRatio-' + i + '-5000.txt','W',[],5000))
        aminoacid_rrm3dpif1m2_crick.append(tRNAGenesWCRatioMeanForStrand('../TestData/rrm3d_pif1m2tRNAGenesRawWCRatio-' + i + '-5000.txt','C',[],5000))
    
    avgratios_cdc9deg_watson = []
    avgratios_cdc9deg_crick = []
    avgratios_rrm3d_watson = []
    avgratios_rrm3d_crick = []
    avgratios_pif1m2_watson = []
    avgratios_pif1m2_crick = []
    avgratios_rrm3dpif1m2_watson = []
    avgratios_rrm3dpif1m2_crick = []
    
    for r in range(10000):
        print r
        ratiosforr = [float(j[r]) for j in aminoacid_cdc9deg_watson]
        avgratios_cdc9deg_watson.append(np.mean(ratiosforr))
        ratiosforr = [float(j[r]) for j in aminoacid_cdc9deg_crick]
        avgratios_cdc9deg_crick.append(np.mean(ratiosforr))
        ratiosforr = [float(j[r]) for j in aminoacid_rrm3d_watson]
        avgratios_rrm3d_watson.append(np.mean(ratiosforr))
        ratiosforr = [float(j[r]) for j in aminoacid_rrm3d_crick]
        avgratios_rrm3d_crick.append(np.mean(ratiosforr))
        ratiosforr = [float(j[r]) for j in aminoacid_pif1m2_watson]
        avgratios_pif1m2_watson.append(np.mean(ratiosforr))
        ratiosforr = [float(j[r]) for j in aminoacid_pif1m2_crick]
        avgratios_pif1m2_crick.append(np.mean(ratiosforr))
        ratiosforr = [float(j[r]) for j in aminoacid_rrm3dpif1m2_watson]
        avgratios_rrm3dpif1m2_watson.append(np.mean(ratiosforr))
        ratiosforr = [float(j[r]) for j in aminoacid_rrm3dpif1m2_crick]
        avgratios_rrm3dpif1m2_crick.append(np.mean(ratiosforr))
        
    fig = plt.figure(1)
    ax = fig.add_subplot(121)
    ax.plot(range(-5000,5000),np.array(avgratios_rrm3d_watson)-np.array(avgratios_cdc9deg_watson),label='rrm3d')
    ax.plot(range(-5000,5000),np.array(avgratios_pif1m2_watson)-np.array(avgratios_cdc9deg_watson),label='pif1m2')
    ax.plot(range(-5000,5000),np.array(avgratios_rrm3dpif1m2_watson)-np.array(avgratios_cdc9deg_watson),label='rrm3dpif1m2')
    plt.legend(loc='best')
    plt.xlabel('Bases')
    plt.xlim(-5000,5000)
    plt.axhline(xmax=5000,color='black')
    plt.axvline(x=0,color='black') 
    plt.ylabel('Log2 Watson-Crick Ratio')
    plt.title('Watson-Crick Ratio Around ' + aminoacidgroup + ' tRNA Watson Genes')

    ax = fig.add_subplot(122)
    ax.plot(range(-5000,5000),np.array(avgratios_rrm3d_crick)-np.array(avgratios_cdc9deg_crick),label='rrm3d')
    ax.plot(range(-5000,5000),np.array(avgratios_pif1m2_crick)-np.array(avgratios_cdc9deg_crick),label='pif1m2')
    ax.plot(range(-5000,5000),np.array(avgratios_rrm3dpif1m2_crick)-np.array(avgratios_cdc9deg_crick),label='rrm3dpif1m2')
    plt.legend(loc='best')
    plt.xlabel('Bases')
    plt.xlim(-5000,5000)
    plt.axhline(xmax=5000,color='black')
    plt.axvline(x=0,color='black') 
    plt.ylabel('Log2 Watson-Crick Ratio')
    plt.title('Watson-Crick Ratio Around ' + aminoacidgroup + ' tRNA Crick Genes')
    
    plt.show()
"""   
convcdc9degW = wcr.tRNAGenesWCRatioMeanForStrand('../TestData/cdc9degConvergentWatsontRNAGenesRawWCRatio-5000.txt','W',[],5000)
convrrm3dW = wcr.tRNAGenesWCRatioMeanForStrand('../TestData/rrm3dConvergentWatsontRNAGenesRawWCRatio-5000.txt','W',[],5000)
convrrm3dpif1m2W = wcr.tRNAGenesWCRatioMeanForStrand('../TestData/rrm3d_pif1m2ConvergentWatsontRNAGenesRawWCRatio-5000.txt','W',[],5000)
convpif1m2W = wcr.tRNAGenesWCRatioMeanForStrand('../TestData/pif1m2ConvergentWatsontRNAGenesRawWCRatio-5000.txt','W',[],5000)
divcdc9degW = wcr.tRNAGenesWCRatioMeanForStrand('../TestData/cdc9degDivergentWatsontRNAGenesRawWCRatio-5000.txt','W',[],5000)
divrrm3dW = wcr.tRNAGenesWCRatioMeanForStrand('../TestData/rrm3dDivergentWatsontRNAGenesRawWCRatio-5000.txt','W',[],5000)
divrrm3dpif1m2W = wcr.tRNAGenesWCRatioMeanForStrand('../TestData/rrm3d_pif1m2DivergentWatsontRNAGenesRawWCRatio-5000.txt','W',[],5000)
divpif1m2W = wcr.tRNAGenesWCRatioMeanForStrand('../TestData/pif1m2DivergentWatsontRNAGenesRawWCRatio-5000.txt','W',[],5000)
convcdc9degC = wcr.tRNAGenesWCRatioMeanForStrand('../TestData/cdc9degConvergentCricktRNAGenesRawWCRatio-5000.txt','C',[],5000)
convrrm3dC = wcr.tRNAGenesWCRatioMeanForStrand('../TestData/rrm3dConvergentCricktRNAGenesRawWCRatio-5000.txt','C',[],5000)
convpif1m2C = wcr.tRNAGenesWCRatioMeanForStrand('../TestData/pif1m2ConvergentCricktRNAGenesRawWCRatio-5000.txt','C',[],5000)
convrrm3dpif1m2C = wcr.tRNAGenesWCRatioMeanForStrand('../TestData/rrm3d_pif1m2ConvergentCricktRNAGenesRawWCRatio-5000.txt','C',[],5000)
divcdc9degC = wcr.tRNAGenesWCRatioMeanForStrand('../TestData/cdc9degDivergentCricktRNAGenesRawWCRatio-5000.txt','C',[],5000)
divrrm3dC = wcr.tRNAGenesWCRatioMeanForStrand('../TestData/rrm3dDivergentCricktRNAGenesRawWCRatio-5000.txt','C',[],5000)
divpif1m2C = wcr.tRNAGenesWCRatioMeanForStrand('../TestData/pif1m2DivergentCricktRNAGenesRawWCRatio-5000.txt','C',[],5000)
divrrm3dpif1m2C = wcr.tRNAGenesWCRatioMeanForStrand('../TestData/rrm3d_pif1m2DivergentCricktRNAGenesRawWCRatio-5000.txt','C',[],5000)

plt.plot(range(-5000,5000),convcdc9degC,label='cdc9deg-Crick',color='blue')
plt.plot(range(-5000,5000),convrrm3dC,label='rrm3d-Crick',color='crimson')
plt.plot(range(-5000,5000),convpif1m2C,label='pif1m2-Crick',color='green')
plt.plot(range(-5000,5000),convrrm3dpif1m2C,label='rrm3d_pif1m2-Crick',color='darkorange')

plt.plot(range(-5000,5000),convcdc9degW,label='cdc9deg-Watson',color='cornflowerblue')
plt.plot(range(-5000,5000),convrrm3dW,label='rrm3d-Watson',color='darkorchid')
plt.plot(range(-5000,5000),convpif1m2W,label='pif1m2-Watson',color='darkseagreen')
plt.plot(range(-5000,5000),convrrm3dpif1m2W,label='rrm3d_pif1m2-Watson',color='tomato')

plt.legend(loc='best')
plt.xlim(-5000,5000)
plt.axhline(xmax=5000,color='black')
plt.axvline(x=0,color='black') 
plt.xlabel('Bases')
plt.ylabel('Watson Crick Log 2 Ratio')
plt.title('Watson-Crick Ratio Around Convergent tRNA Genes')
plt.tight_layout()

plt.plot(range(-5000,5000),divcdc9degC,label='cdc9deg-Crick',color='blue')
plt.plot(range(-5000,5000),divrrm3dC,label='rrm3d-Crick',color='crimson')
plt.plot(range(-5000,5000),divpif1m2C,label='pif1m2-Crick',color='green')
plt.plot(range(-5000,5000),divrrm3dpif1m2C,label='rrm3d_pif1m2-Crick',color='darkorange')

plt.plot(range(-5000,5000),divcdc9degW,label='cdc9deg-Watson',color='cornflowerblue')
plt.plot(range(-5000,5000),divrrm3dW,label='rrm3d-Watson',color='darkorchid')
plt.plot(range(-5000,5000),divpif1m2W,label='pif1m2-Watson',color='darkseagreen')
plt.plot(range(-5000,5000),divrrm3dpif1m2W,label='rrm3d_pif1m2-Watson',color='tomato')

plt.legend(loc='best')
plt.xlim(-5000,5000)
plt.axhline(xmax=5000,color='black')
plt.axvline(x=0,color='black') 
plt.xlabel('Bases')
plt.ylabel('Watson Crick Log 2 Ratio')
plt.title('Watson-Crick Ratio Around Divergent tRNA Genes')
plt.tight_layout()

convcdc9deg_combo = (np.array(convcdc9degC) + np.array(convcdc9degW))/2
convrrm3d_combo = (np.array(convrrm3dC) + np.array(convrrm3dW))/2
convpif1m2_combo = (np.array(convpif1m2C) + np.array(convpif1m2W))/2
convrrm3dpif1m2_combo = (np.array(convrrm3dpif1m2C) + np.array(convrrm3dpif1m2W))/2

plt.plot(range(-5000,5000),convrrm3d_combo - convcdc9deg_combo,label='rrm3d',color='blue')
plt.plot(range(-5000,5000),convpif1m2_combo - convcdc9deg_combo,label='pif1m2',color='crimson')
plt.plot(range(-5000,5000),convrrm3dpif1m2_combo - convcdc9deg_combo,label='rrm3d_pif1m2',color='green')


plt.legend(loc='best')
plt.xlim(-5000,5000)
plt.axhline(xmax=5000,color='black')
plt.axvline(x=0,color='black') 
plt.xlabel('Bases')
plt.ylabel('Watson Crick Log 2 Ratio')
plt.title('Watson-Crick Ratio Around Convergent tRNA Genes Normalized to WT')
plt.tight_layout()

divcdc9deg_combo = (np.array(divcdc9degC) + np.array(divcdc9degW))/2
divrrm3d_combo = (np.array(divrrm3dC) + np.array(divrrm3dW))/2
divpif1m2_combo = (np.array(divpif1m2C) + np.array(divpif1m2W))/2
divrrm3dpif1m2_combo = (np.array(divrrm3dpif1m2C) + np.array(divrrm3dpif1m2W))/2

plt.plot(range(-5000,5000),divrrm3d_combo - divcdc9deg_combo,label='rrm3d',color='blue')
plt.plot(range(-5000,5000),divpif1m2_combo - divcdc9deg_combo,label='pif1m2',color='crimson')
plt.plot(range(-5000,5000),divrrm3dpif1m2_combo - divcdc9deg_combo,label='rrm3d_pif1m2',color='green')


plt.legend(loc='best')
plt.xlim(-5000,5000)
plt.axhline(xmax=5000,color='black')
plt.axvline(x=0,color='black') 
plt.xlabel('Bases')
plt.ylabel('Watson Crick Log 2 Ratio')
plt.title('Watson-Crick Ratio Around Divergent tRNA Genes Normalized to WT')
plt.tight_layout()


plt.plot(range(-5000,5000),convcdc9degC,label='cdc9deg',color='blue')
plt.plot(range(-5000,5000),convrrm3dC,label='rrm3d',color='crimson')
plt.plot(range(-5000,5000),convpif1m2C,label='pif1m2',color='green')
plt.plot(range(-5000,5000),convrrm3dpif1m2C,label='rrm3d_pif1m2',color='darkorange')
plt.plot(range(-5000,5000),divcdc9degC,color='blue')
plt.plot(range(-5000,5000),divrrm3dC,color='crimson')
plt.plot(range(-5000,5000),divpif1m2C,color='green')
plt.plot(range(-5000,5000),divrrm3dpif1m2C,color='darkorange')
plt.legend(loc='best')
plt.xlim(-5000,5000)
plt.axhline(xmax=5000,color='black')
plt.axvline(x=0,color='black') 
plt.xlabel('Bases')
plt.ylabel('Watson Crick Log 2 Ratio')
plt.title('Watson-Crick Ratio Around Crick tRNA Genes - Convergent and Divergent')
plt.tight_layout()

plt.plot(range(-5000,5000),convcdc9degW,label='cdc9deg',color='cornflowerblue')
plt.plot(range(-5000,5000),convrrm3dW,label='rrm3d',color='darkorchid')
plt.plot(range(-5000,5000),convpif1m2W,label='pif1m2',color='darkseagreen')
plt.plot(range(-5000,5000),convrrm3dpif1m2W,label='rrm3d_pif1m2',color='tomato')
plt.plot(range(-5000,5000),divcdc9degW,color='cornflowerblue')
plt.plot(range(-5000,5000),divrrm3dW,color='darkorchid')
plt.plot(range(-5000,5000),divpif1m2W,color='darkseagreen')
plt.plot(range(-5000,5000),divrrm3dpif1m2W,color='tomato')
plt.legend(loc='best')
plt.xlim(-5000,5000)
plt.axhline(xmax=5000,color='black')
plt.axvline(x=0,color='black') 
plt.xlabel('Bases')
plt.ylabel('Watson Crick Log 2 Ratio')
plt.title('Watson-Crick Ratio Around Watson tRNA Genes - Convergent and Divergent')
plt.tight_layout()

fig = plt.figure(1)
ax = fig.add_subplot(211)
ax.plot(range(-5000,5000),np.array(convrrm3dW) - np.array(convcdc9degW),label='rrm3d-Watson',color='blue')
ax.plot(range(-5000,5000),np.array(convpif1m2W) - np.array(convcdc9degW),label='pif1m2-Watson',color='crimson')
ax.plot(range(-5000,5000),np.array(convrrm3dpif1m2W) - np.array(convcdc9degW),label='rrm3d_pif1m2-Watson',color='green')
plt.legend(loc='best')
plt.xlim(-5000,5000)
plt.axhline(xmax=5000,color='black')
plt.axvline(x=0,color='black') 
plt.xlabel('Bases')
plt.ylabel('Watson Crick Log 2 Ratio')

ax2 = fig.add_subplot(212)
plt.plot(range(-5000,5000),np.array(convrrm3dC) - np.array(convcdc9degC),label='rrm3d-Crick',color='blue')
plt.plot(range(-5000,5000),np.array(convpif1m2C) - np.array(convcdc9degC),label='pif1m2-Crick',color='crimson')
plt.plot(range(-5000,5000),np.array(convrrm3dpif1m2C) - np.array(convcdc9degC),label='rrm3d_pif1m2-Crick',color='green')

plt.legend(loc='best')
plt.xlim(-5000,5000)
plt.axhline(xmax=5000,color='black')
plt.axvline(x=0,color='black') 
plt.xlabel('Bases')
plt.ylabel('Watson Crick Log 2 Ratio')
plt.title('Watson-Crick Ratio Around Convergent Watson tRNA Genes Normalized to WT')
plt.tight_layout()  


fig = plt.figure(1)
ax = fig.add_subplot(211)
ax.plot(range(-5000,5000),np.array(divrrm3dW) - np.array(divcdc9degW),label='rrm3d-Watson',color='darkorchid')
ax.plot(range(-5000,5000),np.array(divpif1m2W) - np.array(divcdc9degW),label='pif1m2-Watson',color='darkseagreen')
ax.plot(range(-5000,5000),np.array(divrrm3dpif1m2W) - np.array(divcdc9degW),label='rrm3d_pif1m2-Watson',color='tomato')
plt.legend(loc='best')
plt.xlim(-5000,5000)
plt.axhline(xmax=5000,color='black')
plt.axvline(x=0,color='black') 
plt.xlabel('Bases')
plt.ylabel('Watson Crick Log 2 Ratio')

ax2 = fig.add_subplot(212)
ax2.plot(range(-5000,5000),np.array(divrrm3dC) - np.array(divcdc9degC),label='rrm3d-Crick',color='darkorchid')
ax2.plot(range(-5000,5000),np.array(divpif1m2C) - np.array(divcdc9degC),label='pif1m2-Crick',color='darkseagreen')
ax2.plot(range(-5000,5000),np.array(divrrm3dpif1m2C) - np.array(divcdc9degC),label='rrm3d_pif1m2-Crick',color='tomato')

plt.legend(loc='best')
plt.xlim(-5000,5000)
plt.axhline(xmax=5000,color='black')
plt.axvline(x=0,color='black') 
plt.xlabel('Bases')
plt.ylabel('Watson Crick Log 2 Ratio')
plt.title('Watson-Crick Ratio Around Convergent Watson tRNA Genes Normalized to WT')
plt.tight_layout()  

cdc9degConvC_Crick = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/cdc9degConvergentCricktRNAGenesRawC-5000.txt',[],5000)
cdc9degConvC_Watson = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/cdc9degConvergentCricktRNAGenesRawW-5000.txt',[],5000)
cdc9degConvW_Watson = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/cdc9degConvergentWatsontRNAGenesRawW-5000.txt',[],5000)
cdc9degConvW_Crick = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/cdc9degConvergentWatsontRNAGenesRawC-5000.txt',[],5000)

rrm3dConvC_Crick = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/rrm3dConvergentCricktRNAGenesRawC-5000.txt',[],5000)
rrm3dConvC_Watson = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/rrm3dConvergentCricktRNAGenesRawW-5000.txt',[],5000)
rrm3dConvW_Watson = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/rrm3dConvergentWatsontRNAGenesRawW-5000.txt',[],5000)
rrm3dConvW_Crick = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/rrm3dConvergentWatsontRNAGenesRawC-5000.txt',[],5000)

pif1m2ConvC_Crick = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/pif1m2ConvergentCricktRNAGenesRawC-5000.txt',[],5000)
pif1m2ConvC_Watson = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/pif1m2ConvergentCricktRNAGenesRawW-5000.txt',[],5000)
pif1m2ConvW_Watson = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/pif1m2ConvergentWatsontRNAGenesRawW-5000.txt',[],5000)
pif1m2ConvW_Crick = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/pif1m2ConvergentWatsontRNAGenesRawC-5000.txt',[],5000)

rrm3dpif1m2ConvC_Crick = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/rrm3d_pif1m2ConvergentCricktRNAGenesRawC-5000.txt',[],5000)
rrm3dpif1m2ConvW_Watson = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/rrm3d_pif1m2ConvergentWatsontRNAGenesRawW-5000.txt',[],5000)
rrm3dpif1m2ConvC_Watson = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/rrm3d_pif1m2ConvergentCricktRNAGenesRawW-5000.txt',[],5000)
rrm3dpif1m2ConvW_Crick = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/rrm3d_pif1m2ConvergentWatsontRNAGenesRawC-5000.txt',[],5000)

cdc9degDivC_Crick = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/cdc9degDivergentCricktRNAGenesRawC-5000.txt',[],5000)
cdc9degDivC_Watson = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/cdc9degDivergentCricktRNAGenesRawW-5000.txt',[],5000)
cdc9degDivW_Watson = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/cdc9degDivergentWatsontRNAGenesRawW-5000.txt',[],5000)
cdc9degDivW_Crick = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/cdc9degDivergentWatsontRNAGenesRawC-5000.txt',[],5000)

rrm3dDivC_Crick = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/rrm3dDivergentCricktRNAGenesRawC-5000.txt',[],5000)
rrm3dDivC_Watson = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/rrm3dDivergentCricktRNAGenesRawW-5000.txt',[],5000)
rrm3dDivW_Watson = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/rrm3dDivergentWatsontRNAGenesRawW-5000.txt',[],5000)
rrm3dDivW_Crick = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/rrm3dDivergentWatsontRNAGenesRawC-5000.txt',[],5000)

pif1m2DivC_Crick = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/pif1m2DivergentCricktRNAGenesRawC-5000.txt',[],5000)
pif1m2DivC_Watson = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/pif1m2DivergentCricktRNAGenesRawW-5000.txt',[],5000)
pif1m2DivW_Watson = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/pif1m2DivergentWatsontRNAGenesRawW-5000.txt',[],5000)
pif1m2DivW_Crick = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/pif1m2DivergentWatsontRNAGenesRawC-5000.txt',[],5000)

rrm3dpif1m2DivC_Crick = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/rrm3d_pif1m2DivergentCricktRNAGenesRawC-5000.txt',[],5000)
rrm3dpif1m2DivW_Watson = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/rrm3d_pif1m2DivergentWatsontRNAGenesRawW-5000.txt',[],5000)
rrm3dpif1m2DivC_Watson = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/rrm3d_pif1m2DivergentCricktRNAGenesRawW-5000.txt',[],5000)
rrm3dpif1m2DivW_Crick = wcr.tRNAGenesStrandDensityMeanForStrand('../TestData/rrm3d_pif1m2DivergentWatsontRNAGenesRawC-5000.txt',[],5000)

fig = plt.figure(1)
ax = fig.add_subplot(221)
ax.plot(range(-5000,5000),np.array(rrm3dConvW_Watson)-np.array(cdc9degConvW_Watson),label="rrm3dW_W")
ax.plot(range(-5000,5000),np.array(pif1m2ConvW_Watson)-np.array(cdc9degConvW_Watson),label="pif1m2W_W")
ax.plot(range(-5000,5000),np.array(rrm3dpif1m2ConvW_Watson)-np.array(cdc9degConvW_Watson),label="rrm3dpif1m2W_W")
plt.xlim(-5000,5000)
plt.axhline(xmax=5000,color='black')
plt.axvline(x=0,color='black') 
plt.xlabel('Bases')
plt.ylabel('Watson Density')
plt.legend(loc='best')
ax = fig.add_subplot(222)
ax.plot(range(-5000,5000),np.array(rrm3dConvW_Crick)-np.array(cdc9degConvW_Crick),label="rrm3dW_C")
ax.plot(range(-5000,5000),np.array(pif1m2ConvW_Crick)-np.array(cdc9degConvW_Crick),label="pif1m2W_C")
ax.plot(range(-5000,5000),np.array(rrm3dpif1m2ConvW_Crick)-np.array(cdc9degConvW_Crick),label="rrm3dpif1m2W_C")
plt.xlim(-5000,5000)
plt.axhline(xmax=5000,color='black')
plt.axvline(x=0,color='black') 
plt.xlabel('Bases')
plt.ylabel('Crick Density')
plt.legend(loc='best')

ax = fig.add_subplot(223)
ax.plot(range(-5000,5000),np.array(rrm3dDivC_Watson)-np.array(cdc9degDivC_Watson),label="rrm3dC_W")
ax.plot(range(-5000,5000),np.array(pif1m2DivC_Watson)-np.array(cdc9degDivC_Watson),label="pif1m2C_W")
ax.plot(range(-5000,5000),np.array(rrm3dpif1m2DivC_Watson)-np.array(cdc9degDivC_Watson),label="rrm3dpif1m2C_W")
plt.xlim(-5000,5000)
plt.axhline(xmax=5000,color='black')
plt.axvline(x=0,color='black') 
plt.xlabel('Bases')
plt.ylabel('Watson Density')
plt.legend(loc='best')
ax = fig.add_subplot(224)
ax.plot(range(-5000,5000),np.array(rrm3dDivC_Crick)-np.array(cdc9degDivC_Crick),label="rrm3dC_C")
ax.plot(range(-5000,5000),np.array(pif1m2DivC_Crick)-np.array(cdc9degDivC_Crick),label="pif1m2C_C")
ax.plot(range(-5000,5000),np.array(rrm3dpif1m2DivC_Crick)-np.array(cdc9degDivC_Crick),label="rrm3dpif1m2C_C")
plt.xlim(-5000,5000)
plt.axhline(xmax=5000,color='black')
plt.axvline(x=0,color='black') 
plt.xlabel('Bases')
plt.ylabel('Crick Density')
plt.legend(loc='best')
plt.tight_layout()

fig = plt.figure(1)
ax = fig.add_subplot(221)
ax.plot(range(-5000,5000),np.array(rrm3dDivW_Watson)-np.array(cdc9degDivW_Watson),label="rrm3dW_W")
ax.plot(range(-5000,5000),np.array(pif1m2DivW_Watson)-np.array(cdc9degDivW_Watson),label="pif1m2W_W")
ax.plot(range(-5000,5000),np.array(rrm3dpif1m2DivW_Watson)-np.array(cdc9degDivW_Watson),label="rrm3dpif1m2W_W")
plt.xlim(-5000,5000)
plt.axhline(xmax=5000,color='black')
plt.axvline(x=0,color='black') 
plt.xlabel('Bases')
plt.ylabel('Watson Density')
plt.legend(loc='best')
ax = fig.add_subplot(222)
ax.plot(range(-5000,5000),np.array(rrm3dDivW_Crick)-np.array(cdc9degDivW_Crick),label="rrm3dW_C")
ax.plot(range(-5000,5000),np.array(pif1m2DivW_Crick)-np.array(cdc9degDivW_Crick),label="pif1m2W_C")
ax.plot(range(-5000,5000),np.array(rrm3dpif1m2DivW_Crick)-np.array(cdc9degDivW_Crick),label="rrm3dpif1m2W_C")
plt.xlim(-5000,5000)
plt.axhline(xmax=5000,color='black')
plt.axvline(x=0,color='black') 
plt.xlabel('Bases')
plt.ylabel('Crick Density')
plt.legend(loc='best')

ax = fig.add_subplot(223)
ax.plot(range(-5000,5000),np.array(rrm3dConvC_Watson)-np.array(cdc9degConvC_Watson),label="rrm3dC_W")
ax.plot(range(-5000,5000),np.array(pif1m2ConvC_Watson)-np.array(cdc9degConvC_Watson),label="pif1m2C_W")
ax.plot(range(-5000,5000),np.array(rrm3dpif1m2ConvC_Watson)-np.array(cdc9degConvC_Watson),label="rrm3dpif1m2C_W")
plt.xlim(-5000,5000)
plt.axhline(xmax=5000,color='black')
plt.axvline(x=0,color='black') 
plt.xlabel('Bases')
plt.ylabel('Watson Density')
plt.legend(loc='best')
ax = fig.add_subplot(224)
ax.plot(range(-5000,5000),np.array(rrm3dConvC_Crick)-np.array(cdc9degConvC_Crick),label="rrm3dC_C")
ax.plot(range(-5000,5000),np.array(pif1m2ConvC_Crick)-np.array(cdc9degConvC_Crick),label="pif1m2C_C")
ax.plot(range(-5000,5000),np.array(rrm3dpif1m2ConvC_Crick)-np.array(cdc9degConvC_Crick),label="rrm3dpif1m2C_C")
plt.xlim(-5000,5000)
plt.axhline(xmax=5000,color='black')
plt.axvline(x=0,color='black') 
plt.xlabel('Bases')
plt.ylabel('Crick Density')
plt.legend(loc='best')
plt.tight_layout()

alanine_cdc9deg_watson = wcr.tRNAGenesWCRatioMeanForStrand('../TestData/cdc9degtRNAGenesRawWCRatio-Alanine-5000.txt','W',[],5000)
alanine_cdc9deg_crick = wcr.tRNAGenesWCRatioMeanForStrand('../TestData/cdc9degtRNAGenesRawWCRatio-Alanine-5000.txt','C',[],5000)
alanine_rrm3d_watson = wcr.tRNAGenesWCRatioMeanForStrand('../TestData/rrm3dtRNAGenesRawWCRatio-Alanine-5000.txt','W',[],5000)
alanine_rrm3d_crick = wcr.tRNAGenesWCRatioMeanForStrand('../TestData/rrm3dtRNAGenesRawWCRatio-Alanine-5000.txt','C',[],5000)
alanine_pif1m2_watson = wcr.tRNAGenesWCRatioMeanForStrand('../TestData/pif1m2tRNAGenesRawWCRatio-Alanine-5000.txt','W',[],5000)
alanine_pif1m2_crick = wcr.tRNAGenesWCRatioMeanForStrand('../TestData/pif1m2tRNAGenesRawWCRatio-Alanine-5000.txt','C',[],5000)
alanine_rrm3dpif1m2_watson = wcr.tRNAGenesWCRatioMeanForStrand('../TestData/rrm3d_pif1m2tRNAGenesRawWCRatio-Alanine-5000.txt','W',[],5000)
alanine_rrm3dpif1m2_crick = wcr.tRNAGenesWCRatioMeanForStrand('../TestData/rrm3d_pif1m2tRNAGenesRawWCRatio-Alanine-5000.txt','C',[],5000)

fig = plt.figure(1)
ax = fig.add_subplot(121)
ax.plot(range(-5000,5000),alanine_cdc9deg_watson,label='cdc9deg')
ax.plot(range(-5000,5000),alanine_rrm3d_watson,label='rrm3d')
ax.plot(range(-5000,5000),alanine_pif1m2_watson,label='pif1m2')
ax.plot(range(-5000,5000),alanine_rrm3dpif1m2_watson,label='rrm3dpif1m2')
plt.legend(loc='best')
plt.xlabel('Bases')
plt.xlim(-5000,5000)
plt.axhline(xmax=5000,color='black')
plt.axvline(x=0,color='black') 
plt.ylabel('Log2 Watson-Crick Ratio')
plt.title('Watson-Crick Ratio Around Alanine tRNA Watson Genes')

ax = fig.add_subplot(122)
ax.plot(range(-5000,5000),alanine_cdc9deg_crick,label='cdc9deg')
ax.plot(range(-5000,5000),alanine_rrm3d_crick,label='rrm3d')
ax.plot(range(-5000,5000),alanine_pif1m2_crick,label='pif1m2')
ax.plot(range(-5000,5000),alanine_rrm3dpif1m2_crick,label='rrm3dpif1m2')
plt.legend(loc='best')
plt.xlabel('Bases')
plt.xlim(-5000,5000)
plt.axhline(xmax=5000,color='black')
plt.axvline(x=0,color='black') 
plt.ylabel('Log2 Watson-Crick Ratio')
plt.title('Watson-Crick Ratio Around Alanine tRNA Crick Genes')

"""