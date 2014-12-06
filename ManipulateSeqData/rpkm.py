__author__ = 'jayashreekumar'

#######################################################################################################################
# IMPORT
import numpy as np
import matplotlib.pyplot as plt
import manipulatesgrfiles as sgr
#from math import ceil
import originscomparison as og
from operator import itemgetter
#######################################################################################################################
# FUNCTIONS

# This is a dictionary of yeast chromosome sizes
yeastsize = {'1':230218,'2':813184,'3':316620,'4':1531933,'5':576874,'6':270161,'7':1090940,'8':562643,'9':439888,'10':745751,'11':666816,'12':1078177,'13':924431,'14':784333,'15':1091291,'16':948066,'Total':12157105}

# This function will calculate 'RPKM' for each base in an array of bases for a specific chromosome. The arguments
# needed is a sgrlist for the chromosome and the chromosome number which the array of reads is for (chromosome).
# The function returns an array of RPKMs calculated for each base position and will be the size of the chromosome
def calcRPKMforOneChrom(sgrlistforchrom):
    # Grab the number of reads for every base in sgforchrom which is the third column (index 2) and store in a numpy array
    basereads = np.array([i[2] for i in sgrlistforchrom])    
    # Total number of reads for the chromosome
    totalreads = sum(basereads)
    # Access the chromosome length from the dictionary yeastsize by using key chromosome
    #chrmlen = yeastsize[chromosome]
    #TotalChrmLen = yeastsize['Total']
    # Calculate the RPK for reads
    #rpk = basereads/(float(TotalChrmLen)/1000000)
    # Calculate the RPKM for reads
    #rpkm = rpk/(float(totalreads)/1000000)
    newrpkm = basereads/(float(totalreads)/1000000)
    #return rpkm
    return newrpkm

# This function takes in an sgrfile (sgrfile) which contains the reads per base for every chromosome. It returns a
# dictionary where key is chromosome and value is the array that stores rpkms of every position for the chromosome
def calcRPKMforAllChroms(sgrfileW,sgrfileC,writefile=''):

    # Empty dictionary to store rkpm arrays for every chromosome
    rpkmforallchrm = {}

    # Go through every chromosome 1 through 16
    for c in range(1,17):
        print c
        # Convert c into a string
        chromosome = str(c)
        # Convert sgr file to a sgrlist for a chromosome
        Wsgrlistforchrom = sgr.getSgrList(sgrfileW,chromosome)
        Csgrlistforchrom = sgr.getSgrList(sgrfileC,chromosome)
        # Store RPKM array as value with key chromosome in dictionary rpkmforallchrm
        rpkmforallchrm[chromosome] = calcRPKMforOneChrom(Wsgrlistforchrom) + calcRPKMforOneChrom(Csgrlistforchrom)
    
    if writefile:
        with open(writefile,'w') as fw:
            for c in range(1,17):
                chrm = str(c)
                base = 1
                for i in rpkmforallchrm[chrm]:
                    fw.write(chrm+'\t'+str(base)+'\t'+str(i)+'\n')
                    base += 1
    else:
        return rpkmforallchrm
 

# This function plots RPKM for 2 sgrfiles one for the watson strand (sgrfileW) and one for the crick strand (sgrfileC).
# It also takes in a chromosome value (chromosome) to know which chromosome to plot the RPKMs for. In the top subplot
# it will be the Watson strand and the bottom subplot will be the Crick strand. The strain parm will indicate which strain
# it is for, and displays it in the title of the plot. The parameter sbnum allows you to input the subplot number.
def plotRPKM(sgrlistforchromW,sgrlistforchromC,strain,chromosome,sbnum):

    # Store chromosome length in chrmlen accessing value from dictionary using key chromosome
    chrmlen = yeastsize[chromosome]
        
    # RPKM Array for Watson strand
    rpkmWatson = calcRPKMforOneChrom(sgrlistforchromW)
 
    # RPKM Array for Crick strand. Multiply by -1 to flip it over the x-axis
    rpkmCrick = (calcRPKMforOneChrom(sgrlistforchromC))*-1

    # Create figure
    fig = plt.figure(1)
    # Add subplot on figure
    ax = fig.add_subplot(sbnum)
    # Plot RPKM for both Watson and Crick strand on the same plot
    ax.plot(range(1,chrmlen+1),rpkmWatson,range(1,chrmlen+1),rpkmCrick)
    # Set limits for x and y axis, x axis limits is the number of bases for that chromosome, y axis max is based on max
    # RPKM in watson strand and min is based on min RPKM coverage in crick strand (since negative)
    plt.xlim(xmin=0,xmax=chrmlen)
    plt.ylim(ymin=min(rpkmCrick),ymax=max(rpkmWatson))
    # Y axis ticks are going to be absolute value, since depth coverage for crick strand is also positive but we just
    # want it below the x axis
    ax.set_yticklabels([str(abs(x)) for x in ax.get_yticks()])
    # We want the x-axis line to be black and obvious
    plt.axhline(xmax=chrmlen,color='black')

    # Set x axis label and y axis label and title of plot
    plt.xlabel('Bases')
    plt.ylabel('RKPM')
    plt.title('RPKM for Chromsome ' + chromosome + ' for ' + strain)

# Function finds max RPKM for every 1000 bases and then spits out base for max RPKM and the RPKM value for that base for a particular chromosome
# Across the entire chromosome, Not Local Peaks
def maxRPKMStrandSpecific(sgrlist,window,chromosome,thresholdRPKM,minPeaks):
    # RPKM Array for sgrfile
    rpkm = calcRPKMforOneChrom(sgrlist,chromosome)

    # If the max RPKM is greater than the threshold RPKM than set threshold RPKM to be the max RPKM
    if (max(rpkm) >= thresholdRPKM):
        thresholdRPKM = max(rpkm)
    
    # Dictionary to store the max values
    maxvals = {}
    
    # Keep adjusting the threshold RPKM until you get equal to or more than the minimum number of peaks specified
    while (len(maxvals) < minPeaks):
        # Dictionary to store the max values
        maxvals = {}
        for i in range(0,len(rpkm),window):
            #grab the base at which the max rpkm is occuring
            ind = np.argmax(rpkm[i:i+window]) + i
            if rpkm[ind] >= thresholdRPKM:
                maxvals[ind] = rpkm[ind]

        thresholdRPKM = float(thresholdRPKM)/2

    return maxvals

# Function finds max RPKM for every 1000 bases and then spits out base for max RPKM and the RPKM value for that base for a particular chromosome combines RPKMs from both strands
# Across the entire chromosome, Not Local Peaks
def maxRPKM(sgrlistW,sgrlistC,window,chromosome,thresholdRPKM,minPeaks):
    # RPKM Array for sgrfile
    rpkm = calcRPKMforOneChrom(sgrlistW,chromosome) + calcRPKMforOneChrom(sgrlistC,chromosome)

    if (max(rpkm) >= thresholdRPKM):
        thresholdRPKM = max(rpkm)
    
    # Dictionary to store the max values
    maxvals = {}

    while (len(maxvals)< minPeaks):
        # Empty Dictionary to store the max values
        maxvals = {}
        for i in range(0,len(rpkm),window):
            #print rpkm[i:i+window]
            ind = np.argmax(rpkm[i:i+window]) + i
            if rpkm[ind] >= thresholdRPKM:
                maxvals[ind] = rpkm[ind]

        thresholdRPKM = float(thresholdRPKM)/2

    return maxvals

def maxRPKMLocal(strandoption,sgrliststrand,chromosome,tpmain,tplocal,maxpeaks,originOption=True,sgrlistSecstrand=[],oemfile=''):
    
    # RPKM Array for sgrfile
    if (strandoption == 'Strand'):    
        rpkm = calcRPKMforOneChrom(sgrliststrand)
    else:
        rpkm = calcRPKMforOneChrom(sgrliststrand) + calcRPKMforOneChrom(sgrlistSecstrand)
    
    thresholdrpkm = np.percentile(rpkm,tpmain)
  
    if originOption:
        boundaries = sorted([i[0] for i in og.getPredictedOriginsForOEMList(oemfile,chromosome) if i[1] >= 0.1])
    else:
        boundaries = range(10000,len(rpkm),10000)
    
    # Dictionary to store the max values
    maxvals = []
    
    #Go through every predicted origin, we are looking for peaks between predicted origins
    for i in range(len(boundaries)):
        # For the first predicted origin, we want to take rpkms from 0 to th first predicted origins
        if i==0:
            rpkmlocal = rpkm[0:boundaries[i]]
            baselocs = range(0,boundaries[i])
            rpkmwithBaseslocal = [[baselocs[j],rpkmlocal[j]] for j in range(0,boundaries[i])]
        # For the rest of the predicted origins, we want to take rpkms from the previous predicted origin to the current predicted origin
        else:
            rpkmlocal = rpkm[boundaries[i-1]:boundaries[i]]
            baselocs = range(boundaries[i-1],boundaries[i])
            rpkmwithBaseslocal = [[baselocs[j],rpkmlocal[j]] for j in range(len(rpkmlocal))]
        
        # Only find peaks if the rpkm at percentile tp is greater than the thresholdrpkm
        if len([i for i in rpkmwithBaseslocal if i[1] >= thresholdrpkm])!=0:
            newrpkmWithBaseslocal = rpkmwithBaseslocal
            localthresholdrpkm = thresholdrpkm
            # Continue to find peaks within the local rpkms if the number of peaks found is greater than maxpeaks.
            # The localthresholdrpkm is reduced each time to the rpkm at percentile tp of the new smaller group
            maxrpkmWithBaseslocal = []
            while(len(maxrpkmWithBaseslocal)<= 0 or len(maxrpkmWithBaseslocal)>maxpeaks):
                maxrpkmWithBaseslocal = []
                rpkmsplit = splitListbasedOnParm(newrpkmWithBaseslocal,localthresholdrpkm)
                newrpkmWithBaseslocal = []
                for j in rpkmsplit:
                    if len([i for i in j if i[1] >= localthresholdrpkm])!=0:
                        maxrpkmWithBaseslocal.append(max(j,key=itemgetter(1)))
                        newrpkmWithBaseslocal.extend(j)
                newrpkmlocal = [r[1] for r in newrpkmWithBaseslocal]
                localthresholdrpkm = np.percentile(newrpkmlocal,tplocal)
               
            
            maxvals.extend(maxrpkmWithBaseslocal)             
    realmaxvals = []      
   
    for i in range(len(maxvals)):
        if i==0:
            if maxvals[i][0]+100 < maxvals[i+1][0]:
                realmaxvals.append(maxvals[i])
            else:
                realmaxvals.append(max([maxvals[i],maxvals[i+1]],key=itemgetter(1)))
        else:
            if maxvals[i-1][0]+100 < maxvals[i][0]:
                realmaxvals.append(maxvals[i])
            else:
                if maxvals[i-1] not in realmaxvals and maxvals[i] not in realmaxvals :
                    valuetoputin = max([maxvals[i-1],maxvals[i]],key=itemgetter(1))
                    if maxvals[len(maxvals)-1][0] + 100 < valuetoputin[0]:
                        realmaxvals.append(valuetoputin)
                else:
                    valuetoputin = max([maxvals[i-1],maxvals[i]],key=itemgetter(1))
                    if valuetoputin == maxvals[i-1] and maxvals[i] in realmaxvals:
                        realmaxvals.remove(maxvals[i])
                        if maxvals[len(maxvals)-1][0] + 100 < valuetoputin[0]:
                            realmaxvals.append(valuetoputin)
                    elif valuetoputin == maxvals[i] and maxvals[i-1] in realmaxvals:
                        realmaxvals.remove(maxvals[i-1])
                        if maxvals[len(maxvals)-1][0] + 100 < valuetoputin[0]:
                            realmaxvals.append(valuetoputin)
    #print realmaxvals
    return realmaxvals

def splitListbasedOnParm(list2split,splitval):
    # Empty list to store all the split lists which are dictionaries
    splitlist = []
    # Empty dictionary to store each individual as long as sign does not change between numbers -> used a dictionary
    # because we needed to access the base for each corresponding OEM easily
    splitlistind = []
    # Go through every index in list2split
    for i in range(len(list2split)):
        if i!=len(list2split)-1:
            # If the current value at index and the value at index + 1 is greater than or equal to split val
            if (list2split[i][1] >= splitval and list2split[i+1][1] >=splitval) or (list2split[i][1] < splitval and list2split[i+1][1] < splitval):
                # Let index be the key of dictionary and the value at that index is stored as the paired value
                splitlistind.append(list2split[i])
            else:
                # Add this index and value as key-value pairs to the same dictionary as previously because it would
                # have been greater than split val like the previous value, and the next value is not greater than the splitval
                splitlistind.append(list2split[i])
                #print(splitlistind)
                splitlist.append(splitlistind)
                splitlistind = []
        else:
            splitlistind.append(list2split[i])
            splitlist.append(splitlistind)
    # return this big list of dictionaries
    return splitlist

   
def maxRPKMLocalSgrFile(strandoption,sgrfile1,chromosome,tpmain,tplocal,maxpeaks,originOption=True,sgrfile2='',oemfile=''):
    
    if strandoption == 'Both':
        sgrlist1 = sgr.getSgrList(sgrfile1,chromosome)
        sgrlist2 = sgr.getSgrList(sgrfile2,chromosome)
        return maxRPKMLocal(strandoption,sgrlist1,chromosome,tpmain,tplocal,maxpeaks,originOption,sgrlist2,oemfile)
    else:
        sgrlist1 = sgr.getSgrList(sgrfile1,chromosome)
        return maxRPKMLocal(strandoption,sgrlist1,chromosome,tpmain,tplocal,maxpeaks,originOption,oemfile)
        

def plotmaxRPKMLocalSgrFile(strandoption,sgrfile1,chromosome,tpmain,tplocal,maxpeaks,originOption=True,sgrfile2='',oemfile=''):
    
    peaks = maxRPKMLocalSgrFile(strandoption,sgrfile1,chromosome,tpmain,tplocal,maxpeaks,originOption,sgrfile2,oemfile)

    
    if strandoption == 'Both':
        sgrlist1 = sgr.getSgrList(sgrfile1,chromosome)
        sgrlist2 = sgr.getSgrList(sgrfile2,chromosome)
        totalrpkm = calcRPKMforOneChrom(sgrlist1) + calcRPKMforOneChrom(sgrlist2)
    else:
        sgrlist1 = sgr.getSgrList(sgrfile1,chromosome)
        totalrpkm = calcRPKMforOneChrom(sgrlist1)
        
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.plot(range(len(totalrpkm)),totalrpkm,'green')
    for i in peaks:
        plt.vlines(i[0],0,max(totalrpkm)+1,'red')
    plt.hlines(np.percentile(totalrpkm,tpmain),0,len(totalrpkm),'blue')