__author__ = 'jayashreekumar'

"""
The functions in this script can be used to compare origins for wild type and knockout strains. Many of the functions in this script require an OEM file

"""
################################################################################################
# IMPORT

import oem
import matplotlib.pyplot as plt
from operator import itemgetter
import numpy as np

################################################################################################
# FUNCTIONS

########### THESE BELOW FUNCTIONS DEAL WITH PREDICTING ORIGINS BASED ON OEM ##############################

# This function splits up a list of numbers based on whether they are consecutively positive and negative. For example
# list of numbers like [0,0,2,3,5,-6,-7,0,8,9,-12,-11,0,0], it will split up as [2,3,5], [-6,-7],[8,9],[-12,-11].
def splitListbasedOnSign(list2split):
    # Empty list to store all the split lists which are dictionaries
    splitlist = []
    # Empty dictionary to store each individual as long as sign does not change between numbers -> used a dictionary
    # because we needed to access the base for each corresponding OEM easily
    splitlistind = {}
    # Go through every index in list2split
    for i in range(len(list2split)-1):
        # If the signs of the current value at index and the value at index + 1 is equal and both values are not zero
        if cmp(list2split[i][1],0) == cmp(list2split[i+1][1],0) and cmp(list2split[i][1],0) != 0 and cmp(list2split[i+1][1],0) != 0:
            # Let index be the key of dictionary and the value at that index is stored as the paired value
            splitlistind[list2split[i][0]] = list2split[i][1]
            if i==(len(list2split)-2):
                splitlistind[list2split[i+1][0]] = list2split[i+1][1]
                splitlist.append(splitlistind)
                splitlistind = {}
                
        # If the signs of the current value at index and the value at index + 1 are not equal or either values are zero
        else:
            # If value at current index is not zero
            
            if list2split[i][1] != 0:
                # Add this index and value as key-value pairs to the same dictionary as previously because it would
                # have been equal in sign to previous value, it's not just equal in sign to the next value
                splitlistind[list2split[i][0]] = list2split[i][1]
            # If length of current dictionary is not zero
            if len(splitlistind)!=0:
                # Add this dictionary which contains values of same sign to the big list and reset the dictionary to be
                # empty
                splitlist.append(splitlistind)
                splitlistind = {}
    # return this big list of dictionaries
    return splitlist

# This function returns the key to the maximum value in the dictionary
def findmaxInDict(max2finddict):
    return max(max2finddict.iteritems(), key=itemgetter(1))
# This function returns the key to the minimum value in the dictionary
def findminInDict(min2finddict):
    return min(min2finddict.iteritems(), key=itemgetter(1))
# Get a list of bases which have max OEM in this dictionary as long as OEM is positive
def basesWmaxoem(splitlist):
    # Empty list to store all the bases with max positive OEM values
    maxoembases = []
    # Go through every dictionary in the split list
    for dictoem in splitlist:
        # If first value in dictionary is positive, means all the values in that dictionary are positive since they
        # should all have the same sign
        if dictoem.values()[0] > 0:
            # Multiply each the base with the max OEM value by perbasebox since OEMs were calculated for every perbasebox value
            maxoembases.append(findmaxInDict(dictoem))

    return maxoembases

# This function returns all predicted origins which come within a compWindow of the confirmed origins as well as those
# that don't come within that compWindow
def predictedComparedToOriginsList(predictedOrigins,oriDBorigins,compWindow):
    originsWithinWindow = []
    # For every predicted origin
    for porigin in predictedOrigins:
        # Go through every confirmed origin
        for rorigin in oriDBorigins:
            # If the predicted origin is within a window from either side of confirmed origins coordinates, add to the
            # originsWithinWindow list and break out of for loop for oriDBorigins
            if rorigin[0]-compWindow <= porigin[0] <= rorigin[1]+compWindow:
                originsWithinWindow.append([porigin,rorigin])
                break

    # Grab all the predicted origins that do lie within the compwindow of origins in oriDBorigins
    justpredictedoriginsWithin = [i[0] for i in originsWithinWindow]

    originsNotWithinWindow = [i for i in predictedOrigins if i not in justpredictedoriginsWithin]

    return {'Within':originsWithinWindow,'Not':originsNotWithinWindow}

# This function returns the predicted origins for the oemfile
def getPredictedOriginsForOEMList(oemfile,chromosome):
    
    #Grab the OEM for sgrcomblist
    with open(oemfile) as f:
        oemc = [[int(line.strip().split('\t')[1]),float(line.strip().split('\t')[2])] for line in f.readlines() if line.strip().split('\t')[0]==chromosome]
    
    # Split up the OEM values based on sign, so consecutive positive numbers and negative numbers get grouped
    s = splitListbasedOnSign(oemc)
    # Grab bases which have the max positive OEM value per list
    maxBases = basesWmaxoem(s)

    return maxBases

########### THESE ABOVE FUNCTIONS DEAL WITH PREDICTING ORIGINS BASED ON OEM ##############################

######## THESE BELOW FUNCTIONS DEAL WITH COMPARING PREDICTED ORIGINS TO ORIGINS IN ORIDB BASED ON OEM ###########

# This function gets number of predicted origins in an oem file that are within originwindow bases of confirmed,likely or
# dubious origins. Origins not fitting any of these categories are also recorded
def getnumComparedOriginsWTMUT(oemfile,originwindow):
    # Go through every chromosome in yeast

    # Intialize variables to store total num of predicted matching origins to be zero
    totalConfirmed = 0
    totalLikely = 0
    totalDubious = 0
    totalNoClass = 0

    for i in range(1,17):
        # Need to convert number into string to access it in all lists
        chromosome = str(i)

        # Return predicted origins
        po = getPredictedOriginsForOEMList(oemfile,chromosome)

        # This is a list of confirmed origins from OriDB
        confirmedOrigins = [[i[0],i[1]] for i in oem.readOriginsfile('S-cerevisiae_origins.txt',chromosome,'Confirmed')]
        # This is a list of likely origins from OriDB
        likelyOrigins = [[i[0],i[1]] for i in oem.readOriginsfile('S-cerevisiae_origins.txt',chromosome,'Likely')]
        # This is a list of dubious origins from OriDB
        dubiousOrigins = [[i[0],i[1]] for i in oem.readOriginsfile('S-cerevisiae_origins.txt',chromosome,'Dubious')]

        # Number of origins that are confirmed,likely,dubious and None of those among predicted origins for chromosome
        originsConfirmedInAndOut = predictedComparedToOriginsList(po,confirmedOrigins,originwindow)
        originsLikelyInAndOut = predictedComparedToOriginsList(originsConfirmedInAndOut['Not'],likelyOrigins,originwindow)
        originsDubiousInAndOut = predictedComparedToOriginsList(originsLikelyInAndOut['Not'],dubiousOrigins,originwindow)
        originsNone = originsDubiousInAndOut['Not']

        # Add it to the total number of confirmed, likely, dubious or none of those among predicted origins
        totalConfirmed = totalConfirmed + len(originsConfirmedInAndOut['Within'])
        totalLikely = totalLikely + len(originsLikelyInAndOut['Within'])
        totalDubious = totalDubious + len(originsDubiousInAndOut['Within'])
        totalNoClass = totalNoClass + len(originsNone)



    return {'C':totalConfirmed,'L':totalLikely,'D':totalDubious,'N':totalNoClass}

######## THESE ABOVE FUNCTIONS DEAL WITH COMPARING PREDICTED ORIGINS TO ORIGINS IN ORIDB BASED ON OEM ###########
######## THE FUNCTIONS BELOW DEAL WITH COMPARING PREDICTED ORIGINS OF ONE SPECIES TO ANOTHER ####################

# This function takes in a WT oem file, MUT oem file, a chromosome,spits out a dictionary of lists of base coordinates for predicted origins
#, which need to be above the oem limit.
def predictedOriginsForWTandMUT(oemfileWT,oemfileMUT,chromosome,oemlimit):

    # Get predicted origins for WT oem file
    poWT = getPredictedOriginsForOEMList(oemfileWT,chromosome)

    # Get predicted origins for MUT oem file
    poMUT = getPredictedOriginsForOEMList(oemfileMUT,chromosome)

    # Grab the predicted origins with oem for WT where the oem is greater than the oemlimit
    poWTwithoem = [i for i in poWT if i[1]>oemlimit]

    # Grab the predicted origins with oem for MUT where the oem is greater than the oemlimit
    poMUTwithoem = [i for i in poMUT if i[1]>oemlimit]

    return {'WTbases':[i[0] for i in poWTwithoem],'MUTbases':[i[0] for i in poMUTwithoem],'WTwOEM':poWTwithoem,'MUTwOEM':poMUTwithoem}

# This function takes in a WT predicted origins list, MUT predicted origins list and matchingbasewindow and spits out
# a list of matching origins of base coordinates that come within the matchingbasewindow. The first predicted origin
# within each matching pair of origins is the WT base coordinate and the second predicted origin is MUT base coordinate.
def matchingOriginsBasedOnWindow(predictedoriginsWT,predictedoriginsMUT,matchingbasewindow):
    
    # Store predicted origins for MUT in a temporary variable
    # we need a different pointer to a different list with the same contents hence the use of [:]
    temp_poMUT = predictedoriginsMUT[:]
    
    # matchdbasecoord variable is an empty list to store matched predicted origins
    
    matchedbasecoord = []

    # Go through each base coordinate for predicted origins in WT
    for baseWT in predictedoriginsWT:
        # Store predicted origins for MUT in a temporary variable
        # we need a different pointer to a different list with the same contents hence the use of [:]
        
        # Go through each base coordinate for predicted origins in MUT to compare to base coordinate in WT
        for baseMUT in temp_poMUT:
            # If baseWT is between the matchingbasewindow of baseMUT, proceed
            if baseMUT-matchingbasewindow <= baseWT and baseWT <= baseMUT+matchingbasewindow:
                # Add the matched origins to the matchedbasecoord variable
                matchedbasecoord.append([baseWT,baseMUT])
                # Remove the matched predicted origin MUT from the temporary variable
                temp_poMUT.remove(baseMUT)
                # Break out the loop for baseWT once we have found a matching predicted origin
                break

    return matchedbasecoord

# This function takes in a WT predicted origins plus oem list, MUT predicted origins plus oem list and matchingbasewindow and spits out
# a list of non matching origins of base coordinates that come within the matchingbasewindow.
def nonmatchingOriginsBasedOnWindow(predictedoriginsplusoemWT,predictedoriginsplusoemMUT,matchingbasewindow):

    # Grab the base coordinates for WT where the oem is greater than the oemlimit
    poWTbasecoord = [i[0] for i in predictedoriginsplusoemWT]

    # Grab the base coordinates for MUT where the oem is greater than the oemlimit
    poMUTbasecoord = [i[0] for i in predictedoriginsplusoemMUT]

    # These are the matching origins for WT and MUT
    matchingorigins = matchingOriginsBasedOnWindow(poWTbasecoord,poMUTbasecoord,matchingbasewindow)

    # WT predicted origins that have matches in MUT within the matchingbasewindow
    moWT = [i[0] for i in matchingorigins]

    # MUT predicted origins that have matches in WT within the matchingbasewindow
    moMUT = [i[1] for i in matchingorigins]

    # WT predicted origins that do not have matches in MUT within the matchingbasewindow and grab the oems as well
    nonmoWT = set(poWTbasecoord) - set(moWT)
    nonmoWoemWT = [([i for i in predictedoriginsplusoemWT if i[0]==base]) for base in nonmoWT]

    # MUT predicted origins that do not have matches in WT within the matchingbasewindow and grab the oems as well
    nonmoMUT = set(poMUTbasecoord) - set(moMUT)
    nonmoWoemMUT = [([i for i in predictedoriginsplusoemMUT if i[0]==base]) for base in nonmoMUT]

    #return the non matching predicted origins for WT and MUT
    return {'WT':nonmoWoemWT,'MUT':nonmoWoemMUT}

# Plot fragment density of non matching origins for WT and MUT given a base and sgrcomblists for WT and MUT
def plotFragDensityWTMUTatspecificbase(strain,base,sgrcomblistWT,sgrcomblistMUT,totalNumNMO,indextoplot):
    # If the base-10000 is less than 0, then we only care to grab fragment density starting from 0 to base +10000    
    if base-100000<0:
        print 'Here 1'
        fgdWatsonForbaseWT = [j[2] for j in sgrcomblistWT[0:base+100000]]
        fgdCrickForbaseWT = [-j[3] for j in sgrcomblistWT[0:base+100000]]
        basesWT = [j[1] for j in sgrcomblistWT[0:base+100000]]
        fgdWatsonForbaseMUT = [j[2] for j in sgrcomblistMUT[0:base+100000]]
        fgdCrickForbaseMUT = [-j[3] for j in sgrcomblistMUT[0:base+100000]]
        basesMUT = [j[1] for j in sgrcomblistMUT[0:base+100000]]
    # Else if the base+10000 is greater than the number of bases available in the sgrcomblists, then we only care to grab fragement density from
    # base-10000 to te max number of bases
    elif base+100000>(len(sgrcomblistWT)-1):
        print 'Here 2'
        fgdWatsonForbaseWT = [j[2] for j in sgrcomblistWT[base-100000:len(sgrcomblistWT)]]
        fgdCrickForbaseWT = [-j[3] for j in sgrcomblistWT[base-100000:len(sgrcomblistWT)]]
        basesWT = [j[1] for j in sgrcomblistWT[base-100000:len(sgrcomblistWT)]]
        fgdWatsonForbaseMUT = [j[2] for j in sgrcomblistMUT[base-100000:len(sgrcomblistMUT)]]
        fgdCrickForbaseMUT = [-j[3] for j in sgrcomblistMUT[base-100000:len(sgrcomblistMUT)]]
        basesMUT = [j[1] for j in sgrcomblistMUT[base-100000:len(sgrcomblistMUT)]]
    # Else if none of those conditions are met, we care to grab fragement density from base-10000 to base + 10000   
    else:
        print 'Here 3'
        fgdWatsonForbaseWT = [j[2] for j in sgrcomblistWT[base-100000:base+100000]]
        fgdCrickForbaseWT = [-j[3] for j in sgrcomblistWT[base-100000:base+100000]]
        basesWT = [j[1] for j in sgrcomblistWT[base-100000:base+100000]]
        fgdWatsonForbaseMUT = [j[2] for j in sgrcomblistMUT[base-100000:base+100000]]
        fgdCrickForbaseMUT = [-j[3] for j in sgrcomblistMUT[base-100000:base+100000]]
        basesMUT = [j[1] for j in sgrcomblistMUT[base-100000:base+100000]]
        
        
    fig = plt.figure(1)
    ax = fig.add_subplot(totalNumNMO,2,indextoplot)
    print (len(basesWT))
    print (len(fgdWatsonForbaseWT))
    print (len(fgdCrickForbaseWT))
    ax.plot(basesWT,fgdWatsonForbaseWT,basesWT,fgdCrickForbaseWT)
    print 'Here 4'
    plt.axhline(xmax=len(basesWT),color='black')
    plt.vlines(base,min(fgdCrickForbaseWT),max(fgdWatsonForbaseWT),color='red')
    ax.set_xticklabels([str(x/1000) for x in ax.get_xticks()])
    ax.set_yticklabels([str(abs(x)) for x in ax.get_yticks()])
    plt.title('Fragment Density WT - ' + strain + ' Base ' + str(base))
    ax = fig.add_subplot(totalNumNMO,2,indextoplot+1)
    ax.plot(basesMUT,fgdWatsonForbaseMUT,basesMUT,fgdCrickForbaseMUT)
    plt.axhline(xmax=len(basesMUT),color='black')
    plt.vlines(base,min(fgdCrickForbaseMUT),max(fgdWatsonForbaseMUT),color='red')
    ax.set_xticklabels([str(x/1000) for x in ax.get_xticks()])
    ax.set_yticklabels([str(abs(x)) for x in ax.get_yticks()])
    plt.title('Fragment Density MUT - ' + strain + ' Base ' + str(base))
    indextoplot=indextoplot+2

######## THE FUNCTIONS ABOVE DEAL WITH COMPARING PREDICTED ORIGINS OF ONE SPECIES TO ANOTHER ####################


def writePredictedOrigins(oemfile,predictedoriginsfile):
    # Grab just the OEM values
    with open(predictedoriginsfile,'w') as fw:
        for c in range(1,17):
            chromosome = str(c)
            with open(oemfile) as f:
                oemc = [[int(line.strip().split('\t')[1]),float(line.strip().split('\t')[2])] for line in f.readlines() if line.strip().split('\t')[0]==chromosome]
            #oemb = [float(i[2]) for i in oem if i[0]==chromosome]
            # Split up the OEM values based on sign, so consecutive positive numbers and negative numbers get grouped
            s = splitListbasedOnSign(oemc)
            # Grab bases which have the max positive OEM value per list with their oem value
            maxBases = basesWmaxoem(s)
            for i in maxBases:
                fw.write(str(c) + '\t' + str(i[0]) + '\t' + str(i[1]) + '\n')
                
                
def leastEfficientOrigins(predictedoriginsfile,leastefficientoriginsfile,leastpercent):
    with open(predictedoriginsfile) as f:
        porigins = [[line.strip().split('\t')[0],line.strip().split('\t')[1],float(line.strip().split('\t')[2])] for line in f]
        
    efficiencies = np.array([i[2] for i in porigins])
    
    leastefficientorigins = [i for i in porigins if i[2]<np.percentile(efficiencies,leastpercent)]
    
    leastefficientorigins.sort(key=lambda x: x[2])
    
    with open(leastefficientoriginsfile,'w') as fw:
        for i in leastefficientorigins:
            fw.write(i[0]+'\t'+i[1]+'\t'+str(i[2])+'\n')
    
                
def mostEfficientOrigins(predictedoriginsfile,mostefficientoriginsfile,mostpercent):
    with open(predictedoriginsfile) as f:
        porigins = [[line.strip().split('\t')[0],line.strip().split('\t')[1],float(line.strip().split('\t')[2])] for line in f]
        
    efficiencies = np.array([i[2] for i in porigins])
    
    mostefficientorigins = [i for i in porigins if i[2]>np.percentile(efficiencies,mostpercent)]
    
    mostefficientorigins.sort(key=lambda x: x[2],reverse=True)
    
    with open(mostefficientoriginsfile,'w') as fw:
        for i in mostefficientorigins:
            fw.write(i[0]+'\t'+i[1]+'\t'+str(i[2])+'\n')
