# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 21:15:07 2014

@author: jayashreekumar
"""

#######################################################################################################################
# IMPORT
import originscomparison as oc
import matplotlib.pyplot as plt
import numpy as np 
#######################################################################################################################
# FUNCTIONS

# This function takes in predicted origins which have bases and oems for a chromosome and a replication time file, chromosome and assigns replication times 
# to the predicted origins. It assigns the replication times by looking if the origin falls within a 500b window of the base in the replication time file and then 
# assigns the replication time for that base to the origin.
def assignReplicationTimesForPredictedOrigins(reptime,predictedoriginsforchrom,chromosome):
    
    timesForOrigins = []    
    
    # Grab replication times for the chromosome from reptime file and store in reptimesforchrom    
    with open(reptime) as f:    
        reptimesforchrom = [rt.strip().split('\t') for rt in f.readlines() if rt.strip().split('\t')[0]==chromosome]
    
    # Assign a replication time for predicted origins based on a 500base window
    for i in predictedoriginsforchrom:
        reptimefori = [int(time[2]) for time in reptimesforchrom if int(time[1])-250<=i[0] and i[0]<=int(time[1])+250]
        if len(reptimefori)!=0:        
            timesForOrigins.append([i[0],i[1],reptimefori[0]])
        
    return timesForOrigins

# Separate bases with late replication times and early replication times. The median of replication times is 30, so anything that falls below 30 is early 
# replicating and anything that falls above 30is late replicating.   
def earlyLateReplicatedOrigins(originswithreplicationtimes):
    
    # Assign predicted origins to earlyreplicatedorigins if replication time is less than or equal to 30
    earlyreplicatedorigins = [i for i in originswithreplicationtimes if i[2]<=30]
    
    # Assign predicted origins to latereplicatedorigins if replication time is more than 30
    latereplicatedorigins = [i for i in originswithreplicationtimes if i[2]>30]
    
    # Return early and late replication origins as a dictionary
    return {'Early':earlyreplicatedorigins,'Late':latereplicatedorigins}
    

# This function takes in an oem file and a replication time file, chromosome and assigns replication times to the predicted
# origins found in the sgrfiles for a particular chromosome and then assigns whether the predicted origin falls in the early replicated or later replicated 
# category. A dictionary is returned which contain early replicated origins under the key Early and late replicated origins under the key Late.
def earlyLateReplicatedOriginsForOEMFile(reptime,oemfile,chromosome):
    
    # Grab the predicted origins which include both the base and oem from the sgrcomb list
    predictedorigins = oc.getPredictedOriginsForOEMList(oemfile,chromosome)   
    
    # Assign replication times fr predicted origins from the reptime file
    potimes = assignReplicationTimesForPredictedOrigins(reptime,predictedorigins,chromosome)
    
    # Return early and late replicated origins
    earlylateRO = earlyLateReplicatedOrigins(potimes)
    
    return earlylateRO    

# This function takes in an oem file and a replication time file and goes through all 16 chromosomes and
# separates predicted origins between early and late replicated. A dictionary is returned which contain early replicated origins under 
# the key Early and late replicated origins under the key Late.
def getOEMSForEarlyvsLateReplicatedOriginsAllChromosomes(reptime,oemfile):

    # Store oems for early replicated origins for all chromosomes in this list oemsForEarlyRO
    oemsForEarlyRO = []
    
     # Store oems for late replicated origins for all chromosomes in this list oemsForLateRO
    oemsForLateRO = []    
    
    # Go through every chromosome    
    for c in range(1,17):
        chromosome = str(c)
        print chromosome
        # Grab early and late replicated predicted origins from sgr files for chromosome
        earlylateROforchrom = earlyLateReplicatedOriginsForOEMFile(reptime,oemfile,chromosome)
        # Grab oems for early replicated predicted origins for chromosome and store in oemsforearlyroforchrom
        oemsforearlyroforchrom = [i[1] for i in earlylateROforchrom['Early']]
        # Grab oems for late replicated predicted origins for chromosome and store in oemsforlateroforchrom
        oemsforlateroforchrom = [i[1] for i in earlylateROforchrom['Late']]
        # Add to the oem list for all chromosomes       
        oemsForEarlyRO.extend(oemsforearlyroforchrom)
        oemsForLateRO.extend(oemsforlateroforchrom)
    
    return {'early':oemsForEarlyRO,'late':oemsForLateRO}
    
# This function takes in reptime file and for every chromosome gives the number of bases replicated early and late if plotoption is False and if
# plotoption is True then plot a barplot for all 16 chromosomes showing early and late replication times
def getPlotNumBasesReplicatedEarlyAndLateAllChroms(reptimefile,plotoption):
    
    numearlylateAllChroms = []
    chromosomes = []
    # Go through every chromosome
    for c in range(1,17):
        chromosome = str(c)
        with open(reptimefile) as f:
            reptimesforchrom = [line.strip().split('\t')[1:] for line in f.readlines() if line.strip().split('\t')[0]==chromosome]
        numearlyreptimes = len([int(i[0]) for i in reptimesforchrom if float(i[1])<=30])
        numlatereptimes = len([int(i[0]) for i in reptimesforchrom if float(i[1])>30])
        numearlylateAllChroms.append([c,numearlyreptimes,numlatereptimes])
        chromosomes.append(chromosome)
    if plotoption:
        ind = np.arange(1,17)    # the x locations for the groups
        width = 0.1       # the width of the bars: can also be len(x) sequence
        
        # Grab total numbers for Early and Late replicating bases and store in two different variables early and late
        early = [i[1] for i in numearlylateAllChroms]
        late = [i[2] for i in numearlylateAllChroms]
        props = [float(i[1])/(i[1]+i[2]) for i in numearlylateAllChroms]
        plt1 = plt.bar(ind,early,width,color='g')
        plt2 = plt.bar(ind,late,width,color='r',bottom=early)
        plt.xticks(ind+width/2,chromosomes)
        plt.xlabel('Chromosomes')
        plt.ylabel('Number of Bases')
        plt.legend((plt1,plt2),('Early','Late'))
        for i,rect in enumerate(plt2):
            height = rect.get_height() + plt1[i].get_height()           
            plt.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%s'% float('%.2g' % (props[i])),ha='center', va='bottom')
        
        plt.title('Early vs Late Replicating regions of Chromosomes')
        plt.show()
        
    else:
        return numearlylateAllChroms