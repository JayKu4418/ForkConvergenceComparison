# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 17:24:41 2014

@author: jayashreekumar
"""

"""
This module contains all the functions that will produce graphs and results. All these functions require sgr files
as the inputs
"""
#######################################################################################################################
# IMPORT
import manipulatesgrfiles as ps
import matplotlib.pyplot as plt
import oem
import originscomparison as oc
import numpy as np
import replicationtimes as rt
import rpkm
cerevisiaechrmconvert = {'1':'chrI','2':'chrII','3':'chrIII','4':'chrIV','5':'chrV','6':'chrVI','7':'chrVII','8':'chrVIII','9':'chrIX','10':'chrX','11':'chrXI','12':'chrXII','13':'chrXIII','14':'chrXIV','15':'chrXV','16':'chrXVI'}
#######################################################################################################################
# FUNCTIONS

# This function takes in 4 sgr files representing watson and crick strands for wild type and the mutant strain. It will
# spit out a figure which has two subplots where the top subplot shows the depth coverage for Watson and Crick strands
# for Wild type and the bottom subplot shows the depth coverage for Watson and Crick strands for the mutant strain
def plotSGRWTMUT(WsgrfileWT,CsgrfileWT,WsgrfileMUT,CsgrfileMUT,chromosome):

    # This returns the sgr values from the file WsgrWT as a list and stores it in WatsonWT
    WatsonWT = ps.getSgrList(WsgrfileWT,chromosome)
    # This returns the sgr values from the file CsgrWT as a list and stores it in CrickWT
    CrickWT = ps.getSgrList(CsgrfileWT,chromosome)

    # This returns the sgr values from the file WsgrMUT as a list and stores it in WatsonMUT
    WatsonMUT = ps.getSgrList(WsgrfileMUT,chromosome)
    # This returns the sgr values from the file CsgrMUT as a list and stores it in CrickMUT
    CrickMUT = ps.getSgrList(CsgrfileMUT,chromosome)

    # Plot the fragment density for WT strain in the top subplot
    ps.plotSgrfig(WatsonWT,CrickWT,chromosome,211,'WT')
    # Create some space between the subplots
    plt.subplots_adjust(hspace=0.5)
    # Plot the fragment density for MUT strain in the bottom subplot
    ps.plotSgrfig(WatsonMUT,CrickMUT,chromosome,212,'MUT')
    # Show the plot
    plt.show()


def plotSGR4strains(WsgrfileStrain1,CsgrfileStrain1,Strain1label,WsgrfileStrain2,CsgrfileStrain2,Strain2label,WsgrfileStrain3,CsgrfileStrain3,Strain3label,WsgrfileStrain4,CsgrfileStrain4,Strain4label,chromosome):

    # This returns the sgr values from the file WsgrWT as a list and stores it in WatsonWT
    WatsonStrain1 = ps.getSgrList(WsgrfileStrain1,chromosome)
    # This returns the sgr values from the file CsgrWT as a list and stores it in CrickWT
    CrickStrain1 = ps.getSgrList(CsgrfileStrain1,chromosome)

    # This returns the sgr values from the file WsgrWT as a list and stores it in WatsonWT
    WatsonStrain2 = ps.getSgrList(WsgrfileStrain2,chromosome)
    # This returns the sgr values from the file CsgrWT as a list and stores it in CrickWT
    CrickStrain2 = ps.getSgrList(CsgrfileStrain2,chromosome)
    
    # This returns the sgr values from the file WsgrWT as a list and stores it in WatsonWT
    WatsonStrain3 = ps.getSgrList(WsgrfileStrain3,chromosome)
    # This returns the sgr values from the file CsgrWT as a list and stores it in CrickWT
    CrickStrain3 = ps.getSgrList(CsgrfileStrain3,chromosome)
    
    # This returns the sgr values from the file WsgrWT as a list and stores it in WatsonWT
    WatsonStrain4 = ps.getSgrList(WsgrfileStrain4,chromosome)
    # This returns the sgr values from the file CsgrWT as a list and stores it in CrickWT
    CrickStrain4 = ps.getSgrList(CsgrfileStrain4,chromosome)

    # Plot the fragment density for WT strain in the top subplot
    ps.plotSgrfig(WatsonStrain1,CrickStrain1,chromosome,411,Strain1label,'b','g')
    # Create some space between the subplots
    plt.subplots_adjust(hspace=0.5)
    # Plot the fragment density for MUT strain in the bottom subplot
    ps.plotSgrfig(WatsonStrain2,CrickStrain2,chromosome,412,Strain2label,'r','c')
    
    plt.subplots_adjust(hspace=0.5)
    # Plot the fragment density for MUT strain in the bottom subplot
    ps.plotSgrfig(WatsonStrain3,CrickStrain3,chromosome,413,Strain3label,'m','y')

    plt.subplots_adjust(hspace=0.5)
    # Plot the fragment density for MUT strain in the bottom subplot
    ps.plotSgrfig(WatsonStrain4,CrickStrain4,chromosome,414,Strain4label,'#9B30FF','#FF7F00')
    
    # Show the plot
    plt.show()

# This function takes in 2 oem files for wild type and the mutant strain. It will
# spit out a figure which has two subplots where the top subplot shows the OEM for Watson and Crick strands
# for Wild type and the bottom subplot shows the OEM for Watson and Crick strands for the mutant strain.If the showorigins is True, the
# confirmed and likely origins will be overlaid on the oem plots. If showreptime is true, the replication time will be binned between 
# early and late replicating origins and early, late replication regions will be overlaid on the oem plots
def plotOEMWTMUT(oemfileWT,oemfileMUT,chromosome,showorigins,showreptime,reptimefile=''):
    
    # if showreptime is true, then grab all rep times for chromosome and bin the bases in earlyreptimebases and latereptimebass.
    if showreptime:
        with open(reptimefile) as f:
            reptimeforchrom = [line.strip().split('\t')[1:] for line in f.readlines() if line.strip().split('\t')[0]==chromosome]
        earlyreptimebases = [int(i[0]) for i in reptimeforchrom if float(i[1])<=30]
        latereptimebases = [int(i[0]) for i in reptimeforchrom if float(i[1])>30]    
    
    # Plot OEM for WT with origins if showorigins is TRUE, otherwise plot oem but don't show confirmed or likely orgins
    if showorigins:
        oem.plotOEM(oemfileWT,chromosome,True,211,'WT','green')
    else:
        oem.plotOEM(oemfileWT,chromosome,False,211,'WT','green')
    
    # if showreptime is true, plot vertical lines to indicate where the early and late replicating bases are
    if showreptime:
        for base in earlyreptimebases:
            plt.vlines(base,-1,1,colors='red')
        for base in latereptimebases:
            plt.vlines(base,-1,1,colors='black')
    
    # Create some space between the subplots
    plt.subplots_adjust(hspace=0.5)
    
    # Plot OEM for MUT with origins if showorigins is TRUE, otherwise plot oem but don't show confirmed or likely orgins
    if showorigins:
        oem.plotOEM(oemfileMUT,chromosome,True,212,'MUT','blue')
    else:
        oem.plotOEM(oemfileMUT,chromosome,False,212,'MUT','blue')
    
    # if showreptime is true, then plot vertical lines to indicate where the early and late replicating bases are
    if showreptime:
        for base in earlyreptimebases:
            plt.vlines(base,-1,1,colors='red')
        for base in latereptimebases:
            plt.vlines(base,-1,1,colors='black')    
    
    plt.show()

# This function takes in 2 oem files for wild type and the mutant strain. It will
# spit out a figure which has two subplots where the top subplot shows the OEM for Watson and Crick strands
# for Wild type and the bottom subplot shows the OEM for Watson and Crick strands for the mutant strain.If the showorigins is True, the
# confirmed and likely origins will be overlaid on the oem plots. If showreptime is true, the replication time will be binned between 
# early and late replicating origins and early, late replication regions will be overlaid on the oem plots
def plotOEM4strainsT(oemfilestrain1,strain1label,oemfilestrain2,strain2label,oemfilestrain3,strain3label,oemfilestrain4,strain4label,chromosome,showorigins,showreptime,reptimefile=''):
    
    # if showreptime is true, then grab all rep times for chromosome and bin the bases in earlyreptimebases and latereptimebass.
    if showreptime:
        with open(reptimefile) as f:
            reptimeforchrom = [line.strip().split('\t')[1:] for line in f.readlines() if line.strip().split('\t')[0]==chromosome]
        earlyreptimebases = [int(i[0]) for i in reptimeforchrom if float(i[1])<=30]
        latereptimebases = [int(i[0]) for i in reptimeforchrom if float(i[1])>30]    
    
    # Plot OEM for WT with origins if showorigins is TRUE, otherwise plot oem but don't show confirmed or likely orgins
    if showorigins:
        oem.plotOEM(oemfilestrain1,chromosome,True,411,strain1label,'green')
        plt.subplots_adjust(hspace=0.5)
        oem.plotOEM(oemfilestrain2,chromosome,True,412,strain2label,'blue')
        plt.subplots_adjust(hspace=0.5)
        oem.plotOEM(oemfilestrain3,chromosome,True,413,strain3label,'red')
        plt.subplots_adjust(hspace=0.5)
        oem.plotOEM(oemfilestrain4,chromosome,True,414,strain4label,'#FF7F00')
    else:
        oem.plotOEM(oemfilestrain1,chromosome,False,411,strain1label,'green')
        plt.subplots_adjust(hspace=0.5)
        oem.plotOEM(oemfilestrain2,chromosome,False,412,strain2label,'blue')
        plt.subplots_adjust(hspace=0.5)
        oem.plotOEM(oemfilestrain3,chromosome,False,413,strain3label,'red')
        plt.subplots_adjust(hspace=0.5)
        oem.plotOEM(oemfilestrain4,chromosome,False,414,strain4label,'#FF7F00')
    
    # if showreptime is true, plot vertical lines to indicate where the early and late replicating bases are
    if showreptime:
        for base in earlyreptimebases:
            plt.vlines(base,-1,1,colors='red')
        for base in latereptimebases:
            plt.vlines(base,-1,1,colors='black')  
    
    plt.show()

# Plot a stacked barplot showing total number of predicted origins for WT and MUT within 2.5kB of confirmed, likely,
# dubious origins. Show origins as well that are predicted but not found in any of these classifications.
def plotPredictedOriginsTotals(oemfileWT,oemfileMUT):

    # Get the total number of predicted origins for WT compared against confirmed, likely, dubious origins
    totalPOcomparedWT = oc.getnumComparedOriginsWTMUT(oemfileWT,2500)

    # Get the total number of predicted origins for MUT against confirmed, likely, dubious origins
    totalPOcomparedMUT = oc.getnumComparedOriginsWTMUT(oemfileMUT,2500)

    originsConfirmed = [totalPOcomparedWT['C'],totalPOcomparedMUT['C']]

    originsLikely = [totalPOcomparedWT['L'],totalPOcomparedMUT['L']]

    originsDubious = [totalPOcomparedWT['D'],totalPOcomparedMUT['D']]

    originsNone = [totalPOcomparedWT['N'],totalPOcomparedMUT['N']]

    ind = np.arange(2)    # the x locations for the groups
    width = 0.35       # the width of the bars: can also be len(x) sequence

    plt1 = plt.bar(ind, originsConfirmed,   width, color='g')
    plt2 = plt.bar(ind, originsLikely, width, color='y',bottom=originsConfirmed)
    plt3 = plt.bar(ind, originsDubious, width, color='r',bottom=[originsConfirmed[j]+originsLikely[j] for j in range(len(originsConfirmed))])
    plt4 = plt.bar(ind, originsNone, width, color='grey',bottom=[originsConfirmed[j]+originsLikely[j]+originsDubious[j] for j in range(len(originsConfirmed))])
    plt.xticks(ind+width/2., ('WT', 'KO'))
    plt.legend((plt1,plt2,plt3,plt4),('Confirmed','Likely','Dubious','None'),loc='upper center')
    plt.ylabel('Number of Origins')

    plt.show()

# Plot the percentage of matched origins between predicted origin for WT and MUT as a function of distance between
# the two 'matching' origins
def plotPercentageMatchedOriginsWTMUT(oemfileWT,oemfileMUT,oemlimit,maxdistance):

    windows = range(0,maxdistance+1,50)   
    # Create empty lists to store the total number of matched origins and predicted origins across all chromosomes
    totalMatchedOrigins=[0]*len(windows)
    totalPredictedOrigins = [0]*len(windows)
    # Go through all the chromosomes
    for c in range(1,17):
        chromosome = str(c)
        print('\n')
        print chromosome
        
        # Predicted Origins for WT and MUT
        poWTMUT = oc.predictedOriginsForWTandMUT(oemfileWT,oemfileMUT,chromosome,oemlimit)

        # Base coordinates of predicted origins for WT and MUT
        poWT = poWTMUT['WTbases']
        poMUT = poWTMUT['MUTbases']

        # Go through range of basewindows to find matching predicted origins
        for i in range(0,len(windows)):
            # Add the total number of predicted origins found in WT and MUT - we want a genome wide analysis
            totalPredictedOrigins[i] = totalPredictedOrigins[i] + (len(poWT) + len(poMUT))
            # Find the matching predicted origins for WT and MUT for the base distance in question
            matchedoriginsForchromForDist = oc.matchingOriginsBasedOnWindow(poWT,poMUT,windows[i])
            # Add the number of matching predicted origins for base distance in question to total number of matched
            # predicted origins - we want a genome wide analysis
            totalMatchedOrigins[i] = totalMatchedOrigins[i] + 2*(len(matchedoriginsForchromForDist))

    # Calculate the percentage of matched origins out of predicted origins as the base distance window increases
    percentageOfMatchedOriginsToPredictedOrigins = 100*(np.array(totalMatchedOrigins))/(np.array(totalPredictedOrigins))

    plt.plot(windows,percentageOfMatchedOriginsToPredictedOrigins,'.r')
    plt.xlabel('Distance between matched origins')
    plt.ylabel('Percent of predicted origins Matched')
    plt.show()

# Non matching predicted origins for all chromosomes
def displayNMPOsWTandMUTAllChroms(oemfileWT,oemfileMUT,oemlimit,matchwindow,nonmatchoemlimit,writeFile=''):
    
    allNMPOsWTandMUT = {}    
        
    # Go through all chromosomes 
    for c in range(1,17):
        chromosome = str(c)
        print chromosome

        # Predicted Origins for WT and MUT
        poWTMUT = oc.predictedOriginsForWTandMUT(oemfileWT,oemfileMUT,chromosome,oemlimit)

        # Grab non matching predicted origins for WT and MUT
        nmpoWTMUT = oc.nonmatchingOriginsBasedOnWindow(poWTMUT['WTwOEM'],poWTMUT['MUTwOEM'],matchwindow)

        # Grab non matching predicted origins for WT that have an oem higher than the nonmatchoemlimit
        nmpoWT = [i for i in nmpoWTMUT['WT'] if i[0][1]>nonmatchoemlimit]
        nmpoMUT = [i for i in nmpoWTMUT['MUT'] if i[0][1]>nonmatchoemlimit]
        
        # Add non matching origins as the value and the key is chromosome
        allNMPOsWTandMUT[chromosome] = [nmpoWT,nmpoMUT]
    
    # if writeFile has a value, write non matching predicted origins to writeFile
    if writeFile:
        with open(writeFile,'w') as fw:
            for c in range(1,17):
                chrm = str(c)
                fw.write('Chromosome ' + chrm+'\n')
                if len(allNMPOsWTandMUT[chrm][0])!=0:
                    for i in allNMPOsWTandMUT[chrm][0][0]:
                        fw.write('WT' + '\t'+ str(i[0]) + '\t' + str(i[1]) + '\n')
                if len(allNMPOsWTandMUT[chrm][1])!=0:
                    for i in allNMPOsWTandMUT[chrm][1][0]:
                        fw.write('MUT'+ '\t'+ str(i[0]) + '\t' + str(i[1]) + '\n') 
    else:
        # If writeFile is has no value, then Return non matching predicted origins
        return allNMPOsWTandMUT

# Plot fragment densities for non matching predicted origins in WT and MUT on a per chromosome basis. 
# This function takes in 4 sgr files representing watson and crick strands for wild type and the mutant strain. It also takes in 2 oem files for WT
# and MUT, to calculate where the non matching predicted origins exist for that chromosome. At non matching predicted origins, plot the fragment
# density around the region for both WT and MUT. 
def plotFDforNMPOsWTandMUT(WsgrfileWT,CsgrfileWT,oemfileWT,WsgrfileMUT,CsgrfileMUT,oemfileMUT,chromosome,oemlimit,matchwindow,nonmatchoemlimit):

    # Combine the WsgrWT and CsgrWT files into one combined list where the first column is the base position, second
    # column is the fragment density for watson strand and third column is the fragment density for crick strand
    combsgrWT = oem.combineFRsgrfiles(WsgrfileWT,CsgrfileWT,chromosome)

    # Combine the WsgrMUT and CsgrMUT files into one combined list where the first column is the base position, second
    # column is the fragment density for watson strand and third column is the fragment density for crick strand
    combsgrMUT = oem.combineFRsgrfiles(WsgrfileMUT,CsgrfileMUT,chromosome)

    # Predicted Origins for WT and MUT
    poWTMUT = oc.predictedOriginsForWTandMUT(oemfileWT,oemfileWT,chromosome,oemlimit)

    # Grab non matching predicted origins for WT and MUT
    nmpoWTMUT = oc.nonmatchingOriginsBasedOnWindow(poWTMUT['WTwOEM'],poWTMUT['MUTwOEM'],matchwindow)

    # Grab non matching predicted origins for WT that have an oem higher than the nonmatchoemlimit
    nmpoWT = [i for i in nmpoWTMUT['WT'] if i[0][1]>nonmatchoemlimit]
    nmpoMUT = [i for i in nmpoWTMUT['MUT'] if i[0][1]>nonmatchoemlimit]

    # Grab oem for WT and MUT sgrcomblist
    #oemWT = oem.calculateOEM(combsgrWT,oemwindow)
    #oemMUT = oem.calculateOEM(combsgrMUT,oemwindow)

    # Plot OEM comparisons for non matching predicted origins found in WT
    numNMPOWT = len(nmpoWT)
    numNMPOMUT = len(nmpoMUT)
    totalnumNMPO = numNMPOWT + numNMPOMUT
    
    # Set index to plot to be 1
    k = 1
    # Go through all the non matched predicte origins for WT and plot the fragment density at a 20kB window surrounding the predicted origin for both WT and MUT
    if numNMPOWT > 0:
        print 'Here'
        for i in range(0,numNMPOWT):
            base = nmpoWT[i][0][0]
            oc.plotFragDensityWTMUTatspecificbase('WT',base,combsgrWT,combsgrMUT,totalnumNMPO,k)

   # Go through all the non matched predicte origins for MUT and plot the fragment density at a 20kB window surrounding the predicted origin for both WT and MUT
    if numNMPOMUT > 0:    
        for i in range(0,numNMPOMUT):
            base = nmpoMUT[i][0][0]
            print str(base)
            oc.plotFragDensityWTMUTatspecificbase('MUT',base,combsgrWT,combsgrMUT,totalnumNMPO,k)
    
    plt.show()

# Plot box plot for OEMS of early replicated origins vs late replication origins throughout the genome.This function takes in 2 oem files for WT and MUT 
# and reptime file to assign replication times to origins. It will return box plots of oems for early vs late replicating origins for both WT and MUT strains.
def plotOEMSEarlyLateReplicatedOriginsWTandMUT(oemfilestrain1,oemfilestrain2,oemfilestrain3,oemfilestrain4,reptime):
    
    # Grab OEMS for early and late replicated origins for strain 1
    oemsEarlyLateRO_strain1 = rt.getOEMSForEarlyvsLateReplicatedOriginsAllChromosomes(reptime,oemfilestrain1)
    
    # Grab OEMS for early and late replicated origins for strain 2
    oemsEarlyLateRO_strain2 = rt.getOEMSForEarlyvsLateReplicatedOriginsAllChromosomes(reptime,oemfilestrain2)
    
    # Grab OEMS for early and late replicated origins for strain 3
    oemsEarlyLateRO_strain3 = rt.getOEMSForEarlyvsLateReplicatedOriginsAllChromosomes(reptime,oemfilestrain3)
    
    # Grab OEMS for early and late replicated origin for strain 4
    oemsEarlyLateRO_strain4 = rt.getOEMSForEarlyvsLateReplicatedOriginsAllChromosomes(reptime,oemfilestrain4)
    
    bp = plt.boxplot([oemsEarlyLateRO_strain1['early'],oemsEarlyLateRO_strain1['late'],oemsEarlyLateRO_strain2['early'],oemsEarlyLateRO_strain2['late'],oemsEarlyLateRO_strain3['early'],oemsEarlyLateRO_strain3['late'],oemsEarlyLateRO_strain4['early'],oemsEarlyLateRO_strain4['late']],1)
    
    for line in bp['medians']:
        x1, y = line.get_xydata()[0]
        x2 = line.get_xydata()[1][0]
        x = (x1+x2)/2
        plt.text(x, y, '%.3f' % y,horizontalalignment='center')
        
    plt.xticks(range(1,9),('Early-WT','Late-WT','Early-pif1m2KO','Late-pif1m2KO','Early-rrm3dKO','Late-rrm3dKO','Early-rrm3d_pif1m2KO','Late-rrm3d_pif1m2KO'))

    plt.ylabel('OEM')

    plt.title('OEMs of Early and Late Replicating Origins')
    
    plt.show()

# This function plots RPKM for 4 sgrfiles, 2 files each for two species WT and MUT, one for the watson strand (sgrfileW)
# and one for the crick strand (sgrfileC).
# It also takes in a chromosome value (chromosome) to know which chromosome to plot the RPKMs for. In the top subplot
# it will be the WT strain and the bottom subplot will be the MUT strain
def plotRPKMfor2strains(sgrfileW_WT,sgrfileC_WT,sgrfileW_MUT,sgrfileC_MUT,chromosome,strains):
    
    # Convert Watson sgr file to a watson sgrlist for a chromosome for WT
    sgrlistforchromW_WT = ps.getSgrList(sgrfileW_WT,chromosome)
    
    # Convert Crick sgr file to a crick sgrlist for a chromosome for WT
    sgrlistforchromC_WT = ps.getSgrList(sgrfileC_WT,chromosome)
    
    # Convert Watson sgr file to a watson sgrlist for a chromosome for MUT
    sgrlistforchromW_MUT = ps.getSgrList(sgrfileW_MUT,chromosome)
    
    # Convert Crick sgr file to a crick sgrlist for a chromosome for MUT
    sgrlistforchromC_MUT = ps.getSgrList(sgrfileC_MUT,chromosome)    
    
    # Plot RPKM for Watson and Crick strands for the wild type species
    rpkm.plotRPKM(sgrlistforchromW_WT,sgrlistforchromC_WT,strains[0],chromosome,211)
    
    # Create some space between the subplots
    plt.subplots_adjust(hspace=0.5)    
    
    # Plot RPKM for Watson and Crick strands for the knock out species
    rpkm.plotRPKM(sgrlistforchromW_MUT,sgrlistforchromC_MUT,strains[1],chromosome,212)
    
    # Show plot
    plt.show()

# This function plots base window as the x axis and the proportion of genes of signficance found near windows around peaks
# as the y axis. There will be line for peaks found. There will also be lines to represent proportion of genes of significance 
# found near random bases.  
def plotPropOfSigGenesFound(option,propsForPeaksWT,propsForPeaksMUT1,propsForPeaksMUT2,propsForPeaksMUT3,meanAndStdpropsForRandomBases,maxchromrange):
   
   # grab just the mean props for random bases
   meanPropsForRandomBases = [i[0] for i in meanAndStdpropsForRandomBases]
   
   # grab the max prop for random bases
   maxPropsForRandomBases = [i[0]+i[1] for i in meanAndStdpropsForRandomBases]
   
   # grab the min prop for random bases
   minPropsForRandomBases =[i[0]-i[1] for i in meanAndStdpropsForRandomBases]
   
   basewindows = range(0,maxchromrange+1,100)
   
   # Special case
   #propsForPeaksWT = [propsForPeaksWT[i] for i in range(0,31,2)]
   
   plt.plot(basewindows,meanPropsForRandomBases,'.-r',label='Random Mean')
   plt.plot(basewindows,maxPropsForRandomBases,'.-b',label='Random Max')
   plt.plot(basewindows,minPropsForRandomBases,'.-m',label='Random Min')
   plt.plot(basewindows,propsForPeaksWT,'.-g',label='WT Peaks')
   plt.plot(basewindows,propsForPeaksMUT1,'.-y',label='rrm3dKO Peaks')
   plt.plot(basewindows,propsForPeaksMUT2,'.-c',label='pif1m2KO Peaks')
   plt.plot(basewindows,propsForPeaksMUT3,'.-k',label='doubleKO Peaks')
   plt.legend(loc=2)
   
   plt.xlabel('Window')
   plt.ylabel('Proportion')
   if (option == 'HighlyTxGenes'):
       plt.title('Proportion of Highly Transcribed Genes Found ')
   elif(option=='G4'):
       plt.title('Proportion of G Quadraplexes Found')
   else:
       plt.title('Proportion of tRNA Genes Found')

   plt.show()
   
# This function converts a "sgr like" file to a wig file
def convertToWigFile(readfile,writefile,ifcrick):
    with open(writefile,'w') as fw:
        fw.write('track type=wiggle_0'+'\n')
        for c in range(1,17):
            chrom = str(c)
            print(chrom)
            with open(readfile) as f:
                chrmlines = [line.strip().split('\t') for line in f if line.strip().split('\t')[0]==cerevisiaechrmconvert[chrom]]
            if len(chrmlines)!=0:
                fw.write('variableStep chrom=' + cerevisiaechrmconvert[chrom] + ' span=1' + '\n')
                for i in chrmlines:
                    if ifcrick:
                        if (i[2][0]=='-'):
                            fw.write(i[1] + '\t' +  i[2][1:] + '\n')
                        else:
                            fw.write(i[1] + '\t' +  '-' + i[2] + '\n')
                    else:
                        fw.write(i[1] + '\t' +  i[2] + '\n')
                
   