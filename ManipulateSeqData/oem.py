__author__ = 'jayashreekumar'

"""
The functions in this script will calculate origin efficiency metric (OEM) for given SGR files and plot it optionally.
"""

################################################################################################
# IMPORT
import manipulatesgrfiles as sgr
import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter
################################################################################################
# FUNCTIONS

# for a particular base at index i, calculate either WL for 20000 previous bases or WR for 20000 following bases based
# on option given

# Takes in a file which has confirmed origins and returns a list of origins based on a particular chromosome and status
def readOriginsfile(originfile,chromosome,status):
    with open(originfile) as f:
        originschrom = [line.strip().split('\t') for line in f if line.strip().split('\t')[0] == chromosome and line.strip().split('\t')[5] == status]

    # We only want the base start and end coordinates for origin, but these are strings
    originswpos = [i[1:3] for i in originschrom]
    originswposint = []
    # We convert the base start and end coordinates for origins to ints
    for i in originswpos:
        o = [int(i[0]),int(i[1])]
        originswposint.append(o)

    return originswposint

# Calculate watson/watson+crick for left or right of  every 'perbasebox' base (reference base),
# compare watson and crick depth coverage +/- 'window' number of bases from each reference base
def calculateWLorWR(sgrcomblist,i,option,window):
    if option == 'left':
        # Grab all depth coverage for watson strand left of reference base and crick strand left of reference base
        watsonleft = [j[2] for j in sgrcomblist[i-window:i]]
        crickleft = [j[3] for j in sgrcomblist[i-window:i]]
        # only calculate watsonleft/watsonleft+crickleft is watsonleft+crickleft is not zero
        if (sum(watsonleft) + sum(crickleft)) != 0:
            print(float(sum(watsonleft)))
            print(sum(crickleft))
            WL = float(abs(sum(watsonleft)))/(abs(sum(watsonleft))+abs(sum(crickleft)))
            print(WL)
        else:
            WL = 0
        return WL
    else:
        # Grab all depth coverage for watson strand right of reference base and crick strand right of reference base
        watsonright = [j[2] for j in sgrcomblist[i:i+window]]
        crickright = [j[3] for j in sgrcomblist[i:i+window]]
        if (sum(watsonright)+sum(crickright)) != 0:
            WR = float(abs(sum(watsonright)))/(abs(sum(watsonright))+abs(sum(crickright)))
            print(float(sum(watsonright)))
            print(sum(crickright))
            print(WR)
        else:
            WR = 0
        return WR

# Calculate the origin efficiency metric by take left - minus from values obtained using calculateWLorWR. There is
# a dictionary option if one wanted to return OEM as a dictionary where the base is the key and the corresponding oem it
# the value
def calculateOEM(sgrcomblist,window=10000,skipbasebox=50,dictoption=False):
 # Get number of bases from sgr combination list
    numbases = len(sgrcomblist)
    # Empty list or dictionary to store OEMs for every reference base
    if dictoption:
        oem = {}
    else:
        oem = []
    # We only want to calculate OEM for every perbasebox bases
    for b in range(0,numbases,skipbasebox):
        # we only calculate OEMs for reference bases as long as there are window number of bases left and right of them
        if b >= window and b <= (numbases-window):
            print(b)
            oem4b = calculateWLorWR(sgrcomblist,b,'left',window) - calculateWLorWR(sgrcomblist,b,'right',window)
            print(oem4b)
        else:
            oem4b = 0
        # indoem contains a list which has both reference base and the OEM calculated at that reference base
        if dictoption:
            oem[b] = oem4b
        else:
            indoem = [b,oem4b]
            oem.append(indoem)
    return oem

# Calculate the origin efficiency metric by take left - minus from values obtained using calculateWLorWR. There is
# a dictionary option if one wanted to return OEM as a dictionary where the base is the key and the corresponding oem it
# the value
def calculateOEMAtSpecificBases(baselist,sgrcomblist,window,dictoption=False):
 # Get number of bases from sgr combination list
    numbases = len(sgrcomblist)
    # Empty list or dictionary to store OEMs for every reference base
    if dictoption:
        oem = {}
    else:
        oem = []
    # We only want to calculate OEM for every perbasebox bases
    for b in range(0,len(baselist)):
        # we only calculate OEMs for reference bases as long as there are window number of bases left and right of them
        if baselist[b] >= window and baselist[b] <= (numbases-window):
            oem4b = calculateWLorWR(sgrcomblist,baselist[b],'left',window) - calculateWLorWR(sgrcomblist,baselist[b],'right',window)
        else:
            oem4b = 0
        # indoem contains a list which has both reference base and the OEM calculated at that reference base
        if dictoption:
            oem[baselist[b]] = oem4b
        else:
            indoem = [baselist[b],oem4b]
            oem.append(indoem)
    return oem

# Write oems calculated to a file
def writeOEMForKnownOriginstofile(sgrfilew,sgrfilec,knownoriginsfile,window,writefile):
    with open(knownoriginsfile) as f:
        knownorigins = [line.strip().split('\t') for line in f]
    with open(writefile,'w') as fw:
        for c in range(1,17):
            chromosome = str(c)
            print(chromosome)
            wc = combineFRsgrfiles(sgrfilew,sgrfilec,chromosome)
            specificbasesW = [int(i[1]) for i in knownorigins if i[5]=='+' and i[0]==chromosome]
            #specificbasesC = [int(i[2]) for i in knownorigins if i[5]=='-' and i[0]==chromosome]
            specificbasesC = [int(i[2]) for i in knownorigins if i[5]=='-' and i[0]==chromosome]
            oemcW = calculateOEMAtSpecificBases(specificbasesW,wc,window)
            oemcC = calculateOEMAtSpecificBases(specificbasesC,wc,window)
            for i in oemcW:
                fw.write(chromosome + '\t' + str(i[0]) + '\t' + str(i[1]) + '\t' + '+' + '\n')
            for i in oemcC:
                fw.write(chromosome + '\t' + str(i[0]) + '\t' + str(i[1]) + '\t' + '-' + '\n')


# Write oems calculated to a file
def writeOEMtofile(sgrfilew,sgrfilec,window,writefile):
    with open(writefile,'w') as fw:
        for c in range(1,17):
            chromosome = str(c)
            print(chromosome)
            wc = combineFRsgrfiles(sgrfilew,sgrfilec,chromosome)
            oemc = calculateOEM(wc,window)
            for i in oemc:
                fw.write(chromosome + '\t' + str(i[0]) + '\t' + str(i[1]) + '\n')
                
# Combine 2 sgr watson and crick files into one combined list
def combineFRsgrfiles(sgrfilef,sgrfiler,chromosome):
    f = sgr.getSgrList(sgrfilef,chromosome)
    r = sgr.getSgrList(sgrfiler,chromosome)

    return sgr.combineFRsgrlists(f,r)

def plotOEMDiff2Strains(oemfile1,oemfile2,chromosome,sbnum,color):
    # grab oems for oemfile1
    with open(oemfile1) as f:
        oem_for_Strain1 = [line.strip().split('\t')[1:] for line in f.readlines() if line.strip().split('\t')[0]==chromosome]
        
    # grab oems for oemfile2
    with open(oemfile2) as f:
        oem_for_Strain2 = [line.strip().split('\t')[1:] for line in f.readlines() if line.strip().split('\t')[0]==chromosome]
    
    ind = [int(i[0]) for i in oem_for_Strain1]
    
    oems1 = np.array([float(i[1]) for i in oem_for_Strain1])
    oems2 = np.array([float(i[1]) for i in oem_for_Strain2])
    
    oemsDiff = oems1 - oems2
    
    # Start up a figure
    fig = plt.figure(1)
    # Add a subplot
    ax = fig.add_subplot(sbnum)
    # Plot reference base vs OEM
    ax.plot(ind,oemsDiff,color=color)
    # Set the x limit to be the number of bases for that chromosome
    plt.xlim(xmin=0,xmax=ind[len(ind)-1])
    # Set y limits to be between -1 and 1 since that's max and min of OEM
    plt.ylim(ymin=-1,ymax=1)
    # Set a clear x axis
    plt.axhline(xmax=ind[len(ind)-1],color='black')
    # Fill in the curves so we get filled in graph for WT in black
    plt.fill_between(ind,oemsDiff,color=color)
    plt.ylabel('Origin Efficiency Metric')
    plt.xlabel('Base Position')
    # we want x ticks to be displayed in kilobases
    #ax.set_xticklabels([str(x/1000) for x in ax.get_xticks()])
    # Y axis ticks are going to be absolute value, since OEM for KO is also positive but we just
    # want it below the x axis
    ax.set_yticklabels([str(abs(x)) for x in ax.get_yticks()])
# Plot OEM given a sgrcombfile (which contains sgr for both Watson and Crick strand). The option showorigins if true
# will show confirmed and likely origins from OriDB
def plotOEM(oemfile,chromosome,showorigins,sbnum,strain,color):

    # Calculate OEM for sgrcombfile
    with open(oemfile) as f:
        oem_for_sgr = [line.strip().split('\t')[1:] for line in f.readlines() if line.strip().split('\t')[0]==chromosome]
    # Grab bases at which OEM is calculated for sgrcombfile
    ind = [int(i[0]) for i in oem_for_sgr]
    # Grab OEM for sgrcombfile
    just_oem = [float(i[1]) for i in oem_for_sgr]

    # Start up a figure
    fig = plt.figure(1)
    # Add a subplot
    ax = fig.add_subplot(sbnum)
    # Plot reference base vs OEM
    ax.plot(ind,just_oem,color=color)
    # Set the x limit to be the number of bases for that chromosome
    plt.xlim(xmin=0,xmax=ind[len(ind)-1])
    # Set y limits to be between -1 and 1 since that's max and min of OEM
    plt.ylim(ymin=-1,ymax=1)
    # Set a clear x axis
    plt.axhline(xmax=ind[len(ind)-1],color='black')
    # Fill in the curves so we get filled in graph for WT in black
    plt.fill_between(ind,just_oem,color=color)
    plt.ylabel('Origin Efficiency Metric')
    plt.xlabel('Base Position')
    # we want x ticks to be displayed in kilobases
    #ax.set_xticklabels([str(x/1000) for x in ax.get_xticks()])
    # Y axis ticks are going to be absolute value, since OEM for KO is also positive but we just
    # want it below the x axis
    ax.set_yticklabels([str(abs(x)) for x in ax.get_yticks()])

    # If we want to show confirmed origins and likely origins
    if showorigins:
        origins_confirmed = readOriginsfile('S-cerevisiae_origins.txt',chromosome,'Confirmed')
        origins_likely = readOriginsfile('S-cerevisiae_origins.txt',chromosome,'Likely')
        for i in origins_confirmed:
            # Displays a translucent green box spanning confirmed origin coordinates
            plt.axvspan(i[0], i[1], facecolor='g', alpha=0.5)

        for i in origins_likely:
            # Displays a translucent red box spanning likely origin coordinates
            plt.axvspan(i[0], i[1], facecolor='r', alpha=0.5)
        # Display a legend explaining translucent red and green boxes
        #proxy = [plt.Rectangle((0,0),1,1,fc = 'g',alpha=0.5),plt.Rectangle((0,0),1,1,fc='r',alpha=0.5)]
        #plt.legend(proxy, ["Confirmed Origins","Likely Origins"])

    plt.title(strain + ' - Origin Efficiency Metric of Chromosome ' + chromosome)
    
# split origin files into different bins
def splitOriginsFilesIntoBins(originsfile,numbins,strand):
    with open(originsfile) as f:
        origins = [line.strip().split('\t') for line in f if line.strip().split('\t')[5]==strand]
        
    originsmod = [[i[0],int(i[1]),int(i[2]),float(i[4])] for i in origins]
    
    originsmod = sorted(originsmod, key=itemgetter(3))
    
    numoforigins = len(originsmod)
    
    numperbin = int(numoforigins/numbins)
    
    origindict = {}
    
    for i in range(numbins):
        origindict[i]=originsmod[(i)*numperbin:(i+1)*numperbin]
    
    return origindict
    
def differencebtwn2OEMfiles(oemfile1,oemfile2,writedifffile):
    with open(writedifffile,'w') as fw:
        for c in range(1,17):
            chromosome = str(c)
            print(chromosome)
            with open(oemfile1) as f:
                oem1 = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
            with open(oemfile2) as f:
                oem2 = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
            oemdiff = np.array(oem1)-np.array(oem2)
            for r in range(1,len(oemdiff)+1):
                fw.write(chromosome + '\t' + str(r) + '\t' + str(oemdiff[r-1]) + '\n')
            