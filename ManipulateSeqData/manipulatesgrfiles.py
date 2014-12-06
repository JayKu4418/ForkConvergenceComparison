__author__ = 'jayashreekumar'

"""
The functions in this script will plot 4 sgr files as Bases vs Depth Coverage.
"""
################################################################################################
# IMPORT
import matplotlib.pyplot as plt

################################################################################################
# FUNCTIONS

# Get an sgr list from an sgr file for a particular chromosome
def getSgrList(sgrfile, chromosome):

    # Grab all values from sgrfile and put it into a list sgrlist
    with open(sgrfile) as fileobject:
    	sgrlist = [line.strip().split('\t') for line in fileobject if line.strip().split('\t')[0] == chromosome]

    # Convert entries (base and depth coverage) in sgr list to ints
    sgrlist = [[j[0],int(j[1]),int(j[2])]  for j in sgrlist]

    return sgrlist

# Combine 2 sgr lists (Watson and Crick) into one big list, with 3 columns, first column is base, second column is depth
# coverage for Watson strand, third column is depth coverage for Crick strand
def combineFRsgrlists(sgrFlist, sgrRlist):
    sgrComblist = []
    if len(sgrFlist) == len(sgrRlist):
        for i in range(len(sgrFlist)):
            indrow = [sgrFlist[i][0],sgrFlist[i][1],sgrFlist[i][2],sgrRlist[i][2]]
            sgrComblist.append(indrow)
    return sgrComblist

# Plot Depth Coverage for the Watson and Crick strand for a particular chromosome
def plotSgrfig(sgrlistfor,sgrlistrev,chromosome,sbnum,strain,color1,color2):
    # bases contains a list of bases for chromosome
    bases = [i[1] for i in sgrlistfor]
    # Grab depth coverage for Watson strand
    coveragefor = [i[2] for i in sgrlistfor]
    # Grab depth coverage (multiplied for -1 so that it appears below x axis) for Crick strand
    coveragerev = [-i[2] for i in sgrlistrev]
    # Create figure
    fig = plt.figure(1)
    # Add first subplot on a plot with 1 subplots
    ax = fig.add_subplot(sbnum)
    # Plot both depth coverage for both Watson and Crick strand on the same plot
    print(color1)
    ax.plot(bases,coveragefor,color1)
    #plt.fill_between(bases,coveragefor,color1)
    ax.plot(bases,coveragerev,color2)
    #plt.fill_between(bases,coveragerev,color2)
    # Set limits for x and y axis, x axis limits is the number of bases for that chromosome, y axis max is based on max
    # depth coverage in watson strand and min is based on min depth coverage in crick strand (since negative)
    plt.xlim(xmin=0,xmax=max(bases))
    plt.ylim(ymin=min(coveragerev),ymax=max(coveragefor))
    # Y axis ticks are going to be absolute value, since depth coverage for crick strand is also positive but we just
    # want it below the x axis
    ax.set_yticklabels([str(abs(x)) for x in ax.get_yticks()])
    # X axis ticks shown in kilobases so divide all bases by 1000
    #ax.set_xticklabels([str(x/1000) for x in ax.get_xticks()])
    # We want the x-axis line to be black and obvious
    plt.axhline(xmax=len(bases),color='black')
    plt.ylabel('Fragment Density')
    plt.xlabel('Bases')
    plt.title(strain + '- Chromosome ' + chromosome)