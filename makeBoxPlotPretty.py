# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 10:57:23 2015

@author: jayashreekumar
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from scipy.stats import nanmean
#import wcratiospecificregion as wcs

# Generate some data from five different probability distributions,
# each with different characteristics. We want to play with how an IID
# bootstrap resample of the data preserves the distributional
# properties of the original sample, and a boxplot is one visual tool
# to make this assessment
def createBoxPlot(data,xlabels,title,xlab,ylab,top1,bottom1,top2):
    #randomDists = ['500','1000', '1500', '2000','2500','3000','3500','4000','4500','5000']
#N = 500
#norm = np.random.normal(1,1, N)
#logn = np.random.lognormal(1,1, N)
#expo = np.random.exponential(1, N)
#gumb = np.random.gumbel(6, 4, N)
#tria = np.random.triangular(2, 9, 11, N)

# Generate some random indices that we'll use to resample the original data
# arrays. For code brevity, just use the same random indices for each array
#bootstrapIndices = np.random.random_integers(0, N-1, N)
#normBoot = norm[bootstrapIndices]
#expoBoot = expo[bootstrapIndices]
#gumbBoot = gumb[bootstrapIndices]
#lognBoot = logn[bootstrapIndices]
#triaBoot = tria[bootstrapIndices]

    #data = [norm, normBoot,  logn, lognBoot, expo, expoBoot, gumb, gumbBoot,tria, triaBoot]
    fig, ax1 = plt.subplots(figsize=(10,6))
    fig.canvas.set_window_title('Boxplot')
    plt.tight_layout()    
    #plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

    bp = plt.boxplot(data, notch=0, sym='+', vert=1, whis=1.5)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='red', marker='+')

    # Add a horizontal grid to the plot, but make it very light in color
    # so we can use it for reading data values but not be distracting
    ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
              alpha=0.5)

    # Hide these grid behind plot objects
    ax1.set_axisbelow(True)
    ax1.set_title(title)
    ax1.set_xlabel(xlab)
    ax1.set_ylabel(ylab)

    # Now fill the boxes with desired colors
    boxColors = ['mediumvioletred','darkseagreen']
    numBoxes = len(data)
    medians = range(numBoxes)
    for i in range(numBoxes):
        box = bp['boxes'][i]
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        boxCoords = zip(boxX,boxY)
        # Alternate between Dark Khaki and Royal Blue
        k = i % 2
        boxPolygon = Polygon(boxCoords, facecolor=boxColors[k])
        ax1.add_patch(boxPolygon)
        # Now draw the median lines back over what we just filled in
        med = bp['medians'][i]
        medianX = []
        medianY = []
        for j in range(2):
            medianX.append(med.get_xdata()[j])
            medianY.append(med.get_ydata()[j])
            plt.plot(medianX, medianY, 'k')
            medians[i] = medianY[0]
        # Finally, overplot the sample averages, with horizontal alignment
        # in the center of each box
        plt.plot([nanmean(med.get_xdata())], [nanmean(data[i])],
           color='w', marker='*', markeredgecolor='k')
        
    # Set the axes ranges and axes labels
    ax1.set_xlim(0.5, numBoxes+0.5)
    top = top1
    bottom = bottom1
    ax1.set_ylim(bottom, top)
    xtickNames = plt.setp(ax1, xticklabels=xlabels)
    plt.setp(xtickNames, rotation=45, fontsize=12)

    # Due to the Y-axis scale being different across samples, it can be
    # hard to compare differences in means across the samples. Add upper
    # X-axis tick labels with the sample means to aid in comparison
    # (just use two decimal places of precision)
    pos = np.arange(numBoxes)+1
    upperLabels = [str(np.round(nanmean(s), 2)) for s in data]
    weights = ['bold', 'semibold']
    for tick,label in zip(range(numBoxes),ax1.get_xticklabels()):
        k = tick % 2
        ax1.text(pos[tick], top2, upperLabels[tick],
        horizontalalignment='center', size=12, weight=weights[k],
        color='black')
    #plt.savefig(savefigname)
    plt.show()
"""
# Finally, add a basic legend
plt.figtext(0.80, 0.08,  str(N) + ' Random Numbers' ,
           backgroundcolor=boxColors[0], color='black', weight='roman',
           size='x-small')
plt.figtext(0.80, 0.045, 'IID Bootstrap Resample',
backgroundcolor=boxColors[1],
           color='white', weight='roman', size='x-small')
plt.figtext(0.80, 0.015, '*', color='white', backgroundcolor='silver',
           weight='roman', size='medium')
plt.figtext(0.815, 0.013, ' Average Value', color='black', weight='roman',
           size='x-small')
"""
    