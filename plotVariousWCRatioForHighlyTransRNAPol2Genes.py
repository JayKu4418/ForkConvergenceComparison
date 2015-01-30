# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 15:56:45 2015

@author: jayashreekumar
"""

import watsoncrickRatio as wcr
import matplotlib.pyplot as plt
import numpy as np

################################## Highly Tramscribed RNA Pol 2 Genes######################################

convcdc9degW = wcr.WCRatioMeanForStrand('../TestData/cdc9degConvergentWatsonHighlyTransRNAPol2GenesRawWCRatio-5000.txt','W',[],5000)
convrrm3dW = wcr.WCRatioMeanForStrand('../TestData/rrm3dConvergentWatsonHighlyTransRNAPol2GenesRawWCRatio-5000.txt','W',[],5000)
convrrm3dpif1m2W = wcr.WCRatioMeanForStrand('../TestData/rrm3d_pif1m2ConvergentWatsonHighlyTransRNAPol2GenesRawWCRatio-5000.txt','W',[],5000)
convpif1m2W = wcr.WCRatioMeanForStrand('../TestData/pif1m2ConvergentWatsonHighlyTransRNAPol2GenesRawWCRatio-5000.txt','W',[],5000)
divcdc9degW = wcr.WCRatioMeanForStrand('../TestData/cdc9degDivergentWatsonHighlyTransRNAPol2GenesRawWCRatio-5000.txt','W',[],5000)
divrrm3dW = wcr.WCRatioMeanForStrand('../TestData/rrm3dDivergentWatsonHighlyTransRNAPol2GenesRawWCRatio-5000.txt','W',[],5000)
divrrm3dpif1m2W = wcr.WCRatioMeanForStrand('../TestData/rrm3d_pif1m2DivergentWatsonHighlyTransRNAPol2GenesRawWCRatio-5000.txt','W',[],5000)
divpif1m2W = wcr.WCRatioMeanForStrand('../TestData/pif1m2DivergentWatsonHighlyTransRNAPol2GenesRawWCRatio-5000.txt','W',[],5000)
convcdc9degC = wcr.WCRatioMeanForStrand('../TestData/cdc9degConvergentCrickHighlyTransRNAPol2GenesRawWCRatio-5000.txt','C',[],5000)
convrrm3dC = wcr.WCRatioMeanForStrand('../TestData/rrm3dConvergentCrickHighlyTransRNAPol2GenesRawWCRatio-5000.txt','C',[],5000)
convpif1m2C = wcr.WCRatioMeanForStrand('../TestData/pif1m2ConvergentCrickHighlyTransRNAPol2GenesRawWCRatio-5000.txt','C',[],5000)
convrrm3dpif1m2C = wcr.WCRatioMeanForStrand('../TestData/rrm3d_pif1m2ConvergentCrickHighlyTransRNAPol2GenesRawWCRatio-5000.txt','C',[],5000)
divcdc9degC = wcr.WCRatioMeanForStrand('../TestData/cdc9degDivergentCrickHighlyTransRNAPol2GenesRawWCRatio-5000.txt','C',[],5000)
divrrm3dC = wcr.WCRatioMeanForStrand('../TestData/rrm3dDivergentCrickHighlyTransRNAPol2GenesRawWCRatio-5000.txt','C',[],5000)
divpif1m2C = wcr.WCRatioMeanForStrand('../TestData/pif1m2DivergentCrickHighlyTransRNAPol2GenesRawWCRatio-5000.txt','C',[],5000)
divrrm3dpif1m2C = wcr.WCRatioMeanForStrand('../TestData/rrm3d_pif1m2DivergentCrickHighlyTransRNAPol2GenesRawWCRatio-5000.txt','C',[],5000)

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
plt.plot(range(-5000,5000),np.array(divrrm3dC) - np.array(divcdc9degC),label='rrm3d-Crick',color='blue')
plt.plot(range(-5000,5000),np.array(divpif1m2C) - np.array(divcdc9degC),label='pif1m2-Crick',color='crimson')
plt.plot(range(-5000,5000),np.array(divrrm3dpif1m2C) - np.array(divcdc9degC),label='rrm3d_pif1m2-Crick',color='green')

plt.legend(loc='best')
plt.xlim(-5000,5000)
plt.axhline(xmax=5000,color='black')
plt.axvline(x=0,color='black') 
plt.xlabel('Bases')
plt.ylabel('Watson Crick Log 2 Ratio')
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
plt.plot(range(-5000,5000),np.array(convrrm3dC) - np.array(convcdc9degC),label='rrm3d-Crick',color='darkorchid')
plt.plot(range(-5000,5000),np.array(convpif1m2C) - np.array(convcdc9degC),label='pif1m2-Crick',color='darkseagreen')
plt.plot(range(-5000,5000),np.array(convrrm3dpif1m2C) - np.array(convcdc9degC),label='rrm3d_pif1m2-Crick',color='tomato')

plt.legend(loc='best')
plt.xlim(-5000,5000)
plt.axhline(xmax=5000,color='black')
plt.axvline(x=0,color='black') 
plt.xlabel('Bases')
plt.ylabel('Watson Crick Log 2 Ratio')
plt.tight_layout()  