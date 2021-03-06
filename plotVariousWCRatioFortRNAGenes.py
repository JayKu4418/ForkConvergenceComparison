# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 15:53:38 2015

@author: jayashreekumar
"""

import watsoncrickRatio as wcr
import matplotlib.pyplot as plt
import numpy as np

################################## tRNA Genes######################################
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

cdc9deg_centW = wcr.WCRatioMeanForStrand('../TestData/cdc9degCentromeresRawWCRatio-5000.txt','W',[],5000)
rrm3d_centW = wcr.WCRatioMeanForStrand('../TestData/rrm3dCentromeresRawWCRatio-5000.txt','W',[],5000)
pif1m2_centW = wcr.WCRatioMeanForStrand('../TestData/pif1m2CentromeresRawWCRatio-5000.txt','W',[],5000)
rrm3dpif1m2_centW = wcr.WCRatioMeanForStrand('../TestData/rrm3d_pif1m2CentromeresRawWCRatio-5000.txt','W',[],5000)
plt.plot(range(-5000,5000),cdc9deg_centW,label='cdc9deg')
plt.plot(range(-5000,5000),rrm3d_centW,label='rrm3d')
plt.plot(range(-5000,5000),pif1m2_centW,label='pif1m2')
plt.plot(range(-5000,5000),rrm3dpif1m2_centW,label='rrm3d_pif1m2')
plt.xlim(-5000,5000)
plt.xlabel('Bases')
plt.ylabel('Watson Crick Log 2 Ratio')
plt.axhline(xmax=5000,color='black')
plt.axvline(x=0,color='black') 
plt.title('Watson-crick Ratio around Watson Centromeres')
plt.tight_layout()
plt.legend(loc='best')
cdc9deg_centC = wcr.WCRatioMeanForStrand('../TestData/cdc9degCentromeresRawWCRatio-5000.txt','C',[],5000)
rrm3d_centC = wcr.WCRatioMeanForStrand('../TestData/rrm3dCentromeresRawWCRatio-5000.txt','C',[],5000)
pif1m2_centC = wcr.WCRatioMeanForStrand('../TestData/pif1m2CentromeresRawWCRatio-5000.txt','C',[],5000)
rrm3dpif1m2_centC = wcr.WCRatioMeanForStrand('../TestData/rrm3d_pif1m2CentromeresRawWCRatio-5000.txt','C',[],5000)
plt.plot(range(-5000,5000),cdc9deg_centC,label='cdc9deg')
plt.plot(range(-5000,5000),rrm3d_centC,label='rrm3d')
plt.plot(range(-5000,5000),pif1m2_centC,label='pif1m2')
plt.plot(range(-5000,5000),rrm3dpif1m2_centC,label='rrm3d_pif1m2')
plt.xlim(-5000,5000)
plt.xlabel('Bases')
plt.ylabel('Watson Crick Log 2 Ratio')
plt.title('Watson-crick Ratio around Crick Centromeres')
plt.axhline(xmax=5000,color='black')
plt.axvline(x=0,color='black') 
plt.legend(loc='best')
plt.tight_layout()
reload(wcr)
cdc9deg_cent = wcr.WCRatioMeanForStrand('../TestData/cdc9degCentromeresRawWCRatio-5000.txt','',[],5000)
rrm3d_cent = wcr.WCRatioMeanForStrand('../TestData/rrm3dCentromeresRawWCRatio-5000.txt','',[],5000)
pif1m2_cent = wcr.WCRatioMeanForStrand('../TestData/pif1m2CentromeresRawWCRatio-5000.txt','',[],5000)
rrm3dpif1m2_cent = wcr.WCRatioMeanForStrand('../TestData/rrm3d_pif1m2CentromeresRawWCRatio-5000.txt','',[],5000)
plt.clf()
plt.plot(range(-5000,5000),cdc9deg_cent,label='cdc9deg')
plt.plot(range(-5000,5000),rrm3d_cent,label='rrm3d')
plt.plot(range(-5000,5000),pif1m2_cent,label='pif1m2')
plt.plot(range(-5000,5000),rrm3dpif1m2_cent,label='rrm3d_pif1m2')
plt.legend(loc='best')
plt.xlim(-5000,5000)
plt.xlabel('Bases')
plt.ylabel('Watson Crick Log 2 Ratio')
plt.title('Watson-crick Ratio around Centromeres')
plt.axhline(xmax=5000,color='black')
plt.axvline(x=0,color='black') 
plt.tight_layout()

plt.clf()
plt.plot(range(-5000,5000),np.array(rrm3d_cent)-np.array(cdc9deg_cent),label='rrm3d')
plt.plot(range(-5000,5000),np.array(pif1m2_cent)-np.array(cdc9deg_cent),label='pif1m2')
plt.plot(range(-5000,5000),np.array(rrm3dpif1m2_cent)-np.array(cdc9deg_cent),label='rrm3d_pif1m2')
plt.legend(loc='best')
plt.xlim(-5000,5000)
plt.xlabel('Bases')
plt.ylabel('Watson Crick Log 2 Ratio')
plt.title('Watson-crick Ratio around Centromeres')
plt.axhline(xmax=5000,color='black')
plt.axvline(x=0,color='black') 
plt.tight_layout()

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

