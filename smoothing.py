__author__ = 'jayashreekumar'

# A: finite ordered data set A consisting of pairs of real numbers (x,a) where in our case x is the chromosomal coordinate and a is the depth coverage
# T = # of pairs of points in A
# A+ = ordered set of a values in A
# If L is an ordered list of numbers (NOT pairs), then RotateLeft[L,n] = cycle elements in list L n positions to the left
# Sum[L] is sum of all numbers in list L
#######################################################################################################################
# IMPORT
import matplotlib.pyplot as plt
import numpy as np
import numpy.fft as ft
import pyfftw
from operator import itemgetter
from scipy.stats import pearsonr
import ManipulateSeqData.oem as oem
import ManipulateSeqData.importantfunctions as imptfunc
import ManipulateSeqData.originscomparison as oc
yeastsize = {'1':230218,'2':813184,'3':316620,'4':1531933,'5':576874,'6':270161,'7':1090940,'8':562643,'9':439888,'10':745751,'11':666816,'12':1078177,'13':924431,'14':784333,'15':1091291,'16':948066,'Total':12157105}
#######################################################################################################################
# FUNCTIONS


# This function splits up a list of numbers based on whether they are consecutively positive and negative. For example
# list of numbers like [0,0,2,3,5,-6,-7,0,8,9,-12,-11,0,0], it will split up as [2,3,5], [-6,-7],[8,9],[-12,-11].
def splitListbasedOnSignWithoutIndex(list2split):
    # Empty list to store all the split lists which are dictionaries
    splitlist = []
    # Empty list to store each individual as long as sign does not change between numbers 
    splitlistind = []
    # Go through every index in list2split
    for i in range(len(list2split)-1):
        # If the signs of the current value at index and the value at index + 1 is equal and both values are not zero
        if cmp(list2split[i],0) == cmp(list2split[i+1],0) and cmp(list2split[i],0) != 0 and cmp(list2split[i+1],0) != 0:
            # Let index be the key of dictionary and the value at that index is stored as the paired value
            splitlistind.append(list2split[i])
            if i==(len(list2split)-2):
                splitlistind.append(list2split[i+1])
                splitlist.append(splitlistind)
                splitlistind = []
                
        # If the signs of the current value at index and the value at index + 1 are not equal or either values are zero
        else:
            # If value at current index is not zero
            
            if list2split[i] != 0:
                # Add this index and value as key-value pairs to the same dictionary as previously because it would
                # have been equal in sign to previous value, it's not just equal in sign to the next value
                splitlistind.append(list2split[i])
            # If length of current dictionary is not zero
            if len(splitlistind)!=0:
                # Add this dictionary which contains values of same sign to the big list and reset the dictionary to be
                # empty
                splitlist.append(splitlistind)
                splitlistind = []
    # return this big list of dictionaries
    return splitlist

# Low pass filter -> only lower frequencies go through,high frequencies thrown out
def lowpassusingnp(sgrfile,chromosome,nfreq):
    with open(sgrfile) as f:
        sgrforchrom = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
    
    # Fourier transform on sgrforchrom
    f_sgrforchrom = ft.rfft(sgrforchrom)
    
    # Zero coefficients that from nfreq index and up
    f_sgrforchrom[nfreq:] = 0
    
    # Inverse fourier transform to get back smoothed signal
    sgrforchrom_smooth = ft.irfft(f_sgrforchrom)
    
    fig = plt.figure(1)
    
    ax1 = fig.add_subplot(211)
    
    ax1.plot(sgrforchrom,color='r')
    
    plt.title('Original')

    ax2 = fig.add_subplot(212)
    
    ax2.plot(sgrforchrom_smooth,color='g')
    
    plt.title('Smoothed')
        
    plt.show()
    
# Low pass filter -> only lower frequencies go through,high frequencies thrown out
def lowpassusingfftw(sgrfileW,sgrfileC,chromosome,nfreq1,num):
    with open(sgrfileW) as f:
        sgrforchromW = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
    sgrforchromW = np.array(sgrforchromW)
    with open(sgrfileC) as f:
        sgrforchromC = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
    sgrforchromC = np.array(sgrforchromC)
    
    # Fourier transform on sgrforchrom
    f_sgrforchromW = pyfftw.interfaces.numpy_fft.rfft(sgrforchromW[5000:(len(sgrforchromW)-5000)])
    f_sgrforchromC = pyfftw.interfaces.numpy_fft.rfft(sgrforchromC[5000:(len(sgrforchromC)-5000)])    
    #f_sgrforchrom2 = pyfftw.interfaces.numpy_fft.rfft(sgrforchrom)
    #f_sgrforchrom3 = pyfftw.interfaces.numpy_fft.rfft(sgrforchrom)
    #f_sgrforchrom4 = pyfftw.interfaces.numpy_fft.rfft(sgrforchrom)
    # Zero coefficients that from nfreq index and up
    f_sgrforchromW[nfreq1:] = 0
    f_sgrforchromC[nfreq1:] = 0
    #f_sgrforchrom2[nfreq2:] = 0
    #f_sgrforchrom3[nfreq3:] = 0
    #f_sgrforchrom4[nfreq4:] = 0
    # Inverse fourier transform to get back smoothed signal
    sgrforchrom_smoothW = pyfftw.interfaces.numpy_fft.irfft(f_sgrforchromW)
    sgrforchrom_smoothC = pyfftw.interfaces.numpy_fft.irfft(f_sgrforchromC)
    #sgrforchrom_smooth2 = pyfftw.interfaces.numpy_fft.irfft(f_sgrforchrom2)
    #sgrforchrom_smooth3 = pyfftw.interfaces.numpy_fft.irfft(f_sgrforchrom3)
    #sgrforchrom_smooth4 = pyfftw.interfaces.numpy_fft.irfft(f_sgrforchrom4)
    
    fig = plt.figure(num)
    
    ax1 = fig.add_subplot(211)
    
    ax1.plot(sgrforchromW,color='r')
    ax1.plot(-sgrforchromC,color='g')
    plt.axis('tight')
    #plt.title('Original')

    ax2 = fig.add_subplot(212)
    
    ax2.plot(sgrforchrom_smoothW,color='b')
    ax2.plot(-sgrforchrom_smoothC,color='m')
    plt.axis('tight')    
    #plt.title('Smoothed cutoff= ' + str(nfreq1))
    
    #ax2 = fig.add_subplot(513)
    
    #ax2.plot(sgrforchrom_smooth2,color='m')
    
    #plt.title('Smoothed cutoff= ' + str(nfreq2))
    
    #ax2 = fig.add_subplot(514)
    
    #ax2.plot(sgrforchrom_smooth3,color='b')
    
    #plt.title('Smoothed cutoff= ' + str(nfreq3))
    
    #ax2 = fig.add_subplot(515)
    
    #ax2.plot(sgrforchrom_smooth4,color='y')
    
    #plt.title('Smoothed cutoff= ' + str(nfreq4))
    plt.axis('tight')
    plt.show()

def writefilelowpassusingfftw(sgrfile,chromosome,nfreq,writesgrfileforchrom,writewigfileforchrom,iscrick):
    with open(sgrfile) as f:
        sgrforchrom = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
    
    sgrforchrom = np.array(sgrforchrom)

    f_sgrforchrom = pyfftw.interfaces.numpy_fft.rfft(sgrforchrom)
    f_sgrforchrom[nfreq:] = 0
    sgrforchrom_smooth = pyfftw.interfaces.numpy_fft.irfft(f_sgrforchrom)
    with open(writesgrfileforchrom,'w') as fw:
        for i in range(1,len(sgrforchrom_smooth)+1):
            fw.write(chromosome+'\t'+str(i)+'\t'+str(sgrforchrom_smooth[i-1])+'\n')
    imptfunc.convertToWigFile(writesgrfileforchrom,writewigfileforchrom,iscrick)
            
def splitListbasedOnContinous(list2split):
    # Empty list to store all the split lists which are dictionaries
    splitlist = []
    # Empty dictionary to store each individual as long as sign does not change between numbers -> used a dictionary
    # because we needed to access the base for each corresponding OEM easily
    splitlistind = []
    # Go through every index in list2split
    for i in range(len(list2split)):
        # If the signs of the current value at index and the value at index + 1 is equal and both values are not zero
        if i!=len(list2split)-1 and int(list2split[i][1])+1 == int(list2split[i+1][1]):
            # Let index be the key of dictionary and the value at that index is stored as the paired value
            splitlistind.append(list2split[i])
        # If the signs of the current value at index and the value at index + 1 are not equal or either values are zero
        else:
            # Add this index and value as key-value pairs to the same dictionary as previously because it would
            # have been equal in sign to previous value, it's not just equal in sign to the next value
            splitlistind.append(list2split[i])
            # If length of current dictionary is not zero
            if len(splitlistind)!=0:
                # Add this dictionary which contains values of same sign to the big list and reset the dictionary to be
                # empty
                splitlist.append(splitlistind)
                splitlistind = []
    # return this big list of dictionaries
    return splitlist
    
def linearinterpolatemissingdata(x,sgrlist,chromosome,window,oemfile):
    basesinx = [int(i[1]) for i in x]
    sgrlistofinterest = [i for i in sgrlist if int(i[1]) not in basesinx and (int(i[1])<=(basesinx[len(basesinx)-1]+window) or int(i[1])>=(basesinx[0]-window))]
    valsofsgrlistofinterest = [float(i[2]) for i in sgrlistofinterest]
    basesofsgrlistofinterest = [int(i[1]) for i in sgrlistofinterest]
    
    return np.interp(basesinx,basesofsgrlistofinterest,valsofsgrlistofinterest)

def calculateAndPLotOEMforFFTsmoothedsgr(sgrfileW,sgrfileC,nfreq,chromosome,oemsmoothedwritefile,oemfile):
    with open(sgrfileW) as f:
        sgrforchromW = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
    sgrforchromW = np.array(sgrforchromW)
    print(len(sgrforchromW))
    # Fourier transform on sgrforchrom
    f_sgrforchromW = pyfftw.interfaces.numpy_fft.rfft(sgrforchromW)
    f_sgrforchromW[nfreq:] = 0
    # Inverse fourier transform to get back smoothed signal
    sgrforchrom_smoothW = pyfftw.interfaces.numpy_fft.irfft(f_sgrforchromW)
    
    with open(sgrfileC) as f:
        sgrforchromC = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
    sgrforchromC = np.array(sgrforchromC)
    print(len(sgrforchromC))
    # Fourier transform on sgrforchrom
    f_sgrforchromC = pyfftw.interfaces.numpy_fft.rfft(sgrforchromC)
    f_sgrforchromC[nfreq:] = 0
    # Inverse fourier transform to get back smoothed signal
    sgrforchrom_smoothC = pyfftw.interfaces.numpy_fft.irfft(f_sgrforchromC)
    print(len(sgrforchrom_smoothW))
    print(len(sgrforchrom_smoothC))
    sgrcomblist = []
    for i in range(1,len(sgrforchrom_smoothW)+1):
        sgrcomblist.append([chromosome,i,sgrforchrom_smoothW[i-1],sgrforchrom_smoothC[i-1]])
    print('Before OEM calc')
    x = oem.calculateOEM(sgrcomblist)
    with open(oemsmoothedwritefile,'w') as fw:
        for i in x:
            fw.write(chromosome + '\t' + str(i[0]) + '\t' + str(i[1]) + '\n')
    print('After OEM calc')    
    xoem = [i[1] for i in x]
    print(len(xoem))
    with open(oemfile) as f:
        oemforchrom = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
    print(len(oemforchrom))
    fig = plt.figure(1)
    
    ax1 = fig.add_subplot(211)
    
    ax1.plot(range(1,len(oemforchrom)+1),oemforchrom,color='r')
    
    plt.fill_between(range(1,len(oemforchrom)+1),oemforchrom,color='r')
    
    plt.title('Original')

    ax2 = fig.add_subplot(212)
    
    ax2.plot(range(1,len(xoem)+1),xoem,color='g')
    
    plt.fill_between(range(1,len(xoem)+1),xoem,color='g')
    
    plt.title('Smoothed')
    
    #ax3 = fig.add_subplot(313)
    
    #ax3.plot(np.array(oemforchrom)/np.array(xoem))
    
    plt.show()
    
def plotRatioOEMs(originaloemfile,smoothedoemfile,chromosome,fignum):
    with open(originaloemfile) as f:
        org_oem = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
    with open(smoothedoemfile) as f:
        sm_oem = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
    numoem = min(len(org_oem),len(sm_oem))
    
    ratio = []
    for i in range(numoem):
        ratio.append(sm_oem[i]-org_oem[i])
    plt.figure(fignum)
    plt.clf()
    plt.plot(ratio)
    plt.axis('tight')
    plt.show()
        

def ratioWatsonToCrickForChrom(watsonsgr,cricksgr,chromosome,inverse):
    with open(watsonsgr) as f:
        W = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
    with open(cricksgr) as f:
        C = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
    
    rat = []        
    lenofchrom = len(W)
    for i in range(lenofchrom):
        if inverse:
            if W[i]!=0:
                rat.append(C[i]/W[i])
            else:
                rat.append(C[i])
        else:
            if C[i]!=0:
                rat.append(W[i]/C[i])
            else:
                rat.append(W[i])
            
    return rat
    
    
def smoothTroughs(oemfile,chromosome,startpoint,endpoint,freq,colorg,colsmooth,strain,sbnum):
    with open(oemfile) as f:
        data = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
        
    areaofinterest = data[int(startpoint/50):int(endpoint/50)]
    
    data_fftw = pyfftw.interfaces.numpy_fft.rfft(np.array(areaofinterest))
    
    data_fftw[freq:]=0
    
    data_smooth = pyfftw.interfaces.numpy_fft.irfft(data_fftw)
    
    #plt.plot(areaofinterest,colorg,label='Original-'+strain)
    fig = plt.figure(1)
    # Add a subplot
    ax = fig.add_subplot(sbnum)
    ax.plot(data_smooth,colsmooth,label='Smoothed-'+strain)
    plt.axhline(xmax=len(data_smooth),color='black')
    plt.legend(loc='best')

def smoothTroughsFor2Strains(oemfile1,startpoint1,endpoint1,strain1,oemfile2,startpoint2,endpoint2,strain2,chromosome,freq):
    with open(oemfile1) as f:
        data1 = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
        
    areaofinterest1 = data1[int(startpoint1/50):int(endpoint1/50)]
    
    data_fftw1 = pyfftw.interfaces.numpy_fft.rfft(np.array(areaofinterest1))
    
    data_fftw1[freq:]=0
    
    ind1 = range(startpoint1,endpoint1,50)
    print len(ind1)
    data_smooth1 = pyfftw.interfaces.numpy_fft.irfft(data_fftw1) 
    print len(data_smooth1)
    with open(oemfile2) as f:
        data2 = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
        
    areaofinterest2 = data2[int(startpoint2/50):int(endpoint2/50)]
    
    data_fftw2 = pyfftw.interfaces.numpy_fft.rfft(np.array(areaofinterest2))
    
    data_fftw2[freq:]=0
    ind2 = range(startpoint2,endpoint2,50)
    data_smooth2 = pyfftw.interfaces.numpy_fft.irfft(data_fftw2)
    print len(ind2)
    print len(data_smooth2)
    #plt.plot(areaofinterest,colorg,label='Original-'+strain)
    #fig = plt.figure(1)
    # Add a subplot
    #ax = fig.add_subplot(sbnum)
    plt.plot(ind1[0:len(data_smooth1)],data_smooth1,'.-g',label='Smoothed-'+strain1)
    plt.plot(ind2[0:len(data_smooth2)],data_smooth2,'.-m',label='Smoothed-'+strain2)
    plt.axhline(xmax=max(len(data_smooth2),len(data_smooth1)),color='black')
    plt.axis('tight')
    plt.title('Chromosome ' + chromosome)
    plt.legend(loc='best')

def basesWmaxminoem(splitlist):
    # Empty list to store all the bases with max positive OEM values
    maxminoembases = []
    # Go through every dictionary in the split list
    for dictoem in splitlist:
        # If first value in dictionary is positive, means all the values in that dictionary are positive since they
        # should all have the same sign
        if dictoem.values()[0] > 0:
            # Multiply each the base with the max OEM value by perbasebox since OEMs were calculated for every perbasebox value
            maxminoembases.append(oc.findmaxInDict(dictoem))
        else:
            maxminoembases.append(oc.findminInDict(dictoem))

    return maxminoembases

def propertiesOfTrough(oemfile,startpoint,endpoint,strain,chromosome,freq,col):
    with open(oemfile) as f:
        data = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
        
    areaofinterest = data[int(startpoint/50):int(endpoint/50)]
    
    data_fftw = pyfftw.interfaces.numpy_fft.rfft(np.array(areaofinterest))
    
    data_fftw[freq:]=0
    
    data_smooth = pyfftw.interfaces.numpy_fft.irfft(data_fftw)
    
    ind = range(startpoint,endpoint,50)
    data_with_ind = [[ind[i],data_smooth[i]] for i in range(len(data_smooth))]
    
    splitdata = oc.splitListbasedOnSign(data_with_ind)
    
    maxminoem = basesWmaxminoem(splitdata)    
    
    oems = [abs(i[1]) for i in maxminoem]
    bases = sorted([i[0] for i in maxminoem])
    avgamplitude = float(sum(oems))/3
    midpoint = float(bases[1]-bases[0])/(bases[2]-bases[0])
    
    #data_grad_before = [i[1] for i in data_with_ind if i[0]<bases[0]]
    data_grad = [i[1] for i in data_with_ind if i[0]>=bases[0] and i[0]<=bases[2]]
    
    ind_grad = [i[0] for i in data_with_ind if i[0]>=bases[0] and i[0]<=bases[2]]
    #print len(ind_grad)    
    #data_grad_after = [i[1] for i in data_with_ind if i[0]>bases[2]]
    grad = []
    for i in range(len(data_grad)-1):
        grad.append((data_grad[i+1] - data_grad[i])*100)
    #print(len(grad))
    #print len(range(len(data_grad_before)+1,len(data_grad_before)+len(grad)))
    #plt.plot(range(len(data_grad_before)+1,len(data_grad_before)+len(grad)+1,50),grad,col,label='Gradient-'+strain)
    plt.plot(ind_grad[0:len(grad)],grad,col,label='Gradient-'+strain)
    plt.axhline(xmax=len(data_smooth),color='black')
    
    return [strain+" avgamp:"+str(avgamplitude),strain+" midpoint:"+str(midpoint)]

def gradientOfStrain(oemfile,startpoint,endpoint,chromosome,freq):
    with open(oemfile) as f:
        data = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
        
    areaofinterest = data[int(startpoint/50):int(endpoint/50)]
    
    data_fftw = pyfftw.interfaces.numpy_fft.rfft(np.array(areaofinterest))
    
    data_fftw[freq:]=0
    
    data_smooth = pyfftw.interfaces.numpy_fft.irfft(data_fftw)
    
    ind = range(startpoint,endpoint,50)
    data_with_ind = [[ind[i],data_smooth[i]] for i in range(len(data_smooth))]
    
    splitdata = oc.splitListbasedOnSign(data_with_ind)
    
    maxminoem = basesWmaxminoem(splitdata)    
    
    #oems = [abs(i[1]) for i in maxminoem]
    bases = sorted([i[0] for i in maxminoem])
    #avgamplitude = float(sum(oems))/3
    #midpoint = float(bases[1]-bases[0])/(bases[2]-bases[0])
    
    #data_grad_before = [i[1] for i in data_with_ind if i[0]<bases[0]]
    data_grad = [i[1] for i in data_with_ind if i[0]>=bases[0] and i[0]<=bases[2]]
    #data_grad_after = [i[1] for i in data_with_ind if i[0]>bases[2]]
    grad = []
    for i in range(len(data_grad)-1):
        grad.append((data_grad[i+1] - data_grad[i])*100)
        
    return grad

def numberOfTimesGradientZero(oemfile1,oemfile2,troughfile,troughnumber,freq):
    troughinfo = extractstartandendfromtroughfile(troughfile)
    grad1 = gradientOfStrain(oemfile1,troughinfo[troughnumber-1][1],troughinfo[troughnumber-1][2],troughinfo[troughnumber-1][0],freq)
    grad2 = gradientOfStrain(oemfile2,troughinfo[troughnumber-1][3],troughinfo[troughnumber-1][4],troughinfo[troughnumber-1][0],freq)
    return {'1':len(splitListbasedOnSignWithoutIndex(grad1)),'2':len(splitListbasedOnSignWithoutIndex(grad2))}
    
def numberOfTimesGradientZeroAllTroughs(oemfile1,oemfile2,troughfile,writeGradZeroTroughFile,freq):
    troughinfo = extractstartandendfromtroughfile(troughfile)
    with open(writeGradZeroTroughFile,'w') as fw:
        for i in range(1,len(troughinfo)+1):
            splits = numberOfTimesGradientZero(oemfile1,oemfile2,troughfile,i,freq)
            fw.write(str(splits['1']) + '\t' + str(splits['2']) + '\n')
        
def gradientOf2Strains(oemfile1,strain1,oemfile2,strain2,troughfile,troughnumber,freq):
    troughinfo = extractstartandendfromtroughfile(troughfile)
    grad1 = gradientOfStrain(oemfile1,troughinfo[troughnumber-1][1],troughinfo[troughnumber-1][2],troughinfo[troughnumber-1][0],freq)
    grad2 = gradientOfStrain(oemfile2,troughinfo[troughnumber-1][3],troughinfo[troughnumber-1][4],troughinfo[troughnumber-1][0],freq)
    lenofdatatograb= min(len(grad1),len(grad2))
    #return {"Strain1":grad1[0:lenofdatatograb],"Strain2":grad2[0:lenofdatatograb]}
    return pearsonr(grad1[0:lenofdatatograb],grad2[0:lenofdatatograb])
    #plt.plot(grad1,'.-r',label='Gradient-'+strain1)
    #plt.plot(grad2,'.-b',label='Gradient-'+strain2)
def originsAndTermintionOfStrain(oemfile,startpoint,endpoint,chromosome,freq,col):
    with open(oemfile) as f:
        data = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
        
    areaofinterest = data[int(startpoint/50):int(endpoint/50)]
    
    data_fftw = pyfftw.interfaces.numpy_fft.rfft(np.array(areaofinterest))
    
    data_fftw[freq:]=0
    
    data_smooth = pyfftw.interfaces.numpy_fft.irfft(data_fftw)
    
    ind = range(startpoint,endpoint,50)
    plt.plot(ind,data_smooth,col)
    data_with_ind = [[ind[i],data_smooth[i]] for i in range(len(data_smooth))]
    
    splitdata = oc.splitListbasedOnSign(data_with_ind)
    
    maxminoem = basesWmaxminoem(splitdata)
    print(maxminoem)
    positiveoembases = sorted([i for i in maxminoem if i[1]>0],key=itemgetter(0))
    
    origins = [positiveoembases[0],positiveoembases[len(positiveoembases)-1]]
    
    print(origins)
    negativeoembases = [i for i in maxminoem if i[1]<0]
    if len(negativeoembases) > 1:
        termination = min(negativeoembases.iteritems(), key=itemgetter(1))
    else:
        termination = negativeoembases[0]
    
    return {'o':origins,'t':termination}
    
def originsAndTerminationOfTwoStrains(oemfile1,strain1,oemfile2,strain2,troughfile,troughnumber,freq):
    troughinfo = extractstartandendfromtroughfile(troughfile)
    vals1 = originsAndTermintionOfStrain(oemfile1,troughinfo[troughnumber-1][1],troughinfo[troughnumber-1][2],troughinfo[troughnumber-1][0],freq)
    vals2 = originsAndTermintionOfStrain(oemfile2,troughinfo[troughnumber-1][3],troughinfo[troughnumber-1][4],troughinfo[troughnumber-1][0],freq)
    return [vals1,vals2]
    
def originsAndTerminationOfTwoStrainsForAllTroughs(oemfile1,strain1,oemfile2,strain2,troughfile,freq,writeOrgAndTermOfTroughFile):
    troughinfo = extractstartandendfromtroughfile(troughfile)
    with open(writeOrgAndTermOfTroughFile,'w') as fw:
        for i in range(1,len(troughinfo)+1):
            vals = originsAndTerminationOfTwoStrains(oemfile1,strain1,oemfile2,strain2,troughfile,i,freq)
            fw.write(str(vals[0]["o"][0][0])+'\t'+str(vals[0]["o"][0][1])+'\t'+str(vals[0]["t"][0])+'\t'+str(vals[0]["t"][1])+'\t'+str(vals[0]["o"][1][0])+'\t'+str(vals[0]["o"][1][1])+'\t'+str(vals[1]["o"][0][0])+'\t'+str(vals[1]["o"][0][1])+'\t'+str(vals[1]["t"][0])+'\t'+str(vals[1]["t"][1])+'\t'+str(vals[1]["o"][1][0])+'\t'+str(vals[1]["o"][1][1])+'\n')

def midpointOfStrain(oemfile,startpoint,endpoint,chromosome,freq):
    with open(oemfile) as f:
        data = [float(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome]
        
    areaofinterest = data[int(startpoint/50):int(endpoint/50)]
    
    data_fftw = pyfftw.interfaces.numpy_fft.rfft(np.array(areaofinterest))
    
    data_fftw[freq:]=0
    
    data_smooth = pyfftw.interfaces.numpy_fft.irfft(data_fftw)
    
    ind = range(startpoint,endpoint,50)
    data_with_ind = [[ind[i],data_smooth[i]] for i in range(len(data_smooth))]
    
    splitdata = oc.splitListbasedOnSign(data_with_ind)
    
    maxminoem = basesWmaxminoem(splitdata)    
    
    #oems = [abs(i[1]) for i in maxminoem]
    bases = sorted([i[0] for i in maxminoem])
    #avgamplitude = float(sum(oems))/3
    midpoint = float(bases[1]-bases[0])/(bases[2]-bases[0])
    #return {"avgamplitude":avgamplitude,"midpoint":midpoint}
    return midpoint
def midpointOf2Strains(oemfile1,strain1,oemfile2,strain2,troughfile,troughnumber,freq):
    troughinfo = extractstartandendfromtroughfile(troughfile)
    vals1 = midpointOfStrain(oemfile1,troughinfo[troughnumber-1][1],troughinfo[troughnumber-1][2],troughinfo[troughnumber-1][0],freq)
    vals2 = midpointOfStrain(oemfile2,troughinfo[troughnumber-1][3],troughinfo[troughnumber-1][4],troughinfo[troughnumber-1][0],freq)
    return [vals1,vals2]
    
def midpointOf2StrainsForAllTroughs(oemfile1,strain1,oemfile2,strain2,troughfile,freq,writepropertiesOfTroughFile):
    troughinfo = extractstartandendfromtroughfile(troughfile)
    with open(writepropertiesOfTroughFile,'w') as fw:
        for i in range(1,len(troughinfo)+1):
            vals =  midpointOf2Strains(oemfile1,strain1,oemfile2,strain2,troughfile,i,freq)
            fw.write(str(vals[0]["avgamplitude"])+'\t'+str(vals[1]["avgamplitude"])+'\t'+str(vals[0]["midpoint"])+'\t'+str(vals[1]["midpoint"])+'\n')
def gradientOf2StrainsForAllTroughs(oemfile1,strain1,oemfile2,strain2,troughfile,freq,writepearsoncorrforgradfile):
    troughinfo = extractstartandendfromtroughfile(troughfile)
    with open(writepearsoncorrforgradfile,'w') as fw:
        for i in range(1,len(troughinfo)+1):
            rfortrough = gradientOf2Strains(oemfile1,strain1,oemfile2,strain2,troughfile,i,freq)
            fw.write(str(rfortrough[0])+'\t'+str(rfortrough[1])+'\n')
        
def extractstartandendfromtroughfile(troughfile):
    with open(troughfile) as f:
        alltroughs = [line.strip().split('\t') for line in f]
    
    troughinfo = [[i[1],int(i[4]),int(i[5]),int(i[9]),int(i[10])] for i in alltroughs]
    
    return troughinfo
    
def propertiesOfTroughTwoStrains(oemfile1,strain1,oemfile2,strain2,troughfile,troughnumber,freq):
    
    troughinfo = extractstartandendfromtroughfile(troughfile)
    
    propertiesOfTrough(oemfile1,troughinfo[troughnumber-1][1],troughinfo[troughnumber-1][2],strain1,troughinfo[troughnumber-1][0],freq,'.-r')
    propertiesOfTrough(oemfile2,troughinfo[troughnumber-1][3],troughinfo[troughnumber-1][4],strain2,troughinfo[troughnumber-1][0],freq,'.-b')
    smoothTroughsFor2Strains(oemfile1,troughinfo[troughnumber-1][1],troughinfo[troughnumber-1][2],strain1,oemfile2,troughinfo[troughnumber-1][3],troughinfo[troughnumber-1][4],strain2,troughinfo[troughnumber-1][0],freq)
    plt.legend(loc='best')
    
def differenceOriginPos(originsandtroughfile,writedifferencefile):
    with open(originsandtroughfile) as f:
        alltroughs = [line.strip().split('\t') for line in f]
        
    with open(writedifferencefile,'w') as fw:
        for i in alltroughs:        
            fw.write(str(int(i[6])-int(i[0])) + '\t' + str(int(i[10])-int(i[4])) + '\n')