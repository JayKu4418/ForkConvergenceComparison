# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 10:40:40 2014

@author: jayashreekumar
"""

# The functions in this script will help find CAG repeats based on how many
# repeats you want to find

#######################################################################################################################
# IMPORT
from Bio import SeqIO
import re
# dictionary that converts values from the form 'chrmIII' to '3'
chrmconvert = {'chrI':'1','chrII':'2','chrIII':'3','chrIV':'4','chrV':'5','chrVI':'6','chrVII':'7','chrVIII':'8','chrIX':'9','chrX':'10','chrXI':'11','chrXII':'12','chrXIII':'13','chrXIV':'14','chrXV':'15','chrXVI':'16'}

#######################################################################################################################
# FUNCTIONS

# This function reads in a fasta file and spits out the chromosome sequences
# of interest
def extractChromosomeSeq(fastafile,chromosome):
    for element in chrmconvert.keys():
        if chrmconvert[element] == chromosome:
            chrm = element
    rec = [record for record in SeqIO.parse(open(fastafile,'rU'),'fasta') if record.name == chrm]
    return str(rec[0].seq)
    
# This function finds indexes for where CAG repeats starts. The parameters 
# for this function, seq, the number of CAG repeats to find
#def extractIndexCAGrepeats(seq,numrepeats):
    #regexexpr = "(?<!(CAG))(CAG){" + str(numrepeats) +"}(?!(CAG))"
    #regs = re.finditer(regexexpr,seq)
    #return [[i.start(0),i.end(0)] for i in regs]
    

# This function combines the above two functions were it finds the CAG repeats 
# for a fastafile for a chromosome
#def getCAGRepeatsForChromsome(fastafile,chromosome,numrepeats):
    #seq = extractChromosomeSeq(fastafile,chromosome)
    #return extractIndexCAGrepeats(seq,numrepeats)

# This function finds indexes for any three letter repeats. The parametes 
# for this function, seq, the three letter repeat, number of three letter repeats
def extractIndexrepeats(seq,repeat,numrepeats):
    regexexpr = "(?<!(" + repeat + "))(" + repeat + "){" + str(numrepeats) +"}(?!(" + repeat + "))"
    regs = re.finditer(regexexpr,seq)
    return [[i.start(0),i.end(0)] for i in regs]

# This function combines the above two functions where it finds the trinucleotide repeats 
# for a fastafile for a chromosome
def getRepeatsForChromsome(fastafile,chromosome,repeat,numrepeats):
    seq = extractChromosomeSeq(fastafile,chromosome)
    return extractIndexrepeats(seq,repeat,numrepeats)

# This function gives the total lengths of trinucleotide repeats found for 
# each chromosome
def extractLengthsTrinucleotideRepeats(fastafile,chromosome,repeat,retlength):
    seq = extractChromosomeSeq(fastafile,chromosome)
    regexexpr = "(" + repeat + "){3,}"    
    regs = re.finditer(regexexpr,seq)
    if retlength:
        return [i.end(0)-i.start(0) for i in regs]
    else:
        return [[i.start(0),i.end(0)] for i in regs]

# This function returns the max number of trinucleotide repeats found for all chromosomes
def maxLengthTrinucleotideRepeat(fastafile,repeat):
    maxlengtheachchrom = []
    for c in range(1,17):
        chromosome = str(c)
        lengthsPerChrom = extractLengthsTrinucleotideRepeats(fastafile,chromosome,repeat,True)
        if len(lengthsPerChrom)!=0:
            maxlengtheachchrom.append(max(lengthsPerChrom))
    return maxlengtheachchrom