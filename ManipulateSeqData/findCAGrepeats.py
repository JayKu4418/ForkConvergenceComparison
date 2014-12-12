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
            print chrm
    rec = [record for record in SeqIO.parse(open(fastafile,'rU'),'fasta') if record.name == chrm]
    return str(rec[0].seq)
    
# This function finds indexes for where CAG repeats starts. The parameters 
# for this function, seq, the number of CAG repeats to find
def extractIndexCAGrepeats(seq,numrepeats):
    regexexpr = "(?<!(CAG))(CAG){" + str(numrepeats) +"}(?!(CAG))"
    regs = re.finditer(regexexpr,seq)
    return [[i.start(0),i.end(0)] for i in regs]
    

# This function combines the above two functions were it finds the CAG repeats 
# for a fastafile for a chromosome
def getCAGRepeatsForChromsome(fastafile,chromosome,numrepeats):
    seq = extractChromosomeSeq(fastafile,chromosome)
    return extractIndexCAGrepeats(seq,numrepeats)