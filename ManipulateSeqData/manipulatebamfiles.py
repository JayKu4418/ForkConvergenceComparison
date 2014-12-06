__author__ = 'jayashreekumar'

"""
This script contains functions that will convert BAM files to BED or BEDGRAPH or SGR files.
"""

################################################################################################
# IMPORT
import optparse,subprocess
from operator import itemgetter
yeast_chrm_convert = {'chrI':'1','chrII':'2','chrIII':'3','chrIV':'4','chrV':'5','chrVI':'6','chrVII':'7','chrVIII':'8','chrIX':'9','chrX':'10','chrXI':'11','chrXII':'12','chrXIII':'13','chrXIV':'14','chrXV':'15','chrXVI':'16'}
yeast_chrm_convert_v2 = {'chr1_ref_v2':'1','chr2_ref_v2':'2','chr3_ref_v2':'3','chr4_ref_v2':'4','chr5_ref_v2':'5','chr6_ref_v2':'6','chr7_ref_v2':'7','chr8_ref_v2':'8','chr9_ref_v2':'9','chr10_ref_v2':'10','chr11_ref_v2':'11','chr12_ref_v2':'12','chr13_ref_v2':'13','chr14_ref_v2':'14','chr15_ref_v2':'15','chr16_ref_v2':'16'}
import collections
################################################################################################
# FUNCTIONS

# This function converts a sam file to bam file
def convertsam2bamfile(samfile, bamfile):
    subprocess.call('samtools view -bS '+ samfile + ' > ' + bamfile,shell=True)

# This function extracts only properly paired reads from a bam file
def convert2onlybamfilepairedend(bamfile, pairedendsbamfile):
    subprocess.call('samtools view -bf 0x2 '+ bamfile + ' > ' + pairedendsbamfile,shell=True)

# This function sorts a bamfile into sortedbamfile
def sortbamfile(bamfile,sortedbamfile):
    subprocess.call('samtools sort -n '+ bamfile + ' ' + sortedbamfile,shell=True)
    
# This function creates an index for sortedbamfile
def createbamindexfile(bamfile, indexedbamfile):
    subprocess.call('samtools index '+ bamfile + ' ' + indexedbamfile,shell=True)

# This function converts a bam file to bam file with uniquely mapped reads
def convert2umrbamfile(bamfile, umrbamfile,mapq):
    subprocess.call('samtools view -bq ' + str(mapq) + ' ' + bamfile + ' > ' + umrbamfile,shell=True)

# Unpaired reads This function converts a bam file to a bed file
def convert2bedfile(bamfile,bedfile):
    subprocess.call('bedtools bamtobed -i ' + bamfile + ' > ' + bedfile,shell=True)

# This function removes reads for chromosomes 17 and 18 from bedfile or bedgraph file or sgr file for s.cerevisiae. 17 and 18
# mitochondrial chromosomes
def rmLines(rfile,wfile,option):

    with open(rfile) as f:
        lines = f.readlines()
        linesplit = [i.strip().split('\t') for i in lines]

    with open(wfile,'w') as fw:

        for i in linesplit:
            if (i[0]!='18' and i[0]!='17'):
                fw.write(i[0]+'\t')
                fw.write(i[1]+'\t')
                if option == 'bedgraph':
                    fw.write(i[2]+'\t')
                    fw.write(i[3]+'\n')
                elif option == 'bed':
                    fw.write(i[2]+'\t')
                    fw.write(i[3]+'\t')
                    fw.write(i[4]+'\t')
                    fw.write(i[5]+'\n')
                else:
                    fw.write(i[2]+'\n')

# This function converts a bam file to a bedgraph file for a particular strand. You must specify either the Watson (+) or
# the Crick (-) strand
def convertfrombam2bedgraphfile(bamfile,bedgfile,strand):
    subprocess.call('bedtools genomecov -ibam ' + bamfile + ' -bga -strand ' + strand + ' > ' + bedgfile,shell=True)

# This function converts a bed file to a bedgraph file for a particular strand. You must specify either the Watson (+) or
# the Crick (-) strand. You also need to specify a file which contains the sizes of each chromosome of S. cerevisiae
def convertfrombed2bedgraphfile(bedfile,bedgfile,strand,sizefile):
    subprocess.call('bedtools genomecov -i ' + bedfile + ' -g ' + sizefile + ' -bga -strand ' + strand + ' > ' + bedgfile,shell=True)

# This function converts a bam file to an sgr file which gives depth at each base position. You must specify either the
# Watson (+) strand or the Crick (-) strand
def convertfrombam2sgrfile(bamfile,sgrfile,strand):
    subprocess.call('bedtools genomecov -ibam ' + bamfile + ' -d -strand ' + strand + ' > ' + sgrfile,shell=True)

# This function converts a bed file to an sgr file which gives depth at each base position.You must specify either the
# Watson (+) strand or the Crick (-) strand. You also need to specify a file which contains the sizes of each chromosome of S. cerevisiae
def convertfrombed2sgrfile(bedfile,sgrfile,strand,sizefile):
    subprocess.call('bedtools genomecov -i ' + bedfile + ' -g ' + sizefile + ' -d -strand ' + strand + ' > ' + sgrfile,shell=True)

# This function removes mitochondrial chromosome data from bed files and replaces 'chrI' with '1' and replaces the id with a -
def removeMitochondrialChromosomeAndChrFromBedPeFiles(readbedpefile,writebedpefile,v2option):    
    f = open(readbedpefile)
    fw = open(writebedpefile,'w')
    
    for line in f:
        sl = line.strip().split('\t')
        if sl[0] != 'chrM':
            if v2option:
                fw.write(yeast_chrm_convert_v2[sl[0]] + '\t' + sl[1] + '\t' + sl[2] + '\t' + yeast_chrm_convert_v2[sl[3]] + '\t' +sl[4] + '\t' + sl[5] + '\t' + '-' + '\t' + sl[7] + '\t' + sl[8] + '\t' + sl[9] + '\n')
            else:        
                fw.write(yeast_chrm_convert[sl[0]] + '\t' + sl[1] + '\t' + sl[2] + '\t' + yeast_chrm_convert[sl[3]] + '\t' +sl[4] + '\t' + sl[5] + '\t' + '-' + '\t' + sl[7] + '\t' + sl[8] + '\t' + sl[9] + '\n')
    f.close()
    fw.close()
    
# This function removes mitochondrial chromosome data from sgr files and replaces 'chrI' with '1'
def removeMitochondrialChromosomeAndChrFromSgrFiles(readsgrfile,writesgrfile):
    with open(readsgrfile,'r') as f:
        with open(writesgrfile,'w') as fw:
            for line in f:
                sl = line.strip().split('\t')
                if sl[0] != 'chrM':
                    fw.write(yeast_chrm_convert[sl[0]] + '\t' + sl[1] + '\t' + sl[2] + '\n')
                    
def convertbam2bedpefile(bamfile,bedpefile):
    subprocess.call('bedtools bamtobed -i ' + bamfile + ' -bedpe -mate1 > ' + bedpefile,shell=True)


def convertbedpe2bedfile(bedpefile,bedfile):
    with open(bedpefile) as f:
        with open(bedfile,'w') as fw:
            for line in f:
                sl = line.strip().split('\t')
                if sl[8] == '+' and sl[9] == '-':
                    if int(sl[2]) > int(sl[1]) and int(sl[5]) > int(sl[4]) and int(sl[5]) > int(sl[1]):
                        fw.write(sl[0]+'\t'+ sl[1] + '\t' + sl[5] + '\t' + sl[6] + '\t' + sl[7] + '\t' + sl[8] + '\n')
                if sl[8] == '-' and sl[9] == '+':
                    if int(sl[2]) > int(sl[1]) and int(sl[5]) > int(sl[4]) and int(sl[2]) > int(sl[4]):
                        fw.write(sl[0]+'\t'+ sl[4] + '\t' + sl[2] + '\t' + sl[6] + '\t' + sl[7] + '\t' + sl[8] + '\n')

def convertZeroToOneBedPeFile(bedpereadfile,bedpewritefile,isrDNA):
    with open(bedpereadfile) as f:
        with open(bedpewritefile,'w') as fw:
            for line in f:
                sl = line.strip().split('\t')
                if isrDNA=='True':
                    fw.write(sl[0]+'\t'+ str(int(sl[1])+458797) + '\t' + str(int(sl[2])+458796) + '\t' + sl[3] + '\t' + str(int(sl[4])+458797) + '\t' +  str(int(sl[5])+458796) + '\t' + sl[6] + '\t' + sl[7] + '\t' + sl[8] + '\t' + sl[9] + '\n')
                else:                
                    fw.write(sl[0]+'\t'+ str(int(sl[1])+1) + '\t' + sl[2] + '\t' + sl[3] + '\t' + str(int(sl[4])+1) + '\t' +  sl[5] + '\t' + sl[6] + '\t' + sl[7] + '\t' + sl[8] + '\t' + sl[9] + '\n')

def sortBedFileByChrmAndBase(bedfile,bedfilesorted):
    subprocess.call('sort -k1n,1 -k2n,2 -t $\'\t\' ' + bedfile + ' > ' + bedfilesorted,shell=True)


def getFragmentLength(bedfile,writecountfile):
    with open(writecountfile,'w') as fw:    
        for c in range(1,17):
            chromosome= str(c)
            print chromosome
            with open(bedfile) as f:
                bedlines = [[int(line.strip().split('\t')[1]),int(line.strip().split('\t')[2]),line.strip().split('\t')[5]] for line in f if line.strip().split('\t')[0]==chromosome]
            uniquefrags = [list(i) for i in set(tuple(i) for i in bedlines)]
            uniquefragsSorted = sorted(uniquefrags,key=itemgetter(0,1,2))
            bedlinesCount = collections.Counter([tuple(i) for i in bedlines])
            for i in uniquefragsSorted:
                fw.write(chromosome + '\t' + str(i[0]) + '\t' + str(i[1]) + '\t' + i[2] + '\t' + str(bedlinesCount[tuple(i)]) + '\n')
def limitFragCount(fragcountfile,writefile,lengthlim):
    with open(writefile,'w') as fw:    
        for c in range(1,17):
            chromosome=str(c)
            print(chromosome)
            with open(fragcountfile) as f:
                lines = [line.strip().split('\t') for line in f if line.strip().split('\t')[0]==chromosome]
            for i in lines:
                if int(i[2])-int(i[1])+1 <=lengthlim:
                    fw.write(i[0]+'\t'+i[1]+'\t'+i[2]+'\t'+i[3]+'\t'+i[4]+'\n')
def getFragmentCountByStrand(countfile,writecountfile,strand):
    with open(countfile) as f:
        with open(writecountfile,'w') as fw:
            for line in f:
                sl = line.strip().split('\t')
                if sl[3] == strand:
                    fw.write(sl[0]+'\t'+ sl[1] + '\t' + sl[2] + '\t' + sl[3] + '\t' + sl[4] + '\n')
    
def convertribonucleotidebedpe2bedfile(bedpefile,bedfile):
    with open(bedpefile) as f:
        with open(bedfile,'w') as fw:
            for line in f:
                sl = line.strip().split('\t')
                fw.write(sl[0]+'\t'+ sl[1] + '\t' + sl[2] + '\t' + sl[6] + '\t' + sl[7] + '\t' + sl[8] + '\n')
                
def extract5pendsBedFile(sortedbedfile,strand,end,writefile):
   
    with open(writefile,'w') as fw:
        for c in range(1,17):
            chromosome = str(c)
            print(chromosome)
            with open(sortedbedfile) as f:    
                if (strand=='+' and end=='5') or (strand=='-' and end=='3') :
                    endforchrom = [int(line.strip().split('\t')[1]) for line in f if line.strip().split('\t')[0]==chromosome and line.strip().split('\t')[5]==strand]
                else:
                    endforchrom = [int(line.strip().split('\t')[2]) for line in f if line.strip().split('\t')[0]==chromosome and line.strip().split('\t')[5]==strand]
            if len(endforchrom)!=0:
                uniqueends = sorted(list(set(endforchrom)))
                countForends = collections.Counter()
                for i in endforchrom:
                    countForends[i] += 1
                for i in uniqueends:
                    fw.write(chromosome + '\t' + str(i) + '\t' + str(countForends[i]) + '\n')
                    
def cleanBedPeFile(bedpefile,cleanbedpefile,chromosome,mapqlim):
    with open(bedpefile) as f:
        lines = [line.strip().split('\t') for line in f if line.strip().split('\t')[0]==chromosome and line.strip().split('\t')[0]==line.strip().split('\t')[3] and line.strip().split('\t')[0]!='.' and line.strip().split('\t')[3]!='.']
    
    with open(cleanbedpefile,'w') as fw:
        for i in lines:
            if i[8]!=i[9] and int(i[7])>=mapqlim:
                fw.write(i[0]+'\t'+i[1]+'\t'+i[2]+'\t'+i[3]+'\t'+i[4]+'\t'+i[5]+'\t'+i[6]+'\t'+i[7]+'\t'+i[8]+'\t'+i[9]+'\n')
################################################################################################
# MAIN

def main(indexname,fastqfile1,fastqfile2,samfile,strain,isrDNA):
    
    
    #print fastqfile1
    #print fastqfile2
    #print samfile
    print("Paired end alignment being done to reference genome.......")
    
    subprocess.call('bowtie2 -x ' +  indexname + ' -1 ' + fastqfile1 + ' -2 ' +  fastqfile2 + ' -S ' + samfile,shell=True)   
     
    bamfile = samfile[0:(len(samfile)-4)] + '.bam'
    
    print("Converting SAM file to BAM file....")
    
    convertsam2bamfile(samfile,bamfile)
    
    sortedbamfile = bamfile[0:(len(bamfile)-4)] + 'sorted'    
        
    print("Sorting BAM file.......")
    
    sortbamfile(bamfile,sortedbamfile)
    
    #Create name of umr bam file
    umrbamfile = bamfile[0:(len(bamfile)-4)] + 'sortedumr.bam'

    # Convert bam file to a bam file with only UMRs
    print("Converting BAM file to one with only UMRs.........")
    convert2umrbamfile(sortedbamfile+'.bam',umrbamfile,30)
    
    # Always extract paired reads at the end after all processing has been done and before 
    # you want to convert into a bedpe file!!!!    
    pairedendbamfile = samfile[0:(len(samfile)-4)] + 'pairedreads.bam'
    
    print("Extracting only paired reads from sorted bam file with only UMRs")

    convert2onlybamfilepairedend(umrbamfile,pairedendbamfile)

    #indexbamfile= sortedbamfile + '.bai'
    
    #print ("Creating index for sorted BAM file .....")
    
    #createbamindexfile(sortedbamfile+'.bam',indexbamfile)

    # Create name of umr bedpe file
    umrbedpefile = bamfile[0:(len(bamfile)-4)] + '.bedpe'
    
    print("Converting paired reds BAM file to BEDPE file.........")       
    
    # Convert sorted umr file to bedpe file
    convertbam2bedpefile(pairedendbamfile,umrbedpefile)    
    
    # Create name of umr bed file
    #umrbedfile = bamfile[0:(len(bamfile)-4)] + '.bed'
    
    # Convert bam file with only UMRs to a bed file
    #print("Converting BAM file to BED file.........")
    #convert2bedfile(umrbamfile,umrbedfile)

    # Create name of umr sgr file for forward and reverse strands
    #umrbam2sgrforward = bamfile[0:(len(bamfile)-4)] + 'W.sgr'
    #umrbam2sgrreverse = bamfile[0:(len(bamfile)-4)] + 'C.sgr'

    # Convert bam file with only UMRs to sgr file for forward strand
    #print("Converting BAM file to SGR file for + strand..........")
    #convertfrombam2sgrfile(umrbamfile,umrbam2sgrforward,'+')

    # Convert bam file with only UMRs to sgr file for reverse strand
    #print("Converting BAM file to SGR file for - strand..........")
    #convertfrombam2sgrfile(umrbamfile,umrbam2sgrreverse,'-')

    #print("Removing mitochondiral chromosomes from Sgr files and roman characters")
    #removeMitochondrialChromosomeAndChrFromSgrFiles(umrbam2sgrforward,'new'+umrbam2sgrforward)
    #removeMitochondrialChromosomeAndChrFromSgrFiles(umrbam2sgrreverse,'new'+umrbam2sgrreverse)
    
    
    bedpereadfile = strain + 'CleanZeroBased.bedpe'
    bedpewritefile = strain + 'CleanOneBased.bedpe'
    bedwritefile = strain + '.bed'
    sortedbedfile= strain + 'Sorted.bed'
    sgrWfile = strain + 'W.sgr'
    sgrCfile = strain + 'C.sgr'
    strandcountfile = strain + 'FragCount.txt'
    #watsonstrandCountfile =strain + 'WFragCounttxt'
    #crickstrandCountfile =strain + 'CFragCounttxt'
    
    print(bedpereadfile)
    print(bedwritefile)
    
    print("Removing mitochrondrial chromosomes from Bed file and roman characters")    
    
    removeMitochondrialChromosomeAndChrFromBedPeFiles(umrbedpefile,bedpereadfile)
    
    print('Converting to One-Based BedPe File')
    
    convertZeroToOneBedPeFile(bedpereadfile,bedpewritefile,isrDNA)
    
    print("Converting to Bed file")
    
    convertbedpe2bedfile(bedpewritefile,bedwritefile)
    
    print("Sorting Bed file")
    
    sortBedFileByChrmAndBase(bedwritefile,sortedbedfile)
    
    print("Creating Sgr Watson file")
    
    convertfrombed2sgrfile(sortedbedfile,sgrWfile,'+','yeast.genome')
    
    print("Creating Sgr Crick file")
    
    convertfrombed2sgrfile(sortedbedfile,sgrCfile,'-','yeast.genome')
    
    print('Extracting Fragment Count')
    
    getFragmentLength(sortedbedfile,strandcountfile)
    
    #print("Extracting Fragment Count by Watson strand")
    
    #getFragmentCountByStrand(strandcountfile,watsonstrandCountfile,'+')
    
    #print("Extracting Fragment Count by Crick strand")
    
    #getFragmentCountByStrand(strandcountfile,crickstrandCountfile,'-')
    
    

"""
def main(fastqfile1,fastqfile2,samfile):
    
    print fastqfile1
    print fastqfile2
    print samfile
    #print("Paired end alignment being done to reference genome.......")
    
    #subprocess.call('bowtie2 -x sac-cerevisiae -1 ' + fastqfile1 + ' -2 ' +  fastqfile2 + ' -S ' + samfile,shell=True)   
     
    bamfile = samfile[0:(len(samfile)-4)] + '.bam'
    
    print("Converting SAM file to BAM file....")
    
    convertsam2bamfile(samfile,bamfile)
    
    sortedbamfile = bamfile[0:(len(bamfile)-4)] + 'sorted'    
        
    print("Sorting BAM file.......")
    
    sortbamfile(bamfile,sortedbamfile)
    
    #Create name of umr bam file
    umrbamfile = bamfile[0:(len(bamfile)-4)] + 'sortedumr.bam'

    # Convert bam file to a bam file with only UMRs
    print("Converting BAM file to one with only UMRs.........")
    convert2umrbamfile(sortedbamfile+'.bam',umrbamfile,30)
    
    # Always extract paired reads at the end after all processing has been done and before 
    # you want to convert into a bedpe file!!!!    
    pairedendbamfile = samfile[0:(len(samfile)-4)] + 'pairedreads.bam'
    
    print("Extracting only paired reads from sorted bam file with only UMRs")

    convert2onlybamfilepairedend(umrbamfile,pairedendbamfile)

    #indexbamfile= sortedbamfile + '.bai'
    
    #print ("Creating index for sorted BAM file .....")
    
    #createbamindexfile(sortedbamfile+'.bam',indexbamfile)

    # Create name of umr bedpe file
    umrbedpefile = bamfile[0:(len(bamfile)-4)] + '.bedpe'
    
    print("Converting paired reds BAM file to BEDPE file.........")       
    
    # Convert sorted umr file to bedpe file
    convertbam2bedpefile(pairedendbamfile,umrbedpefile)    
    
    # Create name of umr bed file
    #umrbedfile = bamfile[0:(len(bamfile)-4)] + '.bed'
    
    # Convert bam file with only UMRs to a bed file
    #print("Converting BAM file to BED file.........")
    #convert2bedfile(umrbamfile,umrbedfile)

    # Create name of umr sgr file for forward and reverse strands
    #umrbam2sgrforward = bamfile[0:(len(bamfile)-4)] + 'W.sgr'
    #umrbam2sgrreverse = bamfile[0:(len(bamfile)-4)] + 'C.sgr'

    # Convert bam file with only UMRs to sgr file for forward strand
    #print("Converting BAM file to SGR file for + strand..........")
    #convertfrombam2sgrfile(umrbamfile,umrbam2sgrforward,'+')

    # Convert bam file with only UMRs to sgr file for reverse strand
    #print("Converting BAM file to SGR file for - strand..........")
    #convertfrombam2sgrfile(umrbamfile,umrbam2sgrreverse,'-')

    #print("Removing mitochrondrial chromosomes from Bed file and roman characters")    
    #removeMitochondrialChromosomeAndChrFromBedFiles(umrbedfile,'new'+umrbedfile)

    #print("Removing mitochondiral chromosomes from Sgr files and roman characters")
    #removeMitochondrialChromosomeAndChrFromSgrFiles(umrbam2sgrforward,'new'+umrbam2sgrforward)
    #removeMitochondrialChromosomeAndChrFromSgrFiles(umrbam2sgrreverse,'new'+umrbam2sgrreverse)
"""

 
if __name__ == '__main__':

    parser = optparse.OptionParser()

    parser.add_option('-s','--samfile',dest='samfile',help='provide a sam file with extension .sam')
    
    parser.add_option('-u','--fastqfile1',dest='fastqfile1',help='provide a fastq file with extension .fastq')
    
    parser.add_option('-d','--fastqfile2',dest='fastqfile2',help='provide a fastq file with extension .fastq')
    
    parser.add_option('-e','--strain',dest='strain',help='provide a strain name')

    parser.add_option('-r','--isrdna',dest='isrdna',help='is this rDNA or not')

    parser.add_option('-i','--indexname',dest='indexname',help='provide indexname')    
    
    # load the data
    (options,args) = parser.parse_args()

    # process the inputs
    samfile = options.samfile
    fastqfile1 = options.fastqfile1
    fastqfile2 = options.fastqfile2
    strain = options.strain
    isrdna = options.isrdna
    indexname = options.indexname
    
    main(indexname,fastqfile1,fastqfile2,samfile,strain,isrdna)

"""
if __name__ == '__main__':
    parser = optparse.OptionParser()

    parser.add_option('-s','--strain',dest='strain',help='provide a strain name')
    # load the data
    (options,args) = parser.parse_args()

    # process the inputs
    strain = options.strain
    
    main(strain)
"""