#####################################################################################
##
##  This script calculates linkage disequilibrium (LD) between sites in different alignments (nucleotide or amino acid).
##  It calculates both the Chi-Squared df and D' statistics of linkage disequilibrium.
##  It is meant to accompany the manuscript "Reassortment between influenza B lineages and the emergence of a co-adapted PB1-PB2-HA gene complex"
##  and has only been tested with alignments produced specifically for it.
##  You must have python, numpy and biopython installed on your machine to use this script.
##
##  You also require either protein or nucleotide alignments (in nexus format) with the same number of taxa that have the same names.
##  The sequences used for this project are listed in the acknowledgment table accompanying the manuscript:
##  https://github.com/evogytis/fluB/tree/master/acknowledgement%20tables
##
##  and can be downloaded from GISAID available at:
##  http://platform.gisaid.org
##  
##  open command line, go to folder where the script is located and type in:
##  "python LD_calculator.py" (ignore ")
##
##  To switch between nucleotide and amino acid modes comment out mode='aa' or mode='nt', respectively.
##
##  The following figures in the manuscript were made using this script:
##  Figures 10, S7-S8
##
##
##  Alternatively, instead of running the script the output is available here:
##  https://github.com/evogytis/fluB/tree/master/data
##  
##
##  Gytis Dudas
##  Institute of Evolutionary Biology
##  University of Edinburgh
##  EH9 3JT
##  Edinburgh, UK
##
#####################################################################################

import sys
import numpy as np
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import re
from datetime import datetime as dt
import time
import calendar
from Bio import SeqIO
from Bio import AlignIO
import itertools

def unique(o, idfun=repr):
    """
    Returns a list of unique values in a list.
    """
    seen = {}
    return [seen.setdefault(idfun(e),e) for e in o if idfun(e) not in seen]

def column(data,col):
    """
    Returns a list that is a column in a list of lists.
    """
    return [row [col] for row in data]

def index(data,item):
    """
    Returns a list of indices where an item can be found within a given list.
    """
    return [i for i,x in enumerate(data) if x == item]

def frequency(data):
    """
    Returns a list of counts of each value in a list.
    """
    uvals=unique(data)
    out = dict(zip(uvals, list(0 for u in range(len(uvals)))))
    for i in range(len(uvals)):
        out[uvals[i]]=data.count(uvals[i])
    return out

def decimalDate(sequenceName):
    """
    Returns a yyyy-mm-dd date format as a decimal date.
    """
 #   dateRegex='[A-Z]\/[A-Za-z\/\-\_\.0-9]+\/[0-9]+\_([0-9\-]+)'
    dateRegex='__([0-9\_]+)'
    C=re.search('%s'%(dateRegex),sequenceName)
    if C is not None:
        try:
            yrs = int(C.group(1).split('_')[0])
            mons = int(C.group(1).split('_')[1])
            dds = int(C.group(1).split('_')[2])
            return (float(dt(yrs,mons,dds).strftime("%j"))-1) / 366 + float(dt(yrs,mons,dds).strftime("%Y"))
        except IndexError:
            try:
                u = int(calendar.monthrange(int(C.group(1).split('_')[0]),int(C.group(1).split('_')[1]))[1])
                l = 1
                yrs = int(C.group(1).split('_')[0])
                mons = int(C.group(1).split('_')[1])
                dds=int((u+l)/float(2))
                return (float(dt(yrs,mons,dds).strftime("%j"))-1) / 366 + float(dt(yrs,mons,dds).strftime("%Y"))
            except IndexError:
                yrs = int(C.group(1))
                return float(yrs)+0.5
    else:
        print sequenceName,'no regex match'

def removeItem(someList,itemList):
    """
    Returns a list with an item excluded.
    """
    return [x for x in someList if x not in itemList]

def removeItem2(someList,itemList):
    """
    Returns a list if an item is present.
    """
    return [x for x in someList if x in itemList]

def flatten(someList):
    """
    Reduces a list of lists into a single list.
    """
    return [item for sublist in someList for item in sublist]

import collections
def overlap(a,b):
    """
    Finds the overlap between two lists.
    """
    a_multiset = collections.Counter(a)
    b_multiset = collections.Counter(b)

    overlap = list((a_multiset & b_multiset).elements())
    a_remainder = list((a_multiset - b_multiset).elements())
    b_remainder = list((b_multiset - a_multiset).elements())

    return overlap, a_remainder, b_remainder


## Makes a list of dates if partitioning of LD is desired
# Default value makes a single bin for all sequences isolated between 1900 and 2014
step=40
windowSize=0
timeline=np.arange(1980,2020-windowSize,step)

## Segments to analyze                                                                       

segments=['MERS-CoV.85.GARD.fasta']

## Taxa names to exclude from analysis
# LD analysis in manuscript did not exclude any taxa
#mixedLineage=['B/Alaska/03/1992_1992-06-13', 'B/Alaska/03/1992_1992-06-13', 'B/Alaska/03/1992_1992-06-13', 'B/Cordoba/2979/1991_1991-06-13', 'B/Cordoba/2979/1991_1991-06-13', 'B/Cordoba/2979/1991_1991-06-13', 'B/Bangkok/163/1990_1990', 'B/Bangkok/163/1990_1990', 'B/Bangkok/163/1990_1990', 'B/Norway/1/84_1984', 'B/Norway/1/84_1984', 'B/Norway/1/84_1984', 'B/Houston/1/91_1991', 'B/Houston/1/91_1991', 'B/Houston/1/91_1991', 'B/Memphis/5/93_1993', 'B/Memphis/5/93_1993', 'B/Memphis/5/93_1993', 'B/Bangkok/153/1990_1990', 'B/Bangkok/153/1990_1990', 'B/Bangkok/153/1990_1990', 'B/Johannesburg/06/1994_1994', 'B/Johannesburg/06/1994_1994', 'B/Johannesburg/06/1994_1994', 'B/Lisbon/02/1994_1994-10-26', 'B/Lisbon/02/1994_1994-10-26', 'B/Lisbon/02/1994_1994-10-26', 'B/Connecticut/02/1995_1995-01-05', 'B/Connecticut/02/1995_1995-01-05', 'B/Connecticut/02/1995_1995-01-05', 'B/Texas/14/1991_1991-01-11', 'B/Texas/14/1991_1991-01-11', 'B/Texas/14/1991_1991-01-11', 'B/Connecticut/07/1993_1993-01-26', 'B/Connecticut/07/1993_1993-01-26', 'B/Connecticut/07/1993_1993-01-26', 'B/Ann_Arbor/1994_1994', 'B/Ann_Arbor/1994_1994', 'B/Ann_Arbor/1994_1994', 'B/Oita/15/92_1992-10-16', 'B/Oita/15/92_1992-10-16', 'B/Oita/15/92_1992-10-16', 'B/New_York/24/1993_1993-12-08', 'B/New_York/24/1993_1993-12-08', 'B/New_York/24/1993_1993-12-08', 'B/New_York/39/1991_1991-03-20', 'B/New_York/39/1991_1991-03-20', 'B/New_York/39/1991_1991-03-20', 'B/Nanchang/6/96_1996', 'B/Nanchang/6/96_1996', 'B/Nanchang/6/96_1996', 'B/Nanchang/630/94_1994', 'B/Nanchang/630/94_1994', 'B/Nanchang/630/94_1994']
mixedLineage=[]

## Select mode - nt for nucleotide alignments, aa for amino acid alignments, used later to filter out invalid residues in the alignment
mode='nt'
#mode='aa'


if mode=='nt':
    print 'Running in nucleotide LD mode'
elif mode=='aa':
    print 'Running in amino acid LD mode'


## Iterate over all pairs of segments
for mem1 in range(0,len(segments)):
    print '\n\nFinding associations in %s'%(segments[mem1])

    dictA={}
    dictB={}

    alignmentA=[]

    strains=[]

    segment1=segments[mem1]
    
    ## input and output paths
    path='/Users/admin/Documents/MERS_rec/helper_scripts/'
    outpath='/Users/admin/Documents/MERS_rec/helper_scripts/'

    if mode=='aa':
        ## Define alignment name
        seg1='%s'%(segment1)
    elif mode=='nt':
        ## Define alignment name
        seg1='%s'%(segment1)

    f = open(outpath+'%s_associations.%s.csv'%(segment1,mode),'w')
    header='N associations,N strains,N associations usable,haplotypes,site,number of polymoprhisms at site,number of hosts,minor allele freq,rarest host freq,D,|D\'|,Chi,ChiSqDf'
    print>>f,header

    ## open both alignments
    handle1 = open(path+seg1, "rU")

    seqCheck=[]

    ## Iterate over all sequences, assign to time blocks, exclude some sequences if necessary
    for record in AlignIO.read(handle1, "fasta"):
        seqCheck.append(record)

        
        if record.id not in mixedLineage:
            dictA[record.id]=record.seq
            alignmentA.append(record.seq)
            if record.id not in mixedLineage and record.id not in strains:
                strains.append(record.id)
            if 'Camel' in record.id:
                dictB[record.id]='C'
            else:
                dictB[record.id]='H'
        else:
            pass


    print 'Sequence number check:',len(seqCheck)

    ## iterate through time blocks
    #for tp in range(len(timeline)):
    polymorphicLociA=[]
    polymorphicLociB=[]

    ## polymorphicPairs contains positions of polymorphic loci, reduces loop depth
    polymorphicPairs=[]

    ## assert that numbers of parsed things are the same
    assert len(alignmentA)==len(dictA)==len(strains)==len(dictB)

    ## assert that the time bin has sequences at all
    assert len(dictA)!=0,'Sequences were not captured in a time bin'
    
    ## used for debugging, determines which site LD is calculated from in both alignments, 0 by default
    startFrom=0

    ## only continue if time bin has sequences
    if len(dictA)!=0:
        ## Iterate over all pairs of sites in both alignments, identify polymorphic sites
        for i in range(startFrom,max([len(x) for x in dictA.values()])):
            #print timeline[tp],'site',i+1

            ## ignore invariant sites at locus 1
            if mode=='aa' and len(removeItem(unique(column(alignmentA,i)),['-','?']))==1:
                pass
            elif mode=='nt' and len(removeItem(unique(column(alignmentA,i)),['R','Y','S','W','K','M','B','D','H','H','V','N','-']))==1:
                pass
            ## ignore loci with at least one gap
            #elif '-' in unique(column(alignmentA,i)):
            #    pass
            else:
                ## log polymorphic loci in alignmentA
                polymorphicLociA.append(i)
                

        #print dictB.values()
        for j in range(startFrom,max([len(x) for x in dictB.values()])):
            polymorphicLociB.append(j)
                                                    
        for loc1 in polymorphicLociA:
            for loc2 in polymorphicLociB:
                polymorphicPairs.append((loc1,loc2))

        print '\nNumber of polymorphic pairs of loci on %s %s'%(segments[mem1],len(polymorphicPairs))

        ## seen is a dict that maps haplotype frequencies to their LD statistics, so calculations that were done earlier don't have to be repeated
        seen={}
                        
        ## iterate through potentially polymorphic site pairs
        for z in range(len(polymorphicPairs)):

            ## report progress
            sys.stdout.write('\r')
            sys.stdout.write("[%-100s] %d%%" % ('='*(100*z/len(polymorphicPairs)),100*z/len(polymorphicPairs)))
            sys.stdout.flush()

            ## unpack sites from site pairs
            site1=polymorphicPairs[z][0]
            site2=polymorphicPairs[z][1]
            
            ## reconstruct haplotypes
            checkHaplotypes={strain:(dictA[strain][site1],dictB[strain][site2]) for strain in strains}

            ## filter haplotypes, keep those that do not have invalid residues
            if mode=='aa':
                ## keep haplotypes that don't have gaps or ? as residues
                haplotypes={strain:checkHaplotypes[strain] for strain in checkHaplotypes.keys() if '-' not in checkHaplotypes[strain] and '?' not in checkHaplotypes[strain]}
            elif mode=='nt':
                ## haplotypes are 2 residues long, keep only those that are composed of valid nucleotides
                haplotypes={strain:checkHaplotypes[strain] for strain in checkHaplotypes.keys() if tuple(True for u in checkHaplotypes[strain] if u in ('A','C','T','G','H')).count(True)==2}

            ## some polymorphic loci still have invalid residues
            haps=unique(haplotypes.values())

            ## find alleles at locus 1 and locus 2
            poly1=[haps[q][0] for q in range(len(haps))]
            poly2=[haps[q][1] for q in range(len(haps))]
            poly1=unique(poly1)
            poly2=unique(poly2)

            ## count numbers of each allele
            poly1Count=tuple(column(haplotypes.values(),0).count(x) for x in poly1)
            poly2Count=tuple(column(haplotypes.values(),1).count(x) for x in poly2)

            ## if alleles at locus 1 are associated with invalid alleles at locus 2, things won't work
            if len(poly1)>1 and len(poly2)>1:

                ## make flat contingency table of observed haplotypes
                a=[poly1,poly2]
                allhaps=[x for x in itertools.product(*a)]
                hapCount = [haplotypes.values().count(x) for x in allhaps]
                
                ## assert that haplotype numbers are equal to counts of each allele
                assert sum(poly1Count)==sum(poly2Count)==sum(hapCount)
                total=sum(hapCount)

                ## calculate frequencies of alleles and haplotypes
                poly1Freq=[x/float(total) for x in poly1Count]
                poly2Freq=[x/float(total) for x in poly2Count]
                hapFreq=[x/float(total) for x in hapCount]

##                    # debugging extras
##                    print i,j
##                    print 'alleles at site 1:',poly1,'alleles at site 2:',poly2
##                    print 'all possible haplotypes:',allhaps
##                    print 'observed haplotypes:',hapCount,'observed alleles site 1:',poly1Count,'observed alleles site 2:',poly2Count
##                    print 'haplotype frequencies:',hapFreq,'allele frequencies at site 1:',poly1Freq,'allele frequencies at site 2:',poly2Freq
##                    print 'total number of haplotypes that passed filtering:',total

                ## only keep those haplotypes where all haplotypes have been observed
                #if sum([1 for x in hapCount if x!=0])==len(hapCount):
                if sum([1 for x in hapCount if x!=0])!=0:
                    #print i,j,hapCount

                    ## check whether these haplotype frequencies were seen before
                    if seen.has_key(str(hapCount)):
                        
                        ## give answer computed earlier, but replace site numbers
                        preComp=seen[str(hapCount)].split(',')
                        preComp[4]=str(site1+1)
                        #preComp[5]=str(site2+1)
                        
                        #print '\t'.join(preComp),'(PRECOMPUTED)'
                        print>>f,','.join(preComp)

                    ## calculate everything
                    else:
                        ## calculate D' as well if loci biallelic
                        D='NaN'
                        Dprime='NaN'
                        if len(poly1)==2 and len(poly2)==2:

                            ## calculate D
                            D=hapFreq[0]-(poly1Freq[0]*poly2Freq[0])

                            ## calculate Dmax
                            if D>=0:
                                Dmax=min([poly1Freq[0]*(1-poly2Freq[0]),(1-poly1Freq[0])*poly2Freq[0]])
                            else:
                                Dmax=min([(1-poly1Freq[0])*(1-poly2Freq[0]),poly1Freq[0]*poly2Freq[0]])

                            ## normalize D by Dmax
                            Dprime=np.absolute(D/float(Dmax))
                        
                        ## calculate ChiSq
                        ChiSq=0
                        for q in range(len(poly1)):
                            for w in range(len(poly2)):
                                observed=hapCount[q*len(poly2)+w]
                                expected=poly1Freq[q]*poly2Freq[w]*total
                                ChiSq+=((observed-expected)**2)/float(expected)

                        ## calculate ChiSqdf
                        ChiSqdf=ChiSq/float(total*(len(poly1)-1)*(len(poly2)-1))
                        seen[str(hapCount)]='%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s'%(sum([1 for x in hapCount if x!=0]),len(strains),total,'\t'.join(['%s:%s'%(''.join(allhaps[a]),hapCount[a]) for a in range(len(hapCount))]),site1+1,len(poly1),len(poly2),min(poly1Freq),min(poly2Freq),D,Dprime,ChiSq,ChiSqdf)
                        #print '\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(sum([1 for x in hapCount if x!=0]),len(strains),total,site1+1,len(poly1),len(poly2),min(poly1Freq),min(poly2Freq),D,Dprime,ChiSq,ChiSqdf)
                        print>>f,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s'%(sum([1 for x in hapCount if x!=0]),len(strains),total,'\t'.join(['%s:%s'%(''.join(allhaps[a]),hapCount[a]) for a in range(len(hapCount))]),site1+1,len(poly1),len(poly2),min(poly1Freq),min(poly2Freq),D,Dprime,ChiSq,ChiSqdf)

        ## close output file
        f.close()

        ## close alignment file
        handle1.close()

## done
print '\nDone!'
