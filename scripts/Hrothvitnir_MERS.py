#####################################################################################
##
##  This script analyzes posterior sets of trees produced by BEAST and writes out mutations inferred via robust counting.
##  It is meant to accompany the manuscript "MERS-CoV recombination: implications about the reservoir and potential for adaptation"
##  and has only been tested with trees produced specifically for it.
##  You must have python and numpy installed on your machine to use this script.
##
##  Usage:
##  open command line, go to folder where the script is located and type in:
##  "python Hrothvitnir.py -i <<path to .trees.txt file>> 1><<path to output file>>" (ignore <<, >> and ")
##
##  Gytis Dudas
##  Institute of Evolutionary Biology
##  University of Edinburgh
##  EH9 3JT
##  Edinburgh, UK
##
#####################################################################################

import re
from datetime import datetime as dt
import time
import calendar
import numpy as np
import random
import glob
import os
import sys

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

def makeSliceList(slices,start,end,bins):
    """
    Return a list of dates spaced out regularly within a year.
    """
    for i in np.arange(start,end,1/float(bins)):
        slices.append([i])

def mostRecent(some_list):
    """
    Identify the most recent tip in the tree and its date.
    """
    full=[]
    month=[]
    year=[]
    ## Summarize tip names in regular expressions
    ## Enclose the bit of the sequence name containing the date in brackets
    ## Date should be in format yyyy-mm-dd
##    dateRegex='[A-Z]\/[A-Za-z\/\-\_\.0-9]+\/[0-9]+\_([0-9\-]+)'
    dateRegex='[A-Za-z0-9\_\|\-]+\|([0-9\-]+)$'
    for i in some_list:
        #sys.stderr.write('\n%s'%i)
        C=re.search('%s'%(dateRegex),i[1])
        if C is not None:
            try:
                yrs = int(C.group(1).split('-')[0])
                mons = int(C.group(1).split('-')[1])
                dds = int(C.group(1).split('-')[2])
                full.append([C.group(),(float(dt(yrs,mons,dds).strftime("%j"))-1) / 366 + float(dt(yrs,mons,dds).strftime("%Y"))])
            except IndexError:
                try:
                    u = int(calendar.monthrange(int(C.group(1).split('-')[0]),int(C.group(1).split('-')[1]))[1])
                    l = 1
                    yrs = int(C.group(1).split('-')[0])
                    mons = int(C.group(1).split('-')[1])
                    dds=int((u+l)/float(2))
                    month.append(C.group())
                except IndexError:
                    yrs = int(C.group(1))
                    mons=6
                    dds=15
                    year.append(C.group())
    try:
        assert (len(full)+len(month)+len(year)==len(tips)),'Warning, some dates were not captured by regex!\nReview the \"mostRecent\" function'
        assert (len(full)>1),'No tip has a complete date'
        #sys.stderr.write('\n%s'%full)
        for i in object_list:
            if isinstance(i,leaf)==True:
                #sys.stderr.write('\n%s'%i.name)
                if i.name in column(full,0):
                    #sys.stderr.write('\n%s'%full)
                    #sys.stdout('%s\t%s'%(i.height,int(i.numName)))
                    return i.height,int(i.numName),max(column(full,1))
        
    except ValueError:
        pass
                

def make_tree(data):
    """
    Scans tree string looking for BEAST formatting.
    NOTE - only works with trees drawn from the posterior distribution of trees.
    """
    global i
    global trait_list
    global trait_names
    global partition_name
    stored_i=None
    while i < len(data):
 #       sys.stderr.write('\n%s\t%s'%(i,data[i]))
        assert (stored_i != i),'\nTree string unparseable\nStopped at %s (index %s)\nstring region looks like this: %s'%(data[i],i,data[i-20:i+20])
        stored_i=i
        cerberus=re.match('\[&([A-Za-z0-9\s\.\-\_\,\=\"\']+)\](:|;)',data[i:i+400])
        if cerberus is not None:
            traitComment=cerberus.group(1).split(',')
            for j in traitComment:
                traitName=j.split('=')[0]
                traitValue=j.split('=')[1].strip('"')
                ll.cur_node.traits[traitName]=traitValue
                if traitName not in column(trait_list,0):
                    trait_list.append([traitName,[traitValue]])
                if traitValue not in trait_list[index(column(trait_list,0),traitName)[0]][1]:
                    trait_list[index(column(trait_list,0),traitName)[0]][1].append(traitValue)
                    trait_list[index(column(trait_list,0),traitName)[0]][1]=sorted(trait_list[index(column(trait_list,0),traitName)[0]][1])
            i+=len(cerberus.group())-1

        if data[i] == '(':
            ll.add_node(i)
            i+=1
            cerberus=re.match('([0-9]+)',data[i:i+10])
            if cerberus is not None:
 #               sys.stderr.write('\ntip: %s'%(cerberus.group(1)))
                ll.add_leaf(i,cerberus.group(1))
                i+=len(cerberus.group(1))
            elif data[i]=='(':
                make_tree(data)
                
        if data[i] == ',':
            i+=1
            ll.move_up()
            ll.cur_node.secondChild=True
            cerberus=re.match('([0-9]+)',data[i:i+10])
            if cerberus is not None:
                ll.add_leaf(i,cerberus.group(1))
                i+=len(cerberus.group(1))

        if data[i] == ':':
            i+=1
            #print data[i:i+500]
            comment=re.match('\[&([A-Za-z0-9\s\_\-\.\,\=\{\},]+)\]([0-9\.\-E]+)',data[i:])
            if comment is not None:
                ll.cur_node.branch=float(comment.group(2))
                ll.cur_node.height=float(comment.group(2))
                com_list=[]
                com_list=comment.group(1).split(',')
                for j in com_list:
                    if 'r' in j.split('=')[0] and j.split('=')[0] not in ll.cur_node.traits.keys():
                        if partition_name=='':
                            partition_name=com_list[0].split('=')[0]
                        ll.cur_node.rate=float(com_list[0].split('=')[1])
                    elif 'r' in j.split('=')[0] and j.split('=')[0] in ll.cur_node.traits.keys():
                        ll.cur_node.traitRates[j.split('=')[0]]=float(j.split('=')[1])

                SrobustCounts=re.search('S=\{([0-9\.\{\},ACTG]+)\},N',comment.group(1))
                if SrobustCounts is not None:
                    for rcChunks in SrobustCounts.group(1).split('},{'):
                        codon=re.search('([0-9]+),([0-9\.\-E]+),([ACTG][ACTG][ACTG]),([ACTG][ACTG][ACTG])',rcChunks)
                        if codon is not None:
                            ll.cur_node.S.append((int(codon.group(1)),float(codon.group(2)),codon.group(3),codon.group(4)))




                NrobustCounts=re.search('N=\{([0-9\.\{\},ACTG]+)\}',comment.group(1))
                if NrobustCounts is not None:
                    for rcChunks in NrobustCounts.group(1).split('},{'):
                        codon=re.search('([0-9]+),([0-9\.E\-]+),([ACTG][ACTG][ACTG]),([ACTG][ACTG][ACTG])',rcChunks)
                        if codon is not None:
                            ll.cur_node.N.append((int(codon.group(1)),float(codon.group(2)),codon.group(3),codon.group(4)))



                i+=len(comment.group())

        if data[i] == ')':
            ll.move_up()
            i+=1

        if data[i] == ';':
            assert (len(tips)*2-1 == len(object_list)),'\nTree string has been parsed incorrectly:\nexpected number of objects in tree %s\nobjects found in the tree string: %s'%(len(tips)*2-1,len(object_list))
            if len(trait_names)==0:
                trait_names=ll.cur_node.traits.keys()
            break

################################
## end of string parsing block
################################

## Node, leaf and tree objects used to represent the tree structure.
## Nodes connect to other nodes/leaves (tips) and always have two children

class node:
    def __init__(self):
        self.rate=0.0
        self.branch=None
        self.height=None
        self.absoluteTime=None
        self.parent=None
        self.leftChild=None
        self.rightChild=None
        self.secondChild=False
        self.traits={}
        self.traitRates={}
        self.index=None
        self.childHeight=None
        self.numChildren=0
        self.N=[]
        self.S=[]
        self.buS=None
        self.buN=None
        ## contains references to all tips of this node
        self.leaves=[]

class leaf:
    def __init__(self):
        self.name=None
        self.numName=None
        self.index=None
        self.branch=0
        self.absoluteTime=None
        self.rate=0.0
        self.height=0
        self.parent=None
        self.traits={}
        self.traitRates={}
        self.N=[]
        self.S=[]
        self.buS=None
        self.buN=None


class tree:
    def __init__(self):
        global object_list
        self.cur_node=node()
        self.cur_node.branch=0
        self.cur_node.height=0
        self.root=self.cur_node
        self.cur_node.index='Root'
        self.cur_node.highestTip=0

    def move_up(self):
        """
        Go to current object's parent.
        """
        node=self.cur_node
        self.cur_node=node.parent

    def add_node(self,i):
        """
        Attach a node to current object as a first or second child.
        i refers to the index of the bracket defining the node in the tree string and is a unique identifier.
        """
        if self.cur_node.parent==self.root:
            self.cur_node.branch=0
            self.cur_node.height=0
        new_node=node()
        new_node.index=i
        object_list.append(new_node)
        if self.cur_node.secondChild==False:
            new_node.parent=self.cur_node
            self.cur_node.leftChild=new_node
            self.cur_node=new_node
        else:
            new_node.parent=self.cur_node
            self.cur_node.rightChild=new_node
            self.cur_node=new_node

    def add_leaf(self,i,name):
        """
        Attach a leaf to current object as a first or second child.
        i refers to the index of the tip number in the tree string and is a unique identifier.
        """
        new_leaf=leaf()
        new_leaf.index=i
        object_list.append(new_leaf)
        if self.cur_node.secondChild==False:
            new_leaf.name=name
            new_leaf.numName=name
            new_leaf.parent=self.cur_node
            self.cur_node.leftChild=new_leaf
            self.cur_node=new_leaf
        else:
            new_leaf.name=name
            new_leaf.numName=name
            new_leaf.parent=self.cur_node
            self.cur_node.rightChild=new_leaf
            self.cur_node=new_leaf

    def setAbsoluteTime(self,height,date):
        """
        Sets the absolute time of each object in the tree given the height of the most recent tip and its date.
        """
        for i in object_list:
            i.absoluteTime=date-height+i.height

    def renameTips(self):
        """
        Changes the name attribute of leaf objects to what the sequence name is.
        """
        global tips
        for i in object_list:
            if isinstance(i,leaf)==True:
                i.name=tips[index(column(tips,0),int(i.name))[0]][1]

    def TMRCA(self,objectA,objectB):
        """
        Identifies the most recent common ancestor of two objects.
        Accepts either leaf or node objects.
        """
        cur_nodeA=objectA
        cur_nodeB=objectB
        dictA={}
        dictB={}

        while cur_nodeA.parent!=None:
            cur_nodeA=cur_nodeA.parent
            dictA[cur_nodeA.index]=cur_nodeA

        while cur_nodeB.parent!=None:
            cur_nodeB=cur_nodeB.parent
            dictB[cur_nodeB.index]=cur_nodeB
            if dictA.has_key(cur_nodeB.index):
                key=cur_nodeB.index
                break

        return dictA[key]

    def timeToTMRCA(self,tip1,tip2):
        """
        Finds the total time taken to go from one tip to another.
        Takes two leaf objects.
        """
        cur_nodeA=tip1
        cur_nodeB=tip2
        dictA={}
        dictB={}

        while cur_nodeA.parent!=None:
            cur_nodeA=cur_nodeA.parent
            dictA[cur_nodeA.index]=cur_nodeA

        while cur_nodeB.parent!=None:
            cur_nodeB=cur_nodeB.parent
            dictB[cur_nodeB.index]=cur_nodeB
            if dictA.has_key(cur_nodeB.index):
                key=cur_nodeB.index
                break

        intersect=dictA[key]

        return -2*float(intersect.absoluteTime)+float(tip1.absoluteTime)+float(tip2.absoluteTime)

    def traverse_tree(self):
        """
        Pre-order tree traversal.
        Required to set the height of each object in the tree, to inform each node of its descendent tips.
        """
        cur_node=None
        cur_node=self.root.leftChild

        numChildren=0
        seen=[]
        highestTip=0
        highestTip=cur_node.branch
        height=0
        height=cur_node.branch
        root=False

        while root==False:
            if isinstance(cur_node,node):
                while cur_node.leftChild.index in seen and cur_node.rightChild.index in seen:
                    if cur_node.childHeight <= highestTip:
                        cur_node.childHeight = highestTip
                    elif cur_node.childHeight > highestTip:
                        highestTip = cur_node.childHeight

                    if cur_node.parent.index=='Root':
                        root=True
                        break
                    else:
                        cur_node.parent.numChildren+=cur_node.numChildren
                        cur_node.parent.leaves+=cur_node.leaves
                        cur_node.parent.leaves=map(str,sorted(map(int,cur_node.parent.leaves)))
                        cur_node.height=height
                        height-=cur_node.branch
                        cur_node=cur_node.parent

                    if cur_node.index=='Root':
                        root=True
                        break

                if cur_node.leftChild.index in seen or cur_node.index in seen:
                    height+=cur_node.rightChild.branch
                    cur_node.childHeight=highestTip
                    highestTip=0
                    cur_node=cur_node.rightChild
                    cur_node.height=height
                else:
                    height+=cur_node.leftChild.branch
                    seen.append(cur_node.index)
                    cur_node=cur_node.leftChild
                    cur_node.height=height

            if isinstance(cur_node,leaf):
                cur_node.parent.numChildren+=1
                cur_node.parent.leaves.append(cur_node.name)
                cur_node.parent.leaves=map(str,sorted(map(int,cur_node.parent.leaves)))
                seen.append(cur_node.index)
                cur_node.height=height
                highestTip=height
                height-=cur_node.branch
                cur_node=cur_node.parent

        self.highestTip=highestTip


    def returnTreeString(self):
        """
        Returns a tree string containing the topology of the tree.
        Used as input for RSPR.
        """
        sites=[153,419,442,469,515,725,738,1040,1523,2066,3039,3532,3738,6736,6856,7216,8034,8141,8194,9215]
        outString=[]
        cur_node=self.root.leftChild
        seen=[]
        root=False
 #       allSNPs={x:'' for x in codons}
        while root==False:
            if len(outString)!=0 and outString[-1]==');':
                root=True
            if isinstance(cur_node,node):
                while cur_node.leftChild.index in seen and cur_node.rightChild.index in seen:
                    if cur_node.index=='Root':
                        root=True
                        outString.append(');')
                        break
                    else:
                        outString.append(')')
                        temp_out=['rate=%s,c0=\"X_NNN_NNN\"'%(cur_node.rate)]

                        for st in sites:
                            occupied=''
                            if len(cur_node.N)>0:
                                for Nsnp in cur_node.N:
                                    if Nsnp[0]==st:
                                        occupied=Nsnp[0]
                                        temp_out.append('c%05d=\"N_%s_%s\"'%(Nsnp[0],Nsnp[2],Nsnp[3]))

                            if len(cur_node.S)>0:
                                for Ssnp in cur_node.S:
                                    if Ssnp[0]==st:
                                        occupied=Ssnp[0]
                                        temp_out.append('c%05d=\"S_%s_%s\"'%(Ssnp[0],Ssnp[2],Ssnp[3]))
                            if occupied=='':
                                temp_out.append('c%05d=\"X_NNN_NNN\"'%(st))

                        outString.append(':[&%s]%.4f'%(','.join(temp_out),cur_node.branch))
                                
                        cur_node=cur_node.parent

                    if cur_node.index=='Root':
                        root=True
                        outString.append(');')
                        break

                if cur_node.leftChild.index in seen or cur_node.index in seen:
                    cur_node=cur_node.rightChild
                    if root!=True:
                        outString.append(',')
                else:
                    outString.append('(')
                    seen.append(cur_node.index)
                    cur_node=cur_node.leftChild

            if isinstance(cur_node,leaf):
                outString.append(str(cur_node.numName))

                temp_out=['rate=%s,c0=\"X_NNN_NNN\"'%(cur_node.rate)]


                for st in sites:
                    occupied=''
                    if len(cur_node.N)>0:
                        for Nsnp in cur_node.N:
                            if Nsnp[0]==st:
                                occupied=Nsnp[0]
                                temp_out.append('c%05d=\"N_%s_%s\"'%(Nsnp[0],Nsnp[2],Nsnp[3]))
                    if len(cur_node.S)>0:
                        for Ssnp in cur_node.S:
                            if Ssnp[0]==st:
                                occupied=Ssnp[0]
                                temp_out.append('c%05d=\"S_%s_%s\"'%(Ssnp[0],Ssnp[2],Ssnp[3]))
                    if occupied=='':
                        temp_out.append('c%05d=\"X_NNN_NNN\"'%(st))                          

                outString.append(':[&%s]%.4f'%(','.join(temp_out),cur_node.branch))

                
                seen.append(cur_node.index)
                cur_node=cur_node.parent

        return ''.join(outString)

###################################################

global object_list
object_list=[]
global maxHeight
maxHeight=0
global i
global ll
ll=None

def parseTreeFile(argv):
    """
    Parse mode argument.
    """
    try:
        opts, args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
        print 'Hrothvitnir.py -i <inputfile>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'Hrothvitnir.py -i <inputfile>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            treefile = arg
            sys.stderr.write('tree file is: %s'%(treefile))

    global trait_list
    trait_list=[]
    global trait_slices
    trait_slices=[]
    global trait_names
    trait_names=[]
    global partition_name
    partition_name=''
    global slice_list
    slice_list=[]
    
    global object_list
    object_list=[]
    global maxHeight
    maxHeight=0
    global i
    
    global ll
    ll=None

    global tips
    tips=[]
    tipNum=0

    taxonlist=False
    plate=True
    for line in open(treefile,'r'):
        if plate==True and 'state' not in line.lower():
 #           sys.stdout.write('\n%s'%(line.strip('\n')))
            ## Extract useful information from the bits preceding the actual trees.
            cerberus=re.search('Dimensions ntax\=([0-9]+)\;',line)
            if cerberus is not None:
                tipNum=int(cerberus.group(1))

            if 'Translate' in line:
                taxonlist=True
        
            if taxonlist==True and ';' not in line and 'Translate' not in line:
                try:
                    ## Identifies sequence names given at the beginning of the *.trees.txt file
                    tip_index=int(line.strip('\n').strip('\t').split(' ')[0].strip('\''))
                    tip=line.strip('\n').strip('\r').strip('\t').split(' ')[1].strip(',')
                    if [tip_index,tip] not in tips:
                        tips.append([tip_index,tip])
                except ValueError:
                    pass

        if 'tree STATE_0' in line:
            plate=False
            assert (tipNum == len(tips)),'Expected number of tips: %s\nNumber of tips found: %s'%(tipNum,len(tips))
            global mostrecent

        cerberus=re.match('tree\sSTATE\_([0-9]+).+\[\&R\]\s',line)
        if cerberus is not None:
            #sys.stderr.write('\n%s\t%s'%(len(tips),cerberus.group(1)))
            ## At tree state 0 insert header into output file
            if int(cerberus.group(1))==0:
                i=0
                object_list=[]
                maxHeight=0
                ll=None
                ll=tree()
                start=len(cerberus.group())
                treestring=str(line[start:])
                make_tree(treestring)
                #######################################################
                cods=[]
                for k in object_list:
                    if len(k.S)>0:
                        for x in k.S:
                            if x[0] not in cods:
                                cods.append(x[0])
                            
                    if len(k.N)>0:
                        for y in k.N:
                            if y[0] not in cods:
                                cods.append(y[0])

                sys.stdout.write('states\tN_codon')
                sys.stdout.write('\tN_codon'.join(map(str,range(1,9786))))
                sys.stdout.write('\tS_codon')
                sys.stdout.write('\tS_codon'.join(map(str,range(1,9786))))
                #######################################################

            ## After burnin start processing
            if int(cerberus.group(1)) >= burnin:
                ## i is a global variable used to parse the tree string
                ## object_list contains all the objects in the tree
                ## ll is the tree object
                
                i=0
                object_list=[]
                maxHeight=0
                ll=None
                ll=tree()
                start=len(cerberus.group())
                treestring=str(line[start:])
                make_tree(treestring)

                ## Traverse the tree - sets the height of each object in the tree
                ll.traverse_tree()

                ## Rename tips so their name refers to sequence name
                ll.renameTips()

                ## Calibrate tree - find the height of the most recent tip and its date
                calTipHeight,calTipDate,mostRecentDate=mostRecent(tips)
                ll.setAbsoluteTime(calTipHeight,calTipDate)

                ################################################################################
                highestTip=max([k.height for k in object_list])
                sys.stdout.write('\n%s\t'%cerberus.group(1))

                allMu=[]
                for k in object_list:
                    ## only allow one mutation per site per branch and get rid of double mutations to the same codon
                    # the same site can be hit by a mutation on the same branch stochastically due to CTMC or because of the way BEAST outputs codons with two substitutions
                    allMu+=['N_%s_%s_%s_%s'%(x) for x in k.N if column(k.N,0).count(x[0])==1 and '%s_%s_%s'%(x[0],x[2],x[3]) not in ['%s_%s_%s'%(y[0],y[2],y[3]) for y in k.S]]
                    allMu+=['S_%s_%s_%s_%s'%(x) for x in k.S if column(k.S,0).count(x[0])==1 and '%s_%s_%s'%(x[0],x[2],x[3]) not in ['%s_%s_%s'%(y[0],y[2],y[3]) for y in k.N]]


                sys.stdout.write('\t'.join(allMu))

            if 'End;' in line:
                print '\nEnd of file!'
                sys.stdout.write('End;')

######################################################

global burnin

######################################################

## Define burnin as the number of the state from which to begin analysis.
burnin=10000000

######################################################

import getopt

if __name__ == "__main__":
   parseTreeFile(sys.argv[1:])

######################################################
sys.stderr.write('\n>>>>> Done! <<<<<\n')
