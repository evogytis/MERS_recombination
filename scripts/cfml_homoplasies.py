################################################################################################################################################################
##
##  This script has been written and kindly donated by Jessica Hedge (University of Oxford).
##  It has been modified by Gytis Dudas (University of Edinburgh) such that it reports all mutations, not just homoplasies and reports asymmetric mutations
##
##  16 March, 2015
##
################################################################################################################################################################


from __future__ import print_function
from __future__ import division

import sys
import argparse
import dendropy
import numpy as np
from Bio import AlignIO
from itertools import combinations
from itertools import permutations

"""
Reads in output files from ClonalFrameML and returns a list of
homoplasies, giving:
site number, node, substitution, total #homoplasies at site
Also returns list of branches that harbour substitution at site provided in command line
E.g.:
cfml_homoplasies.py PREFIX --site 20 --substitution AT
"""

def getSeqLen(f):
	"""
	Get alignment length
	"""
	for rec in f:
		seqLen = len(rec.seq)
		return seqLen

def nodeIndex2Name(t):
	"""
	Map node index in array B to node name
	"""
	nodeNames = []
	for node in t.preorder_node_iter():
		if node.level() != 0: # if not root
			nodeNames.append(node.get_node_str())
	nodeNames = np.array(nodeNames)
	return nodeNames

def patternIndex2site(cr):
	"""
	Map pattern index to list of alignment site indices
	"""
	crossDict = {}
	for site, pattern in enumerate(cr):
		patternIndex = pattern - 1
		if patternIndex >= 0:
			if patternIndex in crossDict:
				crossDict[patternIndex].append(site)
			else:
				crossDict[patternIndex] = [site]
	return crossDict

class HomoplasyProcessor(object):
    """
    Class to read the output from clonalFrameML and write two files.
    The first is a list of homoplastic sites and the second is a list
    of branches on which the given substitution occurs.
    """
    def __init__(self, prefix):
        self.prefix = prefix
        labelledTreeFile = prefix + ".labelled_tree.newick" # labelled ML unrooted tree
        MLfastaFile = prefix + ".ML_sequence.fasta" # imputed sequence data
        crossRefFile = prefix + ".position_cross_reference.txt" # cross referencing patterns with sites
        # Read in the files and create the objects we need
        self.tree = dendropy.Tree.get_from_path(labelledTreeFile, schema="newick")
        self.alignment = AlignIO.read(MLfastaFile, format = "fasta")
        self.createNameSequenceMap()
        self.createNodeMap()
        seqLen = getSeqLen(self.alignment)
        #self.subTypes = np.array(["".join(sorted(i)) for i in combinations('ACGT-', 2)]) # substitution combinations
        self.subTypes = np.array(["".join(i) for i in permutations('ACGT-', 2)]) # substitution combinations
        #sys.stderr.write('\nsubTypes:%s\n\n'%self.subTypes)
        # B is a 3D array for multiple substitutions per site (node x pattern x substitution)
        self.B = np.zeros((len(self.alignment)-1, seqLen, self.subTypes.shape[0]), dtype = "int")
        self.fillArray()
        self.nodeNames = nodeIndex2Name(self.tree)
        # read in list for cross-referencing patterns with sites
        with open(crossRefFile, 'r') as f:
            self.cr = [int(s) for s in f.readline().rstrip().rsplit(",")]
        self.crossDict = patternIndex2site(self.cr)

    def createNameSequenceMap(self):
        """
        Map record name to sequence
        """
        self.nameSequenceMap = {}
        for rec in self.alignment:
            self.nameSequenceMap[rec.name] = np.array(rec.seq)

    def createNodeMap(self):
        """
        Map node to sequence
        """
        self.nodeMap = {}
        for node in self.tree.preorder_node_iter():
            node.sequence = self.nameSequenceMap[node.get_node_str()]
            self.nodeMap[node] = node.sequence

    def fillArray(self):
        """
        Fill array B with counts of substitutions across tree
        """
        t = self.tree
        B = self.B
        row = 0
        nodeNames = []
        for node in t.preorder_node_iter():
            if node.level() != 0: # if not root
                nodeNames.append(node.get_node_str())
                dec = node.sequence
                anc = node.parent_node.sequence
                diffs = np.nonzero(anc != dec)
                decBases = dec[diffs]
                ancBases = anc[diffs]
                subs = np.array(["".join([ancBases[i],decBases[i]]) for i in range(len(decBases))])
                #subs = np.array(["".join(sorted([decBases[i], ancBases[i]])) for i in range(len(decBases))])
                #sys.stderr.write('\nsubs:%s'%(subs))
                subIndices = [np.where(self.subTypes == s)[0][0] for s in subs]

                #sys.stderr.write('\nsubIndices:%s'%(subIndices))
                B[row, diffs, subIndices] +=1
                row += 1

    def writeSubstitutions(self, siteNumber, substitution):
        """
        Writes list of branch with substitution at given site
        """
        formattedSub = "".join(sorted([i for i in substitution]))
        outFilename = self.prefix + ".site" + str(siteNumber) + ".branchSubs.txt"
        out = open(outFilename, "wt")
        header = "branch\tsubstitution\tsite\n"
        out.write(header)
        patternIndex = self.cr[siteNumber-1]-1
        for node in range(self.B.shape[0]):
            subsPresentList = np.nonzero(self.B[node, patternIndex, :] >0)[0]
            if len(subsPresentList) > 0:
                nodeName = self.nodeNames[node]
                subType = self.subTypes[subsPresentList]
                if subType == substitution:
                    out.write(nodeName)
                    out.write("\t")
                    out.write(subType)
                    out.write("\t")
                    out.write(str(siteNumber))
                    out.write("\n")
        out.close()
        print("Distribution of substitution " + substitution + " at site "
                + str(siteNumber) + " across tree saved in " + outFilename)

    def writeHomoplasies(self):
        """
        Writes all homoplasies to file
        """
        hpCount = 0
        outFilenameSites = self.prefix + ".homoplasies.txt"
        outSites = open(outFilenameSites, "wt")
        outSites.write("site\tnode\tsubstitution\ttotal_#nodes_with_substitution\n")
        for pattern in range(self.B.shape[1]):
            sites = self.crossDict[pattern]
            for sub in range(self.B.shape[2]):
                substitution = self.subTypes[sub]
                v = self.B[:, pattern, sub]
                if sum(v) >0:
                    nodes = np.nonzero(self.B[:, pattern, sub] >0)[0]
                    hpCount += len(sites)
                    for node in nodes:
                        nodeName = self.nodeNames[node]
                        for site in sites:
                            row = str(site + 1) + "\t" + str(nodeName) + "\t" + str(substitution) + "\t" + str(sum(v)) + "\n"
                            outSites.write(row)
        outSites.close()

        print("\nThere are " + str(hpCount) + " homoplastic sites in the alignment")
        print("These are saved in " + outFilenameSites)


def main():
    parser = argparse.ArgumentParser("Process ClonalFrameML files")
    parser.add_argument("prefix", metavar="PREFIX", nargs=1,
            help="The input file prefix")
    parser.add_argument("--site", metavar="SITE", type=int, default=-1)
    parser.add_argument("--substitution", metavar="SUBSTITUTION", default=None)
    args = parser.parse_args()
    hp = HomoplasyProcessor(args.prefix[0])
    hp.writeHomoplasies()
    if args.site != -1:
        hp.writeSubstitutions(args.site, args.substitution)

if __name__ == "__main__":
    main()

