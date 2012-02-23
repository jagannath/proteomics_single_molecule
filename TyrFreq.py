#! usr/bin/python

# This programme tries to capture the frequency of Tyrosines and Histidines in the Human proteome. Maybe a histogram of the frequency for proteins and then for tryptic peptides. 
from __future__ import division
import os
import sys
import pdb
import pickle
import shelve
import numpy as np
import pylab as p
import re
import matplotlib


class Protein:
    def __init__(self,protID,seq):
	self.protID = protID
	self.seq = seq
	self.seqLength = len(seq)
    def proteinFrequency(self,aa):
	return ((self.seq.count(aa)/self.seqLength)*100)
    def cutProtein(self,cut1,cut2):
	seq = self.seq
	cut1Pept = seq.split(cut1)
	allcutPept = []
	for pept in cut1Pept:
	    cut2Pept = pept.split(cut2)
	    allcutPept.extend(cut2Pept)
	return allcutPept
    def peptideFrequency(self,cut1='E',cut2='X',aa='Y'):
	allPeptideFreq = []
	peptides = self.cutProtein(cut1,cut2)
	allPeptideFreq = [(item.count(aa)/len(item))*100 for item in peptides if len(item) > 20]
	#print peptides
	return allPeptideFreq
	
def parse_seq(inseq):
    """
    @param inseq: This is a sequence information from the non_redundant proteome data set; This must be parsed
    @function: Parse the sequence to return back only the sequence information
    """
    inseq_lines = inseq.split('\n')
    sequence = ''
    for lines in inseq_lines:
	if not lines.startswith('SQ'): 
	    if not '//' in lines:
		sequence+= lines
    sequence = sequence.replace(' ','')
    return sequence

def getHumanIDSeq(fname):
    ifile = open(fname)
    lines = ifile.readlines()
    ifile.close()
    flag = False
    sequence = ''
    gi_seq_dictionary = {}
    inseq = ''
    for line in lines:
	if line.startswith('ID'):
	    gi = line.split('  ')[1]
	    gi = gi.replace(' ','')
	if line.startswith('SQ'):
	    flag = True
	if (flag):
	    inseq += line
	if line.startswith('//'):
	    flag = False
	    sequence = parse_seq(inseq)
	    gi_seq_dictionary[gi] = sequence
	    inseq = ''
    return gi_seq_dictionary

def loadRelevantFiles(org = 'human'):
    if org == 'human':
	fname = '/project/marcotte/jagannath/projectfiles/Downloads/human_proteome/uniprot_sprot_human.dat'
	idSeq_dict = getHumanIDSeq(fname)
    
    return idSeq_dict

def main(argument):
    allProteinFreq = []
    allPeptideFreq = []
    for protID, seq in idSeq_dict.items():
	prot = Protein(protID,seq)
	if prot.proteinFrequency('H')>25:
	    print protID
	    
	allProteinFreq.append(prot.proteinFrequency('Y'))
	allPeptideFreq.extend(prot.peptideFrequency(cut1='K',cut2='R',aa='Y'))

    #Drawing a Histogram
    p.rcParams.update({'font.size':14})
    fig = p.figure(figsize=(10,10))
    p.hist(allProteinFreq,100)
    p.xlabel('Percentage Histidines in Proteins')
    p.ylabel('Number of Proteins')
    #p.savefig('HisFrequencyHumanProteome_histogram.png')
    p.show()

    print "Protein frequency mean (His) = %f"%(np.mean(allProteinFreq))
    print "Peptide frequency mean = %f"%(np.mean(allPeptideFreq))
    print "Protein Frequency max, min : %f , %f "%(max(allProteinFreq), min(allProteinFreq))
if __name__ == '__main__':
    
    argument = sys.argv[1:]
    
    print "Processing %s ..."%(argument)
    idSeq_dict = loadRelevantFiles('human')
    
    main(argument)
    
    import time
    print "Script - %s \t Completed \t %s"%(sys.argv, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))