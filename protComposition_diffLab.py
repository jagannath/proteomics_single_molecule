#!/usr/bin/python

# This programme simulates the (1) Cut the proteome into peptides (2) Transform peptides into #labels for (D&E;Y;K;S&T) Identify Unique compositions and obtain the proteome uniquely covered. 

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

# Opening Classes
class ProteinID():
    def __init__(self,protID):
	self.protID = protID
	self.seq = idSeq_dict[protID]
    def sequence(self):
	return self.seq
    def cut(self,cutAA='E'):
	self.listPeptides = self.seq.split(cutAA)
	self.listPeptides = [item for item in self.listPeptides if not item == '']
    def anchor(self,anchorAA='C'):
	self.listPeptides = [item for item in self.listPeptides if not item.count(anchorAA) == 0]
	self.anchorAA = anchorAA	#Will use this to assert that the anchor is present
	return len(self.listPeptides)
    def getCounts(self,pept,code):
	if code == 'd': counts = pept.count('D') #End should be an E, which is removed while splitting peptides; So no practical difference; Acidic
	elif code == 'y': counts = pept.count('Y') #Tyrosine label
	elif code == 'k': counts = pept.count('K') #Basic label
	elif code == 's': counts = pept.count('S') + pept.count('T') # Alcohol label Serine and Threonine
	else: print "Improper code"
	return counts
    def transform(self,labelCodes):
	assert len(self.listPeptides)>0
	labelList = list(labelCodes)
	trPeptList = []
	for pept in self.listPeptides:
	    assert pept.count(self.anchorAA)>0 
	    trPept = ''
	    for code in labelList:
		counts = self.getCounts(pept,code)
		trPept += str(counts)+code+'-'
	    trPeptList.append(trPept)
	return trPeptList
    



def getYeastIDSeq(fname):
    ifile = open(fname)
    lines = ifile.readlines()
    ifile.close()
    seq = ''
    count = 0
    protIDSeq_dict = {}
    for line in lines:
	if line.startswith('>'):
	    count +=1 
	    if count>1:
		seq = seq.replace('\n','')
		seq = seq.upper()[:-1]	#Ensuring all sequences are in upper case and removing an annoying trailing *
		protIDSeq_dict[protID] = seq
		seq = ''
	    protID = line.split(' ').pop(0)[1:]	    
	    flag = True
	else:
	    seq+=line
    return protIDSeq_dict

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

def getCDR3(fname):
    ifile = open(fname)
    lines = ifile.readlines()
    ifile.close()
    count = 0
    cdr3_seq_dict = {}
    for line in lines:
	count+=1 
	index = 'cdr3_'+str(count)
	seq = line.split('\t')[3]
	cdr3_seq_dict[index] = seq

    return cdr3_seq_dict

def pickle_file(file_name, contents, path_dir=os.getcwd()):
    """
    @param file_name: This is the file name to which the contents must be dumped
    @param contents: This is the contents to be pickled
    @function: The contents passed are pickled to the file in the pkl_directory
    """
    #pkl_dir = os.getcwd() + '/pkl_files/'
    pkl_dir = path_dir + '/pkl_files/'
    pkl_file = pkl_dir + file_name
    ofile = open(pkl_file,'wb')
    pickle.dump(contents, ofile)
    ofile.close()
    
    return True

def unpickle_file(file_name, path_dir=os.getcwd()):
    """
    @param file_name: This is the file name to which the contents must be dumped
    @param contents: This is the contents to be pickled
    @function: The contents passed are pickled to the file in the pkl_directory
    """
    pkl_dir = path_dir + '/pkl_files/'
    pkl_file = pkl_dir + file_name
    ifile = open(pkl_file)
    contents = pickle.load(ifile)
    ifile.close()
    
    return contents

def loadRelevantFiles(org='yeast'):
    if org == 'yeast':
	fname = '/project/marcotte/jagannath/projectfiles/Downloads/yeast_proteins_sgd/orf_trans.fasta'
	# Get ID: Sequence dictionary
	idSeq_dict = getYeastIDSeq(fname)
	print len(idSeq_dict)
	pickle_file('yeastIDSeq_dict.pkl',idSeq_dict)
    if org == 'human':
	fname = '/project/marcotte/jagannath/projectfiles/Downloads/human_proteome/uniprot_sprot_human.dat'
	idSeq_dict = getHumanIDSeq(fname)
    if org == 'cdr3':
	fname = '/project/marcotte/jagannath/projectfiles/Downloads/human_proteome/CDR3withFR4_rabbit.csv'
	idSeq_dict = getCDR3(fname)

    return idSeq_dict

def testIT():
    testProtID = 'YMR065W'
    test = ProteinID(testProtID)
    print testProtID, test.sequence()
    test.cut()
    print test.anchor()
    print test.transform('dkys')

def main(argument):
    trPeptIDs_dict = {}
    alltransformedPeptides = []
    counter = 0
    for protID in idSeq_dict.keys():
	counter += 1
	protein = ProteinID(protID)
	#print protID, protein.sequence()
	protein.cut('Z')
	if protein.anchor('G'):
	    transformedPeptides = protein.transform('dyks')	#Each letter code in small represents a label attaching to a defined AA
	    #transformedPeptides = protein.fiterTransformed(False,nbrLabels=10)
	    for trPept in transformedPeptides:
		if trPept in alltransformedPeptides:
		    protIDList = trPeptIDs_dict[trPept]
		    protIDList.append(protID)
		    trPeptIDs_dict[trPept] = protIDList
		else:
		    trPeptIDs_dict[trPept] = [protID]
		    alltransformedPeptides.append(trPept)
    
    print len(alltransformedPeptides)
    print len(idSeq_dict)
    uniqProteins = []
    allSizeProtID = []
    for trPept, protID in trPeptIDs_dict.items():
	#print trPept, protID, len(set(protID))
	if len(set(protID)) == 1: uniqProteins.extend(protID)
	allSizeProtID.append(len(set(protID)))
    
    print len(set(uniqProteins))
    histNbrs = []
    sumTotal = 0
    
    for i in xrange(1,15):
	histNbrs.append(allSizeProtID.count(i))
	sumTotal+= allSizeProtID.count(i)

    histNbrs.append(len(allSizeProtID)-sumTotal)
    
    print histNbrs

    p.rcParams.update({'font.size':14})
    fig = p.figure(figsize=(10,10))
    p.plot(histNbrs)
    p.xlabel('Proteins Identifiable')
    p.ylabel('Number of peptides - transformed by compositional code')
    #p.savefig('HumanPeptideCompositional_histogram.png')
    p.show()

if __name__ == '__main__':
    
    argument = sys.argv[1:]
    
    print "Processing %s ..."%(argument)
    idSeq_dict = loadRelevantFiles('cdr3')

    main(argument)
    
    import time
    print "Script - %s \t Completed \t %s"%(sys.argv, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))