#! usr/bin/python

from __future__ import division
import numpy as np
import pylab as p
import matplotlib
import os
import sys
import pickle



def main(argument):
    
    ## (1) Open the pickle file to get the dictionary of peptide : peptide count
    #pkl_file = open('pept_lib_yeast_dict_counts_cyslab_glu.pkl','rb')
    #pept_count_dict = pickle.load(pkl_file)
    # Initializing
    all_lengths = []
    
    # (1) Open the txt file that contains the peptide: peptide count 
    file_name = 'bin_pept_lib_yeast_split.glu'
    ifile = open(file_name,'r')
    lines = ifile.readlines()
    ifile.close()
    
    # (2) Obtain the peptide length for all entries where the count was 1
    for line in lines:
	if (line.split('\t')[1])[:-1] == '1':
	    peptide = line.split('\t')[0]
	    len_peptide = len(peptide)
	    all_lengths.append(len_peptide)

    print max(all_lengths)
    print min(all_lengths)
    print sum(all_lengths)
    print len(all_lengths)
    print sum(all_lengths)/len(all_lengths)
    
    # (3) Make a frequency histogram of the length vs occurrences
    bins = 500
    fig = p.figure(figsize=(15,15))
    ax1 = fig.add_subplot(1,1,1)
    n,bins,patches = p.hist(all_lengths,bins,color='b')
    p.xlabel('# Length of peptide ')
    p.ylabel('# Occurrences')
    p.title('Frequency histogram of length of peptide when count=1; GluC cut - Lys,Arg labeled')
    fig_dir = os.getcwd() + '/figures/'
    fig_name = fig_dir + 'len_distr_count_1_glucut.png'
    p.savefig(fig_name,format = "png", orientation='landscape')
    



if __name__ == '__main__':
    
    argument = sys.argv[1:]
    
    print "Processing %s ..."%(argument)
    main(argument)
    
    import time
    print "Script - len_distr.py %s \t Completed \t %s"%(argument, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
    

