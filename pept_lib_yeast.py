#! usr/bin/python

# Script to parse the yeast sgd transl file (in fasta format) and obtain the peptide library. 

import os
import sys
import re
import numpy as np
import pylab as p
import matplotlib
import pickle


def open_file(name_file, open_status = 'r'):
    """ This function just opens the file for reading and returns back the file handle for the file. If the file cannot be opened it just exits! It can open and read any type of files. It can also return back the file handle for writing. The default status for opening is to read the file. Note that the file handle is and must be closed in the place from where this function was called """

    #Opening and reading/writing the passed file_name """
    try:
	file = open(name_file,open_status)
    except IOError,err:	# If the file cannot be opened i am exiting
	print "File %s cannot be opened : %s "%(name_file,err.strerror)
	sys.exit(0)

    return file

def get_gi_seq_dict(lines):
    """
    @param lines: The lines in the non redundant proteome data base (for yeast this time). 
    @function: Generates a key:value dictionary pair; gi id is the key and the value is the sequence
    """
    #Initializing
    flag = False
    all_sequences = []
    sequence = ''
    gi_seq_dictionary = {}
    
    # (1) Iterate lines and get gi
    for line in lines:
	if line.startswith('>'):
	    # (2) Replace the sequence \n with nothing
	    sequence = sequence.replace('\n','')
	    # (3) Obtain gi
	    gi = line.split('|')[1]
	    # (4) Making the key:value dictionary
	    sequence = sequence.upper()		#Ensuring the sequence is in uppercase
	    gi_seq_dictionary[gi] = sequence
	    sequence = ''
	if not line.startswith('>'):	# This is to add all lines following the first line of the fasta file delimiter '>'
	    sequence += line

    return gi_seq_dictionary

def convert_gi_seq_dict_to_peptides(gi_seq_dict,cut_site):
    """
    @param gi_seq_dict: This is the key:value pair of protein_gi and sequence
    @param cut_site: It tells the site at which the sequence should be cut or split
    @function: Breaks the sequence according to the cut site. Relabells K and R to 1 and all others to 0 for each of the peptide fragment. 
    """
    # Initializing
    bin_peptide_lib = []
    gi_bin_pept_dict = {}
    
    # (1) Iterate through each key:value pair
    for gi in gi_seq_dict.keys():
	# Initializing
	gi_pept = ''
	# (2) Obtain sequence
	sequence = gi_seq_dict[gi]
	# (3) Break sequence according to cut_site
	assert len(cut_site) == 1	# Ensure that only one amino acid letter is provided
	peptides_list = sequence.split(cut_site)
	# (4) For each member in the peptide_list, convert K and R into 1. Rest all to 0
	for peptide in peptides_list:
	    peptide = re.sub('[K,R]','1',peptide)
	    peptide = re.sub('C','2',peptide)
	    peptide = re.sub('[^1,^2]','0',peptide)
	    # (5) Adding to the peptide library
	    bin_peptide_lib.append(peptide)
	    # (6) Adding to gi_pept string
	    gi_pept += peptide + '|'

	    gi_bin_pept_dict[gi] = gi_pept

    return gi_bin_pept_dict, bin_peptide_lib
    
def main(argument):
    
    # Initializing
    all_counts = []
    peptide_count_dictionary = {}
    
    # (1) Open and read lines of fasta file containing the translated sequence
    file_name = 'yeast_nrpep.fasta'
    ifile = open_file(file_name)
    lines = ifile.readlines()
    ifile.close()
    
    # (2) Parse the file to retrieve the gi:sequence key-value pair
    gi_sequence_dictionary = get_gi_seq_dict(lines)
    
    # (3) Convert the sequence of each protein gi into shorter peptides and translate to 1s and 0s
    cut_site = 'P'	# Can cut at the site as needed
    converted_gi_peptide_dictionary, bin_peptide_lib = convert_gi_seq_dict_to_peptides(gi_sequence_dictionary,cut_site)
    
    # (4) Make a non redundant set of all the entries in the binary peptide library
    nr_set_bin_pept_lib = set(bin_peptide_lib)
    print len(nr_set_bin_pept_lib)
    print len(bin_peptide_lib)
    print bin_peptide_lib[1]
    
    # (5) Writing output of peptide - count; Checking if it exists.
    out_file = 'bin_pept_lib_yeast_split_cys_lab.pro'
    if os.path.exists(out_file):
	os.remove(out_file)
    
    # (6) Obtain counts of each peptide
    for peptide in nr_set_bin_pept_lib:
	count = bin_peptide_lib.count(peptide)
	ofile = open_file(out_file,'a')
	ofile.write(peptide + '\t' + str(count) + '\n')
	all_counts.append(count)
	# (7) Making a key dictionary pair for peptide:counts
	peptide_count_dictionary[peptide] = count
    
    ofile.close()
    
    # (8) Pickling the dictionary
    out_pkl_file = open('pept_lib_yeast_dict_counts_cyslab_pro.pkl','wb')
    pickle.dump(peptide_count_dictionary,out_pkl_file)
    out_pkl_file.close()
	
    # (9) Make a frequency histogram of the bin_peptide_library
    bins = 500
    fig = p.figure(figsize=(15,15))
    ax1 = fig.add_subplot(1,1,1)
    n,bins,patches = p.hist(all_counts,bins,color='g')
    p.xlabel('# Counts ')
    p.ylabel('# Peptides')
    p.title('Frequency histogram of number of counts of each peptide in the converted peptide library ')
    fig_dir = os.getcwd() + '/figures/'
    fig_name = fig_dir + 'yeast_converted_peptide_freq.png'
    p.show()
    p.savefig(fig_name,format = "png", orientation='landscape')

if __name__ == '__main__':
    
    argument = sys.argv[1:]
    
    print "Processing %s ..."%(argument)
    main(argument)
    
    import time
    print "Script - pept_lib_yeast.py %s \t Completed \t %s"%(argument, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
    
