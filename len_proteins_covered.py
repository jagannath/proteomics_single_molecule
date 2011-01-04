#! usr/bin/python

# Script obtains a list of all unique peptides. Makes a cut off based on their length and gets a filtered list of all peptides less than the given length. I then make a dictionary of binary peptide: gi (by comparing the whether the unique bin peptide was in the gi). Then I can obtain the number of gis considered based on the length cut off. 

import os
import sys
import pickle
import re


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
    @function: Breaks the sequence according to the cut site. Relabells K and R to 1 and all others to 0 for each of the peptide fragment. Makes a gi:peptide(converted binary) dictionary
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


def get_unique_bin_pept_list(bin_peptide_count_dict):
    """
    @param bin_peptide_count_dict: A dictionary with key: value = binarypeptide: count
    @function: Runs through the key value pair to generate a list of all binary peptides with count = 1 (i.e. unique peptides)
    """
    # Initializing
    unique_bin_peptides_list = []
    
    # (1) Iterate through the dictionary
    for key, values in bin_peptide_count_dict.items():
	# (2) Check if value = 1. Then append to list
	if values == 1:
	    unique_bin_peptides_list.append(key)
    
    return unique_bin_peptides_list

def get_unique_bin_peptide_gi_dict(unique_bin_peptide_list,gi_bin_peptide_dict):
    """
    @ param unique_bin_peptide_list: A list of unique peptides
    @ function: Iterates through the list - for each peptide it checks which of the gi it occurs in. Then makes a unique peptide: gi dictionary pair
    """
    # Initializing
    unique_bin_peptide_gi_dictionary = {}
    
    for peptide in unique_bin_peptide_list:
	for key,value in gi_bin_peptide_dict.items():
	    bin_pept_list_for_a_gi = value.split('|')
	    if peptide in bin_pept_list_for_a_gi:
		print peptide, key
		unique_bin_peptide_gi_dictionary[peptide] = key
    
    return unique_bin_peptide_gi_dictionary
    

def get_gi_list_length_cutoff(unique_bin_peptide_gi_dict, cutoff):
    """
    @ param unique_bin_peptide_gi_dict: A dictionary unique_bin_peptide: gi
    @ param cutoff: This is the cut off of the string length. Only bin_peptides less than this length are considered. 
    @ function: Iterates through the unique peptide_list and retrieves only those peptides that have a length of less than equal 20. It obtains the gi and makes a list of the gis. 
    """
    
    # Initializing
    list_gis = []
    
    # (1) Iterate the dictionary to get all the unique_peptides
    for key,value in unique_bin_peptide_gi_dict.items():
	# (2) Check length of the peptide. Append the gi to the list
	if len(key) <= cutoff:
	    list_gis.append(value)
    
    return list_gis

def main(argument):
    
    # Initializing
    [argv] = argument
    assert len(argument) == 1
    
    if argv == '1':
	# (1) Make a key dictionary pair gi:binary peptide list. 
	# (1a) Open and read lines of fasta file containing the translated sequence
	file_name = 'yeast_nrpep.fasta'
	ifile = open_file(file_name)
	lines = ifile.readlines()
	ifile.close()
	
	# (1b) Parse the file to retrieve the gi:sequence key-value pair
	gi_sequence_dictionary = get_gi_seq_dict(lines)
	
	# (1c) Convert the sequence of each protein gi into shorter peptides and translate to 1s and 0s
	cut_site = 'E'	# Can cut at the site as needed
	converted_gi_peptide_dictionary, bin_peptide_lib = convert_gi_seq_dict_to_peptides(gi_sequence_dictionary,cut_site)
	print "binary peptide dictionary created"
	
	# (1d) Pickle the dictionary
	pkl_file_name = 'gi_binpeptide_dict_glucut.pkl'
	ofile = open(pkl_file_name,'wb')
	pickle.dump(converted_gi_peptide_dictionary,ofile)
	ofile.close()
	print "pickled"
	
	# (2) Unpickle the dictionary file containing the binarypeptide:count
	pkl_file_name = 'pept_lib_yeast_dict_counts_cyslab_glu.pkl'
	ifile = open(pkl_file_name)
	bin_peptide_count_dict = pickle.load(ifile)
	ifile.close()
	
	# (3) Obtain a list of all unique binary peptides (i.e with count 1)
	unique_bin_peptides_list = get_unique_bin_pept_list(bin_peptide_count_dict)
	
	# (4) Generate a dictionary pair unique_bin_peptide:gi
	unique_bin_peptide_gi_dictionary = get_unique_bin_peptide_gi_dict(unique_bin_peptides_list, converted_gi_peptide_dictionary)
	
	# (5) Pickle the dictionary
	pkl_file_name = 'unique_bin_peptide_gi_dictionary_glucut.pkl'
	ofile = open(pkl_file_name,'wb')
	pickle.dump(unique_bin_peptide_gi_dictionary, ofile)
	ofile.close()
	print "Pickle of unique bin peptide : gi dictionary done! "

    if argv == '2':
	# (6) Unpickle the unique_bin_peptide_gi_dictionary_glucut
	pkl_file_name = 'unique_bin_peptide_gi_dictionary_glucut.pkl'
	ifile = open(pkl_file_name)
	unique_bin_peptide_gi_dictionary = pickle.load(ifile)
	ifile.close()
	
	# (7) Generate list of all gis (redundant list) for a determined cut off length of sequence
	length_cutoff = 20
	filtered_gi_list = get_gi_list_length_cutoff(unique_bin_peptide_gi_dictionary,length_cutoff)
	print "Number of redundant gis %d"%(len(filtered_gi_list))
	print "Number of non-redundant gis covered %d"%(len(set(filtered_gi_list)))


if __name__ == '__main__':
    
    argument = sys.argv[1:]
    
    print "Processing %s ..."%(argument)
    main(argument)
    
    import time
    print "Script - len_proteins_covered.py %s \t Completed \t %s"%(argument, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
    