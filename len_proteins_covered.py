#! usr/bin/python

# Script obtains a list of all unique peptides. Makes a cut off based on their length and gets a filtered list of all peptides less than the given length. I then make a dictionary of binary peptide: gi (by comparing the whether the unique bin peptide was in the gi). Then I can obtain the number of gis considered based on the length cut off. 
from __future__ import division
import os
import sys
import pickle
import re
import numpy as np
import pylab as p
import matplotlib



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
	    gi = (line.split(' ')[0])[1:]
	    # (4) Making the key:value dictionary
	    sequence = sequence.upper()		#Ensuring the sequence is in uppercase
	    sequence = sequence[:-1]		# The end has a trailing *
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
	    #peptide = re.sub('[^1]','0',peptide)
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

def get_perc_non_red_gis(file_name):
    """
    @param file_name: pickle file containing dictionary of unique_bin_peptide: gi
    @function: Opens the dictionary and filters the unique peptide list based on the length of the cut off. It then calculates the number of redundant gis and finally the no of non-redundant gis from the filtered peptide list.
    """
    # Initializing
    perc_non_red_gis_list = []
    
    # (1) Unpickle the unique_bin_peptide_gi_dictionary_glucut
    pkl_file_name = file_name
    ifile = open(pkl_file_name)
    unique_bin_peptide_gi_dictionary = pickle.load(ifile)
    ifile.close()
    
    # (2) Iterate through a range of length_cutoff from 5 to 50. 
    for length_cutoff in xrange(5,51,5):
	# (3) Generate list of all gis (redundant list) for a determined cut off length of sequence
	filtered_gi_list = get_gi_list_length_cutoff(unique_bin_peptide_gi_dictionary,length_cutoff)
	print "Number of redundant gis %d"%(len(filtered_gi_list))
	total_yeast_gis = 11081
	nbr_non_red = len(set(filtered_gi_list))
	perc_gis_covered = (nbr_non_red/total_yeast_gis)*100
	print "Number of non-redundant gis covered %d"%(nbr_non_red)
	perc_non_red_gis_list.append(perc_gis_covered)
	
    return perc_non_red_gis_list

def get_truc_gis_list(file_name):
    """
    @param file_name: Pickle file - dictionary of unique binary peptide : gi
    @function: Truncates each binary peptide to first length_cutoff AA. Retains the gi. Generates a gi list that can be made into a non redundant set
    """
    all_perc_coverage = []
    
    # (1) Unpickling the file and streaming the dictionary
    pkl_file_name = file_name
    ifile = open(pkl_file_name)
    unique_bin_peptide_gi_dictionary = pickle.load(ifile)
    ifile.close()
    
    for cut_off in xrange(5,51,5):
	#Initializing 
	trunc_bin_pept_gi_list = []
	trunc_peptide_list = []
	filtered_trunc_bin_peptide_gi_list = []
	filtered_trunc_peptide_list = []
	all_gi_list = []


	# (2) Iterating the unique peptide list keys and values
	for bin_peptide, gi in unique_bin_peptide_gi_dictionary.items():
	    trunc_peptide = bin_peptide[0:cut_off]
	    trunc_bin_pept_gi_list.append([trunc_peptide, gi])
	    trunc_peptide_list.append(trunc_peptide)
	print "Cut off : %d"%(cut_off)
	print "Number of Truncated peptide library : %d" %(len(trunc_peptide_list))
	    
	# (3) Filter the trunc_peptide_list to contain only unique peptides
	set_trunc_pept_list = set(trunc_peptide_list)
	for peptide in set_trunc_pept_list:
	    if trunc_peptide_list.count(peptide) == 1:
		filtered_trunc_peptide_list.append(peptide)
	
	print "Number of Unique truncated peptides : %d" %(len(filtered_trunc_peptide_list))
	
	# (4) Filter the trunc_bin_pept_gi_list to contain only those peptides in the filtered_trunc_peptide_list
	for trunc_bin_peptide, gi in trunc_bin_pept_gi_list:
	    if trunc_bin_peptide in filtered_trunc_peptide_list:
		filtered_trunc_bin_peptide_gi_list.append([trunc_bin_peptide, gi])
		all_gi_list.append(gi)
    
	# (5) Number of proteins covered
	nbr_proteins_covered = len(set(all_gi_list))
	perc_coverage = (nbr_proteins_covered / 5885)*100
	print "Percentage proteome covered : %f"%(perc_coverage)
	all_perc_coverage.append(perc_coverage)

    return all_perc_coverage

   



def main(argument):
    
    # Initializing
    [argv] = argument
    assert len(argument) == 1
    perc_non_red_gis_list = []
    
    if argv == '1':
	# (1) Make a key dictionary pair gi:binary peptide list. 
	# (1a) Open and read lines of fasta file containing the translated sequence
	file_name = 'orf_trans.fasta'
	ifile = open_file(file_name)
	lines = ifile.readlines()
	ifile.close()
	
	# (1b) Parse the file to retrieve the gi:sequence key-value pair
	gi_sequence_dictionary = get_gi_seq_dict(lines)
	print len(gi_sequence_dictionary)
	
	# (1c) Convert the sequence of each protein gi into shorter peptides and translate to 1s and 0s
	cut_site = 'E'	# Can cut at the site as needed
	#cut_site = 'P'	# Can cut at the site as needed
	converted_gi_peptide_dictionary, bin_peptide_lib = convert_gi_seq_dict_to_peptides(gi_sequence_dictionary,cut_site)
	print "binary peptide dictionary created"
	
	# (1d) Pickle the dictionary
	pkl_file_name = 'gi_binpeptide_dict_glucut_nocys.ver2.pkl'
	ofile = open(pkl_file_name,'wb')
	pickle.dump(converted_gi_peptide_dictionary,ofile)
	ofile.close()
	print "pickled - glucut nocys"
	
	# (2) Unpickle the dictionary file containing the binarypeptide:count
	pkl_file_name = 'pept_lib_yeast_dict_counts_glucut_nocys.ver2.pkl'
	ifile = open(pkl_file_name)
	bin_peptide_count_dict = pickle.load(ifile)
	ifile.close()
	
	# (3) Obtain a list of all unique binary peptides (i.e with count 1)
	unique_bin_peptides_list = get_unique_bin_pept_list(bin_peptide_count_dict)
	
	# (4) Generate a dictionary pair unique_bin_peptide:gi
	unique_bin_peptide_gi_dictionary = get_unique_bin_peptide_gi_dict(unique_bin_peptides_list, converted_gi_peptide_dictionary)
	
	# (5) Pickle the dictionary
	pkl_file_name = 'unique_bin_peptide_gi_dictionary_glucut_nocys.ver2.pkl'
	ofile = open(pkl_file_name,'wb')
	pickle.dump(unique_bin_peptide_gi_dictionary, ofile)
	ofile.close()
	print "Pickle of unique bin peptide : gi dictionary done! "

    if argv == '2':
	
	# (1) Obtain the non-redundant list of gis coverage for different length cutoffs
	file_name1 = 'unique_bin_peptide_gi_dictionary_glucut_labcys.ver2.pkl'
	perc_non_red_gis_list1 = get_perc_non_red_gis(file_name1)
	print "Glucut Nocys :"
	print perc_non_red_gis_list1
	
	file_name2 = 'unique_bin_peptide_gi_dictionary_glucut_nocys.ver2.pkl'
	perc_non_red_gis_list2 = get_perc_non_red_gis(file_name2)
	print "Glucut cyslab :"
	print perc_non_red_gis_list2
	
	file_name3 = 'unique_bin_peptide_gi_dictionary_procut_labcys.ver2.pkl'
	perc_non_red_gis_list3 = get_perc_non_red_gis(file_name3)
	print "Procut Nocys :"
	print perc_non_red_gis_list3
		
	file_name4 = 'unique_bin_peptide_gi_dictionary_procut_nocys.ver2.pkl'
	perc_non_red_gis_list4 = get_perc_non_red_gis(file_name4)
	print "Procut cyslab :"
	print perc_non_red_gis_list4
    
	# (2) Plot graph 
	xlist = [i for i in xrange(5,51,5)]
	
	fig = p.figure(figsize=(15,15))
	ax1 = fig.add_subplot(1,1,1)
	p.xlabel('# Length of peptide ')
	p.ylabel('# Percentage of proteome covered')
	p.title('Frequency histogram of number of counts of each peptide in the converted peptide library ')
	line1 = p.plot(xlist, perc_non_red_gis_list1, 'g-o',linewidth=2.0)
	line2 = p.plot(xlist, perc_non_red_gis_list2, 'b-s',linewidth=2.0)
	line3 = p.plot(xlist, perc_non_red_gis_list3, 'r-^',linewidth=2.0)
	line4 = p.plot(xlist, perc_non_red_gis_list4, 'k-*',linewidth=2.0)
	p.legend((line1,line2,line3, line4),('E cut; K,R labeled', 'E cut; K, R, C labeled', 'P cut; K, R labeled','P cut; K, R, C labeled'), loc=4)
	fig_dir = os.getcwd() + '/figures/'
	fig_name = fig_dir + 'proteome_coverage_length_all.png'
	p.savefig(fig_name,format = "png", orientation='landscape')
	#p.show()

    if argv == '3':
	# This considers the unique set of binary peptides. It then parses it to obtain a truncated list of the 1st 20 code. Then checks if they are unique and among the unique codes - checks how much of the proteome is covered. 
	file_name1 = 'unique_bin_peptide_gi_dictionary_glucut_labcys.ver2.pkl'
	trunc_perc_nr_gis_list1 = get_truc_gis_list(file_name1)
	print trunc_perc_nr_gis_list1
	
	file_name2 = 'unique_bin_peptide_gi_dictionary_glucut_nocys.ver2.pkl'
	trunc_perc_nr_gis_list2 = get_truc_gis_list(file_name2)
	print trunc_perc_nr_gis_list2
	
	file_name3 = 'unique_bin_peptide_gi_dictionary_procut_labcys.ver2.pkl'
	trunc_perc_nr_gis_list3 = get_truc_gis_list(file_name3)
	print trunc_perc_nr_gis_list3
	
	file_name4 = 'unique_bin_peptide_gi_dictionary_procut_nocys.ver2.pkl'
	trunc_perc_nr_gis_list4 = get_truc_gis_list(file_name4)
	print trunc_perc_nr_gis_list4
	
	# (2) Plot graph
	xlist = [i for i in xrange(5,51,5)]
	fig = p.figure(figsize=(15,15))
	ax1 = fig.add_subplot(1,1,1)
	p.xlabel('# Length of peptide ')
	p.ylabel('# Percentage of proteome covered')
	p.title('Proteome coverage as a function of sequenced length ')
	line1 = p.plot(xlist, trunc_perc_nr_gis_list1, 'g-o',linewidth=2.0)
	line2 = p.plot(xlist, trunc_perc_nr_gis_list2, 'b-s',linewidth=2.0)
	line3 = p.plot(xlist, trunc_perc_nr_gis_list3, 'r-^',linewidth=2.0)
	line4 = p.plot(xlist, trunc_perc_nr_gis_list4, 'k-*',linewidth=2.0)
	p.legend((line1,line2,line3, line4),('E cut; K,R labeled', 'E cut; K, R, C labeled', 'P cut; K, R C labeled','P cut; K, R labeled'), loc=4)
	fig_dir = os.getcwd() + '/figures/'
	fig_name = fig_dir + 'proteome_coverage_length_all_ver2.png'
	p.savefig(fig_name,format = "png", orientation='landscape')
	
    
if __name__ == '__main__':
    
    argument = sys.argv[1:]
    
    print "Processing %s ..."%(argument)
    main(argument)
    
    import time
    print "Script - len_proteins_covered.py %s \t Completed \t %s"%(argument, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
    