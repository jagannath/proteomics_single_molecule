#! usr/bin/python

# The previous script len_proteins_covered made a mistake. It considered truncating only the unique peptides which were filtered. In reality all the peptides will be truncated irrespective of the counts. I will then have to filter out the ones that still remain unique. And then estimate the coverage of proteome. 

from __future__ import division
import os
import sys
import numpy as np
import pylab as p
import pickle
import len_proteins_covered as old
import re


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
	    #peptide = re.sub('C','2',peptide)
	    #peptide = re.sub('[^1,^2]','0',peptide)
	    peptide = re.sub('[^1]','0',peptide)
	    # (5) Adding to the peptide library
	    bin_peptide_lib.append(peptide)
	    # (6) Adding to gi_pept string
	    gi_pept += peptide + '|'

	    gi_bin_pept_dict[gi] = gi_pept

    return gi_bin_pept_dict, bin_peptide_lib


def get_trunc_gi_binpept_dict(gi_binpeptide_dict,cutoff):
    """
    @param gi_binpeptide_dict: This is dictionary gi:binpeptide1|binpeptide2...
    @param cutoff: This is the length at which the peptide will be cutoff or the limit of the edman sequencing
    @function: Truncate the peptide accoring to the cutoff. Make a resulting gi:truncpeptide1|truncpeptide2...
    """
    #Initializing
    trunc_gi_binpeptide_dictionary = {}
    all_truncated_peptides_list = []

    
    
    # (1) Iterate through the dictionary of gi:binpeptide
    for gi, all_peptides in gi_binpeptide_dict.items():
	all_trunc_peptides = ''
	# (2) Split the peptide list into the individual binary peptides
	peptide_list = all_peptides.split('|')
	# (3) Iterate through the peptides
	for peptide in peptide_list:
	    # (4) Truncate the peptide to 1st cutoff
	    trunc_peptide = peptide[0:cutoff]
	    # (5) Add this truncate peptide to the others
	    all_trunc_peptides += trunc_peptide + '|'
	    all_truncated_peptides_list.append(trunc_peptide)
	# (6) Update the truncated dictionary
	trunc_gi_binpeptide_dictionary[gi] = all_trunc_peptides

    return trunc_gi_binpeptide_dictionary, all_truncated_peptides_list

def get_pept_count_dict(all_peptides,set_peptides):
    """
    @param all_peptides: This is the list of all peptides (binary and may/maynot be truncated. truncated in this case). 
    @param set_peptides: This is the non-redundant set of all peptides (set(all_peptides))
    @function: Generates the count of each of the peptide and makes a peptide:count dictionary
    """
    # Initializing
    peptide_count_dictionary = {}
    
    print "Creating peptide:count dictionary .."
    for peptide in set_peptides:
	# (1) Count the occurances in the all_peptides list
	count = all_peptides.count(peptide)
	# (2) Make the dictionary
	peptide_count_dictionary[peptide] = count
	
    return peptide_count_dictionary

def filter_peptide_count_dict(peptide_count_dictionary):
    """
    @param peptide_count_dictionary: This is the truncated_peptide:count dictionary
    @function: Iterate through all the key:value and get a list of all peptides which had a count of only one
    """
    # Initializing
    unique_trunc_pept_list = []
    
    # (1) Iterate through the peptide_count_dictionary
    for peptide, count in peptide_count_dictionary.items():
	# (2) Check if count == 1 and add peptide to the list then
	if count == 1:
	    unique_trunc_pept_list.append(peptide)
	    
    return unique_trunc_pept_list

def get_gi_list(trunc_gi_binpeptide_dictionary, unique_peptides_list,cutoff):
    """
    @param trunc_gi_binpeptide_dictionary: Dictionary of gi: truncated peptides list
    @param unique_peptides_list: A list of all peptides (truncated) that had only one count
    """
    # Initializing
    all_gi_list = []
    print "Generating unique truncated peptide : gi dictionary ..."
    unique_trunc_pept_gi_dictionary = old.get_unique_bin_peptide_gi_dict(unique_peptides_list,trunc_gi_binpeptide_dictionary)   
    #Pickle the dictionary
    pkl_file_name = 'unique_trunc_peptide_gi_dict_procut_cys_anc.ver2.'+ str(cutoff) +'.pkl'
    print "pickling file - %s"%(pkl_file_name)
    ofile = open(pkl_file_name,'wb')
    pickle.dump(unique_trunc_pept_gi_dictionary,ofile)
    ofile.close()
    # Unpickle the dictionary
    #pkl_file_name = 'unique_trunc_peptide_gi_dict_glucut_nocys.ver2.pkl'
    #ifile = open(pkl_file_name)
    #unique_trunc_pept_gi_dictionary = pickle.load(ifile)
    #ifile.close()
    
    # (2) Get list of all gis by iterating through this dictionary
    for gi in unique_trunc_pept_gi_dictionary.values():
	all_gi_list.append(gi)
    
    return all_gi_list

def filter_proteins_wcys(gi_sequence_dictionary):
    """
    @param gi_sequence_dictionary: This is the gi:sequence dictionary. 
    @function: Run through the dictionary and retrieve only sequences that contain atleast a cysteine. 
    """
    
    gi_protein_wcys_dictionary = {}
    # (1) Iterate through the dictionary
    for gi, proteins in gi_sequence_dictionary.items():
	# (2) Check if there is cys in the proteins
	if proteins.count('C') > 0:
	    gi_protein_wcys_dictionary[gi] = proteins

    return gi_protein_wcys_dictionary
    
def convert_gi_seq_dict_to_peptides(gi_seq_dict,cut_site):
    """
    @param gi_seq_dict: This is the key:value pair of protein_gi and sequence
    @param cut_site: It tells the site at which the sequence should be cut or split
    @function: Breaks the sequence according to the cut site. Relabells K and R to 1 and all others to 0 for each of the peptide fragment. Makes a gi:peptide(converted binary) dictionary; 
    #**********************
    # Note - Unlike what was used in other programs this function has an added filtering step. It totally neglects all sequences where cysteine is not present. 
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
	    # Before adding check if the peptide contains cysteine. This is critical as the peptide is anchored only via cysteine
	    if peptide.count('2') > 0:
		bin_peptide_lib.append(peptide)
		# (6) Adding to gi_pept string
		gi_pept += peptide + '|'

		gi_bin_pept_dict[gi] = gi_pept

    return gi_bin_pept_dict, bin_peptide_lib
    
def get_trunc_gi_binpept_dict_wcys_anc(gi_bin_pept_wcys_anc_dict, cutoff):
    """
    @param gi_protein_wcys_dictionary: This is the gi:sequence (containing cysteine)
    @function: Then based on the cutoff truncate the peptide list to contain the first cutoff or first cysteine. Generates a list and dictionary of bin_peptide:count 
    """
    # Initializing
    all_trunc_bin_pept_wcys_cutoff_list = []
    trunc_gi_binpeptide_cys_anc_dictionary = {}
    
    # (1) Iterate through the gi_peptide
    for gi, all_peptides in gi_bin_pept_wcys_anc_dict.items():
	all_trunc_peptides = ''
	# (2) Make a list of all peptides
	peptide_list = all_peptides.split('|')
	# (3) Iterate through each peptide in the list and cut it (a) at first cys and then (b) at the cutoff. It really wouldnt matter as if the length is short it [x:y] will still return back all the characters present
	for peptide in peptide_list:
	    # (3a) Find first cysteine position and truncate
	    cpos = peptide.find('2')
	    peptide = peptide[0:cpos]
	    # (3b) Truncate according to the cutoff argument passed
	    peptide = peptide[0:cutoff]
	    # (4) Add the peptide to the others in the list
	    all_trunc_peptides += peptide + '|'
	    # (5) Update the truncated_peptide_with cys anchored list
	    all_trunc_bin_pept_wcys_cutoff_list.append(peptide)
	    
	# (5) Update the gi_binpeptide_cys_anc_dictionary
	trunc_gi_binpeptide_cys_anc_dictionary[gi] = all_trunc_peptides
    
    return trunc_gi_binpeptide_cys_anc_dictionary, all_trunc_bin_pept_wcys_cutoff_list
        
    

def main(argument):
    
    [argv] = argument
    assert len(argument) == 1

    if argv == '1':
	    
	# (1) Make a gi_sequence dictionary
	# (1a) Open and read lines of fasta file containing the translated sequence
	file_name = 'orf_trans.fasta'
	ifile = old.open_file(file_name)
	lines = ifile.readlines()
	ifile.close()
	
	# (1b) Parse the file to retrieve the gi:sequence key-value pair
	gi_sequence_dictionary = old.get_gi_seq_dict(lines)
	print len(gi_sequence_dictionary)
	
	# (2) Unpickle dictionary - # These are gi - cut by E or P and labelled at K,R and/not C. Already made this in the previous programme. So will simply use the picke file. All starts with gi_binpeptide_dict_*.ver2.pkl
	# I am going to focus only on glucut_nocys
	pkl_file_name = 'gi_binpeptide_dict_glucut_nocys.ver2.pkl'
	ifile = open(pkl_file_name)
	gi_binpeptide_dict = pickle.load(ifile)
	ifile.close()
	
	# (3) Truncate all peptides according to a cutoff. And obtain a new trunc_gi_binpeptide_dictionary with gi:truncated peptides| and a list of all truncated peptides
	for cutoff in xrange(5,31,5):
	    print "Truncating peptide library for %d"%(cutoff)
	    trunc_gi_binpeptide_dictionary, all_truncated_peptides_list = get_trunc_gi_binpept_dict(gi_binpeptide_dict, cutoff)	
	    # Pickle the trunc_gi_binpeptide_dictionary
	    pkl_file_name = 'trunc_gi_binpeptide_dict_glucut_nocys.ver2.' + str(cutoff) + '.pkl'
	    ofile = open(pkl_file_name,'wb')
	    pickle.dump(trunc_gi_binpeptide_dictionary,ofile)
	    ofile.close()
	    print len(all_truncated_peptides_list)
	    set_trunc_peptides_list = set(all_truncated_peptides_list)
	    print len(set_trunc_peptides_list)
	    
	    # (4) Make a dictionary of peptide: count
	    peptide_count_dictionary = get_pept_count_dict(all_truncated_peptides_list,set_trunc_peptides_list)
	    # Pickling the dictionary
	    pkl_file_name = 'trunc_peptide_count_dic_glucut_nocys.ver2.'+ str(cutoff) + '.pkl'
	    ofile = open(pkl_file_name,'wb')
	    pickle.dump(peptide_count_dictionary,ofile)
	    ofile.close()
	
    if argv == '2':
	# Initializing
	proteome_coverage_list = []
	# This part generates the proteome coverage for different length cutoffs
	for cutoff in xrange(5,31,5):
	    print "Processing for %d "%(cutoff)
	    # (1) Unpickle the trunc_peptide:count dictionary
	    pkl_file_name = 'trunc_peptide_count_dic_glucut_nocys.ver2.' + str(cutoff) + '.pkl'
	    ifile = open(pkl_file_name)
	    peptide_count_dictionary = pickle.load(ifile)
	    ifile.close()
	    
	    # (2) Filter out from this dictionary a list of unique peptides
	    unique_peptides_list = filter_peptide_count_dict(peptide_count_dictionary)
 
	    # (3) For each of the unique peptide check in which gi it is present
	    # (3a) Unpickle the trunc_gi_binpeptide_dictionary
	    pkl_file_name = 'trunc_gi_binpeptide_dict_glucut_nocys.ver2.'+ str(cutoff) + '.pkl'
	    ifile = open(pkl_file_name)
	    trunc_gi_binpeptide_dictionary = pickle.load(ifile)
	    ifile.close()
	    
	    # (3b) Get list of all gis for which the unique peptides are a member of
	    all_gi_list = get_gi_list(trunc_gi_binpeptide_dictionary, unique_peptides_list,cutoff)
	    print len(all_gi_list)
	    nr_gi_list = set(all_gi_list)
	    nbr_nr_gis = len(nr_gi_list)
	    proteome_covered = (nbr_nr_gis/5885)*100
	    proteome_coverage_list.append(proteome_covered)

	# (4) Draw graph
	xlist = [i for i in xrange(5,51,5)]
	fig = p.figure(figsize=(15,15))
	ax1 = fig.add_subplot(1,1,1)
	p.xlabel('# Length of peptide ')
	p.ylabel('# Percentage of proteome covered')
	p.title('Proteome coverage as a function of sequenced length ')
	line1 = p.plot(xlist, proteome_coverage_list, 'g-o',linewidth=2.0)
	#line2 = p.plot(xlist, trunc_perc_nr_gis_list2, 'b-s',linewidth=2.0)
	#line3 = p.plot(xlist, trunc_perc_nr_gis_list3, 'r-^',linewidth=2.0)
	#line4 = p.plot(xlist, trunc_perc_nr_gis_list4, 'k-*',linewidth=2.0)
	#p.legend((line1,line2,line3, line4),('E cut; K,R labeled', 'E cut; K, R, C labeled', 'P cut; K, R C labeled','P cut; K, R labeled'), loc=4)
	fig_dir = os.getcwd() + '/figures/'
	fig_name = fig_dir + 'proteome_coverage_length_glucut_nocys_allpeptides_considered_ver2.png'
	p.savefig(fig_name,format = "png", orientation='landscape')
	#p.show()

    if argv == '3':
	# This is the module where the Cys containing proteins are first filtered. Then the truncation of the peptide is so done such that it is cut off at the first cysteine or the cutoff that is passed. All the rest is identical. 
	
	# (1) Obtain the key dictionary pair of gi_protein_dictionary
	# (1a) Open and read lines of fasta file containing the translated sequence
	file_name = 'orf_trans.fasta'
	ifile = old.open_file(file_name)
	lines = ifile.readlines()
	ifile.close()
	
	# (1b) Parse the file to retrieve the gi:sequence key-value pair
	gi_sequence_dictionary = old.get_gi_seq_dict(lines)
	print len(gi_sequence_dictionary)
	
	# (2) Filter the dictionary to contain only cys. Get a count of proteins with cysteine
	gi_protein_wcys_dictionary = filter_proteins_wcys(gi_sequence_dictionary)
	print "Number of proteins containing Cysteine : %d or proteome covered %f" %(len(gi_protein_wcys_dictionary), (len(gi_protein_wcys_dictionary)/5885)*100)

	# (3) Convert the sequence into binary peptides according to glu cut and cys labelled
	# (3a) Generate the list of all the binary peptide library (with cysteine included) and the gi:binary peptide dictionary
	cut_site = 'P'	# Glu cut
	gi_bin_pept_wcys_anc_dict, bin_pept_wcys_anc_list = convert_gi_seq_dict_to_peptides(gi_protein_wcys_dictionary,cut_site)
	print "Size of the peptide library with cysteine %d"%(len(bin_pept_wcys_anc_list))
	
	# (4) Truncate all peptides according to a cutoff. And obtain a new trunc_gi_binpeptide_dictionary with gi:truncated peptides| and a list of all truncated peptides
	# ************* Note : This function is modified. It will cut also at first C (or 2) it sees apart from the defined cutoff argument passed. 
	for cutoff in xrange(5,51,5):
	    print "Truncating peptide (with cys anchored) library for %d"%(cutoff)
	    trunc_gi_binpeptide_cys_anc_dictionary, all_truncated_peptides_list = get_trunc_gi_binpept_dict_wcys_anc(gi_bin_pept_wcys_anc_dict, cutoff)	
	    # Pickle the trunc_gi_binpeptide_dictionary
	    pkl_file_name = 'trunc_gi_binpeptide_dict_procut_cys_anc.ver2.' + str(cutoff) + '.pkl'
	    ofile = open(pkl_file_name,'wb')
	    pickle.dump(trunc_gi_binpeptide_cys_anc_dictionary,ofile)
	    ofile.close()
	    print "Size of all truncated peptides (with cysteine anchored) : %d"%(len(all_truncated_peptides_list))
	    set_trunc_peptides_list = set(all_truncated_peptides_list)
	    print "Nonredundant set of all truncated peptides (with cysteine anchored) : %d"%(len(set_trunc_peptides_list))
	    
	    # (5) Make a dictionary of peptide: count
	    peptide_count_dictionary = get_pept_count_dict(all_truncated_peptides_list,set_trunc_peptides_list)
	    # Pickling the dictionary
	    pkl_file_name = 'trunc_peptide_count_dic_procut_cys_anc.ver2.'+ str(cutoff) + '.pkl'
	    ofile = open(pkl_file_name,'wb')
	    pickle.dump(peptide_count_dictionary,ofile)
	    ofile.close()
    
    if argv == '4':
	# This computes the proteome coverages for Human proteome or E.coli proteome
	



    if argv == '5':
	# This will now compare the unique peptides in the truncated peptided library (reminder = truncated at first C or first cutoff) and map it to the gi it is present in. A dictionary of trunc_gi_binpeptide_cys_anc_dictionary is already created
	# Initializing
	proteome_coverage_list = []
	
	for cutoff in xrange(5,51,5):

	    # (1) Unpickle the dictionary of trunc peptide:count to obtain a list of only unique peptides
	    pkl_file_name = 'trunc_peptide_count_dic_procut_cys_anc.ver2.' + str(cutoff) + '.pkl'
	    ifile = open(pkl_file_name)
	    peptide_count_dictionary = pickle.load(ifile)
	    ifile.close()
	    
	    # (2) Filter out from this dictionary a list of unique peptides
	    unique_peptides_list = filter_peptide_count_dict(peptide_count_dictionary)
	    print "Number of unique peptides in the library truncated by (a) cutoff %d and (b) first Cysteine observed : %d"%(cutoff,len(unique_peptides_list))
	    
	    # (3) For each of the unique peptide check in which gi it is present
	    # (3a) Unpickle the trunc_gi_binpeptide_dictionary
	    pkl_file_name = 'trunc_gi_binpeptide_dict_procut_cys_anc.ver2.'+ str(cutoff) + '.pkl'
	    ifile = open(pkl_file_name)
	    trunc_gi_binpeptide_cys_anc_dictionary = pickle.load(ifile)
	    ifile.close()

	    # (3b) Get list of all gis for which the unique peptides are a member of
	    all_gi_list = get_gi_list(trunc_gi_binpeptide_cys_anc_dictionary, unique_peptides_list,cutoff)
	    nr_gi_list = set(all_gi_list)
	    nbr_nr_gis = len(nr_gi_list)
	    print "Number of Proteins uniquely called : %d"%(nbr_nr_gis)
	    proteome_covered = (nbr_nr_gis/5885)*100
	    print "Proteome covered with cutoff %d : %f "%(cutoff,proteome_covered)
	    proteome_coverage_list.append(proteome_covered)
	    
	# (4) Draw graph
	xlist = [i for i in xrange(5,51,5)]
	fig = p.figure(figsize=(15,15))
	ax1 = fig.add_subplot(1,1,1)
	p.xlabel('# Length of peptide ')
	p.ylabel('# Percentage of proteome covered')
	p.title('Proteome coverage as a function of sequenced length : Procut - Cys anchored ')
	#line1 = p.plot(xlist, proteome_coverage_list, 'g-o',linewidth=2.0)
	line2 = p.plot(xlist, proteome_coverage_list, 'b-s',linewidth=2.0)
	#line3 = p.plot(xlist, trunc_perc_nr_gis_list3, 'r-^',linewidth=2.0)
	#line4 = p.plot(xlist, trunc_perc_nr_gis_list4, 'k-*',linewidth=2.0)
	#p.legend((line1,line2,line3, line4),('E cut; K,R labeled', 'E cut; K, R, C labeled', 'P cut; K, R C labeled','P cut; K, R labeled'), loc=4)
	fig_dir = os.getcwd() + '/figures/'
	fig_name = fig_dir + 'proteome_coverage_length_procut_cys_anc_ver2.png'
	p.savefig(fig_name,format = "png", orientation='landscape')
	#p.show()

    
    


if __name__ == '__main__':
    
    argument = sys.argv[1:]
    
    print "Processing %s ..."%(argument)
    main(argument)
    
    import time
    print "Script - proteome_coverage_all_peptides.py %s \t Completed \t %s"%(argument, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
    