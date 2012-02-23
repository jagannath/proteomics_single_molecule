#! usr/bin/python 

# This is an improvised version of the earlier files for cutting and analysing the proteome coverages. Will be importing a lot of functions taken from proteome _coverage_all_peptides
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

def get_gi_seq_dict_human(lines):
    """
    @param lines: The lines in the non redundant proteome data base (for yeast this time). 
    @function: Generates a key:value dictionary pair; gi id is the key and the value is the sequence
    """
    #Initializing
    flag = False
    all_sequences = []
    sequence = ''
    gi_seq_dictionary = {}
    inseq = ''
    i =0
    # (1) Iterate lines and get gi
    for line in lines:
	i+=1
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
    
def get_gi_seq_dict_yeast(lines):
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
	    gi = line.split(' ')[0][1:]
	    # (4) Making the key:value dictionary
	    sequence = sequence.upper()		#Ensuring the sequence is in uppercase
	    sequence = sequence[:-1]		# The end has a trailing *
	    gi_seq_dictionary[gi] = sequence
	    sequence = ''
	if not line.startswith('>'):	# This is to add all lines following the first line of the fasta file delimiter '>'
	    sequence += line

    return gi_seq_dictionary    

def get_gi_seq_dict_ecoli(lines):
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
	    gi = line.split('|')[1][:-1]
	if not line.startswith('>'):	# This is to add all lines following the first line of the fasta file delimiter '>'
	    sequence += line
	    # (4) Making the key:value dictionary
	    sequence = sequence.upper()		#Ensuring the sequence is in uppercase
	    sequence = sequence		# The end has a trailing *
	    gi_seq_dictionary[gi] = sequence
	    sequence = ''

    return gi_seq_dictionary    


def convert_gi_seq_dict_to_peptides(gi_seq_dict,cut_site,fluor_limit):
    """
    @param gi_seq_dict: This is the key:value pair of protein_gi and sequence
    @param cut_site: It tells the site at which the sequence should be cut or split
    @param fluor_limit: This is the number of fluorescent K's you can ideally have to be in the detectable range
    @function: Breaks the sequence according to the cut site. Relabells K and R to 1 and all others to 0 for each of the peptide fragment. 
    """
    # Initializing
    bin_peptide_lib = []
    gi_bin_pept_dict = {}
    k_count_list = []
    
    # (1) Iterate through each key:value pair
    for gi, sequence in gi_seq_dict.items():
	# Initializing
	gi_pept = []
	# (2) Break sequence according to cut_site
	assert len(cut_site) == 1	# Ensure that only one amino acid letter is provided
	peptides_list = sequence.split(cut_site)
	# (4) For each member in the peptide_list, convert K and R into 1. Rest all to 0
	for peptide in peptides_list:
	    peptide = re.sub('[K]','1',peptide)
	    peptide = re.sub('C','2',peptide)
	    #peptide = re.sub('[E,D]','3',peptide)	#Note this is only for Threonine labelling
	    #peptide = re.sub('Y','4',peptide)
	    peptide = re.sub('[^1,^2,^3,^4]','0',peptide)
	    # (5) Adding to the peptide library only those peptides containing cysteine; Truncate it to the first Cys position
	    # Note = this is  a change from earlier use of the function
	    if peptide.find('2') > -1:
		k_count_list.append(peptide.count('1'))
		if peptide.count('1') < fluor_limit:	#This is for ensuring you dont have too much fluorescence
		    cys_pos = peptide.find('2')
		    peptide = peptide[0:cys_pos]
		    bin_peptide_lib.append(peptide)
		    # (6) Adding to gi_pept string
		    gi_pept.append(peptide)
		    gi_bin_pept_dict[gi] = gi_pept

    return gi_bin_pept_dict, bin_peptide_lib, k_count_list

def pickle_file(file_name, contents):
    """
    @param file_name: This is the file name to which the contents must be dumped
    @param contents: This is the contents to be pickled
    @function: The contents passed are pickled to the file in the pkl_directory
    """
    pkl_dir = os.getcwd() + '/pkl_files/'
    pkl_file = pkl_dir + file_name
    ofile = open(pkl_file,'wb')
    pickle.dump(contents, ofile)
    ofile.close()
    
    return True
    
def unpickle_file(file_name):
    """
    @param file_name: This is the file name to which the contents must be dumped
    @param contents: This is the contents to be pickled
    @function: The contents passed are pickled to the file in the pkl_directory
    """
    pkl_dir = os.getcwd() + '/pkl_files/'
    pkl_file = pkl_dir + file_name
    ifile = open(pkl_file)
    contents = pickle.load(ifile)
    ifile.close()
    
    return contents
    
def make_bin_pept_gis_dict(gi_bin_pept_dict, bin_pept_lib, org, mod):
    """
    @param gi_bin_pept_dict: This is gi: binary peptide dictionary
    @param bin_pept_lib: This is a list of all the binary peptides in the library (it is a redundant set)
    @function: For each of the binary peptide make a gi list
    """
    pept_gilist_dict = {}
    shelve_file = os.getcwd() + '/shelve_files/' + org + 'prot_bin_pept'+ mod + ':gi_list'
    s = shelve.open(shelve_file)
    i = 0
    
    # (1) Iterate through the bin peptide set. For each run through the gis_pept dict and ask if the bin peptide was in the gis. 
    for pept in set(bin_pept_lib):
	i+=1
	print "Processing %d : %s "%(i,str(pept))
	if (pept):	# There is a nonsensical empty peptide
	    gi_list = []
	    # (2) Iterate through the gi_bin_pept_dict and ask if this pept was in the list
	    for gi, pept_list in gi_bin_pept_dict.items():
		if pept in pept_list:
		    gi_list.append(gi)
	    
	    # (3) Make the key: value dictionary pair and shelve it
	    pept_gilist_dict[pept] = gi_list
	    s[pept] = gi_list
    
    s.close()
    
    return pept_gilist_dict

def cut_all_peptides(pept_gilist_dict, cutoff):
    """
    @param pept_gilist_dict: This is peptide: gi_list dictionary;
    @function: Iterate through the list. Make another non_unique_list; Truncate the peptide according to the cutoff and put the extra ones in the non unique list; Also generate a total peptide list (truncated). From these two list - unique truncated peptides will be obtained that can be mapped to the gi list to obtain list of all gis
    """
    all_cut_peptides = []
    non_uniq_peptides = []
    trunc_pept_gilist_dict = {}
    
    for pept, gi_list in pept_gilist_dict.items():
	trunc_pept = pept[0:cutoff]
	if not trunc_pept in all_cut_peptides:
	    all_cut_peptides.append(trunc_pept)
	    trunc_pept_gilist_dict[trunc_pept] = gi_list
	else:
	    non_uniq_peptides.append(trunc_pept)
    
    return all_cut_peptides, non_uniq_peptides, trunc_pept_gilist_dict

def get_all_gi_list(uniq_cut_peptides, trunc_pept_gilist_dict):
    """
    @param uniq_cut_peptides: This is a unique list of cut peptides
    @param trunc_pept_gilist_dict: This is truncated peptide: gilist dictionary; 
    @function: Obtain a list of all gis in all these unique cut peptides list
    """
    all_gis = []
    for pept in uniq_cut_peptides:
	gi_list = trunc_pept_gilist_dict[pept]
	if len(gi_list) == 1:
	    for item in gi_list:
		all_gis.append(item)
		
    return all_gis

def main(argument):
    
    [argv] = argument
    org = 'human'
    mod = '_K_cutE'
    ver = '6'
    cut_site = 'E'
    fluor_limit = 100
    req_hist = True
    #pdb.set_trace()
    if argv == '0':
	# Just makes gi: sequence dictionary and gi: bin_peptide dictionary; The bin_peptide_lib is generated by cutting with GluC and then filtering out peptides with Cys and truncating peptides to the first cysteine
	# (1) Open the database - human.protein.faa in Downloads/human_proteome directory. Contains a lot of genes -34287; So this is a large dataset to play with
	#
	print org, mod
	if org == 'ecoli':
	    file_name = 'NC_000913.fasta'
	if org == 'yeast':
	    file_name = '/project/marcotte/jagannath/projectfiles/Downloads/yeast_proteins_sgd/orf_trans.fasta'
	if org == 'human':
	    file_name = '/project/marcotte/jagannath/projectfiles/Downloads/human_proteome/uniprot_sprot_human.dat'
	    
	ifile = open(file_name,'r')
	lines = ifile.readlines()
	ifile.close()
	# (1b) Obtain the dictionary of gi: sequence; 
	if org == 'ecoli':
	    gi_sequence_dictionary = get_gi_seq_dict_ecoli(lines)
	if org == 'yeast':
	    gi_sequence_dictionary = get_gi_seq_dict_yeast(lines)
	if org == 'human':
	    gi_sequence_dictionary = get_gi_seq_dict_human(lines)
	    
	
	print len(gi_sequence_dictionary)
	# (1c) Pickle this
	file_name = org+'prot_gi:sequence' + mod + '.pkl'
	pickle_file(file_name, gi_sequence_dictionary)

	# (2) Cut the proteome with a defined cut site - like E for GluC. Generate a list of peptides - Peptides need to contain cysteine and truncated to the first cys position
	#cut_site = 'E'
	gi_bin_pept_dict, bin_pept_lib,k_count_list = convert_gi_seq_dict_to_peptides(gi_sequence_dictionary,cut_site,fluor_limit)
	nbr_peptides = len(k_count_list)
	k_counts = [((k_count_list.count(item)/nbr_peptides)*100) for item in xrange(0,25)]
	# Making a histogram if needed and terminating the programme
	if req_hist:
	    # Make a frequency histogram of the #Peptides with Different number of fluorescent K's  vs Counts
	    fig = p.figure(figsize=(15,15))
	    ax1 = fig.add_subplot(1,1,1)
	    xaxis = np.arange(25)
	    width = 0.35
	    rect = ax1.bar(xaxis, k_counts,width, color = 'r')
	    p.xlabel('# Number of Lysines in peptide ')
	    p.ylabel('# Frequency of occurrence (in Perc of total peptides in library)')
	    title_name = org + 'Frequency histogram of Number of Lysines in peptides observed (with Cys anchored)' + mod
	    p.title(title_name)
	    fig_dir = os.getcwd() + '/figures/'
	    fig_name = fig_dir + org + 'freq_nbr_labeled_K' + mod + '.png'
	    #p.show()
	    p.savefig(fig_name,format = "png", orientation='landscape')
	    sys.exit(1)
	
	# (2b) Pickle the dictionary and the binary peptide library
	file_name = org + 'prot_gi:bin_peptide' + mod + '.pkl'
	pickle_file(file_name, gi_bin_pept_dict)
	file_name = org + 'prot_bin_pept_lib' + mod + '.pkl'
	pickle_file(file_name,bin_pept_lib)
	print len(bin_pept_lib)
	print len(set(bin_pept_lib))
	# (2c) Creating a unique peptide list
	unique_peptides_in_lib = [item for item in set(bin_pept_lib) if bin_pept_lib.count(item) == 1]
	print len(unique_peptides_in_lib)
	file_name = org + 'prot_uniq_pept_lib' + mod + '.pkl'
	pickle_file(file_name, unique_peptides_in_lib)
    
    if argv == '1':
    
	print org, mod
	# (1) Unpickle gi:bin_pept dictionary and the binary peptide library
	file_name = org + 'prot_gi:bin_peptide' + mod + '.pkl'
	gi_bin_pept_dict = unpickle_file(file_name)
	file_name = org + 'prot_bin_pept_lib' + mod + '.pkl'
	bin_pept_lib = unpickle_file(file_name)

	# (2) Make a dictionary of bin_pept: [gis] for each of the ones in the library
	pept_gilist_dict = make_bin_pept_gis_dict(gi_bin_pept_dict, bin_pept_lib,org, mod)
	
	# (3) Pickling in case
	file_name = org + 'prot_bin_pept:gi_list' + mod + '.pkl'
	pickle_file(file_name, pept_gilist_dict)
	
    if argv == '2':
	# (3) cutoff = 5
	# Iterate through all the peptides in the pept_gilist_dict
	# (1) Making cuts. So each cut get a dictionary of cut_peptide (according to cutoff): [gi_list]; This arises simply by truncating the peptides in the peptide: gi list; As you run through, will make a redundant set - i.e. peptides that were observed before as well. So that among the list of all the truncated peptides, unique peptides may be identifiable. Then get a list of all gis and then get number of unique gis
	coverage_list = []
	print org, mod
	# (1) Unpickle the bin_pept:gi_list dictionary
	file_name = org + 'prot_bin_pept:gi_list' + mod + '.pkl'
	pept_gilist_dict = unpickle_file(file_name)
	
	for cutoff in xrange(0,51,5):
	    # (2) Get a list of (a) all cut peptides (b) list of non unique peptides (c) cut_peptide: gi_list
	    all_cut_peptides, non_uniq_peptides, trunc_pept_gilist_dict = cut_all_peptides(pept_gilist_dict, cutoff)
	    print len(all_cut_peptides)
	    print len(non_uniq_peptides)
	    # (3) Generate a set of unique peptides from the two  list
	    uniq_cut_peptides = [item for item in all_cut_peptides if non_uniq_peptides.count(item) == 0]
	    
	    print len(uniq_cut_peptides)
	    # (4) For these unique peptides get a list of all uniquely occurring gis
	    all_gis = get_all_gi_list(uniq_cut_peptides, trunc_pept_gilist_dict)
	    print len(all_gis)
	    print len(set(all_gis))
	    ## (5) Get a list of unique gis that corresponds to the proteome covered
	    # ********This a bug **********#
	    #uniq_gis = [item for item in set(all_gis) if all_gis.count(item) == 1]	
	    if org == 'ecoli':nbr_proteins = 4317
	    if org == 'yeast':nbr_proteins = 5885
	    if org == 'human':nbr_proteins = 20252
	    prot_covered = (len(set(all_gis))/nbr_proteins) * 100
	    
	    coverage_list.append(prot_covered)
	    print "Proteome coverage for %d cycles : %f"%(cutoff, prot_covered)

	print "Coverage list "
	print coverage_list
	# (6) Draw graph
	xlist = [i for i in xrange(0,51,5)]
	fig = p.figure(figsize=(15,15))
	ax1 = fig.add_subplot(1,1,1)
	p.xlabel('# Number of Edman Cycles ')
	p.ylabel('# Percentage of proteome covered')
	title_name = org + ' proteome coverage as a function of sequenced length : Glucut - Cys anchored and labelled at' + mod
	p.title(title_name)
	line = p.plot(xlist, coverage_list, 'b-s',linewidth=2.0)
	fig_dir = os.getcwd() + '/figures/'
	fig_name = fig_dir + org + 'proteome_coverage_' + mod + ver + '.png'
	#p.savefig(fig_name,format = "png", orientation='landscape')
	p.show()

    if argv == '3':
	
	# When constraining the number of K's that can be observed. i.e if #Ks > 10 it is hard to distinguish. But if #Ks <5 (limit) it is easier to detect
	# No Limit 
	# If we labeled only K ; cut E
	cov1 = [0.0, 0.0049377839225755484, 1.0517479755085917, 5.3130555006912896, 12.137072881690697, 18.427809599051944, 23.256962275330832, 26.817104483507805, 29.281058660873001, 31.053723089077621, 32.115346632431361]
	
	# Limit = 2
	cov2 = [0.0, 0.0, 1.1011258147343472, 5.8364605964842973, 13.070314043057477, 21.346039897294094, 27.666403318190795, 31.80920402923168, 34.712620975706102, 36.79636579103299, 38.208571992889588]
	# Limit = 4
	cov3 =[0.0, 0.22220027651589969, 8.4979261307525178, 27.335571795378232, 43.279676081374681, 51.945486865494757, 56.809204029231687, 59.139838040687344, 60.305155046415173, 60.724866679834086, 60.976693659885441]
	# Limit = 6
	cov4 = [0.0, 0.20244914082559748, 8.5917440252814536, 28.979853841595894, 46.336164329448941, 55.421686746987952, 60.532293106853643, 62.917242741457635, 63.840608334979265, 64.299822239778791, 64.507209164526955]
	# Limit = 8
	cov5 = [0.0, 0.05925340707090658, 16.146553426822042, 45.546118901836856, 60.581670946079399, 66.067548884060841, 68.299427217064974, 69.326486272960693, 69.825202449140818, 70.106656132727636, 70.294291921785501]
	# Limit = 10
	cov6 = [0.0, 0.05097706032285472, 2.4978759558198811, 10.620220900594731, 21.546304163126592, 30.416312659303312, 36.805437553101108, 40.832625318606624, 42.990654205607477, 43.959218351741718, 44.706881903143589]
	# Limit = 12
	cov7 = [0.0, 0.05097706032285472, 2.7017841971113001, 10.994052676295667, 21.988105352591333, 30.926083262531861, 37.349192863211556, 41.308411214953274, 43.483432455395068, 44.468988954970264, 45.250637213254038]
	
	# (6) Draw graph
	xlist = [i for i in xrange(0,51,5)]
	fig = p.figure(figsize=(10,10))
	ax1 = fig.add_subplot(1,1,1)
	p.xlabel('# Number of Edman Cycles ')
	p.ylabel('# Percentage of proteome covered')
	p.title('Human proteome coverage as a function of sequenced length : Different Labels; Cut - E or R')
	line1 = p.plot(xlist, cov1, 'b-s',linewidth=2.0)
	line2 = p.plot(xlist, cov2, 'g-o',linewidth=2.0)
	line3 = p.plot(xlist, cov3, 'r-^',linewidth=2.0)
	line4 = p.plot(xlist, cov4, 'c-*',linewidth=2.0)
	line5 = p.plot(xlist, cov5, 'k-+',linewidth=2.0)
	#line6 = p.plot(xlist, cov6, 'm-D',linewidth=2.0)
	#line7 = p.plot(xlist, cov7, 'y-H',linewidth=2.0)
	#line5 = p.plot(xlist, cov5, 'm-p',linewidth=2.0)
	p.legend((line1,line2,line3, line4, line5),('Label-K;Cut-E','Label-K;Cut-R','Label-K,T;Cut-E','Label-K,T;Cut-R','Label-K,[E,D];Cut-R'), loc=4)
	fig_dir = os.getcwd() + '/figures/'
	fig_name = fig_dir + org + 'proteome_coverage_diffLabels.eps'
	p.savefig(fig_name,format = "eps", orientation='landscape')
	p.show()
    
    if argv == '4':
	uniq_bin_pept = '0101100011'
	# (1) Unpickle the bin_pept:gi_list dictionary
	file_name = org + 'prot_bin_pept:gi_list' + mod + '.pkl'
	bin_pept_gilist_dict = unpickle_file(file_name)
	# (2) Unpickle gi:pept list
	file_name = org+'prot_gi:sequence' + mod + '.pkl'
	gi_seq_dict = unpickle_file(file_name)
	
	cutoff_long_pept = [(pept,gi) for pept, gi in bin_pept_gilist_dict.items() if len(pept) == 10 and len(gi) == 1]
	print cutoff_long_pept
	
	file_name = org + 'prot_gi:bin_peptide' + mod + '.pkl'
	gi_bin_pept_dict = unpickle_file(file_name)
	
	seq = gi_seq_dict['YFR031C-A']
	split_seq = [item for item in seq.split(cut_site) if not item.count('C') == 0]
	print split_seq
	print gi_bin_pept_dict['YFR031C-A']
	
    if argv == '5':
	# This just computes the frequency of each of the amino acid 
	aa_list = ['R','H','K', 'D','E','S','T','N','Q', 'C', 'G','P', 'A','I','L','M', 'F','W','Y','V']
	all_sequence = ''
	amino_acid_freq_dictionary = {}
	all_freq = 0
	
	print org, mod
	if org == 'ecoli':
	    file_name = 'NC_000913.fasta'
	if org == 'yeast':
	    file_name = '/project/marcotte/jagannath/projectfiles/Downloads/yeast_proteins_sgd/orf_trans.fasta'
	if org == 'human':
	    file_name = '/project/marcotte/jagannath/projectfiles/Downloads/human_proteome/uniprot_sprot_human.dat'
	    
	ifile = open(file_name,'r')
	lines = ifile.readlines()
	ifile.close()
	# (1b) Obtain the dictionary of gi: sequence; 
	if org == 'ecoli':
	    gi_sequence_dictionary = get_gi_seq_dict_ecoli(lines)
	if org == 'yeast':
	    gi_sequence_dictionary = get_gi_seq_dict_yeast(lines)
	if org == 'human':
	    gi_sequence_dictionary = get_gi_seq_dict_human(lines)
	    
	# Concatenate all the Amino acid sequences to get a gigantic sequence
	for sequence in gi_sequence_dictionary.values():
	    all_sequence += sequence
	print all_sequence[0:500]
	# Iterate through each of the amino acid in the list and get its counts. Make a dictionary pair - Aminoacid: frequency (counts/length of the all_sequence)
	for aa in aa_list:
	    count = all_sequence.count(aa)
	    freq = (count/len(all_sequence))*100
	    all_freq += freq
	    amino_acid_freq_dictionary[aa] = freq
	    
	print len(all_sequence)
	print all_freq
	print amino_acid_freq_dictionary
    

if __name__ == '__main__':
    
    argument = sys.argv[1:]
    
    print "Processing %s ..."%(argument)
    main(argument)
    
    import time
    print "Script - proteome_coverage_single_mol_edman.py %s \t Completed \t %s"%(argument, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
    