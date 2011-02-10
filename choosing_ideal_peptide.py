#! usr/bin/python

# Script is aimed to choose an ideal length peptide that has ideal properties like - (a) Length after GluC cut : approx 25 - 30; (b) Has a cysteine within the sequence and probably spaced > 20 AA from N-terminus (c) Have atleast more than 2 (Lys and Arg combined) in the first 15 bases (d) Sequencing 15 rounds should provide it with some uniqueness. 
import pickle
import os
import sys
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
	    gi = (line.split(' ')[0])[1:]
	    # (4) Making the key:value dictionary
	    sequence = sequence.upper()		#Ensuring the sequence is in uppercase
	    sequence = sequence[:-1]		# The end has a trailing *
	    gi_seq_dictionary[gi] = sequence
	    sequence = ''
	if not line.startswith('>'):	# This is to add all lines following the first line of the fasta file delimiter '>'
	    sequence += line

    return gi_seq_dictionary


def filter_sequences_ideal(gi_sequence_dict):
    """
    @param gi_sequence_dict: This is the gi: sequence dictionary 
    @function: Iterates through every sequence and (1) breaks them into peptides depending on a GluC cut (2) Chooses peptides containing cysteine (3) Filters these cysteine containing peptides for cysteine present > 20AA from N-terminus (4) From these chooses those peptides containing atleast 2 K, R (combined) (5) Eliminates those peptides that are greater than length 30. 
    """
    # Initializing
    cys_peptides = []
    ideal_peptide_list = []
    # (1) Iterate through the dictionary
    for gi, sequence in gi_sequence_dict.items():
	# (2) Make a GluC cut
	cut_site = 'E'
	peptide_list = sequence.split('E')
	# (3) Iterate through the list and choose those containing cysteine
	for peptide in peptide_list:
	    if peptide.find('C') > -1:
		# (4) Determine the position of C. Get only sequences where C position is > 20
		cys_pos = peptide.find('C')
		if cys_pos > 10 and cys_pos < 15:
		    # (5) Only one cysteine
		    if peptide.count('C') == 1:
			cys_peptides.append(peptide)
		    # (6) Choose peptides with atleast 2 K and R combined before the first C occurrence
			if (peptide[0:(peptide.find('C'))].count('K') + peptide[0:(peptide.find('C'))].count('R')) > 3:
			    # (7) Choose those peptides whose length is < 30
			    if len(peptide) > 15 and len(peptide) < 30:
				ideal_peptide_list.append(peptide)
				#print " %s : %s : \t Cysteine position = %d \t length = %d \t Nbr(K,R) = %d"%(gi, peptide, peptide.find('C'), len(peptide),(peptide.count('K')+ peptide.count('R')))

    for peptide in ideal_peptide_list:
	# (1) Convert K, R to 1 and rest to 0. Cys remains C
	bin_peptide = re.sub('[K,R]','1',peptide)
	bin_peptide = re.sub('[^1,^C]','0',bin_peptide)
	print " %s : %s : \t Cysteine position = %d \t length = %d \t Nbr(K,R) = %d \t %s"%(gi, peptide, peptide.find('C'), len(peptide),(peptide.count('K')+ peptide.count('R')),bin_peptide)
    
    return True
		
    
    





def main(argument):
    
    #[argv] = argument
    
    # (1) Make a key dictionary pair gi:binary peptide list. 
    # (1a) Open and read lines of fasta file containing the translated sequence
    file_name = 'orf_trans.fasta'
    ifile = open_file(file_name)
    lines = ifile.readlines()
    ifile.close()
    
    # (1b) Parse the file to retrieve the gi:sequence key-value pair
    gi_sequence_dictionary = get_gi_seq_dict(lines)
    print len(gi_sequence_dictionary)
    
    # (2) Filter sequences according to the various criterias 
    gi_ideal_seqs_dictionary = filter_sequences_ideal(gi_sequence_dictionary)




if __name__ == '__main__':
    
    argument = sys.argv[1:]
    
    print "Processing %s ..."%(argument)
    main(argument)
    
    import time
    print "Script - choosing_ideal_peptide.py %s \t Completed \t %s"%(argument, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
    