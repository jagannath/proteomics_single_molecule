#!usr/bin/python

# Last modified: 10th Feb 2010
# Script for retrieving a next gene neighborhood upstream and next gene downstream for all the genes in all organisms. 

import os
import sys
import pickle







def get_pkl(path_dir,pattern):
    """
    @param path_dir: This is the directory that must be searched (os.walk) to retrieve a list of all files matching the pattern
    @param pattern: This is the pattern that must be looked for in the file
    @function: This function walks through the path_dir looking for the file with the matching pattern and returning a list
    """
    file_list = []
    for root, dirs,files in os.walk(path_dir):
	if files.find(pattern)>-1:
	    path = root + dirs + files
	    print path
	    sys.exit(1)
	    file_list.append(root+dirs+files)
	   
	    
    
    





def main(argument):
    [argv] = argument
    
    # Initializing
    base_dir = os.getcwd()
    path_dir = base_dir + '/parsed_gbk_files_genome'
    locus_cds_pkl = []
    
    pattern = 'locus_cds_information_dictionary.pkl'
    locus_cds_pkl = get_files(path_dir,pattern)
    
    





if __name__ == '__main__':
    
    argument = sys.argv[1:]
    
    print "Processing %s ..."%(argument)
    main(argument)
    
    import time
    print "Script - next_prev_gene.py %s \t Completed \t %s"%(argument, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
    


