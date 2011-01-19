#! usr/bin/python
from __future__ import division
import sys
import os
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


def main(argument):
    
    ylist = []
    
    # (1) Open file : bin_pept_lib_yeast_split.glu
    file_name = 'bin_pept_lib_yeast_split_procut_nocys.ver2.txt'
    ifile = open_file(file_name)
    lines = ifile.readlines()
    ifile.close()

    
    # (2) Obtain list of all counts
    count_list = [int(line.split('\t')[1]) for line in lines]
    for i in xrange(1,21):
	ylist.append(count_list.count(i))
	print "Number of peptides that occurrenced %d : %d"%(i,count_list.count(i))

    print ylist
    

    ## (3) Make a frequency histogram of the bin_peptide_library
    #bins = 500
    #fig = p.figure(figsize=(15,15))
    #ax1 = fig.add_subplot(1,1,1)
    #n,bins,patches = p.hist(count_list,bins,color='g')
    #p.xlabel('# Counts ')
    #p.ylabel('# Peptides')
    #p.title('Frequency histogram of number of counts of each peptide in the converted peptide library ')
    #fig_dir = os.getcwd() + '/figures/'
    #fig_name = fig_dir + 'yeast_converted_peptide_freq.png'
    #p.show()
    #p.savefig(fig_name,format = "png", orientation='landscape')
    
    return True





if __name__ == '__main__':
    
    argument = sys.argv[1:]
    
    print "Processing %s ..."%(argument)
    main(argument)
    
    import time
    print "Script - pept_lib_yeast.py %s \t Completed \t %s"%(argument, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
    