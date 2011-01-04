#! usr/bin/python

from __future__ import division
import numpy as np
import pylab as p
import sys
import os



def main(argument):
    
    xlist = np.arange(1,21)
    # E cut (K, R labeled)
    list1 = [26638, 11904, 4201, 1590, 850, 510, 391, 326, 240, 177, 146, 116, 115, 79, 67, 67, 56, 49, 24, 37]
    norm_list1 = [(item*100/48520) for item in list1]
    # E cut (K, R, C labeled)
    list2 = [31768, 13952, 4851, 1790, 883, 531, 425, 325, 220, 169, 158, 132, 96,  87, 73, 59, 46, 40, 47, 41]
    norm_list2 = [(item*100/48520) for item in list2]
    # P cut (K, R, C labeled)
    list3 = [32876, 14171, 4688, 1581, 735, 437, 354, 272, 167, 136, 112, 106, 63,  59, 43, 51, 46, 41, 26, 34]
    norm_list3 = [(item*100/48520) for item in list3]
    print norm_list1
    print norm_list2
    print norm_list3
    
    fig = p.figure(figsize=(15,15))
    ax1 = fig.add_subplot(1,1,1)
    p.xlabel('# Counts ')
    p.ylabel('# Peptides')
    p.title('Frequency histogram of number of counts of each peptide in the converted peptide library ')
    line1 = p.plot(xlist, norm_list1, linewidth=2.0)
    line2 = p.plot(xlist, norm_list2, linewidth=2.0)
    line3 = p.plot(xlist, norm_list3, linewidth=2.0)
    p.legend((line1,line2,line3),('E cut; K,R labeled', 'E cut; K, R, C labeled', 'P cut; K, R, C labeled'), loc=4)
    fig_dir = os.getcwd() + '/figures/'
    fig_name = fig_dir + 'yeast_converted_peptide_freq.png'
    p.savefig(fig_name,format = "png", orientation='landscape')
    
    
	
    
    





if __name__ == '__main__':
    
    argument = sys.argv[1:]
    
    print "Processing %s ..."%(argument)
    main(argument)
    
    import time
    print "Script - pept_lib_yeast.py %s \t Completed \t %s"%(argument, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
    