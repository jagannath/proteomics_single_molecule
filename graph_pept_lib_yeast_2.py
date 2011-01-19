#! usr/bin/python

from __future__ import division
import numpy as np
import pylab as p
import sys
import os



def main(argument):
    
    xlist = np.arange(1,21)
    # E cut (K, R labeled) i.e glucut_nocys.ver2
    #list1 = [26638, 11904, 4201, 1590, 850, 510, 391, 326, 240, 177, 146, 116, 115, 79, 67, 67, 56, 49, 24, 37]
    list1 = [37930, 3223, 1082, 609, 358, 244, 184, 135, 106, 96, 73, 65, 63, 37, 39, 48, 33, 37, 32, 32]
    norm_list1 = [(item*100/45008) for item in list1]
    # E cut (K, R, C labeled) i.e glucut cyslab.ver2	# Nrset of peptides = 52,297
    list2 = [45019, 3419, 1124, 605, 363, 242, 178, 144, 108, 89, 71, 67, 38, 47, 48, 39, 38, 39, 29, 27]
    #list2 = [31768, 13952, 4851, 1790, 883, 531, 425, 325, 220, 169, 158, 132, 96,  87, 73, 59, 46, 40, 47, 41]
    norm_list2 = [(item*100/52297) for item in list2]
    # P cut (K, R, C labeled) i.e procut_cyslab.ver2
    list3 = [46034, 2892, 865, 447, 255, 151, 132, 96, 89, 84, 50, 46, 42, 35, 27, 22, 22, 21, 21, 17]
    #list3 = [32876, 14171, 4688, 1581, 735, 437, 354, 272, 167, 136, 112, 106, 63,  59, 43, 51, 46, 41, 26, 34]
    norm_list3 = [(item*100/51731) for item in list3]
    # Pro cut No labelled cys
    list4 = [40565, 2820, 848, 494, 241, 189, 131, 109, 87, 79, 69, 42, 48, 39, 37, 23, 27, 28, 23, 10]
    norm_list4 = [(item*100/46307) for item in list4]
    print norm_list1
    print norm_list2
    print norm_list3
    print norm_list4
    
    fig = p.figure(figsize=(15,15))
    ax1 = fig.add_subplot(1,1,1)
    p.xlabel('# Counts ')
    p.ylabel('# Peptides')
    p.title('Frequency histogram of number of counts of each peptide in the converted peptide library ')
    line1 = p.plot(xlist, norm_list1, linewidth=2.0, color='g',marker='s')
    line2 = p.plot(xlist, norm_list2, linewidth=2.0, color='b',marker='^')
    line3 = p.plot(xlist, norm_list3, linewidth=2.0, color='r',marker='o')
    line4 = p.plot(xlist, norm_list4, linewidth=2.0, color='k',marker='*')
    p.legend((line1,line2,line3, line4),('E cut; K,R labeled', 'E cut; K, R, C labeled', 'P cut; K, R, C labeled', 'P cut; K, R labeled'), loc=4)
    fig_dir = os.getcwd() + '/figures/'
    fig_name = fig_dir + 'yeast_bin_peptide_freq_lineplot.png'
    p.savefig(fig_name,format = "png", orientation='landscape')
    
    
	
    
    





if __name__ == '__main__':
    
    argument = sys.argv[1:]
    
    print "Processing %s ..."%(argument)
    main(argument)
    
    import time
    print "Script - pept_lib_yeast.py %s \t Completed \t %s"%(argument, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))
    