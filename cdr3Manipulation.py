#! usr/bin/python


if __name__ == '__main__':
    
    argument = sys.argv[1:]
    
    print "Processing %s ..."%(argument)
    idSeq_dict = loadRelevantFiles('cdr3')
    
    main(argument)
    
    import time
    print "Script - %s \t Completed \t %s"%(sys.argv, time.strftime("%d %b %Y %H:%M:%S",time.localtime()))