#!/bin/env python

from sys import *
from getopt import getopt

def printUsageAndExit(programName):
    print >> stderr,programName,"[options] matrix > normalizedMatrix"
    print >> stderr,"options:"
    print >> stderr,"--factor f. multiply normalized values by this factor [1.0]"
    print >> stderr,"--col c. the column number of the value to be normalized [8]"
    print >> stderr,"--start-row r. The starting row for normalization in order to skip header in some cases. [1]"
    print >> stderr,"--normalizer filename:col. specify normalizer filename and column [default uses the matrix and col specified]"
    exit(1)


if __name__=='__main__':
    programName=argv[0]
    opts,args=getopt(argv[1:],'',['factor=','col=','start-row=','normalizer='])
    
    factor=1.0
    col=7
    normalizerCol=-1
    startRow=1
    normalizerFilename=None
    
    
    
    for o,v in opts:
        if o=='--factor':
            factor=float(v)
        elif o=='--col':
            col=int(v)-1 #convert to 0-based
        elif o=='--start-row':
            startRow=int(v)
        elif o=='--normalizer':
            normalizerFilename,normalizerCol=v.split(":")
            normalizerCol=int(normalizerCol)-1




    try:
        filename,=args
    except:
        printUsageAndExit(programName)

    if normalizerCol<0:
        normalizerCol=col
    
    if not normalizerFilename:
        #use same file for the normalizer
        normalizerFilename=filename

    #An naive approach:
    #first pass: read all values and find sum
    #second pass replace values with value/sum*factor

    sum=0.0
    
    
    
    #first pass
    lino=0
    fil=open(normalizerFilename)
    for lin in fil:
        lino+=1
        if lino<startRow:
            continue
        
        lin=lin.rstrip()
        fields=lin.split("\t")
        value=float(fields[normalizerCol])
        sum+=value
        if lino==startRow:
            minV=value
            maxV=value
        else:
            minV=min(minV,value)
            maxV=max(maxV,value)

    fil.close()
    
    print >> stderr,"first pass: lino=%d sum=%f min=%f max=%f" %(lino,sum,minV,maxV)
    
    #second pass
    lino=0
    normalizedSum=0.0
    
    fil=open(filename)
    for lin in fil:
        lino+=1
        lin=lin.rstrip()
        if lino<startRow: #directly output this line
            print >> stdout,lin
            continue
        
        fields=lin.split("\t")
        value=float(fields[col])
        normalizedValue=value/sum*factor
        #now replace
        fields[col]=str(normalizedValue)
        #now output the fields joined by tab
        print >> stdout,"\t".join(fields)
        normalizedSum+=normalizedValue
        if lino==startRow:
            minV=normalizedValue
            maxV=normalizedValue
        else:
            minV=min(minV,normalizedValue)
            maxV=max(maxV,normalizedValue)

    fil.close()
    
    print >> stderr,"second pass normalized: lino=%d factor=%f sum=%f min=%f max=%f" %(lino,factor,normalizedSum,minV,maxV)

