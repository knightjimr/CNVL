
import os,sys,csv
import re
import argparse
import subprocess
import numpy as np
from numpy import genfromtxt


def generate_breakpoint(bpfile, challfile):
    cov = []
    with open(challfile) as fp:
        cr = csv.reader(fp, delimiter='\t')
        for row in cr:
            cov.append(row)

    dict={}
    for ch in range(1,23)+[24]:
        dict.update({ch:[]})
    for i in range(1,len(cov)):
        dict[int(cov[i][0])].append([cov[i][1],int(cov[i][3]),int(cov[i][4]),float(cov[i][8])])

    breaks = []
    with open(bpfile) as fp:
        br = csv.reader(fp, delimiter='\t')
        for row in br:
            breaks.append(row)

    boundlist=[]
    for k in range(1,len(breaks)):
        chr=int(breaks[k][1])
        start=int(breaks[k][2])
        end=int(breaks[k][3])
        num=int(breaks[k][4])
        meanlog=breaks[k][5]
        sdlog=np.std([m[3] for m in dict[chr][start-1:end]])
        startpos=dict[chr][start-1][1]
        endpos=dict[chr][end-1][2]
        length=endpos-startpos
        if length>0:
            boundlist.append([chr,startpos,endpos,length,num,meanlog,sdlog])
    return boundlist

def main():
    parser = argparse.ArgumentParser(description='CBS segmentation')
    parser.add_argument('bpfile',help='breakpoint file from segmentation.R')
    parser.add_argument('challfile',help='chall file from segmentation.R')
    parser.add_argument('samplename')
    args = parser.parse_args()
    bpfile=args.bpfile
    challfile=args.challfile
    samplename=args.samplename
    boundlist=generate_breakpoint(bpfile, challfile)
    fp = open(samplename+'_CBS_calling.txt','w')
    outputbound=csv.writer(fp,delimiter='\t')
    outputbound.writerow(['#chr','start','end','length','#logRmarkers','meanlogR','sdlogR'])
    for line in boundlist:
        outputbound.writerow([str(m) for m in line])
    fp.close()
            
if __name__=='__main__':
    main()
