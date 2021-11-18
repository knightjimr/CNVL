'''
Created on Feb 11, 2013

@author: siming
'''
import os,sys,csv
import re
import argparse
import rpy2.robjects as robjects
import numpy as np
from numpy import genfromtxt


def generate_breakpoint(covfile):
    robjects.r('''
    library("DNAcopy")
    cov <- read.table("%s", sep="\t", header=F)
    ch.all <- {}
    for (ch in c(1:22,24)) {
        ch.temp <- cov[cov[,1]==paste("chr",ch,sep=""),]
        ch.temp <- cbind(seq(1:dim(ch.temp)[1]),ch.temp)
        ch.temp <- cbind(rep(ch,dim(ch.temp)[1]),ch.temp)
        ch.all <- rbind(ch.all,ch.temp)
        }
    CNA.object <- CNA(ch.all[,9], ch.all[,1], ch.all[,2])
    smoothed.CNA.object <- smooth.CNA(CNA.object)
    breakpoint <- segment(smoothed.CNA.object,min.width=5,alpha=0.00001,undo.splits="sdundo",undo.SD=1.0,verbose=2)
    '''%(covfile))
    cov=robjects.r["ch.all"]
    breaks=robjects.r["breakpoint"]
    breaks=breaks[1]
    dict={}
    for ch in range(1,23)+[24]:
        dict.update({ch:[]})
    for i in range(len(cov[1])):
            dict[cov[0][i]].append([cov[1][i],cov[3][i],cov[4][i],cov[8][i]])
    boundlist=[]
    for k in range(len(breaks[1])):
        chr=int(breaks[1][k])
        start=breaks[2][k]
        end=breaks[3][k]
        num=int(breaks[4][k])
        meanlog=breaks[5][k]
        sdlog=np.std([m[3] for m in dict[chr][start-1:end]])
        startpos=dict[chr][start-1][1]
        endpos=dict[chr][end-1][2]
        length=endpos-startpos
        if length>0:
            boundlist.append([chr,startpos,endpos,length,num,meanlog,sdlog])
    return boundlist

def main():
    parser = argparse.ArgumentParser(description='CBS segmentation')
    parser.add_argument('covfile',help='cov.ratio file from Murim p01TumorCNV')
    parser.add_argument('samplename')
    args = parser.parse_args()
    covfile=args.covfile
    samplename=args.samplename
    boundlist=generate_breakpoint(covfile)
    fp = open(samplename+'_CBS_calling.txt','w')
    outputbound=csv.writer(fp,delimiter='\t')
    outputbound.writerow(['#chr','start','end','length','#logRmarkers','meanlogR','sdlogR'])
    for line in boundlist:
        outputbound.writerow([str(m) for m in line])
    fp.close()
            
if __name__=='__main__':
    main()
