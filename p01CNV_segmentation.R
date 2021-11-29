
library("DNAcopy")

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
   stop("Usage: p01CNV_segmentation.R filePrefix\n")
}

prefix <- args[1]
ratio.file <- paste(prefix, "_covRatio.txt", sep="")
bp.file <- paste(prefix, "_covRatio.breakpoint.txt", sep="")
chall.file <- paste(prefix, "_covRatio.chall.txt", sep="")

cov <- read.table(ratio.file, sep="\t", header=F)
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

write.table(breakpoint[2], file=bp.file, sep="\t", row.names=F, quote=F)

write.table(ch.all, file=chall.file, sep="\t", row.names=F, quote=F)

