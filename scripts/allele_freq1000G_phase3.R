#!/usr/bin/env Rscript
#Cesare de Filippo, MPI-EVA
## Date: 03/08/2017
## Task: Generate ALternetive allele frequency per populations from vcf files of the 1000G phase 3.
#Modified by Barbara Bitarello May/2019

## The script is called after piping:
## zgrep -v ^# VCF > test.vcf
##./allele_freq1000G_phase3.R test.vcf
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (a name for this run).n", call.=FALSE)
}
library("optparse")
info <- read.table("/project/mathilab/data/1kg/20130502_phase3_final/integrated_call_samples_v3.20130502.ALL.panel", as.is=T, header=T)
popN <- table(info[,2])
nt=c("A", "C", "G", "T")
## An example of "info" or *.panel file
##  sample pop super_pop gender
## HG00096 GBR       EUR   male
## HG00097 GBR       EUR female
## HG00099 GBR       EUR female
## HG00100 GBR       EUR female
## HG00101 GBR       EUR   male
## HG00102 GBR       EUR female

## Allele separator in the GT field
SEPARATOR="\\|"

options(warn=2)
cat(paste(c("#chr", "pos", "ref", "alt", names(popN)),collapse="\t",sep=""),"\n")
f <- file(args[1])
open(f)
while(length(line <- readLines(f,n=1)) > 0) {
    x <- unlist(strsplit(line,split="\t"))
    g <- x[10:length(x)]
    ref <- x[4]; alt=unlist(strsplit(x[5], split=","))
    if(sum(unique(nchar(alt)) > 1) == 0) { ## disregard indels
        g <- gsub("0",ref, g); for (i in 1:length(alt)) { g <- gsub(i, alt[i], g)}
        res <- vector("list", length(popN)); names(res) <- names(popN)
        for (i in names(popN)) {
            ids <- which(info[,2] == i)
            a <- unlist(strsplit(g[ids], split=SEPARATOR))
            res[[i]] <- round(sum(a %in% alt)/sum(a %in% c(ref,alt)),4)
        }
        res <- unlist(res)
        cat(paste(c(x[c(1:2,4,5)], res),collapse="\t"),"\n",sep="")
    }
}

