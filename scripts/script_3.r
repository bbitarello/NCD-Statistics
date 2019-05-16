################################################################################
#	Barbara D Bitarello
#
#	Last modified: 17.10.2016
#
#	Read in bin simulations, etc
##############################################################################
library(parallel)
library(SOAR)
library(ggplot2)
library(vioplot)
library(dplyr)
library(reshape)
library(reshape2)
Sys.setenv(R_LOCAL_CACHE="estsession")
library(data.table)

#first, load the scan data
##########################skip the next block, as it has already been saved in the R object ############################
pops<-c("AWS","LWK","YRI","CEU", "FIN","GBR","TSI", "CHB","CHS" ,"JPT","MXL", "CLM","PUR")

#I copied the sims from cee's directory: /mnt/scratch/cee/bs_genomescan/simulations/SuSt/

system.time(list.MSMS.rec.1e_09<-mclapply(c(1:6), function(x) setDT(read.table(paste0("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/cesare_sims/neutral_n100_mu1e-07_rho1e-09_3000bp.downsampled.allstats",x,".tsv.gz"), header=T))))
#1937.807
list.MSMS<-vector('list', 3)

do.call('rbind', list.MSMS.rec.1e_09)-> list.MSMS[[1]]

remove(list.MSMS.rec.1e_09)
#
list.MSMS.rec.1e_08<-mclapply(c(1:6), function(x) setDT(read.table(paste0("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/cesare_sims/neutral_n100_mu1e-07_rho1e-08_3000bp.downsampled.allstats",x,".tsv.gz"), header=T)))

do.call('rbind', list.MSMS.rec.1e_08)-> list.MSMS[[2]]

remove(list.MSMS.rec.1e_08)

list.MSMS.rec.1e_07<-mclapply(c(1:6), function(x) setDT(read.table(paste0("/mnt/sequencedb/PopGen/barbara/simulations/msms_sims/cesare_sims/neutral_n100_mu1e-07_rho1e-07_3000bp.downsampled.allstats",x,".tsv.gz"), header=T)))

do.call('rbind', list.MSMS.rec.1e_07)-> list.MSMS[[3]]

remove(list.MSMS.rec.1e_07)

#
#mclapply(list.MSMS, function(x) cbind(x, Nr.IS=x$S+x$FD))-> list.MSMS #not enough memory...

cbind(list.MSMS[[1]], Nr.IS=list.MSMS[[1]]$S+list.MSMS[[1]]$FD)-> list.MSMS1[[1]]
cbind(list.MSMS[[2]], Nr.IS=list.MSMS[[2]]$S+list.MSMS[[2]]$FD)-> list.MSMS1[[2]]
cbind(list.MSMS[[3]], Nr.IS=list.MSMS[[3]]$S+list.MSMS[[3]]$FD)-> list.MSMS1[[3]]


list.MSMS1->list.MSMS

remove(list.MSMS1)

#first separate AF< EU, AS
names(list.MSMS)<-c('rec.1e_09', 'rec.1e_08', 'rec.1e_07')

AFRICA<-vector('list', 3)
names(AFRICA)<-c('rec.1e_09', 'rec.1e_08', 'rec.1e_07')
EUROPE<-vector('list', 3)
names(EUROPE)<-c('rec.1e_09', 'rec.1e_0.8', 'rec.1e_07')
ASIA<-vector('list', 3)
names(ASIA)<-c('rec.1e_09', 'rec.1e_0.8', 'rec.1e_07')
#58,000 simulations for AFRICA for each rec rate.
AFRICA[[1]]<-subset(list.MSMS[[1]], pop==0)
AFRICA[[2]]<-subset(list.MSMS[[2]], pop==0)
AFRICA[[3]]<-subset(list.MSMS[[3]], pop==0)

EUROPE[[1]]<-subset(list.MSMS[[1]], pop==1)
EUROPE[[2]]<-subset(list.MSMS[[2]], pop==1)
EUROPE[[3]]<-subset(list.MSMS[[3]], pop==1)

ASIA[[1]]<-subset(list.MSMS[[1]], pop==2)
ASIA[[2]]<-subset(list.MSMS[[2]], pop==2)
ASIA[[3]]<-subset(list.MSMS[[3]], pop==2)
#separate sims in bins of Nr.Inf. SItes

Store(list.MSMS)
Store(ASIA)
Store(EUROPE)
Store(AFRICA)

###
