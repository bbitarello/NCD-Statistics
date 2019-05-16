#Author: Barbara Bitarello
#
# Read in scan data (NCD-with-FD and NCD-no-FD)
#
# Last modified: 21.09.2016
######################################################################################

library(parallel)
library(SOAR)  #speed up workspace loading.
pops<-c("AWS","LWK","YRI","CEU", "FIN","GBR","TSI", "CHB","CHS" ,"JPT","MXL", "CLM","PUR")

######################################################################################
#Part I: read in scan results and save in .RData file for easy manipulation later.####
######################################################################################

CHR<-seq(1:22)

PATH.1<-paste('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/tmp/chr',CHR, '/', sep='')

All.Results=vector("list",length(pops)); names(All.Results)=pops

All.Results.Final=vector("list",length(pops)); names(All.Results.Final)=pops




for (j in 1:length(pops)){

res=vector("list",22); names(res) = paste0("CHR",1:22)

for ( i in 1:22){
        setwd(PATH.1[i])
        TMP.1<-data.frame( Beg.Win=NA, End.Win=NA,Initial_seg_sites=NA, Initial_fds_sites=NA, NCDf1=NA, NCDf2=NA, NCDf3=NA, NCDf4=NA, NCDf5=NA, Nr.SNPs=NA, Nr.FDs=NA)

        as.numeric(system('ls |grep bin -c', intern=T))->nBINS

        badbin<-NA

                for (k in 1:nBINS){

                        try(load(paste(PATH.1[i], 'bin', k,'/', 'res__',CHR[i],'_',k,'_scan_',pops[j],'.RData', sep='')))->tmp

                                if(inherits(tmp, "try-error"))

                                        badbin<-c(badbin,k)
                                                next
                                }
        a<-seq(1:nBINS)

        NAMES.A1<-paste('res__',CHR[i], '_',a,'_scan_',pops[j],sep='')

        badbin<-badbin[-1]
        if(length(badbin)>0){
        a1<-a[-which(a %in% badbin)]
        }
        if(length(badbin)==0){a1<-a
        }

for (w in a1){

TMP.1<-rbind(TMP.1, get(NAMES.A1[w]))
}
res<- TMP.1[-1,]
All.Results[[j]][[i]] <- res
}
}
remove(list=ls(pattern='res__'))

#join results from all chromosomes
for (j in 1:length(pops)){

        for (i in 1:22){

                All.Results[[j]][[i]]<-cbind(data.frame(Chr=rep(CHR[i]),All.Results[[j]][[i]]))
}       }

#set directory to save this workspace and save the workspace so it can easily be loaded later.

setwd('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/read_scan_data')
for (j in 1:length(pops)){

        do.call(rbind, All.Results[[j]])->All.Results.Final[[j]]
}
##
#*******************************************************************************************
#Sept 2016: this part is obsolete and we don't have an update version of this window coverage bedfile and actually we don't need it anymore.
#*******************************************************************************************

#include window coverage of the inputs in the output

#cov.win<-read.table('windows_coordinates_cov.bed.gz')

#names(cov.win)<-c('CHR', 'Beg.Win', 'End.Win', 'Nr.Map.Seg', 'Total.Cov.Leng', 'Total.Win.Leng', 'Proportion.Covered')

#add coverage to dataframes

#cov.win[order(cov.win$CHR, cov.win$Beg.Win),]-> cov.win2   #the coverage values are not in oder (the windows)

#now windows are sorted in NCV output. 

#remove(cov.win)

#mclapply(All.Results.Final, function(x) cbind(x, Proportion.Covered=cov.win2$Proportion.Covered))->All.Results.Final

objectName<-'All.Results.Final'

save(list=objectName, file= 'All.Results.Final.RData')   #this workspace is already saved in the directory.

#remove(All.Results.Final)

#now the workspace with NCV results for all pops is saved. 
##################################################################
##################################################################

######################################################################################
#I will commnet this entire section because it refers to NCD-no-FD, and we don't actually go into it.

#CHR<-seq(1:22)

#PATH.1<-paste('/mnt/scratch/barbara/ncv_allpops_no_FD/',CHR, '/', sep='')

#All.Results=vector("list",length(pops)); names(All.Results)=pops

#All.Results.Final.no.FD=vector("list", length(pops)); names(All.Results.Final.no.FD)=pops

#for (j in 1:length(pops)){

#res=vector("list",22); names(res) = paste0("CHR",1:22)

#for ( i in 1:22){
#	setwd(PATH.1[i])
#	TMP.1<--data.frame( Beg.Win=NA, End.Win=NA,Initial_seg_sites=NA, Initial_fds_sites=NA, NCVf1=NA, NCVf2=NA, NCVf3=NA, NCVf4=NA, NCVf5=NA, Nr.SNPs=NA, Nr.FDs=NA)
#	as.numeric(system('ls |grep bin -c', intern=T))->nBINS

#	badbin<-NA

#		for (k in 1:nBINS){

#			try(load(paste(PATH.1[i], 'bin', k,'/', 'res__',CHR[i],'_',k,'scanv7_',pops[j],'.RData', sep='')))->tmp

#				if(inherits(tmp, "try-error"))


#					badbin<-c(badbin,k)
#						next
#				}
#	a<-seq(1:nBINS)

#	NAMES.A1<-paste('res__',CHR[i], '_',a,'scanv7_',pops[j],sep='')
	
#	badbin<-badbin[-1]
#	if(length(badbin)>0){
#	a1<-a[-which(a %in% badbin)]
#	}
#	if(length(badbin)==0){a1<-a
#	}

#for (w in a1){

#TMP.1<-rbind(TMP.1, get(NAMES.A1[w]))
#}
#res<- TMP.1[-1,]
#All.Results[[j]][[i]] <- res
#}
#}
#remove(list=ls(pattern='res__'))

#for (j in 1:length(pops)){

#	for (i in 1:22){

#		All.Results[[j]][[i]]<-cbind(data.frame(Chr=rep(CHR[i]),All.Results[[j]][[i]]))
#}	}

#set directory to save this workspace.
#setwd('/mnt/sequencedb/PopGen/barbara/scan_may_2014')
#for (j in 1:length(pops)){

#	do.call(rbind, All.Results[[j]])->All.Results.Final.no.FD[[j]]
#}


#include window coverage of the inputs in the output

#mclapply(All.Results.Final.no.FD, function(x) cbind(x, Proportion.Covered=cov.win2$Proportion.Covered))->All.Results.Final.no.FD


#objectName<-'All.Results.Final.no.FD'

#save(list=objectName, file= 'All.Results.Final.no.FD.RData')   #this workspace is already saved in the directory.

#now the workspace with NCV results for all pops is saved. 
##################################################################
##################################################################
##################################################################

remove(All.Results)
