##################################
#	Bárbara D.Bitarello
#
#	Last modified: 15.12.2016
##################################

library(parallel)
library(SOAR)
library(ggplot2)
library(vioplot)
library(dplyr)
library(reshape)
library(reshape2)
Sys.setenv(R_LOCAL_CACHE="estsession")
library(data.table)
library(qqman)

pops<-c("AWS","LWK","YRI","CEU", "FIN","GBR","TSI", "CHB","CHS" ,"JPT","MXL", "CLM","PUR")
read.table('/mnt/sequencedb/PopGen/cesare/hg19/bedfiles/ensembl_genes_hg19.bed.gz')->hg19.coding.coords.bed
names(hg19.coding.coords.bed)<-c('chr', 'beg', 'end','name', 'type')


my.function.improved<-function(B, E, df=LWK.chr, chr=6){
as.numeric(gsub("chr", "", chr))-> chr
rbind(filter(df[[chr]], Chr==chr,End.Win > B,End.Win < E), filter(df[[chr]], Chr==chr,Beg.Win > B,Beg.Win < E), filter(df[[chr]], Chr==chr,Beg.Win<B, End.Win>E), filter(df[[chr]], Chr==chr,Beg.Win>B ,End.Win<E))->res
setDT(res)
setkey(res, Win.ID)
#df[rownames(res[!duplicated(res),]),]-> res2
unique(res)-> res2
return(res2)
}

# find.gene - that's literally what it does. But not just protein coding genes.
#
source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/find_gene.R')
#
source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/assign_tf_window_original.R')
####
####
Objects()


#lapply(1:22, function(x) subset(prd.bed, chr==x))-> prot.cod.bed.list
mclapply(1:22, function(x) setDT(subset(hg19.coding.coords.bed, chr==x)))-> coding.per.chr.list 

Store(coding.per.chr.list)

mclapply(1:22, function(x) within(coding.per.chr.list[[x]], comp.names<-paste0(name, ":", type))$comp.names)->names.all.coding

Store(names.all.coding)

list.SCAN[[2]]-> LWK;list.SCAN[[3]]-> YRI

list.SCAN[[6]]-> GBR;list.SCAN[[7]]-> TSI

setDT(LWK); setDT(YRI); setDT(GBR); setDT(TSI)

mclapply(1:22, function(x) setDT(filter(LWK, Chr==x)))-> LWK.chr

mclapply(1:22, function(x) setDT(filter(YRI, Chr==x)))-> YRI.chr

mclapply(1:22, function(x) setDT(filter(GBR, Chr==x)))-> GBR.chr

mclapply(1:22, function(x) setDT(filter(TSI, Chr==x)))-> TSI.chr

#test for chr 22:

system.time(mclapply(1:nrow(coding.per.chr.list[[22]]), function(x) my.function.improved(B=coding.per.chr.list[[22]]$beg[x], E=coding.per.chr.list[[22]]$end[x], chr=paste0("chr", 22), df=LWK.chr))-> test)  #9 seconds!!
#

LWK.win<-vector('list', 22)
YRI.win<-vector('list', 22)
GBR.win<-vector('list', 22)
TSI.win<-vector('list', 22)


system.time(for ( i in 1:22){
mclapply(1:nrow(coding.per.chr.list[[i]]), function(x) my.function.improved(B=coding.per.chr.list[[i]]$beg[x], E=coding.per.chr.list[[i]]$end[x], chr=paste0("chr", i), df=LWK.chr))-> LWK.win[[i]]
gc()
mclapply(1:nrow(coding.per.chr.list[[i]]), function(x) my.function.improved(B=coding.per.chr.list[[i]]$beg[x], E=coding.per.chr.list[[i]]$end[x], chr=paste0("chr", i), df=YRI.chr))-> YRI.win[[i]]
gc()
mclapply(1:nrow(coding.per.chr.list[[i]]), function(x) my.function.improved(B=coding.per.chr.list[[i]]$beg[x], E=coding.per.chr.list[[i]]$end[x], chr=paste0("chr", i), df=GBR.chr))-> GBR.win[[i]]
gc()
mclapply(1:nrow(coding.per.chr.list[[i]]), function(x) my.function.improved(B=coding.per.chr.list[[i]]$beg[x], E=coding.per.chr.list[[i]]$end[x], chr=paste0("chr", i), df=TSI.chr))-> TSI.win[[i]]
gc()
cat ('chromosome ',i, ' done\n')
})

gc() #2143.112  great!

#stopped here

Objects()

Store(YRI.win, LWK.win, GBR.win, TSI.win)

Objects()

for(i in 1:22){

names.all.coding[[i]]-> names(LWK.win[[i]])
names.all.coding[[i]]-> names(YRI.win[[i]])
names.all.coding[[i]]-> names(GBR.win[[i]])
names.all.coding[[i]]-> names(TSI.win[[i]])}


Store(YRI.win, LWK.win, GBR.win, TSI.win)



#playground


system.time(mclapply(1:length(as.character(filter(hg19.coding.coords.bed, type=='protein_coding')$name)), function(x) try(find.gene(name1=as.character(filter(hg19.coding.coords.bed, type=='protein_coding')$name)[x])))-> test)

mclapply(test, function(x) try(assign.tf(x$query_subset)))-> test.2

#########################################################################################################################
#manhattan plots
source('../scripts/my.man.R')
#mclapply(Union.top0.5_0.4_0.3, function(x) select(x, Chr:End.Win, Dist.NCD.f0.5, Z.f0.5.P.val))->tes.manhattan.union
mclapply(top829f0.5, function(x) select(x, Chr:End.Win, Dist.NCD.f0.5, Z.f0.5.P.val))-> tes.manha.f0.5.top
mclapply(CANDf0.5, function(x) select(x, Chr:End.Win, Dist.NCD.f0.5, Z.f0.5.P.val))->  tes.manha.f0.5.cand
#mclapply(tes.manhattan.union, function(x) cbind(x, BP=(x$Beg.Win+x$End.Win)/2))->tes.manhattan.2.union
mclapply(tes.manha.f0.5.top, function(x) cbind(x, BP=(x$Beg.Win+x$End.Win)/2))->tes.manha.f0.5.2.top
mclapply(tes.manha.f0.5.cand, function(x) cbind(x, BP=(x$Beg.Win+x$End.Win)/2))->tes.manha.f0.5.2.cand


#for(i in 1:7){
#as.numeric(gsub("chr","", tes.manhattan.2.union[[i]]$Chr))->tes.manhattan.2.union[[i]]$Chr

#colnames(tes.manhattan.2.union[[i]])[c(1,4,5)]<-c('CHR', 'SNP', 'P')}

for(i in 1:7){
colnames(tes.manha.f0.5.2.top[[i]])[c(1,4,5)]<-c('CHR', 'SNP', 'P')
colnames(tes.manha.f0.5.2.cand[[i]])[c(1,4,5)]<-c('CHR', 'SNP', 'P')
}

mclapply(tes.manha.f0.5.2.top, function(x) setDT(arrange(x, P)))-> tes.manha.f0.5.top
mclapply(tes.manha.f0.5.2.cand, function(x) setDT(arrange(x, P)))-> tes.manha.f0.5.cand

#mclapply(tes.manhattan.2.union, function(x) setDT(arrange(x, P)))->tes.manhattan.union

#union.man.top829<-mclapply(tes.manhattan.union,function(x) head(x,829))

#mclapply(union.man.top829, function(x) arrange(x, CHR, Beg.Win))->sort.union.man.top829 

mclapply(tes.manha.f0.5.top, function(x) arrange(x, CHR, Beg.Win))->sort.man.top829 
mclapply(tes.manha.f0.5.cand, function(x) arrange(x, CHR, Beg.Win))->sort.man.cand
#mclapply(tes.manhattan.union, function(x) arrange(x, CHR, Beg.Win))->sort.tes.manhattan.union



#now background:

select(list.SCAN[[2]],Chr:End.Win, Dist.NCD.f0.5, Z.f0.5.P.val)-> LWK.bg
setDT(LWK.bg)

with(LWK.bg, cbind(LWK.bg, BP=(Beg.Win+End.Win)/2))-> LWK.bg

colnames(LWK.bg)[c(1,4,5)]<-c('CHR', 'SNP', 'P')

arrange(LWK.bg, CHR, Beg.Win)-> LWK.bg2
remove(LWK.bg)

#pdf('bedfiles/Union.top829.my.manhattan.test.pdf') #or 
png("figures/manhattan.png", width=3000, height=1000, pointsize=18)
as.numeric(LWK.bg2$CHR)->LWK.bg2$CHR
my.manhattan(LWK.bg2,highlight=as.character(sort.man.top829[[2]]$SNP),highlight2=as.character(sort.man.cand[[2]]$SNP), suggestiveline=F,genomewideline=F, cex.axis=1.2, cex.lab=1.2)
#legend('topright', c("0.05% cutoff","sim-based cutoff"),col=c('violetred1', 'darkorange'), lty=1) #legend does not work
dev.off()


png("figures/aida_manhattan.png", width=3000, height=1000, pointsize=18)


as.numeric(LWK.bg2$CHR)->LWK.bg2$CHR
my.manhattan(LWK.bg2,highlight2=as.character(sort.man.cand[[2]]$SNP), suggestiveline=F,genomewideline=F, cex.axis=1.2, cex.lab=1.2)
dev.off()


###
