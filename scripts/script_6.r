###############################################
#
#	Barbara D Bitarello
#
#	Last modified: 20.10.2016
#
#	Make bedfiles for Joao
#################################################


library(parallel)
library(SOAR)
library(ggplot2)
library(vioplot)
library(dplyr)
library(reshape)
library(reshape2)
Sys.setenv(R_LOCAL_CACHE="estsession")
library(data.table)

pops<-c("AWS","LWK","YRI","CEU", "FIN","GBR","TSI", "CHB","CHS" ,"JPT","MXL", "CLM","PUR")

BED.PATH<-'/mnt/sequencedb/PopGen/barbara/NCV_dir_package/read_scan_data/bedfiles/'

#generate a background bed file (all scanned genes)

###MAKE UNION OF TF FOR EACH POP ##########
mclapply(1:7, function(x) setDT(rbind(CANDf0.5[[x]], CANDf0.4[[x]], CANDf0.3[[x]])[-(which(duplicated(rbind(CANDf0.5[[x]], CANDf0.4[[x]], CANDf0.3[[x]])))),]))-> Union.CANDf0.5_0.4_0.3

mclapply(1:7, function(x) setDT(rbind(top829f0.5[[x]], top829f0.4[[x]], top829f0.3[[x]])[-(which(duplicated(rbind(top829f0.5[[x]], top829f0.4[[x]], top829f0.3[[x]])))),]))-> Union.top0.5_0.4_0.3

Store(Union.CANDf0.5_0.4_0.3)
Store(Union.top0.5_0.4_0.3)

Objects()

#######################
#a few sanity checks

####
#test
for(i in 1:7){
nrow(list.SCAN[[i]])-> n1;nrow(Union.CANDf0.5_0.4_0.3[[i]])-> n2;nrow(Union.top0.5_0.4_0.3[[i]])-> n3

cbind(rbind(select(list.SCAN[[i]], Nr.IS, Nr.FDs, Nr.SNPs, PtoD), select(Union.CANDf0.5_0.4_0.3[[i]], Nr.IS, Nr.FDs, Nr.SNPs, PtoD), select(Union.top0.5_0.4_0.3[[i]], Nr.IS, Nr.FDs, Nr.SNPs, PtoD)), Type=c(rep('genomic', n1), rep('Significant', n2), rep('Outliers', n3)))-> temp

tmpname<-paste0('Nr.IS.',pops[[i]], '.pdf')
tmpname1<-paste0('Nr.FDs.',pops[[i]], '.pdf')
tmpname2<-paste0('Nr.SNPs.',pops[[i]], '.pdf')
tmpname3<-paste0('PtoD.',pops[[i]], '.pdf')

ggplot(temp) + geom_density(aes(x = Nr.IS, colour = Type))
ggsave(paste0('figures/',tmpname))

ggplot(temp) + geom_density(aes(x = Nr.FDs, colour = Type))
ggsave(paste0('figures/',tmpname1))

ggplot(temp) + geom_density(aes(x = Nr.SNPs, colour = Type))
ggsave(paste0('figures/',tmpname2))

ggplot(temp) + geom_density(aes(x = PtoD, colour = Type))
ggsave(paste0('figures/',tmpname3))

cat('Finished', pops[i], '\n')
}


#

#plots SFS

source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/SFS_script.r')

#tf=0.5

#unlist(mclapply(1:nrow(top829f0.5[[2]]), function(x) SFS.function(CHR=top829f0.5[[2]][x,]$Chr, BEG=top829f0.5[[2]][x,]$Beg.Win, END=top829f0.5[[2]][x,]$End.Win, POP=2)))-> LWK.top.f0.5
#LWK.top.f0.5/100-> LWK.top.f0.5

unlist(mclapply(1:nrow(merge.top829f0.5[[2]]), function(x) SFS.function(CHR=merge.top829f0.5[[2]][x,]$V1, BEG=merge.top829f0.5[[2]][x,]$V2, END=merge.top829f0.5[[2]][x,]$V3, POP=2)))-> LWK.top.f0.5
LWK.top.f0.5/100-> LWK.top.f0.5
unlist(mclapply(1:nrow(merge.top829f0.5[[3]]), function(x) SFS.function(CHR=merge.top829f0.5[[3]][x,]$V1, BEG=merge.top829f0.5[[3]][x,]$V2, END=merge.top829f0.5[[3]][x,]$V3, POP=3)))->YRI.top.f0.5
YRI.top.f0.5/100-> YRI.top.f0.5
unlist(mclapply(1:nrow(merge.top829f0.5[[6]]), function(x) SFS.function(CHR=merge.top829f0.5[[6]][x,]$V1, BEG=merge.top829f0.5[[6]][x,]$V2, END=merge.top829f0.5[[6]][x,]$V3, POP=6)))->GBR.top.f0.5
GBR.top.f0.5/100-> GBR.top.f0.5
unlist(mclapply(1:nrow(merge.top829f0.5[[7]]), function(x) SFS.function(CHR=merge.top829f0.5[[7]][x,]$V1, BEG=merge.top829f0.5[[7]][x,]$V2, END=merge.top829f0.5[[7]][x,]$V3, POP=7)))->TSI.top.f0.5
TSI.top.f0.5/100-> TSI.top.f0.5
#
unlist(mclapply(1:nrow(merge.CANDf0.5[[2]]), function(x) SFS.function(CHR=merge.CANDf0.5[[2]][x,]$V1, BEG=merge.CANDf0.5[[2]][x,]$V2, END=merge.CANDf0.5[[2]][x,]$V3, POP=2)))->LWK.cand.f0.5
LWK.cand.f0.5/100-> LWK.cand.f0.5

unlist(mclapply(1:nrow(merge.CANDf0.5[[3]]), function(x) SFS.function(CHR=merge.CANDf0.5[[3]][x,]$V1, BEG=merge.CANDf0.5[[3]][x,]$V2, END=merge.CANDf0.5[[3]][x,]$V3, POP=3)))->YRI.cand.f0.5
YRI.cand.f0.5/100-> YRI.cand.f0.5

unlist(mclapply(1:nrow(merge.CANDf0.5[[6]]), function(x) SFS.function(CHR=merge.CANDf0.5[[6]][x,]$V1, BEG=merge.CANDf0.5[[6]][x,]$V2, END=merge.CANDf0.5[[6]][x,]$V3, POP=6)))->GBR.cand.f0.5
GBR.cand.f0.5/100-> GBR.cand.f0.5

unlist(mclapply(1:nrow(merge.CANDf0.5[[7]]), function(x) SFS.function(CHR=merge.CANDf0.5[[7]][x,]$V1, BEG=merge.CANDf0.5[[7]][x,]$V2, END=merge.CANDf0.5[[7]][x,]$V3, POP=7)))->TSI.cand.f0.5
TSI.cand.f0.5/100-> TSI.cand.f0.5

#tf=0.3


#unlist(mclapply(1:nrow(top829f0.3[[2]]), function(x) SFS.function(CHR=top829f0.3[[2]][x,]$Chr, BEG=top829f0.3[[2]][x,]$Beg.Win, END=top829f0.3[[2]][x,]$End.Win, POP=2)))-> LWK.top.f0.3
#LWK.top.f0.3/100->LWK.top.f0.3
#unlist(mclapply(1:nrow(top829f0.3[[3]]), function(x) SFS.function(CHR=top829f0.3[[3]][x,]$Chr, BEG=top829f0.3[[3]][x,]$Beg.Win, END=top829f0.3[[3]][x,]$End.Win, POP=3)))-> YRI.top.f0.3
#YRI.top.f0.3/100->YRI.top.f0.3
#unlist(mclapply(1:nrow(top829f0.3[[6]]), function(x) SFS.function(CHR=top829f0.3[[6]][x,]$Chr, BEG=top829f0.3[[6]][x,]$Beg.Win, END=top829f0.3[[6]][x,]$End.Win, POP=6)))-> GBR.top.f0.3
#GBR.top.f0.3/100->GBR.top.f0.3
#unlist(mclapply(1:nrow(top829f0.3[[7]]), function(x) SFS.function(CHR=top829f0.3[[7]][x,]$Chr, BEG=top829f0.3[[7]][x,]$Beg.Win, END=top829f0.3[[7]][x,]$End.Win, POP=7)))-> TSI.top.f0.3
#TSI.top.f0.3/100->TSI.top.f0.3

unlist(mclapply(1:nrow(merge.top829f0.3[[2]]), function(x) SFS.function(CHR=merge.top829f0.3[[2]][x,]$V1, BEG=merge.top829f0.3[[2]][x,]$V2, END=merge.top829f0.3[[2]][x,]$V3, POP=2)))-> LWK.top.f0.3
LWK.top.f0.3/100->LWK.top.f0.3
unlist(mclapply(1:nrow(merge.top829f0.3[[3]]), function(x) SFS.function(CHR=merge.top829f0.3[[3]][x,]$V1, BEG=merge.top829f0.3[[3]][x,]$V2, END=merge.top829f0.3[[3]][x,]$V3, POP=3)))-> YRI.top.f0.3
YRI.top.f0.3/100->YRI.top.f0.3
unlist(mclapply(1:nrow(merge.top829f0.3[[6]]), function(x) SFS.function(CHR=merge.top829f0.3[[6]][x,]$V1, BEG=merge.top829f0.3[[6]][x,]$V2, END=merge.top829f0.3[[6]][x,]$V3, POP=6)))-> GBR.top.f0.3
GBR.top.f0.3/100->GBR.top.f0.3
unlist(mclapply(1:nrow(merge.top829f0.3[[7]]), function(x) SFS.function(CHR=merge.top829f0.3[[7]][x,]$V1, BEG=merge.top829f0.3[[7]][x,]$V2, END=merge.top829f0.3[[7]][x,]$V3, POP=7)))-> TSI.top.f0.3
TSI.top.f0.3/100->TSI.top.f0.3
#
#unlist(mclapply(1:nrow(CANDf0.3[[2]]), function(x) SFS.function(CHR=CANDf0.3[[2]][x,]$Chr, BEG=CANDf0.3[[2]][x,]$Beg.Win, END=CANDf0.3[[2]][x,]$End.Win, POP=2)))-> LWK.cand.f0.3
#LWK.cand.f0.3/100->LWK.cand.f0.3
#unlist(mclapply(1:nrow(CANDf0.3[[3]]), function(x) SFS.function(CHR=CANDf0.3[[3]][x,]$Chr, BEG=CANDf0.3[[3]][x,]$Beg.Win, END=CANDf0.3[[3]][x,]$End.Win, POP=3)))-> YRI.cand.f0.3
#YRI.cand.f0.3/100->YRI.cand.f0.3
#unlist(mclapply(1:nrow(CANDf0.3[[6]]), function(x) SFS.function(CHR=CANDf0.3[[6]][x,]$Chr, BEG=CANDf0.3[[6]][x,]$Beg.Win, END=CANDf0.3[[6]][x,]$End.Win, POP=6)))-> GBR.cand.f0.3
#GBR.cand.f0.3/100->GBR.cand.f0.3
#unlist(mclapply(1:nrow(CANDf0.3[[7]]), function(x) SFS.function(CHR=CANDf0.3[[7]][x,]$Chr, BEG=CANDf0.3[[7]][x,]$Beg.Win, END=CANDf0.3[[7]][x,]$End.Win, POP=7)))-> TSI.cand.f0.3
#TSI.cand.f0.3/100->TSI.cand.f0.3


unlist(mclapply(1:nrow(merge.CANDf0.3[[2]]), function(x) SFS.function(CHR=merge.CANDf0.3[[2]][x,]$V1, BEG=merge.CANDf0.3[[2]][x,]$V2, END=merge.CANDf0.3[[2]][x,]$V3, POP=2)))-> LWK.cand.f0.3
LWK.cand.f0.3/100->LWK.cand.f0.3
unlist(mclapply(1:nrow(merge.CANDf0.3[[3]]), function(x) SFS.function(CHR=merge.CANDf0.3[[3]][x,]$V1, BEG=merge.CANDf0.3[[3]][x,]$V2, END=merge.CANDf0.3[[3]][x,]$V3, POP=3)))-> YRI.cand.f0.3
YRI.cand.f0.3/100->YRI.cand.f0.3
unlist(mclapply(1:nrow(merge.CANDf0.3[[6]]), function(x) SFS.function(CHR=merge.CANDf0.3[[6]][x,]$V1, BEG=merge.CANDf0.3[[6]][x,]$V2, END=merge.CANDf0.3[[6]][x,]$V3, POP=6)))-> GBR.cand.f0.3
GBR.cand.f0.3/100->GBR.cand.f0.3
unlist(mclapply(1:nrow(merge.CANDf0.3[[7]]), function(x) SFS.function(CHR=merge.CANDf0.3[[7]][x,]$V1, BEG=merge.CANDf0.3[[7]][x,]$V2, END=merge.CANDf0.3[[7]][x,]$V3, POP=7)))-> TSI.cand.f0.3
TSI.cand.f0.3/100->TSI.cand.f0.3

#tf=0.4

#unlist(mclapply(1:nrow(top829f0.4[[2]]), function(x) SFS.function(CHR=top829f0.4[[2]][x,]$Chr, BEG=top829f0.4[[2]][x,]$Beg.Win, END=top829f0.4[[2]][x,]$End.Win, POP=2)))-> LWK.top.f0.4
#LWK.top.f0.4/100->LWK.top.f0.4
#unlist(mclapply(1:nrow(top829f0.4[[3]]), function(x) SFS.function(CHR=top829f0.4[[3]][x,]$Chr, BEG=top829f0.4[[3]][x,]$Beg.Win, END=top829f0.4[[3]][x,]$End.Win, POP=3)))-> YRI.top.f0.4
#YRI.top.f0.4/100->YRI.top.f0.4
#unlist(mclapply(1:nrow(top829f0.4[[6]]), function(x) SFS.function(CHR=top829f0.4[[6]][x,]$Chr, BEG=top829f0.4[[6]][x,]$Beg.Win, END=top829f0.4[[6]][x,]$End.Win, POP=6)))-> GBR.top.f0.4
#GBR.top.f0.4/100->GBR.top.f0.4
#unlist(mclapply(1:nrow(top829f0.4[[7]]), function(x) SFS.function(CHR=top829f0.4[[7]][x,]$Chr, BEG=top829f0.4[[7]][x,]$Beg.Win, END=top829f0.4[[7]][x,]$End.Win, POP=7)))-> TSI.top.f0.4
#TSI.top.f0.4/100->TSI.top.f0.4
#

unlist(mclapply(1:nrow(merge.top829f0.4[[2]]), function(x) SFS.function(CHR=merge.top829f0.4[[2]][x,]$V1, BEG=merge.top829f0.4[[2]][x,]$V2, END=merge.top829f0.4[[2]][x,]$V3, POP=2)))-> LWK.top.f0.4
LWK.top.f0.4/100->LWK.top.f0.4
unlist(mclapply(1:nrow(merge.top829f0.4[[3]]), function(x) SFS.function(CHR=merge.top829f0.4[[3]][x,]$V1, BEG=merge.top829f0.4[[3]][x,]$V2, END=merge.top829f0.4[[3]][x,]$V3, POP=3)))-> YRI.top.f0.4
YRI.top.f0.4/100->YRI.top.f0.4
unlist(mclapply(1:nrow(merge.top829f0.4[[6]]), function(x) SFS.function(CHR=merge.top829f0.4[[6]][x,]$V1, BEG=merge.top829f0.4[[6]][x,]$V2, END=merge.top829f0.4[[6]][x,]$V3, POP=6)))-> GBR.top.f0.4
GBR.top.f0.4/100->GBR.top.f0.4
unlist(mclapply(1:nrow(merge.top829f0.4[[7]]), function(x) SFS.function(CHR=merge.top829f0.4[[7]][x,]$V1, BEG=merge.top829f0.4[[7]][x,]$V2, END=merge.top829f0.4[[7]][x,]$V3, POP=7)))-> TSI.top.f0.4
TSI.top.f0.4/100->TSI.top.f0.4


unlist(mclapply(1:nrow(merge.CANDf0.4[[2]]), function(x) SFS.function(CHR=merge.CANDf0.4[[2]][x,]$V1, BEG=merge.CANDf0.4[[2]][x,]$V2, END=merge.CANDf0.4[[2]][x,]$V3, POP=2)))-> LWK.cand.f0.4
LWK.cand.f0.4/100->LWK.cand.f0.4
unlist(mclapply(1:nrow(merge.CANDf0.4[[3]]), function(x) SFS.function(CHR=merge.CANDf0.4[[3]][x,]$V1, BEG=merge.CANDf0.4[[3]][x,]$V2, END=merge.CANDf0.4[[3]][x,]$V3, POP=3)))-> YRI.cand.f0.4
YRI.cand.f0.4/100->YRI.cand.f0.4
unlist(mclapply(1:nrow(merge.CANDf0.4[[6]]), function(x) SFS.function(CHR=merge.CANDf0.4[[6]][x,]$V1, BEG=merge.CANDf0.4[[6]][x,]$V2, END=merge.CANDf0.4[[6]][x,]$V3, POP=6)))-> GBR.cand.f0.4
GBR.cand.f0.4/100->GBR.cand.f0.4
unlist(mclapply(1:nrow(merge.CANDf0.4[[7]]), function(x) SFS.function(CHR=merge.CANDf0.4[[7]][x,]$V1, BEG=merge.CANDf0.4[[7]][x,]$V2, END=merge.CANDf0.4[[7]][x,]$V3, POP=7)))-> TSI.cand.f0.4
TSI.cand.f0.4/100->TSI.cand.f0.4

#now neutral
list.SCAN[[2]]-> all.chrs.LWK
setDT(all.chrs.LWK)

arrange(select(all.chrs.LWK, Chr:End.Win, Win.ID), Chr, Beg.Win)-> all.chrs.LWK

setDT(all.chrs.LWK)

with(all.chrs.LWK, paste0('chr',Chr))->all.chrs.LWK$Chr

bedTools.merge(bed1=all.chrs.LWK)-> merge.all.chrs.LWK

system.time(unlist(mclapply(1:nrow(merge.all.chrs.LWK), function(x) try(SFS.function(CHR=merge.all.chrs.LWK$V1[x], BEG=merge.all.chrs.LWK$V2[x], END=merge.all.chrs.LWK$V3[x], POP=2))))-> genomicSFS.LWK) # 445

genomicSFS.LWK[-(which(genomicSFS.LWK=="Error in read.table(out, header = F) : no lines available in input\n"))]-> test

as.numeric(test)-> genomic.SFS.LWK

genomic.SFS.LWK/100-> genomic.SFS.LWK

Store(genomic.SFS.LWK);Store(LWK.top.f0.5, LWK.top.f0.3, LWK.top.f0.4)
Store(YRI.top.f0.5, YRI.top.f0.3, YRI.top.f0.4, GBR.top.f0.5, GBR.top.f0.4, GBR.top.f0.3, TSI.top.f0.5, TSI.top.f0.4, TSI.top.f0.3)
Store(YRI.cand.f0.5, YRI.cand.f0.3, YRI.cand.f0.4, GBR.cand.f0.5, GBR.cand.f0.4, GBR.cand.f0.3, TSI.cand.f0.5, TSI.cand.f0.4, TSI.cand.f0.3, LWK.cand.f0.5, LWK.cand.f0.4, LWK.cand.f0.3)
 

#TO DO: SFS FOR THE OTHER # POPS.

list.SCAN[[3]]-> all.chrs.YRI

setDT(all.chrs.YRI)

arrange(select(all.chrs.YRI, Chr:End.Win, Win.ID), Chr, Beg.Win)-> all.chrs.YRI

setDT(all.chrs.YRI)

with(all.chrs.YRI, paste0('chr',Chr))->all.chrs.YRI$Chr

bedTools.merge(bed1=all.chrs.YRI)-> merge.all.chrs.YRI

system.time(unlist(mclapply(1:nrow(merge.all.chrs.YRI), function(x) try(SFS.function(CHR=merge.all.chrs.YRI$V1[x], BEG=merge.all.chrs.YRI$V2[x], END=merge.all.chrs.YRI$V3[x], POP=2))))-> genomicSFS.YRI) # 624

genomicSFS.YRI[-(which(genomicSFS.YRI=="Error in read.table(out, header = F) : no lines available in input\n"))]-> test

as.numeric(test)-> genomic.SFS.YRI

genomic.SFS.YRI/100-> genomic.SFS.YRI

remove(all.chrs.YRI)

#

list.SCAN[[6]]-> all.chrs.GBR
setDT(all.chrs.GBR)

arrange(select(all.chrs.GBR, Chr:End.Win, Win.ID), Chr, Beg.Win)-> all.chrs.GBR

setDT(all.chrs.GBR)

with(all.chrs.GBR, paste0('chr',Chr))->all.chrs.GBR$Chr

bedTools.merge(bed1=all.chrs.GBR)-> merge.all.chrs.GBR

system.time(unlist(mclapply(1:nrow(merge.all.chrs.GBR), function(x) try(SFS.function(CHR=merge.all.chrs.GBR$V1[x], BEG=merge.all.chrs.GBR$V2[x], END=merge.all.chrs.GBR$V3[x], POP=2))))-> genomicSFS.GBR) # 654

genomicSFS.GBR[-(which(genomicSFS.GBR=="Error in read.table(out, header = F) : no lines available in input\n"))]-> test

as.numeric(test)-> genomic.SFS.GBR

genomic.SFS.GBR/100-> genomic.SFS.GBR

remove(all.chrs.GBR)
#
list.SCAN[[7]]-> all.chrs.TSI
setDT(all.chrs.TSI)

arrange(select(all.chrs.TSI, Chr:End.Win, Win.ID), Chr, Beg.Win)-> all.chrs.TSI

setDT(all.chrs.TSI)

with(all.chrs.TSI, paste0('chr',Chr))->all.chrs.TSI$Chr

bedTools.merge(bed1=all.chrs.TSI)-> merge.all.chrs.TSI

system.time(unlist(mclapply(1:nrow(merge.all.chrs.TSI), function(x) try(SFS.function(CHR=merge.all.chrs.TSI$V1[x], BEG=merge.all.chrs.TSI$V2[x], END=merge.all.chrs.TSI$V3[x], POP=2))))-> genomicSFS.TSI) # 1258.333

genomicSFS.TSI[-(which(genomicSFS.TSI=="Error in read.table(out, header = F) : no lines available in input\n"))]-> test

as.numeric(test)-> genomic.SFS.TSI
genomic.SFS.TSI/100-> genomic.SFS.TSI
#
#actual plots
library(lattice)

pdf('figures/SFS.LWK.pdf')
par(mfrow=c(1,7))

histogram(genomic.SFS.LWK[genomic.SFS.LWK !=0 & genomic.SFS.LWK !=1], col='lightgray', main='Neutral', xlab=list(label='DAF', cex=1.2), ylab=list(label='Relative Frequency (%)',cex=1.2),breaks=seq(from=0,to=1,by=0.05), scales=list(x=list(cex=1.2),y=list(cex=1.2)), auto.key = list(lines=F),ylim=c(0, 50))

histogram(LWK.cand.f0.5[LWK.cand.f0.5 != 0 & LWK.cand.f0.5 !=1], col='cornflowerblue', main='Significant tf=0.5', xlab=list(label='DAF', cex=1.2), ylab=list(label='Relative Frequency (%)',cex=1.2), breaks=seq(from=0,to=1,by=0.05), scales=list(x=list(cex=1.2),y=list(cex=1.2)), auto.key=list(lines=F),ylim=c(0, 25))

histogram(LWK.top.f0.5[LWK.top.f0.5 != 0 & LWK.top.f0.5 !=1], col='cornflowerblue', main='Outlier tf=0.5', xlab=list(label='DAF', cex=1.2), ylab=list(label='Relative Frequency (%)',cex=1.2),breaks=seq(from=0,to=1,by=0.05),scales=list(x=list(cex=1.2),y=list(cex=1.2)), auto.key=list(lines=F),ylim=c(0, 25))

histogram(LWK.cand.f0.4[LWK.cand.f0.4 != 0 & LWK.cand.f0.4 !=1], col='sienna1', main='Signficant tf=0.4', xlab=list(label='DAF', cex=1.2), ylab=list(label='Relative Frequency (%)',cex=1.2), breaks=seq(from=0,to=1,by=0.05),scales=list(x=list(cex=1.2), y=list(cex=1.2)),auto.key=list(lines=F),ylim=c(0, 25))

histogram(LWK.top.f0.4[LWK.top.f0.4 != 0 & LWK.top.f0.4 !=1], col='sienna1', main='Outlier tf=0.4', xlab=list(label='DAF', cex=1.2), ylab=list(label='Relative Frequency (%)',cex=1.2), breaks=seq(from=0,to=1,by=0.05), scales=list(x=list(cex=1.2), y=list(cex=1.2)),auto.key=list(lines=F),ylim=c(0, 25))

histogram(LWK.cand.f0.3[LWK.cand.f0.3 != 0 & LWK.cand.f0.3 !=1], col='violetred1', main='Signficant tf=0.3',  xlab=list(label='DAF', cex=1.2), ylab=list(label='Relative Frequency (%)',cex=1.2),  breaks=seq(from=0,to=1,by=0.05),scales=list(x=list(cex=1.2), y=list(cex=1.2)),auto.key=list(lines=F),ylim=c(0, 25))

histogram(LWK.top.f0.3[LWK.top.f0.3 != 0 & LWK.top.f0.3 !=1], col='violetred1', main='Outlier tf=0.3',  xlab=list(label='DAF', cex=1.2), ylab=list(label='Relative Frequency (%)',cex=1.2), breaks=seq(from=0,to=1,by=0.05), scales=list(x=list(cex=1.2), y=list(cex=1.2)),auto.key=list(lines=F),ylim=c(0, 25))

dev.off()

#
pdf('figures/SFS.YRI.pdf')
par(mfrow=c(1,7))

histogram(genomic.SFS.YRI[genomic.SFS.YRI != 0 & genomic.SFS.YRI !=1], col='lightgray', main='Neutral', xlab=list(label='DAF', cex=1.2), ylab=list(label='Relative Frequency (%)',cex=1.2), breaks=seq(from=0,to=1,by=0.05), scales=list(x=list(cex=1.2), y=list(cex=1.2)),auto.key=list(lines=F),ylim=c(0, 50))

histogram(YRI.cand.f0.5[YRI.cand.f0.5 != 0 & YRI.cand.f0.5 !=1], col='cornflowerblue', main='Significant tf=0.5', xlab=list(label='DAF', cex=1.2), ylab=list(label='Relative Frequency (%)',cex=1.2),  breaks=seq(from=0,to=1,by=0.05),scales=list(x=list(cex=1.2), y=list(cex=1.2)),auto.key=list(lines=F),ylim=c(0, 25))

histogram(YRI.top.f0.5[YRI.top.f0.5 != 0 & YRI.top.f0.5 !=1], col='cornflowerblue', main='Outlier tf=0.5',xlab=list(label='DAF', cex=1.2), ylab=list(label='Relative Frequency (%)',cex=1.2), breaks=seq(from=0,to=1,by=0.05),scales=list(x=list(cex=1.2), y=list(cex=1.2)),auto.key=list(lines=F),ylim=c(0, 25))

histogram(YRI.cand.f0.4[YRI.cand.f0.4 != 0 & YRI.cand.f0.4 !=1], col='sienna1', main='Significant tf=0.4', xlab=list(label='DAF', cex=1.2), ylab=list(label='Relative Frequency (%)',cex=1.2),  breaks=seq(from=0,to=1,by=0.05),scales=list(x=list(cex=1.2), y=list(cex=1.2)),auto.key=list(lines=F),ylim=c(0, 25))

histogram(YRI.top.f0.4[YRI.top.f0.4 != 0 & YRI.top.f0.4 !=1], col='sienna1', main='Outlier tf=0.4', xlab=list(label='DAF', cex=1.2), ylab=list(label='Relative Frequency (%)',cex=1.2), breaks=seq(from=0,to=1,by=0.05),scales=list(x=list(cex=1.2), y=list(cex=1.2)),auto.key=list(lines=F),ylim=c(0, 25))

histogram(YRI.cand.f0.3[YRI.cand.f0.3 != 0 & YRI.cand.f0.3 !=1], col='violetred1', main='Significant tf=0.3', xlab=list(label='DAF', cex=1.2), ylab=list(label='Relative Frequency (%)',cex=1.2),  breaks=seq(from=0,to=1,by=0.05),scales=list(x=list(cex=1.2), y=list(cex=1.2)),auto.key=list(lines=F),ylim=c(0, 25))

histogram(YRI.top.f0.3[YRI.top.f0.3 != 0 & YRI.top.f0.3 !=1], col='violetred1', main='Outlier tf=0.3', xlab=list(label='DAF', cex=1.2), ylab=list(label='Relative Frequency (%)',cex=1.2), breaks=seq(from=0,to=1,by=0.05),scales=list(x=list(cex=1.2), y=list(cex=1.2)),auto.key=list(lines=F),ylim=c(0, 25))

dev.off()

#

pdf('figures/SFS.GBR.pdf')
par(mfrow=c(1,7))
histogram(genomic.SFS.GBR[genomic.SFS.GBR != 0 & genomic.SFS.GBR !=1], col='lightgray', main='Neutral', xlab=list(label='DAF', cex=1.2),ylab=list(label='Relative Frequency (%)',cex=1.2), breaks=seq(from=0,to=1,by=0.05),scales=list(x=list(cex=1.2), y=list(cex=1.2)),auto.key=list(lines=F), ylim=c(0,50)) #GBR

histogram(GBR.cand.f0.5[GBR.cand.f0.5 != 0 & GBR.cand.f0.5 !=1], col='cornflowerblue', main='Significant tf=0.5', xlab=list(label='DAF', cex=1.2), ylab=list(label='Relative Frequency (%)',cex=1.2), breaks=seq(from=0,to=1,by=0.05),scales=list(x=list(cex=1.2), y=list(cex=1.2)),auto.key=list(lines=F), ylim=c(0, 25))

histogram(GBR.top.f0.5[GBR.top.f0.5 != 0 & GBR.top.f0.5 !=1], col='cornflowerblue', main='Outlier tf=0.5', xlab=list(label='DAF', cex=1.2),ylab=list(label='Relative Frequency (%)',cex=1.2), breaks=seq(from=0,to=1,by=0.05),scales=list(x=list(cex=1.2), y=list(cex=1.2)),auto.key=list(lines=F), ylim=c(0, 25))

histogram(GBR.cand.f0.4[GBR.cand.f0.4 != 0 & GBR.cand.f0.4 !=1], col='sienna1', main='Significant tf=0.4', xlab=list(label='DAF', cex=1.2), ylab=list(label='Relative Frequency (%)',cex=1.2), breaks=seq(from=0,to=1,by=0.05),scales=list(x=list(cex=1.2), y=list(cex=1.2)),auto.key=list(lines=F), ylim=c(0,25))

histogram(GBR.top.f0.4[GBR.top.f0.4 != 0 & GBR.top.f0.4 !=1], col='sienna1', main='Outlier tf=0.4', xlab=list(label='DAF', cex=1.2), ylab=list(label='Relative Frequency (%)',cex=1.2),breaks=seq(from=0,to=1,by=0.05),scales=list(x=list(cex=1.2), y=list(cex=1.2)),auto.key=list(lines=F), ylim=c(0,25))

histogram(GBR.cand.f0.3[GBR.cand.f0.3 != 0 & GBR.cand.f0.3 !=1], col='violetred1', main='Significant tf=0.3', xlab=list(label='DAF', cex=1.2),  ylab=list(label='Relative Frequency (%)',cex=1.2),breaks=seq(from=0,to=1,by=0.05),scales=list(x=list(cex=1.2), y=list(cex=1.2)),auto.key=list(lines=F), ylim=c(0,25))

histogram(GBR.top.f0.3[GBR.top.f0.3 != 0 & GBR.top.f0.3 !=1], col='violetred1', main='Outlier tf=0.3', xlab=list(label='DAF', cex=1.2), ylab=list(label='Relative Frequency (%)',cex=1.2),breaks=seq(from=0,to=1,by=0.05),scales=list(x=list(cex=1.2), y=list(cex=1.2)),auto.key=list(lines=F), ylim=c(0,25))

dev.off()



pdf('figures/SFS.TSI.pdf')
par(mfrow=c(1,7))
histogram(genomic.SFS.TSI[genomic.SFS.TSI != 0 & genomic.SFS.TSI !=1], col='lightgray', main='Neutral',xlab=list(label='DAF', cex=1.2), ylab=list(label='Relative Frequency (%)',cex=1.2), breaks=seq(from=0,to=1,by=0.05),scales=list(x=list(cex=1.2), y=list(cex=1.2)),auto.key=list(lines=F),ylim=c(0, 50)) #LWK

histogram(TSI.cand.f0.5[TSI.cand.f0.5 != 0 & TSI.cand.f0.5 !=1], col='cornflowerblue', main='Significant tf=0.5', xlab=list(label='DAF', cex=1.2), ylab=list(label='Relative Frequency (%)',cex=1.2),  breaks=seq(from=0,to=1,by=0.05),scales=list(x=list(cex=1.2), y=list(cex=1.2)),auto.key=list(lines=F),ylim=c(0, 25))

histogram(TSI.top.f0.5[TSI.top.f0.5 != 0 & TSI.top.f0.5 !=1], col='cornflowerblue', main='Outlier tf=0.5', xlab=list(label='DAF', cex=1.2), ylab=list(label='Relative Frequency (%)',cex=1.2), breaks=seq(from=0,to=1,by=0.05),scales=list(x=list(cex=1.2), y=list(cex=1.2)),auto.key=list(lines=F),ylim=c(0, 25))

histogram(TSI.cand.f0.4[TSI.cand.f0.4 != 0 & TSI.cand.f0.4 !=1], col='sienna1', main='Significant tf=0.4', xlab=list(label='DAF', cex=1.2), ylab=list(label='Relative Frequency (%)',cex=1.2),  breaks=seq(from=0,to=1,by=0.05),scales=list(x=list(cex=1.2), y=list(cex=1.2)),auto.key=list(lines=F),ylim=c(0, 25))

histogram(TSI.top.f0.4[TSI.top.f0.4 != 0 & TSI.top.f0.4 !=1], col='sienna1', main='Outlier tf=0.4', xlab=list(label='DAF', cex=1.2), ylab=list(label='Relative Frequency (%)',cex=1.2), breaks=seq(from=0,to=1,by=0.05),scales=list(x=list(cex=1.2), y=list(cex=1.2)),auto.key=list(lines=F),ylim=c(0, 25))

histogram(TSI.cand.f0.3[TSI.cand.f0.3 != 0 & TSI.cand.f0.3 !=1], col='violetred1', main='Significant tf=0.3', xlab=list(label='DAF', cex=1.2), ylab=list(label='Relative Frequency (%)',cex=1.2),  breaks=seq(from=0,to=1,by=0.05),scales=list(x=list(cex=1.2), y=list(cex=1.2)),auto.key=list(lines=F),ylim=c(0, 25))

histogram(TSI.top.f0.3[TSI.top.f0.3 != 0 & TSI.top.f0.3 !=1], col='violetred1', main='Outlier tf=0.3', xlab=list(label='DAF', cex=1.2), ylab=list(label='Relative Frequency (%)',cex=1.2), breaks=seq(from=0,to=1,by=0.05),scales=list(x=list(cex=1.2), y=list(cex=1.2)),auto.key=list(lines=F),ylim=c(0, 25))

dev.off()


Store(genomic.SFS, genomic.SFS.GBR, genomic.SFS.TSI, genomic.SFS.YRI)

##############
#write bed files

BED.PATH<-'/mnt/sequencedb/PopGen/barbara/NCV_dir_package/read_scan_data/bedfiles/'


write.table(select(Union.CANDf0.5_0.4_0.3[[3]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,"Union.CANDf0.5_0.4_0.3_YRI.bed"), quote=F, sep="\t", col.names=F, row.names=F) 

write.table(select(Union.CANDf0.5_0.4_0.3[[2]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,"Union.CANDf0.5_0.4_0.3_LWK.bed"), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(Union.CANDf0.5_0.4_0.3[[6]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,"Union.CANDf0.5_0.4_0.3_GBR.bed"), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(Union.CANDf0.5_0.4_0.3[[7]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,"Union.CANDf0.5_0.4_0.3_TSI.bed"), quote=F, sep="\t", col.names=F, row.names=F)




write.table(select(Union.top0.5_0.4_0.3[[3]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,"Union.top829.0.5_0.4_0.3_YRI.bed"), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(Union.top0.5_0.4_0.3[[2]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,"Union.top829.0.5_0.4_0.3_LWK.bed"), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(Union.top0.5_0.4_0.3[[6]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,"Union.top829.0.5_0.4_0.3_GBR.bed"), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(Union.top0.5_0.4_0.3[[7]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1),  file=paste0(BED.PATH,"Union.top829.0.5_0.4_0.3_TSI.bed"), quote=F, sep="\t", col.names=F, row.names=F)


write.table(select(top829f0.5[[2]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'top829f0.5.', pops[2], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(top829f0.5[[3]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'top829f0.5.', pops[3], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(top829f0.5[[6]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'top829f0.5.', pops[6], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(top829f0.5[[7]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'top829f0.5.', pops[7], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)


write.table(select(top829f0.4[[2]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'top829f0.4.', pops[2], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(top829f0.4[[3]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'top829f0.4.', pops[3], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(top829f0.4[[6]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'top829f0.4.', pops[6], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(top829f0.4[[7]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'top829f0.4.', pops[7], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)


write.table(select(top829f0.3[[2]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'top829f0.3.', pops[2], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(top829f0.3[[3]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'top829f0.3.', pops[3], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(top829f0.3[[6]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'top829f0.3.', pops[6], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(top829f0.3[[7]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'top829f0.3.', pops[7], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

#

write.table(select(CANDf0.5[[2]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'CANDf0.5.', pops[2], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(CANDf0.5[[3]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'CANDf0.5.', pops[3], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(CANDf0.5[[6]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'CANDf0.5.', pops[6], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(CANDf0.5[[7]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'CANDf0.5.', pops[7], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)


write.table(select(CANDf0.4[[2]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'CANDf0.4.', pops[2], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(CANDf0.4[[3]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'CANDf0.4.', pops[3], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(CANDf0.4[[6]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'CANDf0.4.', pops[6], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(CANDf0.4[[7]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'CANDf0.4.', pops[7], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)


write.table(select(CANDf0.3[[2]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'CANDf0.3.', pops[2], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(CANDf0.3[[3]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'CANDf0.3.', pops[3], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(CANDf0.3[[6]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'CANDf0.3.', pops[6], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

write.table(select(CANDf0.3[[7]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'CANDf0.3.', pops[7], '.bed'), quote=F, sep="\t", col.names=F, row.names=F)

#background windows

write.table(select(list.SCAN[[2]], Chr, Beg.Win, End.Win, Win.ID), options(scipen=1), file=paste0(BED.PATH,'background_windows.bed'), quote=F, sep="\t", col.names=F, row.names=F)


####Scanned windows annotation ####


#sort ensembl bed file

read.table('/mnt/sequencedb/PopGen/barbara/scan_may_2014/final_encode.bed')-> bed2

gsub("chr", "", bed2$V1)-> bed2$V1

arrange(bed2, V1, V2, V3)-> bed2

paste0("chr", bed2$V1)-> bed2$V1

#background windows

gsub("chr", "", list.SCAN[[2]]$Chr)-> list.SCAN[[2]]$Chr

arrange(list.SCAN[[2]], Chr, Beg.Win, End.Win)-> bed1

paste0("chr", bed1$Chr)-> bed1$Chr

select(bed1, Chr, Beg.Win, End.Win, Win.ID)-> bed1

setDT(bed1)


#this session below needs refinement...i mention proportion fo windows overlapping something, but these 'windows' are actually merged windows from the merge and then intersect sessions...need to find a way to make this clearer.


bedTools.merge(bed1=bed1)-> merge.scanned.windows

setDT(merge.scanned.windows)

#with(merge.scanned.windows, cbind(merge.scanned.windows, Win.ID=paste0(V1,"|", V2, "|", V3)))-> merge.scanned.windows

length(unique(sort(merge.scanned.windows$V4)))  #14,195 merged windows in the scan

system.time(bedTools.2in(bed1=merge.scanned.windows, bed2=bed2)-> intersect.scanned.windows) #1699.756 

setDT(intersect.scanned.windows)

Store(intersect.scanned.windows)

#stopped here 27.10
length(unique(sort(intersect.scanned.windows$V8))) #48,254 number of 'coding elements' scanned

length(unique(sort(filter(intersect.scanned.windows, V11=="protein_coding")$V8))) #18,633 genes scanned

length(unique(intersect.scanned.windows$V4)) #number of merged windows overlapping coding elements 1,670, i,em 11670/14195=82% of the merged windows

length(unique(sort(filter(intersect.scanned.windows, V11=="protein_coding")$V4)))  #8514 number of merged windows scanning a prot. coding gene, i.e, 8514/14195=60%

length(unique(filter(intersect.scanned.windows, V11=="protein_coding")$V4))/length(unique(merge.scanned.windows$V4)) # 60%  of (merged) background windows overlap protein_coding genes

mclapply(c(2,3,6,7), function(x) length(unique(filter(intersect.top829f0.5[[x]], V11=="protein_coding")$V4))/length(unique(merge.top829f0.5[[x]]$V4))) #for the top windows, around 37% overap of merged windows overlap prot coding genes

mclapply(c(2,3,6,7), function(x) length(unique(filter(intersect.CANDf0.5[[x]], V11=="protein_coding")$V4))/length(unique(intersect.CANDf0.5[[x]]$V4))) #and around 75% of the candidte windows: why this difference??


#End
