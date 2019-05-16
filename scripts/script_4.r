################################################################################
#       Barbara D Bitarello
#
#       Last modified: 17.10.2016
#
#       Calculate normalized NCD, etc.
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


###
#join results from test and test2. finally, collapse bins >250.

#MAKE A TABLE LIKE THE ONE ABOVE, BUT WITH BINS 1:18, AND THEN (19,20), (21,22)...250 AND THEN 250+

Objects()

setwd('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/read_scan_data/')

X<-AFRICA[[2]]

#3 vectors for the bins
bin.vec1<-seq(from=10, to=229) #1        #all windows with at least 10 IS
bin.vec2<-230
nsims<-10000
system.time(mclapply(bin.vec1, function(x) subset(X, Nr.IS==x))->list.bin.vec1) #266
system.time(mclapply(bin.vec2, function(x) subset(X, Nr.IS>=x))-> list.bin.vec2)
system.time(mclapply(list.bin.vec1, function(x) x[sample(seq(1:dim(x)[1]), nsims),])-> l.bin.vec1)
system.time(mclapply(bin.vec2, function(x) subset(X, Nr.IS>=x))-> list.bin.vec2)
system.time(mclapply(list.bin.vec2, function(x) x[sample(seq(1:dim(x)[1]),nsims),])-> l.bin.vec2)

#
mclapply(All.Res2, function(x) setDT(cbind(x, P.val.NCDf0.5=rep(NA, dim(x)[1]), P.val.NCDf0.4=rep(NA, dim(x)[1]), P.val.NCDf0.3=rep(NA, dim(x)[1]), P.val.NCDf0.2=rep(NA, dim(x)[1]), P.val.NCDf0.1=rep(NA, dim(x)[1]), Dist.NCD.f0.5=rep(NA, dim(x)[1]), Dist.NCD.f0.4=rep(NA, dim(x)[1]),  Dist.NCD.f0.3=rep(NA, dim(x)[1]),  Dist.NCD.f0.2=rep(NA, dim(x)[1]),  Dist.NCD.f0.1=rep(NA, dim(x)[1]))))-> All.Res3

gc()
YRI.2<-All.Res3[[3]]


nsims<-10000

mclapply(bin.vec1, function(x) (which(YRI.2$Nr.IS==x)))->temp.YRI
mclapply(temp.YRI, function(x) length(x))-> temp3.YRI
length(which(YRI.2$Nr.IS>=bin.vec2[[1]]))->temp4.YRI

which(YRI.2$Nr.IS>=bin.vec2[[1]])->temp2.YRI

test<-vector('list' ,length(l.bin.vec1))
for ( i in 1:length(l.bin.vec1)){
if(temp3.YRI[[i]]<=nsims){l.bin.vec1[[i]][sample(seq(1:nsims), temp3.YRI[[i]]),]->test[[i]]}
if(temp3.YRI[[i]]>nsims){l.bin.vec1[[i]]->test[[i]]}} #now do vioplots with this and the real data.

pdf('figures/oct_2016.vioplots.YRI.subsampled.pdf')
par(mfrow=c(4,4))
sapply(1:16, function(x) vioplot(test[[x]]$ncvFD_f0.5, YRI.2[temp.YRI[[x]],]$NCDf5, col='cornflowerblue', border='gray', rectCol=' white', colMed='black',names=c('sims', 'data')))
dev.off()

Store(All.Res3);gc();remove(All.Res2);gc()
#
system.time(for (i in 1: length(temp.YRI)){
I<-temp.YRI[[i]]
unlist(mclapply(I, function(x) (sum(YRI.2$NCDf5[x]>=l.bin.vec1[[i]]$ncvFD_f0.5)/nsims)))->YRI.2$P.val.NCDf0.5[I]
unlist(mclapply(I, function(x) (sum(YRI.2$NCDf4[x]>=l.bin.vec1[[i]]$ncvFD_f0.4)/nsims)))->YRI.2$P.val.NCDf0.4[I]
unlist(mclapply(I, function(x) (sum(YRI.2$NCDf3[x]>=l.bin.vec1[[i]]$ncvFD_f0.3)/nsims)))->YRI.2$P.val.NCDf0.3[I]
unlist(mclapply(I, function(x) (sum(YRI.2$NCDf2[x]>=l.bin.vec1[[i]]$ncvFD_f0.2)/nsims)))->YRI.2$P.val.NCDf0.2[I]

unlist(mclapply(I, function(x) (sum(YRI.2$NCDf1[x]>=l.bin.vec1[[i]]$ncvFD_f0.1)/nsims)))->YRI.2$P.val.NCDf0.1[I]
cat(i, 'done\n')
}) #10269.308

#stopped here 07.10.2016

unlist(mclapply(temp2.YRI, function(x) (sum(YRI.2$NCDf5[x]>=l.bin.vec2[[1]]$ncvFD_f0.5)/nsims)))->YRI.2$P.val.NCDf0.5[temp2.YRI]
unlist(mclapply(temp2.YRI, function(x) (sum(YRI.2$NCDf4[x]>=l.bin.vec2[[1]]$ncvFD_f0.4)/nsims)))->YRI.2$P.val.NCDf0.4[temp2.YRI]
unlist(mclapply(temp2.YRI, function(x) (sum(YRI.2$NCDf3[x]>=l.bin.vec2[[1]]$ncvFD_f0.3)/nsims)))->YRI.2$P.val.NCDf0.3[temp2.YRI]
unlist(mclapply(temp2.YRI, function(x) (sum(YRI.2$NCDf2[x]>=l.bin.vec2[[1]]$ncvFD_f0.2)/nsims)))->YRI.2$P.val.NCDf0.2[temp2.YRI]
unlist(mclapply(temp2.YRI, function(x) (sum(YRI.2$NCDf1[x]>=l.bin.vec2[[1]]$ncvFD_f0.1)/nsims)))->YRI.2$P.val.NCDf0.1[temp2.YRI]

pdf('figures/october.2016.vioplots.YRI.pdf')
par(mfrow=c(4,4))
sapply(1:16, function(x) vioplot(l.bin.vec1[[x]]$ncvFD_f0.5, YRI.2[temp.YRI[[x]],]$NCDf5, col='cornflowerblue', border='gray', rectCol=' white', colMed='black',names=c('sims', 'data')))
dev.off()


pdf('figures/october.2016.YRI.vs.neutral.sims.p0.5.pdf')
unlist(mclapply(l.bin.vec1, function(x) quantile(x$ncvFD_f0.5, prob=0.5)))-> sim
unlist(mclapply(temp.YRI, function(x) quantile(YRI.2[x,]$NCVf5, prob=0.5)))-> data
plot(data~sim, col='cornflowerblue', ylim=c(min(sort(c(data,sim))), 0.5), xlim=c(min(sort(c(data,sim))), 0.5), type= 'n')
points(data[1:10]~sim[1:10], col='cornflowerblue', ylim=c(min(sort(c(data,sim))), 0.5), xlim=c(min(sort(c(data,sim))), 0.5))
points(data[11:20]~sim[11:20], col='blue', ylim=c(min(sort(c(data,sim))), 0.5), xlim=c(min(sort(c(data,sim))), 0.5))
points(data[21:30]~sim[21:30], col='purple', ylim=c(min(sort(c(data,sim))), 0.5), xlim=c(min(sort(c(data,sim))), 0.5))
points(data[31:40]~sim[31:40], col='magenta', ylim=c(min(sort(c(data,sim))), 0.5), xlim=c(min(sort(c(data,sim))), 0.5))
points(data[41:226]~sim[41:226], col='violetred1')
lines(y=seq(from=0.43, to=0.45, by=0.01), x=seq(from=0.43, to=0.45, by=0.01),lty=2, col= 'red')

legend('topleft',c('4:13 I.S', '14:23 I.S','24:33 I.S', '33:43 I.S', '43+ I.S'), col=c('cornflowerblue','blue','purple','magenta','violetred1'), pch=19)
dev.off()

pdf('figures/october.2016.vioplots.4_54.IS.YRI.pdf')
par(mfrow=c(2,3))
vioplot(l.bin.vec1[[1]]$ncvFD_f0.5, YRI.2[temp.YRI[[1]],]$NCDf5)
vioplot(l.bin.vec1[[11]]$ncvFD_f0.5, YRI.2[temp.YRI[[11]],]$NCDf5)
vioplot(l.bin.vec1[[21]]$ncvFD_f0.5, YRI.2[temp.YRI[[21]],]$NCDf5)
vioplot(l.bin.vec1[[31]]$ncvFD_f0.5, YRI.2[temp.YRI[[31]],]$NCDf5)
vioplot(l.bin.vec1[[41]]$ncvFD_f0.5, YRI.2[temp.YRI[[41]],]$NCDf5)
vioplot(l.bin.vec1[[51]]$ncvFD_f0.5, YRI.2[temp.YRI[[51]],]$NCDf5)
dev.off()

objName<-'YRI.2'

save(list=objName, file= 'YRI.2.RData')

Store(YRI.2)
############
#now for LWK
Objects()

LWK.2<-All.Res3[[2]]

mclapply(bin.vec1, function(x) (which(LWK.2$Nr.IS==x)))->temp.LWK

which(LWK.2$Nr.IS>=bin.vec2[[1]])->temp2.LWK

system.time(for(i in 1:length(temp.LWK)){
I<-temp.LWK[[i]]
unlist(mclapply(I, function(x) (sum(LWK.2$NCDf5[x]>=l.bin.vec1[[i]]$ncvFD_f0.5)/nsims)))->LWK.2$P.val.NCDf0.5[I]
unlist(mclapply(I, function(x) (sum(LWK.2$NCDf4[x]>=l.bin.vec1[[i]]$ncvFD_f0.4)/nsims)))->LWK.2$P.val.NCDf0.4[I]
unlist(mclapply(I, function(x) (sum(LWK.2$NCDf3[x]>=l.bin.vec1[[i]]$ncvFD_f0.3)/nsims)))->LWK.2$P.val.NCDf0.3[I]
unlist(mclapply(I, function(x) (sum(LWK.2$NCDf2[x]>=l.bin.vec1[[i]]$ncvFD_f0.2)/nsims)))->LWK.2$P.val.NCDf0.2[I]
unlist(mclapply(I, function(x) (sum(LWK.2$NCDf1[x]>=l.bin.vec1[[i]]$ncvFD_f0.1)/nsims)))->LWK.2$P.val.NCDf0.1[I]
cat(i, 'done\n') #10450.778
})

unlist(mclapply(temp2.LWK, function(x) (sum(LWK.2$NCDf5[x]>=l.bin.vec2[[1]]$ncvFD_f0.5)/nsims)))->LWK.2$P.val.NCDf0.5[temp2.LWK]
unlist(mclapply(temp2.LWK, function(x) (sum(LWK.2$NCDf4[x]>=l.bin.vec2[[1]]$ncvFD_f0.4)/nsims)))->LWK.2$P.val.NCDf0.4[temp2.LWK]
unlist(mclapply(temp2.LWK, function(x) (sum(LWK.2$NCDf3[x]>=l.bin.vec2[[1]]$ncvFD_f0.3)/nsims)))->LWK.2$P.val.NCDf0.3[temp2.LWK]
unlist(mclapply(temp2.LWK, function(x) (sum(LWK.2$NCDf2[x]>=l.bin.vec2[[1]]$ncvFD_f0.2)/nsims)))->LWK.2$P.val.NCDf0.2[temp2.LWK]
unlist(mclapply(temp2.LWK, function(x) (sum(LWK.2$NCDf1[x]>=l.bin.vec2[[1]]$ncvFD_f0.1)/nsims)))->LWK.2$P.val.NCDf0.1[temp2.LWK]

pdf('figures/october.2016.vioplots.LWK.pdf')
par(mfrow=c(4,4))
sapply(1:16, function(x) vioplot(l.bin.vec1[[x]]$ncvFD_f0.5, LWK.2[temp.LWK[[x]],]$NCDf5, col='cornflowerblue', border='gray', rectCol=' white', colMed='black',names=c('sims', 'data')))
dev.off()


pdf('figures/october.2016.AFRICA.NrIS.NCV.neutral.p0.01.pdf')
plot(c(seq(from=0.12, to=0.45, by=0.01), rep(0.3,196))~seq(1:230), type='n', ylab='NCD', xlab='Number of Informative Sites')
points(c(unlist(lapply(l.bin.vec1, function(x) quantile(x$ncvFD_f0.5, prob=0.01))),quantile(l.bin.vec2[[1]]$ncvFD_f0.5, prob=0.01))~c(bin.vec1,bin.vec2), col='cornflowerblue', pch=20)
lines(c(unlist(lapply(l.bin.vec1, function(x) quantile(x$ncvFD_f0.5, prob=0.01))),quantile(l.bin.vec2[[1]]$ncvFD_f0.5, prob=0.01))~c(bin.vec1,bin.vec2), col= 'lightgray', cex=0.2)
points(c(unlist(lapply(l.bin.vec1, function(x) quantile(x$ncvFD_f0.4, prob=0.01))),quantile(l.bin.vec2[[1]]$ncvFD_f0.4, prob=0.01))~c(bin.vec1,bin.vec2), col='sienna1', pch=20)
lines(c(unlist(lapply(l.bin.vec1, function(x) quantile(x$ncvFD_f0.4, prob=0.01))),quantile(l.bin.vec2[[1]]$ncvFD_f0.4, prob=0.01))~c(bin.vec1,bin.vec2), col='lightgray', cex=0.2)
points(c(unlist(lapply(l.bin.vec1, function(x) quantile(x$ncvFD_f0.3, prob=0.01))),quantile(l.bin.vec2[[1]]$ncvFD_f0.3, prob=0.01))~c(bin.vec1,bin.vec2), col='violetred1', pch=20)
lines(c(unlist(lapply(l.bin.vec1, function(x) quantile(x$ncvFD_f0.3, prob=0.01))),quantile(l.bin.vec2[[1]]$ncvFD_f0.3, prob=0.01))~c(bin.vec1,bin.vec2), col='lightgray', cex=0.2)
#points(c(unlist(lapply(l.bin.vec1, function(x) quantile(x$ncvFD_f0.2, prob=0.01))),quantile(l.bin.vec2[[1]]$ncvFD_f0.2, prob=0.01))~c(bin.vec1,bin.vec2), col='darkolivegreen', pch=20)
#lines(c(unlist(lapply(l.bin.vec1, function(x) quantile(x$ncvFD_f0.2, prob=0.01))),quantile(l.bin.vec2[[1]]$ncvFD_f0.2, prob=0.01))~c(bin.vec1,bin.vec2), col='lightgray', cex=0.2)
#legend('bottomright', c('tf=feq=0.5', 'tf=feq=0.4','tf=feq=0.3', 'tf=feq=0.2'), col=c('cornflowerblue', 'sienna1', 'violetred1', 'darkolivegreen'), pch=20, bty='n')
legend('bottomright', c('tf=feq=0.5', 'tf=feq=0.4','tf=feq=0.3'), col=c('cornflowerblue', 'sienna1', 'violetred1'), pch=20, bty='n')
dev.off()


pdf('figures/october.2016.LWK.vs.neutral.sims.p0.5.pdf')
unlist(mclapply(l.bin.vec1, function(x) quantile(x$ncvFD_f0.5, prob=0.5)))-> sim
unlist(mclapply(temp.LWK, function(x) quantile(LWK.2[x,]$NCDf5, prob=0.5)))-> data
plot(data~sim, col='cornflowerblue', ylim=c(min(sort(c(data,sim))), 0.5), xlim=c(min(sort(c(data,sim))), 0.5), type= 'n')
points(data[1:10]~sim[1:10], col='cornflowerblue', ylim=c(min(sort(c(data,sim))), 0.5), xlim=c(min(sort(c(data,sim))), 0.5))
points(data[11:20]~sim[11:20], col='blue', ylim=c(min(sort(c(data,sim))), 0.5), xlim=c(min(sort(c(data,sim))), 0.5))
points(data[21:30]~sim[21:30], col='purple', ylim=c(min(sort(c(data,sim))), 0.5), xlim=c(min(sort(c(data,sim))), 0.5))
points(data[31:40]~sim[31:40], col='magenta', ylim=c(min(sort(c(data,sim))), 0.5), xlim=c(min(sort(c(data,sim))), 0.5))
points(data[41:211]~sim[41:211], col='violetred1')
lines(y=seq(from=0.43, to=0.45, by=0.01), x=seq(from=0.43, to=0.45, by=0.01),lty=2, col= 'red')
legend('topleft',c('19:28 I.S', '29:38 I.S','39:48 I.S', '49:58 I.S', '58+ I.S'), col=c('cornflowerblue','blue','purple','magenta','violetred1'), pch=19)

dev.off()

pdf('figures/october.2016.vioplots.4_54.IS.LWK.pdf')
par(mfrow=c(2,3))
vioplot(l.bin.vec1[[1]]$ncvFD_f0.5, LWK.2[temp.LWK[[1]],]$NCDf5)
vioplot(l.bin.vec1[[11]]$ncvFD_f0.5, LWK.2[temp.LWK[[11]],]$NCDf5)
vioplot(l.bin.vec1[[21]]$ncvFD_f0.5, LWK.2[temp.LWK[[21]],]$NCDf5)
vioplot(l.bin.vec1[[31]]$ncvFD_f0.5, LWK.2[temp.LWK[[31]],]$NCDf5)
vioplot(l.bin.vec1[[41]]$ncvFD_f0.5, LWK.2[temp.LWK[[41]],]$NCDf5)
vioplot(l.bin.vec1[[51]]$ncvFD_f0.5, LWK.2[temp.LWK[[51]],]$NCDf5)
dev.off()

objName<-'LWK.2'

save(list=objName, file= 'LWK.2.RData')

Store(LWK.2)
#############
#############

AWS.2<-All.Res3[[1]]

mclapply(bin.vec1, function(x) (which(AWS.2$Nr.IS==x)))->temp.AWS

which(AWS.2$Nr.IS>=bin.vec2[[1]])->temp2.AWS

system.time(for(i in 1:length(temp.AWS)) {
I<-temp.AWS[[i]]
unlist(mclapply(I, function(x) (sum(AWS.2$NCDf5[x]>=l.bin.vec1[[i]]$ncvFD_f0.5)/nsims)))->AWS.2$P.val.NCDf0.5[I]
unlist(mclapply(I, function(x) (sum(AWS.2$NCDf4[x]>=l.bin.vec1[[i]]$ncvFD_f0.4)/nsims)))->AWS.2$P.val.NCDf0.4[I]
unlist(mclapply(I, function(x) (sum(AWS.2$NCDf3[x]>=l.bin.vec1[[i]]$ncvFD_f0.3)/nsims)))->AWS.2$P.val.NCDf0.3[I]
unlist(mclapply(I, function(x) (sum(AWS.2$NCDf2[x]>=l.bin.vec1[[i]]$ncvFD_f0.2)/nsims)))->AWS.2$P.val.NCDf0.2[I]
unlist(mclapply(I, function(x) (sum(AWS.2$NCDf1[x]>=l.bin.vec1[[i]]$ncvFD_f0.1)/nsims)))->AWS.2$P.val.NCDf0.1[I]
cat(i, 'done\n')
})

unlist(mclapply(temp2.AWS, function(x) (sum(AWS.2$NCDf5[x]>=l.bin.vec2[[1]]$ncvFD_f0.5)/nsims)))->AWS.2$P.val.NCDf0.5[temp2.AWS]
unlist(mclapply(temp2.AWS, function(x) (sum(AWS.2$NCDf4[x]>=l.bin.vec2[[1]]$ncvFD_f0.4)/nsims)))->AWS.2$P.val.NCDf0.4[temp2.AWS]
unlist(mclapply(temp2.AWS, function(x) (sum(AWS.2$NCDf3[x]>=l.bin.vec2[[1]]$ncvFD_f0.3)/nsims)))->AWS.2$P.val.NCDf0.3[temp2.AWS]
unlist(mclapply(temp2.AWS, function(x) (sum(AWS.2$NCDf2[x]>=l.bin.vec2[[1]]$ncvFD_f0.2)/nsims)))->AWS.2$P.val.NCDf0.2[temp2.AWS]
unlist(mclapply(temp2.AWS, function(x) (sum(AWS.2$NCDf1[x]>=l.bin.vec2[[1]]$ncvFD_f0.1)/nsims)))->AWS.2$P.val.NCDf0.1[temp2.AWS]

pdf('/mnt/sequencedb/PopGen/barbara/scan_may_2014/figures/october.2016.vioplots.AWS.pdf')
par(mfrow=c(4,4))
sapply(1:16, function(x) vioplot(l.bin.vec1[[x]]$ncvFD_f0.5, AWS.2[temp.AWS[[x]],]$NCDf5, col='cornflowerblue', border='gray', rectCol=' white', colMed='black',names=c('sims', 'data')))
dev.off()

objName<-'AWS.2'

save(list=objName, file= 'AWS.2.RData')

Store(AWS.2)
########################################################################################################################################
########################################################################################################################################
X<-EUROPE[[2]]


bin.vec1.eu<-seq(from=10, to=207) #1        by 1 bins
bin.vec2.eu<-208
nsims<-10000

system.time(mclapply(bin.vec1.eu, function(x) subset(X, Nr.IS==x))->list.bin.vec1.eu)
system.time(mclapply(list.bin.vec1.eu, function(x) x[sample(seq(1:dim(x)[1]), nsims),])-> l.bin.vec1.eu)
system.time(mclapply(bin.vec2.eu, function(x) subset(X, Nr.IS>=x))-> list.bin.vec2.eu)
system.time(mclapply(list.bin.vec2.eu, function(x) x[sample(seq(1:dim(x)[1]), nsims),])-> l.bin.vec2.eu)

CEU.2<-All.Res3[[4]]
mclapply(bin.vec1.eu, function(x) (which(CEU.2$Nr.IS==x)))->temp.CEU
which(CEU.2$Nr.IS>=bin.vec2.eu[[1]])->temp2.CEU

pdf('figures/october.2016.vioplots.CEU.pdf')
par(mfrow=c(4,4))
sapply(1:16, function(x) vioplot(l.bin.vec1.eu[[x]]$ncvFD_f0.5, CEU.2[temp.CEU[[x]],]$NCDf5, col='cornflowerblue', border='gray', rectCol=' white', colMed='black',names=c('sims', 'data')))
dev.off()

system.time(for (i in 1:length(temp.CEU)) {
I<-temp.CEU[[i]]
unlist(mclapply(I, function(x) (sum(CEU.2$NCDf5[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.5)/nsims)))->CEU.2$P.val.NCDf0.5[I]
unlist(mclapply(I, function(x) (sum(CEU.2$NCDf4[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.4)/nsims)))->CEU.2$P.val.NCDf0.4[I]
unlist(mclapply(I, function(x) (sum(CEU.2$NCDf3[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.3)/nsims)))->CEU.2$P.val.NCDf0.3[I]
unlist(mclapply(I, function(x) (sum(CEU.2$NCDf2[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.2)/nsims)))->CEU.2$P.val.NCDf0.2[I]
unlist(mclapply(I, function(x) (sum(CEU.2$NCDf1[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.1)/nsims)))->CEU.2$P.val.NCDf0.1[I]
}) #1955.988
                         
unlist(mclapply(temp2.CEU, function(x) (sum(CEU.2$NCDf5[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.5)/nsims)))->CEU.2$P.val.NCDf0.5[temp2.CEU]
unlist(mclapply(temp2.CEU, function(x) (sum(CEU.2$NCDf4[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.4)/nsims)))->CEU.2$P.val.NCDf0.4[temp2.CEU]
unlist(mclapply(temp2.CEU, function(x) (sum(CEU.2$NCDf3[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.3)/nsims)))->CEU.2$P.val.NCDf0.3[temp2.CEU]
unlist(mclapply(temp2.CEU, function(x) (sum(CEU.2$NCDf2[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.2)/nsims)))->CEU.2$P.val.NCDf0.2[temp2.CEU]
unlist(mclapply(temp2.CEU, function(x) (sum(CEU.2$NCDf1[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.1)/nsims)))->CEU.2$P.val.NCDf0.1[temp2.CEU]



gc()

pdf('figures/october.2016.vioplots.4_54.IS.CEU.pdf')
par(mfrow=c(2,3))
vioplot(l.bin.vec1.eu[[1]]$ncvFD_f0.5, CEU.2[temp.CEU[[1]],]$NCDf5,names=c('sims', 'data'))
vioplot(l.bin.vec1.eu[[11]]$ncvFD_f0.5, CEU.2[temp.CEU[[11]],]$NCDf5,names=c('sims', 'data'))
vioplot(l.bin.vec1.eu[[21]]$ncvFD_f0.5, CEU.2[temp.CEU[[21]],]$NCDf5,names=c('sims', 'data'))
vioplot(l.bin.vec1.eu[[31]]$ncvFD_f0.5, CEU.2[temp.CEU[[31]],]$NCDf5,names=c('sims', 'data'))
vioplot(l.bin.vec1.eu[[41]]$ncvFD_f0.5, CEU.2[temp.CEU[[41]],]$NCDf5,names=c('sims', 'data'))
vioplot(l.bin.vec1.eu[[51]]$ncvFD_f0.5, CEU.2[temp.CEU[[51]],]$NCDf5,names=c('sims', 'data'))
dev.off()


pdf('figures/october.2016.vioplots.CEU.pdf')
par(mfrow=c(4,4))
sapply(1:16, function(x) vioplot(l.bin.vec1.eu[[x]]$ncvFD_f0.5, CEU.2[temp.CEU[[x]],]$NCDf5, col='cornflowerblue', border='gray', rectCol=' white', colMed='black',names=c('sims', 'data')))
dev.off()

pdf('figures/october.2016.EUROPE.NrIS.NCV.neutral.p0.01.pdf')
plot(c(seq(from=0.12, to=0.45, by=0.01), rep(0.3,196))~seq(1:230), type='n', ylab='NCV', xlab='Number of Informative Sites')
points(c(unlist(lapply(l.bin.vec1.eu, function(x) quantile(x$ncvFD_f0.5, prob=0.01))),quantile(l.bin.vec2.eu[[1]]$ncvFD_f0.5, prob=0.01))~c(bin.vec1.eu,bin.vec2.eu), col='cornflowerblue', pch=20)
lines(c(unlist(lapply(l.bin.vec1.eu, function(x) quantile(x$ncvFD_f0.5, prob=0.01))),quantile(l.bin.vec2.eu[[1]]$ncvFD_f0.5, prob=0.01))~c(bin.vec1.eu,bin.vec2.eu), col= 'lightgray', cex=0.2)
points(c(unlist(lapply(l.bin.vec1.eu, function(x) quantile(x$ncvFD_f0.4, prob=0.01))),quantile(l.bin.vec2.eu[[1]]$ncvFD_f0.4, prob=0.01))~c(bin.vec1.eu,bin.vec2.eu), col='sienna1', pch=20)
lines(c(unlist(lapply(l.bin.vec1.eu, function(x) quantile(x$ncvFD_f0.4, prob=0.01))),quantile(l.bin.vec2.eu[[1]]$ncvFD_f0.4, prob=0.01))~c(bin.vec1.eu,bin.vec2.eu), col='lightgray', cex=0.2)
points(c(unlist(lapply(l.bin.vec1.eu, function(x) quantile(x$ncvFD_f0.3, prob=0.01))),quantile(l.bin.vec2.eu[[1]]$ncvFD_f0.3, prob=0.01))~c(bin.vec1.eu,bin.vec2.eu), col='violetred1', pch=20)
lines(c(unlist(lapply(l.bin.vec1.eu, function(x) quantile(x$ncvFD_f0.3, prob=0.01))),quantile(l.bin.vec2.eu[[1]]$ncvFD_f0.3, prob=0.01))~c(bin.vec1.eu,bin.vec2.eu), col='lightgray', cex=0.2)
#points(c(unlist(lapply(l.bin.vec1.eu, function(x) quantile(x$ncvFD_f0.2, prob=0.01))),quantile(l.bin.vec2.eu[[1]]$ncvFD_f0.2, prob=0.01))~c(bin.vec1.eu,bin.vec2.eu), col='darkolivegreen', pch=20)
#lines(c(unlist(lapply(l.bin.vec1.eu, function(x) quantile(x$ncvFD_f0.2, prob=0.01))),quantile(l.bin.vec2.eu[[1]]$ncvFD_f0.2, prob=0.01))~c(bin.vec1.eu,bin.vec2.eu), col='lightgray', cex=0.2)
#legend('bottomright', c('feq=0.5', 'feq=0.4','feq=0.3', 'feq=0.2'), col=c('cornflowerblue', 'sienna1', 'violetred1', 'darkolivegreen'), pch=20, bty='n')
legend('bottomright', c('feq=0.5', 'feq=0.4','feq=0.3'), col=c('cornflowerblue', 'sienna1', 'violetred1'), pch=20, bty='n')
dev.off()

objName<-'CEU.2'

save(list=objName, file= 'CEU.2.RData')

Store(CEU.2)
################
#now for the remaining European pops.
FIN.2<-All.Res3[[5]]
GBR.2<-All.Res3[[6]]
TSI.2<-All.Res3[[7]]
mclapply(bin.vec1.eu, function(x) (which(GBR.2$Nr.IS==x)))->temp.GBR
mclapply(bin.vec1.eu, function(x) (which(TSI.2$Nr.IS==x)))->temp.TSI
mclapply(bin.vec1.eu, function(x) (which(FIN.2$Nr.IS==x)))->temp.FIN

which(GBR.2$Nr.IS>=bin.vec2.eu[[1]])->temp2.GBR
which(TSI.2$Nr.IS>=bin.vec2.eu[[1]])->temp2.TSI
which(FIN.2$Nr.IS>=bin.vec2.eu[[1]])->temp2.FIN


system.time(for (i in 1:length(temp.GBR)){
I<-temp.GBR[[i]]
unlist(mclapply(I, function(x) (sum(GBR.2$NCDf5[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.5)/nsims)))->GBR.2$P.val.NCDf0.5[I]
unlist(mclapply(I, function(x) (sum(GBR.2$NCDf4[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.4)/nsims)))->GBR.2$P.val.NCDf0.4[I]
unlist(mclapply(I, function(x) (sum(GBR.2$NCDf3[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.3)/nsims)))->GBR.2$P.val.NCDf0.3[I]
unlist(mclapply(I, function(x) (sum(GBR.2$NCDf2[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.2)/nsims)))->GBR.2$P.val.NCDf0.2[I]
unlist(mclapply(I, function(x) (sum(GBR.2$NCDf1[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.1)/nsims)))->GBR.2$P.val.NCDf0.1[I]
cat(i, 'done\n') #1936.337 
})


unlist(mclapply(temp2.GBR, function(x) (sum(GBR.2$NCDf5[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.5)/nsims)))->GBR.2$P.val.NCDf0.5[temp2.GBR]
unlist(mclapply(temp2.GBR, function(x) (sum(GBR.2$NCDf4[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.4)/nsims)))->GBR.2$P.val.NCDf0.4[temp2.GBR]
unlist(mclapply(temp2.GBR, function(x) (sum(GBR.2$NCDf3[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.3)/nsims)))->GBR.2$P.val.NCDf0.3[temp2.GBR]
unlist(mclapply(temp2.GBR, function(x) (sum(GBR.2$NCDf2[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.2)/nsims)))->GBR.2$P.val.NCDf0.2[temp2.GBR]
unlist(mclapply(temp2.GBR, function(x) (sum(GBR.2$NCDf1[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.1)/nsims)))->GBR.2$P.val.NCDf0.1[temp2.GBR]

system.time(for (i in 1:length(temp.FIN)){
I<-temp.FIN[[i]]
unlist(mclapply(I, function(x) (sum(FIN.2$NCDf5[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.5)/nsims)))->FIN.2$P.val.NCDf0.5[I]
unlist(mclapply(I, function(x) (sum(FIN.2$NCDf4[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.4)/nsims)))->FIN.2$P.val.NCDf0.4[I]
unlist(mclapply(I, function(x) (sum(FIN.2$NCDf3[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.3)/nsims)))->FIN.2$P.val.NCDf0.3[I]
unlist(mclapply(I, function(x) (sum(FIN.2$NCDf2[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.2)/nsims)))->FIN.2$P.val.NCDf0.2[I]
unlist(mclapply(I, function(x) (sum(FIN.2$NCDf1[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.1)/nsims)))->FIN.2$P.val.NCDf0.1[I]
cat(i, 'done\n')  #1897.888 
})

unlist(mclapply(temp2.FIN, function(x) (sum(FIN.2$NCDf5[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.5)/nsims)))->FIN.2$P.val.NCDf0.5[temp2.FIN]
unlist(mclapply(temp2.FIN, function(x) (sum(FIN.2$NCDf4[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.4)/nsims)))->FIN.2$P.val.NCDf0.4[temp2.FIN]
unlist(mclapply(temp2.FIN, function(x) (sum(FIN.2$NCDf3[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.3)/nsims)))->FIN.2$P.val.NCDf0.3[temp2.FIN]
unlist(mclapply(temp2.FIN, function(x) (sum(FIN.2$NCDf2[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.2)/nsims)))->FIN.2$P.val.NCDf0.2[temp2.FIN]
unlist(mclapply(temp2.FIN, function(x) (sum(FIN.2$NCDf1[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.1)/nsims)))->FIN.2$P.val.NCDf0.1[temp2.FIN]

system.time(for (i in 1:length(temp.TSI)){
I<-temp.TSI[[i]]
unlist(mclapply(I, function(x) (sum(TSI.2$NCDf5[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.5)/nsims)))->TSI.2$P.val.NCDf0.5[I]
unlist(mclapply(I, function(x) (sum(TSI.2$NCDf4[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.4)/nsims)))->TSI.2$P.val.NCDf0.4[I]
unlist(mclapply(I, function(x) (sum(TSI.2$NCDf3[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.3)/nsims)))->TSI.2$P.val.NCDf0.3[I]
unlist(mclapply(I, function(x) (sum(TSI.2$NCDf2[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.2)/nsims)))->TSI.2$P.val.NCDf0.2[I]
unlist(mclapply(I, function(x) (sum(TSI.2$NCDf1[x]>=l.bin.vec1.eu[[i]]$ncvFD_f0.1)/nsims)))->TSI.2$P.val.NCDf0.1[I]
cat(i, 'done\n')
})

unlist(mclapply(temp2.TSI, function(x) (sum(TSI.2$NCDf5[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.5)/nsims)))->TSI.2$P.val.NCDf0.5[temp2.TSI]
unlist(mclapply(temp2.TSI, function(x) (sum(TSI.2$NCDf4[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.4)/nsims)))->TSI.2$P.val.NCDf0.4[temp2.TSI]
unlist(mclapply(temp2.TSI, function(x) (sum(TSI.2$NCDf3[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.3)/nsims)))->TSI.2$P.val.NCDf0.3[temp2.TSI]
unlist(mclapply(temp2.TSI, function(x) (sum(TSI.2$NCDf2[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.2)/nsims)))->TSI.2$P.val.NCDf0.2[temp2.TSI]
unlist(mclapply(temp2.TSI, function(x) (sum(TSI.2$NCDf1[x]>=l.bin.vec2.eu[[1]]$ncvFD_f0.1)/nsims)))->TSI.2$P.val.NCDf0.1[temp2.TSI]

pdf('figures/october.2016.vioplots.GBR.pdf')
par(mfrow=c(4,4))
sapply(1:16, function(x) vioplot(l.bin.vec1.eu[[x]]$ncvFD_f0.5, GBR.2[temp.GBR[[x]],]$NCDf5, col='cornflowerblue', border='gray', rectCol=' white', colMed='black',names=c('sims', 'data')))
dev.off()


pdf('figures/october.2016.vioplots.FIN.pdf')
par(mfrow=c(4,4))
sapply(1:16, function(x) vioplot(l.bin.vec1.eu[[x]]$ncvFD_f0.5, FIN.2[temp.FIN[[x]],]$NCDf5, col='cornflowerblue', border='gray', rectCol=' white', colMed='black',names=c('sims', 'data')))
dev.off()

pdf('figures/october.2016.vioplots.TSI.pdf')
par(mfrow=c(4,4))
sapply(1:16, function(x) vioplot(l.bin.vec1.eu[[x]]$ncvFD_f0.5, TSI.2[temp.TSI[[x]],]$NCDf5, col='cornflowerblue', border='gray', rectCol=' white', colMed='black',names=c('sims', 'data')))
dev.off()

objName1<-'GBR.2';objName2<-'TSI.2';objName3<-'FIN.2'

save(list=objName1, file='GBR.2.RData')
save(list=objName2, file='TSI.2.RData')
save(list=objName3, file='FIN.2.RData')

Store(FIN.2);Store(GBR.2);Store(TSI.2)
###############################################################################################
####################################################################################################################################################################################
Objects()
list.SCAN<-vector('list', 7)
names(list.SCAN)<-c("AWS","LWK", "YRI", "CEU", "FIN", "GBR", "TSI")

list.SCAN[[1]]<-AWS.2;list.SCAN[[2]]<-LWK.2;list.SCAN[[3]]<-YRI.2;list.SCAN[[4]]<-CEU.2;list.SCAN[[5]]<-FIN.2;list.SCAN[[6]]<-GBR.2;list.SCAN[[7]]<-TSI.2

Store(list.SCAN)

Store(AWS.2);Store(LWK.2);Store(YRI.2);Store(CEU.2);Store(FIN.2);Store(GBR.2);Store(TSI.2)
######################################################################################################################################################################################

Store(bin.vec1.eu);Store(bin.vec1);Store(bin.vec2.eu)
Store(bin.vec2);Store(l.bin.vec1)
Store(l.bin.vec1.eu);Store(l.bin.vec2.eu)
Store(l.bin.vec2);Store(list.bin.vec1.eu);Store(list.bin.vec1)
Store(list.bin.vec2);Store(list.bin.vec2.eu)

gc()
####################################################################################################

#adiitional plot
X<-AFRICA[[2]]

#3 vectors for the bins
bin.vec1<-seq(from=1, to=229) #1        #all windows with at least 10 IS
bin.vec2<-230
nsims<-10000
system.time(mclapply(bin.vec1, function(x) subset(X, Nr.IS==x))->list.bin.vec1) #266
system.time(mclapply(bin.vec2, function(x) subset(X, Nr.IS>=x))-> list.bin.vec2)
#system.time(mclapply(list.bin.vec1, function(x) x[sample(seq(1:dim(x)[1]), nsims),])-> l.bin.vec1)
#system.time(mclapply(list.bin.vec2, function(x) x[sample(seq(1:dim(x)[1]),nsims),])-> l.bin.vec2)




pdf('figures/october.2016.AFRICA.NrIS.NCV.neutral.p0.01.ALLBINS.pdf')
plot(c(seq(from=0.01, to=0.5, by=0.01), rep(0.3,180))~seq(1:230), type='n', ylab='NCD2', xlab='Number of Informative Sites', main='Simulations for Africa (1% quantile)')
points(c(unlist(lapply(list.bin.vec1, function(x) quantile(x$ncvFD_f0.5, prob=0.01))),quantile(list.bin.vec2[[1]]$ncvFD_f0.5, prob=0.01))~c(bin.vec1,bin.vec2), col='cornflowerblue', pch=20)
lines(c(unlist(lapply(list.bin.vec1, function(x) quantile(x$ncvFD_f0.5, prob=0.01))),quantile(list.bin.vec2[[1]]$ncvFD_f0.5, prob=0.01))~c(bin.vec1,bin.vec2), col= 'lightgray', cex=0.2)
points(c(unlist(lapply(list.bin.vec1, function(x) quantile(x$ncvFD_f0.4, prob=0.01))),quantile(list.bin.vec2[[1]]$ncvFD_f0.4, prob=0.01))~c(bin.vec1,bin.vec2), col='sienna1', pch=20)
lines(c(unlist(lapply(list.bin.vec1, function(x) quantile(x$ncvFD_f0.4, prob=0.01))),quantile(list.bin.vec2[[1]]$ncvFD_f0.4, prob=0.01))~c(bin.vec1,bin.vec2), col='lightgray', cex=0.2)
points(c(unlist(lapply(list.bin.vec1, function(x) quantile(x$ncvFD_f0.3, prob=0.01))),quantile(list.bin.vec2[[1]]$ncvFD_f0.3, prob=0.01))~c(bin.vec1,bin.vec2), col='violetred1', pch=20)
lines(c(unlist(lapply(list.bin.vec1, function(x) quantile(x$ncvFD_f0.3, prob=0.01))),quantile(list.bin.vec2[[1]]$ncvFD_f0.3, prob=0.01))~c(bin.vec1,bin.vec2), col='lightgray', cex=0.2)
#points(c(unlist(lapply(l.bin.vec1, function(x) quantile(x$ncvFD_f0.2, prob=0.01))),quantile(l.bin.vec2[[1]]$ncvFD_f0.2, prob=0.01))~c(bin.vec1,bin.vec2), col='darkolivegreen', pch=20)
#lines(c(unlist(lapply(l.bin.vec1, function(x) quantile(x$ncvFD_f0.2, prob=0.01))),quantile(l.bin.vec2[[1]]$ncvFD_f0.2, prob=0.01))~c(bin.vec1,bin.vec2), col='lightgray', cex=0.2)
abline(v=10, col='gray', lty=2)
legend('bottomright', c('tf=feq=0.5', 'tf=feq=0.4','tf=feq=0.3'), col=c('cornflowerblue', 'sienna1', 'violetred1'), pch=20, bty='n')
dev.off()




