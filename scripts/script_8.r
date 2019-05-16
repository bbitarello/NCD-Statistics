##################################
#	Make tables for paper
#
#	Barbara Bitarello
#	Last modified: 13.12.2016	
##################################

library(SOAR)
Sys.setenv(R_LOCAL_CACHE="estsession")
library(parallel)
source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/mclapply2.R')
source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/assign_tf_gene.R')
source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/find_gene.R')
source('/mnt/sequencedb/PopGen/barbara/NCV_dir_package/scripts/assign_tf_window_original.R')
library(data.table)
library(dplyr)
Objects()

as.character(read.table('bedfiles/African.Genes.bed')$V1)-> Afr
mclapply2(Afr, function(y) funcA(name2=y))-> res.afr
setDT(do.call('rbind',res.afr))-> res.afr2
as.character(res.afr2$Acronym)-> res.afr2$Acronym
as.numeric(res.afr2$Chr)-> res.afr2$Chr
setDT(arrange(res.afr2, p.LWK, p.YRI, p.GBR, p.TSI)-> res.afr2)
write.table(res.afr2, file='bedfiles/Table_afr_manuscript.txt', row.names=F) #check
Store(res.afr2)
#
as.character(read.table('bedfiles/European.Genes.bed')$V1)-> Eur
mclapply2(Eur, function(y) funcA(name2=y))-> res.eur
setDT(do.call('rbind',res.eur))-> res.eur2
as.character(res.eur2$Acronym)-> res.eur2$Acronym
as.numeric(res.eur2$Chr)-> res.eur2$Chr
setDT(arrange(res.eur2, p.GBR, p.TSI, p.LWK, p.YRI)-> res.eur2)
write.table(res.eur2, file='bedfiles/Table_eur_manuscript.txt', row.names=F)
Store(res.eur2)


#significant shared genes  #gettin error here.
as.character(read.table('bedfiles/cand.African.Genes.bed')$V1)-> Afr.cand
mclapply2(Afr.cand, function(y) funcA(name2=y))-> res.afr.cand
setDT(do.call('rbind',res.afr.cand))-> res.afr.cand2
as.character(res.afr.cand2$Acronym)->res.afr.cand2$Acronym
as.numeric(res.afr.cand2$Chr)-> res.afr.cand2$Chr
setkey(res.afr.cand2, Acronym); unique(res.afr.cand2)-> res.afr.cand2 #remove duplicated gene entries.
setDT(arrange(res.afr.cand2, p.LWK, p.YRI, p.GBR, p.TSI)-> res.afr.cand2)
write.table(res.afr.cand2, file='bedfiles/Table_afr_cand_manuscript.txt', row.names=F)
Store(res.afr.cand2)
#
as.character(read.table('bedfiles/cand.just.Afr.genes.bed')$V1)-> just.Afr.cand
filter(res.afr.cand2, Acronym %in% just.Afr.cand)-> res.just.Afr.cand
as.character(read.table('bedfiles/cand.just.Eur.genes.bed')$V1)-> just.Eur.cand
filter(res.eur.cand2, Acronym %in% just.Eur.cand)-> res.just.Eur.cand
as.character(read.table('bedfiles/cand.afrANDeur.genes.bed')$V1)-> Afr.and.Eur.cand
filter(res.eur.cand2, Acronym %in% Afr.and.Eur.cand)-> res.afrandeur.cand

setDT(arrange(res.just.Afr.cand, p.LWK, p.YRI, p.GBR, p.TSI)-> res.just.Afr.cand2)

setDT(arrange(res.just.Eur.cand, p.GBR, p.TSI, p.LWK, p.YRI)-> res.just.Eur.cand2)

setDT(arrange(res.afrandeur.cand, p.LWK, p.YRI, p.GBR, p.TSI))-> res.afrandeur.cand2

write.table(res.just.Eur.cand2, file='bedfiles/Table_justeur_cand_manuscript.txt', row.names=F)

write.table(res.just.Afr.cand2, file='bedfiles/Table_justafr_cand_manuscript.txt', row.names=F)

write.table(res.afrandeur.cand2, file='bedfiles/Table_afrandeur_cand_manuscript.txt', row.names=F)

Store(res.just.Eur.cand2); Store(res.just.Afr.cand2)
#
as.character(read.table('bedfiles/cand.European.Genes.bed')$V1)-> Eur.cand
mclapply2(Eur.cand, function(y) funcA(y))-> res.eur.cand
setDT(do.call('rbind',res.eur.cand))-> res.eur.cand2
setkey(res.eur.cand2, Acronym); unique(res.eur.cand2)-> res.eur.cand2
as.character(res.eur.cand2$Acronym)-> res.eur.cand2$Acronym
as.numeric(res.eur.cand2$Chr)-> res.eur.cand2$Chr
setDT(arrange(res.eur.cand2, p.GBR, p.TSI, p.LWK, p.YRI)-> res.eur.cand2);
write.table(res.eur.cand2, file='bedfiles/Table_eur_cand_manuscript.txt', row.names=F)
Store(res.eur.cand2)



#playground

table((res.afr2 %>% filter(Acronym %in% read.table('bedfiles/afrANDeur.genes.bed')[,1]) %>% group_by(Acronym) %>% summarise(Question=sum(tf.LWK==tf.YRI)))$Question)

#0  1
#31  71



table(((res.afr2 %>% filter(Acronym %in% read.table('bedfiles/afrANDeur.genes.bed')[,1]) %>% group_by(Acronym) %>% summarise(Question=sum(tf.LWK==tf.YRI), Question2=sum(tf.GBR==tf.TSI), Third=Question+Question2)))$Third)


rbind(res.afr2, res.eur2)-> temp
setkey(temo, Acronym)
unique(temp)-> temp

rbind(res.afr.cand2, res.eur.cand2)-> temp2
setkey(temp2, Acronym)
unique(temp2)-> temp2

