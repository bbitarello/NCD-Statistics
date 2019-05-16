library(data.table)
library(dplyr)
load("/mnt/expressions/miguel/gwas_expr/expr.rdat") #Data provided by Michael Dannemmann
fread('~/NCD-Statistics/read_scan_data/bedfiles/all_scanned_genes.txt', header=F)-> scanned

fread('~/NCD-Statistics/read_scan_data/bedfiles/extreme.genes.bed')-> extreme_genes
fread('~/NCD-Statistics/read_scan_data/bedfiles/cand.extreme.genes.bed')-> cand_extreme_genes

extreme_genes[which(extreme_genes$V2 %in% ex.t.gc.dt.scanned[[1]]$Gene_name)]-> edited_extreme_genes #233
extreme_genes[which(cand_extreme_genes$V2 %in% ex.t.gc.dt.scanned[[1]]$Gene_name)]-> edited_cand_extreme_genes #1445
library(parallel)

mclapply(1:53, function(x) setDT(as.data.frame(ex.t.gc[[x]]), keep.rownames=T))-> ex.t.gc.dt
for(i in 1: 53){
colnames(ex.t.gc.dt[[i]])[1]<-'Gene_name'
}

ex.t.gc.dt.scanned<-vector('list', 53)

for (i in 1:53){
ex.t.gc.dt[[i]] %>% dplyr::filter(Gene_name %in% scanned$V1) %>% as.data.table-> ex.t.gc.dt.scanned[[i]]
}
remove(ex.t.gc)
remove(ex.t.gc.dt)


#arbitrary cutofss
#these measures are not normalized by gene length, but only by all idnividuals. So some genes have huge values because they are larger.
#if more than 5 individuals out of 361 have exrpession above an also arbitrary cutoff (100) we say the gene is expressed.
#also, this dataset only contains genes that gave >0 for at least one individual, so non-expressed genes not considered.

for(i in 1:53){
	cbind(ex.t.gc.dt.scanned[[i]], data.frame(RS0=rowSums(ex.t.gc.dt.scanned[[i]][,-1]>0), 
	RS1=rowSums(ex.t.gc.dt.scanned[[i]][,-1]>1), RS10=rowSums(ex.t.gc.dt.scanned[[i]][,-1]>10), 
	RS100=rowSums(ex.t.gc.dt.scanned[[i]][,-1]>100), RS500=rowSums(ex.t.gc.dt.scanned[[i]][,-1]>500))) %>% as.data.table %>% mutate(Expr_Status=ifelse(RS100>5, "TRUE", "FALSE")) %>% as.data.table -> ex.t.gc.dt.scanned[[i]]
}


mclapply(1:265, function(x) mclapply(1:53,function(i) ex.t.gc.dt.scanned[[i]] %>% dplyr::filter(Gene_name %in% extreme_genes$V2[x]) %>% select(Expr_Status)))-> top
unlist(mclapply(1:265, function(x) sum(unlist(top[[x]])=="TRUE")))
summary(unlist(mclapply(1:265, function(x) sum(unlist(top[[x]])=="TRUE"))))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   0.00   14.00   45.00   33.39   52.00   52.00 

system.time(mclapply(1:1594, function(x) mclapply(1:53,function(i) ex.t.gc.dt.scanned[[i]] %>% dplyr::filter(Gene_name %in% cand_extreme_genes$V2[x]) %>% select(Expr_Status)))-> cand) #700s 


all<-vector('list', 53)


n<-sapply(1:53, function(i) nrow(ex.t.gc.dt.scanned[[i]]))

for(j in 1:53){
	mclapply(1:n[j], function(x) ex.t.gc.dt.scanned[[i]] %>% dplyr::filter(Gene_name %in% ex.t.gc.dt.scanned[[i]]$Gene_name[x]) %>% select(Expr_Status))-> all[[i]]
	print(i)
}

