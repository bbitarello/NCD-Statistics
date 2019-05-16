## Cesare de Filippo
## 15-01-2014
## Script to generate positions where the human and chimpanzee reference genomes differ.
## OUTPUTS:
## - Map50_100.TRF.SDs.hg19_pantro2.[CHROM].bed
##   file where human and chimp have sequences. The positions where chimp is deleted are not considered.  
## - fds.hg19_pantro2.[CHROM].tsv
##   file with the positions where human and chimp differ. #rd and 4th columns are the hg19 and pantro2 genomes.
## INPUTS:
## - ~/ucsc/goldenPath/hg19/vsPanTro2/axtNet/mafnet/chr[CHROM].hg19.panTro2.net.maf
## - Map50_100.TRF.SDs.bed.gz
#Modified by Bárbara Bitarello (May 2019)



options(scipen=1);
library(multicore)
Seq <- function(START, END, bed=TRUE) {
  if(length(START) == length(END)) {
    if (bed == TRUE) {
      n = 1
    } else {
      n=0
    }
    return(unlist(sapply(1:length(START), function(x) (START[x]+n):END[x])))
  }
}

for(k in c(1:22,"X")) {
  BED <- paste("tmp",k,".bed",sep="")
  system(paste("tabix -p bed Map50_100.TRF.SDs.bed.gz",k,">",BED))
  x <- read.table(BED)
  n <- nrow(x)
  w <- 1000
  b <- seq(1,n,w)
  bins <- cbind(b, c(b[-1]-1, n))
  k=unique(x[,1])
  START <- as.numeric(x[,2])
  END <- as.numeric(x[,3])
  out <- unlist(mclapply(1:nrow(bins), function(y) Seq(START[bins[y,1]:bins[y,2]],END[bins[y,1]:bins[y,2]],bed=TRUE)))
  write.table(out, file="tmp.pos", col.names=F,quote=F,sep="\t",row.names=F)
  POS=paste("pos",k,sep=".")
  POS.retrieved <- paste(POS,".tsv",sep="")
  FDS <- paste("fds",k,"tsv",sep=".")
  FDS.output <- paste("fds.hg19_pantro2.",k,".tsv",sep="")
  system(paste("sed 's/^/chr",k,"\t/' tmp.pos > ", POS, sep=""))
  system("rm tmp.pos")
  CHAINFILE=paste("~/ucsc/goldenPath/hg19/vsPanTro2/axtNet/mafnet/chr",k,".hg19.panTro2.net.maf",sep="")
  system(paste("cat", POS, " | /home/pruefer/src/BamSNPTool/BamSNPAddMaf",CHAINFILE, "hg19 pantro2 | sed 's/chr//' > ", POS.retrieved)) #this software was written by K.Pruefer
  system(paste("rm",POS))
  system(paste("cat", POS.retrieved, "| awk ' tolower($3) != tolower($4) {print $0}' > ", FDS))
  system(paste("grep '-' ",FDS," | cut -f 1,2 > tmp.pos")) # position where pantro2 is deleted. 
  system(paste("grep -v '-' ",FDS," > ", FDS.output))
  system("pos2bed.sh tmp.pos > pos.rm.bed")
  system(paste("subtractBed -a tmp",k,".bed -b pos.rm.bed > Map50_100.TRF.SDs.hg19_pantro2.",k,".bed",sep=""))
  system(paste("rm *tmp* pos.rm.bed",POS.retrieved, FDS))
  cat(k,"done\n")
}

