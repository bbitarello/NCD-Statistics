##########################################################################################
#       Use Phase 1 1000G data to verify allele frequencies within inversion in chr4
#       Bárbara Bitarello
#       Created: 21.04.2017
#       Last modified: 25.04.2017
##########################################################################################

#do: system.time(source('run_NCD2.r')) 808.963 seconds for four populations!
#preamble
library(pegas);library(dplyr)
library(plyr);library(data.table)
library(parallel);library(lattice)
library(SOAR);Sys.setenv(R_LOCAL_CACHE="inversions")
library(ggplot2);library(splitstackshape)
library(pryr)
library(doMC)
registerDoMC(22)
#library(checkpoint)
#checkpoint("2017-04-25")
Objects()
##################################################################################

#test:

#victor's code with modifications by me:o
#victor did NCD1 and I am doing NCD1
NCD2 <- function(X, Y,  W = 100000, S = 50000) {
#X[, c("AF1","AF2","AF3","nAL","N_chr","AF"):=NULL]  
X[order(as.numeric(POS))]-> X 
  windows_dt <- 
    data.table(POS = seq(X[1, POS], X[nrow(X), POS], S))[
      , POS2 := POS + W][
        -length(POS)]; #create dt with POS, then add POS2, then remove last row.

    print (paste0('Finished setting up coordinates for chr ', unique(X$CHR)))
  
  setkey(windows_dt, POS, POS2) #set two keys makes processing faster.
  X[, POS2 := POS] # := creates new column. Why is he doing this again? I think it's because last clolum needs to be pos for the foroverlaps function
  
  Y[, POS2 := POS]
  X_windows <-
    foverlaps(X, windows_dt, type = "within", nomatch = 0L)[ #this is not ideal because counts end position but it's okay as long as i document this well.
      , window := .GRP, by = .(POS, POS2)][
        order(window, i.POS)][
          , .(Win.ID = paste(CHR, POS, POS2, sep = "|"), MAF)][,N_Raw:=.N,by = Win.ID][,N_SNPs_cor:= sum(MAF != 0 & MAF != 1),by = Win.ID]
  print (paste0('Finished setting up SNPs for chr ', unique(X$CHR)))
   Y_windows<-
   foverlaps(Y, windows_dt, type = "within", nomatch = 0L)[ #this is not ideal because counts end position but it's okay as long as i document this well.
      , window := .GRP, by = .(POS, POS2)][
        order(window, i.POS)][
          , .(Win.ID = paste(CHR, POS, POS2, sep = "|"))][,N_FDs_Raw:=.N,by = Win.ID]
	setkey(X_windows, Win.ID)
	setkey(Y_windows, Win.ID)

	setkey(X, ID, POS2)
	setkey(Y, ID, POS2)
	Y[!(ID %in% Y[X[MAF!=0 & MAF!=1], nomatch=0][ALT==Chimp_REF][,ID])]-> Y2
	setkey(Y2, POS, POS2)
	Y2[order(POS)]-> Y2 #corrected set of FDs

    Y2_windows<-
   foverlaps(Y2, windows_dt, type = "within", nomatch = 0L)[ #this is not ideal because counts end position but it's okay as long as i document this well.
      , window := .GRP, by = .(POS, POS2)][
        order(window, i.POS)][
          , .(Win.ID = paste(CHR, POS, POS2, sep = "|"))][,N_FDs_cor:=.N,by = Win.ID]
	
	setkey(Y2_windows, Win.ID)
	unique(X_windows)-> X2_windows
	unique(Y2_windows)->Y3_windows
	unique(Y_windows)-> Y_u_windows
	setkey(X2_windows, Win.ID)
	setkey(Y3_windows, Win.ID)
	setkey(Y_u_windows, Win.ID)
  print (paste0('Finished selecting SNPs and FDs  per window for chr ', unique(X$CHR)))
	#final dt
	X2_windows[Y_u_windows][Y3_windows][,PtoD:= N_SNPs_cor/(N_FDs_cor+1)][, MAF:=NULL][, IS:=N_SNPs_cor+N_FDs_cor]-> Z

  X_NCD <-
    X_windows[Z][MAF!=0 & MAF!=1][
      , .(N_Raw= N_Raw,
	  N_SNPs_cor=N_SNPs_cor,
		N_FDs_Raw=N_FDs_Raw,
		N_FDs_cor=N_FDs_cor,
		PtoD=PtoD,
		IS=IS,	
          NCD2_tf0.5 = sqrt(sum((c(rep(0,unique(N_FDs_cor)),MAF)-0.5)^2)/IS),  
          NCD2_tf0.4 = sqrt(sum((c(rep(0,unique(N_FDs_cor)),MAF)-0.4)^2)/IS),
          NCD2_tf0.3 = sqrt(sum((c(rep(0,unique(N_FDs_cor)),MAF)-0.3)^2)/IS)), 
      by = Win.ID]
  
  print (paste0('NCD2 calculations done for chr ', unique(X$CHR)))
  
	setkey(X_NCD, Win.ID)
	unique(X_NCD)
}

#it appears to be working fine. Need to do some manual calculatioins and check results for both NCD1 and NCD2. 

system.time(LWK.run2<-foreach(x=1:22, .combine="rbind", .packages=c("data.table")) %dopar%
         NCD2(X=LWK[[x]], Y=FD_list[[x]],  W=3000, S=1500)); # 157 seconds!

na.omit(LWK.run2)-> LWK_NCD2_gen;

LWK_NCD2_gen[IS>=10]-> LWK_NCD2_gen_IS; #filter min informative sitesa

LWK_NCD2_gen_IS[, c("Chr", "POS1","POS2") := tstrsplit(Win.ID, "|", fixed=TRUE)];
LWK_NCD2_gen_IS[order(as.numeric(Chr), as.numeric(POS1))]-> LWK_NCD2_gen_IS;

#

system.time(YRI.run2<-foreach(x=1:22, .combine="rbind", .packages=c("data.table")) %dopar%
         NCD2(X=YRI[[x]],FD_list[[x]], W=3000, S=1500)); # 58

na.omit(YRI.run2)-> YRI_NCD2_gen;

YRI_NCD2_gen[IS>=10]-> YRI_NCD2_gen_IS; #filter min informative sitesa

YRI_NCD2_gen_IS[, c("Chr", "POS1","POS2") := tstrsplit(Win.ID, "|", fixed=TRUE)];
YRI_NCD2_gen_IS[order(as.numeric(Chr), as.numeric(POS1))]-> YRI_NCD2_gen_IS;

#
system.time(GBR.run2<-foreach(x=1:22, .combine="rbind", .packages=c("data.table")) %dopar%
         NCD2(X=GBR[[x]],FD_list[[x]], W=3000, S=1500)); # 

na.omit(GBR.run2)-> GBR_NCD2_gen;

GBR_NCD2_gen[IS>=10]-> GBR_NCD2_gen_IS; #filter min informative sitesa

GBR_NCD2_gen_IS[, c("Chr", "POS1","POS2") := tstrsplit(Win.ID, "|", fixed=TRUE)];
GBR_NCD2_gen_IS[order(as.numeric(Chr), as.numeric(POS1))]-> GBR_NCD2_gen_IS;

#

system.time(TSI.run2<-foreach(x=1:22, .combine="rbind", .packages=c("data.table")) %dopar%
         NCD2(X=TSI[[x]],Y=FD_list[[x]], W=3000, S=1500)); # 

na.omit(TSI.run2)-> TSI_NCD2_gen;

TSI_NCD2_gen[IS>=10]-> TSI_NCD2_gen_IS; #filter min informative sitesa

TSI_NCD2_gen_IS[, c("Chr", "POS1","POS2") := tstrsplit(Win.ID, "|", fixed=TRUE)];
TSI_NCD2_gen_IS[order(as.numeric(Chr), as.numeric(POS1))]-> TSI_NCD2_gen_IS;


NCD2_res_pops<-vector('list', 4)

NCD2_res_pops[[1]]<-LWK_NCD2_gen_IS
NCD2_res_pops[[2]]<-YRI_NCD2_gen_IS
NCD2_res_pops[[3]]<-GBR_NCD2_gen_IS
NCD2_res_pops[[4]]<-TSI_NCD2_gen_IS

names(NCD2_res_pops)<-c('LWK','YRI','GBR','TSI')






