require(data.table)
library(doMC)
registerDoMC(11)
library(dplyr)

NCD1 <- function(X, W = 100000, S = 50000) {

  windows_dt <-
    data.table(POS = seq(X[1, POS], X[nrow(X), POS], S))[
      , POS2 := POS + W][
        -length(POS)] #

    print (paste0('Finished setting up coordinates for chr ', unique(X$CHR)))

    setkey(windows_dt, POS, POS2) #
    X[, POS2 := POS] #
    setkey(X, POS, POS2)
  
   X_windows <-
    foverlaps(X, windows_dt, type = "within", nomatch = 0L)[ #this is not ideal but for now it's fine
      , window := .GRP, by = .(POS, POS2)][
        order(window, i.POS)][
          , .(Win.ID = paste(CHR, POS, POS2, sep = "|"), MAF)][,N_Raw:=.N,by = Win.ID]

  print (paste0('Finished selecting SNPs per window for chr ', unique(X$CHR)))

  X_NCD <-
    X_windows[MAF != 0 & MAF != 1][
      , .(N_Raw= N_Raw,
          N_SNPs_cor = .N,
          NCD1_tf0.5 = sqrt(sum((MAF-0.5)^2)/.N),  #.N is the equivalent of n() in dplyr.
          NCD1_tf0.4 = sqrt(sum((MAF-0.4)^2)/.N),
          NCD1_tf0.3=sqrt(sum((MAF-0.3)^2)/.N)),
      by = Win.ID]

  print (paste0('NCD1 calculations done for chr ', unique(X$CHR)))

  unique_windows <-

    X_windows[, .(N_SNPs_cor = sum(MAF != 0 & MAF != 1)), by = Win.ID]

  setkey(unique_windows, Win.ID, N_SNPs_cor)
  setkey(X_NCD, Win.ID, N_SNPs_cor)

  unique(X_NCD[unique_windows])
}


#

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
  setkey(X, POS,POS2)
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
	if(nrow(Y2_windows)>0){
        X2_windows[Y_u_windows][Y3_windows][,PtoD:= N_SNPs_cor/(N_FDs_cor+1)][, MAF:=NULL][, IS:=N_SNPs_cor+N_FDs_cor]-> Z


  X_NCD <-
    X_windows[Z, allow.cartesian=TRUE][MAF!=0 & MAF!=1][
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
if(nrow(Y2_windows)==0){
      X2_NCD <-
        X_windows[!(Win.ID %in% unique(Y2_windows$Win.ID))][MAF!=0 & MAF!=1][
                 , .(N_Raw= N_Raw,
                N_SNPs_cor=N_SNPs_cor,
                N_FDs_Raw=0,
                N_FDs_cor=0,
                PtoD=N_SNPs_cor/(0+1),
                IS=N_SNPs_cor,
          NCD2_tf0.5 = sqrt(sum((MAF-0.5)^2)/N_SNPs_cor),
          NCD2_tf0.4 = sqrt(sum((MAF-0.4)^2)/N_SNPs_cor),
          NCD2_tf0.3 = sqrt(sum((MAF-0.3)^2)/N_SNPs_cor)),
      by = Win.ID]

        setkey(X2_NCD, Win.ID)
        unique(X2_NCD)

     #   rbind(X_NCD, X2_NCD)-> X3_NCD

  print (paste0('NCD2 calculations done for chr ', unique(X$CHR)))

        setkey(X2_NCD, Win.ID)
        unique(X2_NCD)

}

}
