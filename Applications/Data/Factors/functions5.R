##############################################
###########    Functions         #############
#********************************************#


#################################################
###### data transformations and factor estimation function: 


datafac.fun <- function(datapanel, iu_measure, polinst, nfac, beobachtbare=c("ffr","CPIAUCSL"), inf=FALSE,  
                        fakmethod="naive", stand = FALSE, period=c("1977-01-01", "2015-12-01")){
  ###Input: 
  # 1. main panel
  # 2. iu measure
  # 3. period
  # 4. number of factors, factor extraction method, observable factors, p, 
  ### output: 
  # -> save image? or list 
  pcfac <- function(data, anzahlfac){
    
    ## PC factor extraction
    ## data must be (TxN)  
    nvars <- dim(data)[2]
    
    xx <- t(data)%*%data# X'X is proportional to covariance matrix
    eigen_x <- eigen(xx) 
    
    # applying the factor restricition
    lam <- sqrt(nvars)*eigen_x[[2]][,c(1:anzahlfac)] # lam'lam/n = I. loadings are, where n = number of time series in x
    factors <- data%*%lam/nvars     
    
    return(list(factors, lam))
  }
  olsfac <- function(y, z){ 
    
    # Input: Y and Z matrices
    # Output: Esimtators
    
    # for this the input must be data(KxT)
    result <- solve(crossprod(z),crossprod(z,y))
    return(result)
  }
  
  
  ##### first data adjustment  #####
  # decide on IU measure:
  iu <- iu_measure[[1]][which(iu_measure[[2]]==period[1]):which(iu_measure[[2]]==period[2])]
  #length(iu)
  
  # policy instrument: 
  ffr_raw <- polinst[[1]][which(polinst[[2]]==period[1]):which(polinst[[2]]==period[2])]
  #length(ffr_raw)
  
  
  ### observable factors:  
  
  if (length(beobachtbare)>1){
  xindex <- which(all_names%in%beobachtbare[-1])  # index for obserable factors in data matrix
  ### adjust data for observable factors: (i.e. exclude), iff more than the ffr is included
  data <- datapanel[[1]][,-xindex]
  } else {
    
    data <- datapanel[[1]]  
    
  }
  

  dataorg <- datapanel[[1]][-c(1,2),]   ## prob: time refers to the rows without rows 1 and 2...
  yindex <- which(all_names%in%beobachtbare)  # index for obserable factors in factorest algo/function
  adj_speed <- datapanel[[1]][2, -1]
  # code for needed transformation of each time series
  tcode <- data[1,-1] 
  #data <- data[which(datapanel[[2]]==period[1]):which(datapanel[[2]]==period[2]),]
  ##
  x_da <- data[-c(1,2), -1]
  x_da <- x_da[which(datapanel[[2]]==period[1]):which(datapanel[[2]]==period[2]),]
  #### Values for Model ####
  
  #length(all_names)
  
  
  k_numb <- nfac # number of factors
  o_numb <- length(beobachtbare)  # number of additional observed variables
  var_numb <- k_numb  + o_numb  # number of explanatory variables
  
  ##### Preprocess data ####
  # Make time series stationary #
  
  # therefore use the tcodes
  # (1) no transformation
  # (2) first difference
  # (3) second difference
  # (4) log
  # (5) first difference and log
  # (6) second difference and log
  # (7) first difference & ((x_t/x_{t-1})-1)
  
  # adjust dimension for differencing
  t_data <- data.matrix(x_da)
  t_tcode <- tcode
  
  
  if((is.element(3, t_tcode) || is.element(6, t_tcode) || is.element(7, t_tcode))){
    x_dim <- dim(t_data)[1] - 2 # here we loose two observations
  }else{
    if((is.element(2, t_tcode) || is.element(5, t_tcode) )){
      x_dim <- dim(t_data)[1] -1 # we loose one observation
    }else{
      x_dim <- dim(t_data)[1] # for no transformation and log, number of 
    }                         # observation remains the same
  }
  
  begin <- (dim(t_data)[1] - x_dim) +1
  ending <- dim(t_data)[1]
  
  # apply transformation to time series
  xxx <- matrix(0, nrow = x_dim, ncol = dim(t_data)[2])
  
  for(i in 1:dim(t_data)[2]){
    if(t_tcode[i] ==1){
      xxx[,i] <- t(t_data[c(begin:ending),i]) # no transformation
      
    }else{
      if(t_tcode[i] == 2){
        xxx[,i] <- t(diff(t_data[,i]))[c((begin-1):(ending-1))] # first difference
        
      } else {
        if(t_tcode[i] == 3){
          xxx[,i] <- t(diff(diff(t_data[,i]))) # second difference
          
        } else {
          if(t_tcode[i] == 4){
            xxx[,i] <- t(log(t_data[c(begin:ending),i])) # log
            
          } else {
            if(t_tcode[i] == 5){
              xxx[,i] <- t(diff(log(t_data[,i])))[c((begin-1):(ending-1))] # log first difference
              
            } else {
              if(t_tcode[i] == 6){
                xxx[,i] <- t(diff(diff(log(t_data[,i])))) # log second difference
                
              } else {
                if(t_tcode[i] == 7){
                  xxx[,i] <- t(diff((t_data[-1,i]/t_data[-dim(t_data)[1],i])-1))
                }
              }
            }
          }
        }
      }
    }
  }
  
  x_data <- xxx
  #which(x_data==NaN)
  nob <- dim(x_data)[1] # number of observations
  xvars <- dim(x_data)[2] # number of time series used for factor estimation
  
  ##### apply data transformations to observable factors (eg. first diffs)
  if (length(beobachtbare)>1){
  if (inf==TRUE){  
    observables <- cbind(ffr_raw, c(0, diff(log(dataorg[which(datapanel[[2]]==period[1]):which(datapanel[[2]]==period[2]),xindex] )))*100) # --> inflation
    #diff(log(dataorg[which(datapanel[[2]]==period[1]):which(datapanel[[2]]==period[2]),111] )) # --> inflation
    
  } else {
    observables <- cbind(ffr_raw, dataorg[which(datapanel[[2]]==period[1]):which(datapanel[[2]]==period[2]),xindex] )  # --> CPI index
    }  
   
    # adjust length of observable variables
    ffr <- ffr_raw[(length(ffr_raw) - x_dim+1): length(ffr_raw)]     ### DO we need this? 
    nobs <- length(ffr)                                              ### DO we need this?
    obser <- observables[(length(ffr_raw) - x_dim+1): length(ffr_raw),]   
     
  } else {
    
  observables <- ffr_raw 
      # adjust length of observable variables
  ffr <- ffr_raw[(length(ffr_raw) - x_dim+1): length(ffr_raw)]     ### DO we need this? 
  nobs <- length(ffr)                                              ### DO we need this?
  obser <- observables[(length(ffr_raw) - x_dim+1): length(ffr_raw)]
  }  
    
  # Standardize data 
  x_means <- apply(x_data, 2, mean)
  x_sds <- apply(x_data, 2, sd)
  mean_mat <- matrix(rep(x_means, nob), ncol = xvars, byrow = TRUE)
  std_mat <- matrix(rep(x_sds, nob), ncol = xvars, byrow = T)
  x_sd <- (x_data -  mean_mat)/std_mat
  
  
  ##### BEGIN FACTOR EXTRACTION PART #####
  
  
  if(fakmethod == "BGM"){
    
    ### BGM factors
    facest <- fac.fun.bgh(data = x_sd, beob = obser, nfac = k_numb, niter=20, ctol = 10^(-8), standard = FALSE, yindex=yindex) 
    fac_est <- facest$factors
    loadings <- facest$loadings 
    Loadings <- facest$Loadings
  } else {
    
    if (fakmethod == "2step"){
      
     facest <- twostepfac.fun(data = x_sd, beob = obser, numfac = k_numb, standard = stand, yindex=yindex) 
     fac_est <- facest$factors
     loadings <- facest$loadings 
     Loadings <- facest$Loadings
      
    } else {
  
    ### naive factors:
    
    unobfac <- pcfac(x_sd, nfac)
    ufac <- unobfac[[1]]
    unobload <- unobfac[[2]]
    
    loadings <- t(olsfac(x_sd, cbind(obser, ufac)))
    
    facnaive <- cbind(obser, ufac)
    fac_est  <- facnaive
    Loadings <- rbind(diag(var_numb)[1:o_numb,], loadings)
    }
  }
  
  return(list(fac_est = fac_est, loadings = loadings, Loadings = Loadings, 
              datasd = x_sd, iu = iu, K=var_numb, N = xvars, stdx = x_sds, tcode = tcode, period = period ))
  
} ### for function

## fakmethod = c("naive", "BGM", "2step")
## 



###########################################################
###### Sign restriction functions:           ##############

irf.favar <- function(armat, B, loads, neg=FALSE, hor, normalize=F, shockpos, transcode, impnames, stdX){  
  nlags <- ncol(armat)/nrow(armat)
  K <- nrow(armat)
  nshocks <- nrow(B)

  
  Kp <- K*nlags
  if(nlags==1){
    A <- armat
    shockcomp <- B
  } else {
    A <- rbind(armat, diag(Kp)[1:(Kp-K),])  
    shockcomp <- rbind(B, matrix(0,nrow=(K * (nlags-1)), ncol=K))
  }
  #shock <- diag(nshocks)
  #shockcomp <- rbind(CholOmegaNL, matrix(0,nrow=(var_numb * (nlags-1)), ncol=ncol(CholOmegaNL)))
  #CholOmegaNL[,shockpos]/CholOmegaNL[shockpos, shockpos]  for 100 bps shock
  
  if (is.numeric(normalize) == TRUE){  # normalization is now done on all columns=shocks, owing to sign restrictions
    
    #for (i in 1:K){
    
    shockcomp[, shockpos] <- normalize*shockcomp[,shockpos]/shockcomp[1,shockpos]  
    #}
  }
  if(neg == TRUE){
    shockcomp <- -1*shockcomp
  } 
  
  # irf with a shock
  irf_lin <- shockcomp
  irf_imp <- irf_lin
  for(i in 2:hor){
    irf_lin <- A %*% irf_lin
    irf_imp <- cbind(irf_imp, irf_lin) # matrix with kp rows and K^2*h columns, only first K rows are interesting 
  }                                    # each K columns correspond to one h, ie. h=1,2,3 etc.; columns are shocks, 
  #rows are responses of variables
  # point estimate of the ir
  irf_Dimo <- loads %*% irf_imp[1:K,]
  rownames(irf_Dimo) <- impnames
  ##### Preprocess data ####
  # Make time series stationary #
  
  # therefore use the tcodes
  # (1) no transformation
  # (2) first difference
  # (3) second difference
  # (4) log
  # (5) first difference and log
  # (6) second difference and log
  # (7) first difference & ((x_t/x_{t-1})-1)
  select <- which(rownames(irf_Dimo) %in% impnames)
  
  #c("ffr", "INDPRO", "UNRATE", "AMBSL", "M1SL", "M2SL", "TOTRESNS", "TB3MS", 
  #                       "TB6MS", "GS10", "T5YFFM", "T10YFFM", "CPIAUCSL", "PCEPI", "CUSR0000SAC","DJINDUS","VXOCLSx", "CSUSHPISA"))
  #c("ffr", "CPIAUCSL","INDPRO", "TB3MS", "M1SL","CUSR0000SAC")
  
  for(n in select ){
    
    temp <- vector("numeric", hor)
    
    for(s in 1:K){
      
      for(j in 1:hor){
        
        ### get IRFs for one shock over all horizons in 
        temp[j]  <- irf_Dimo[n,(s+((j-1)*K))]
        
      } #for hor
      
      
      if(transcode[n]==1){
        temp2 <- temp*stdX[n]
        
      } else{
        
        if(transcode[n]==2){
          temp <- temp*stdX[n]
          temp2   <- cumsum(temp)
          
        } else {
          
          if(transcode[n]==3){
            temp <- temp*stdX[n]
            temp2   <- cumsum(cumsum(temp))
            
          } else {
            if(transcode[n]==4){
              temp <- temp*stdX[n]
              temp2   <- (exp(temp) - rep(1, (hor)))
              
            } else {
              if(transcode[n]==5){
                if(n == 2) {
                  temp <- temp*stdX[n]
                  temp2 <- (exp(cumsum(temp)) - rep(1, hor))
                  } else {
                temp <- temp*stdX[n]
                temp2 <- (exp(cumsum(temp)) - rep(1, hor))*100
                  }
              } else {
                if(transcode[n]==6){
                  temp <- temp*stdX[n]
                  temp2 <- (exp(cumsum(temp)) - rep(1, hor))*100
                  #temp2 <- (cumsum(exp(cumsum(temp)) - rep(1, hor)))
                }
                
              }
            }
          }
        }
      }
      
      
      for(j in 1:hor){
        
        irf_Dimo[n,(s+((j-1)*K))]  <- temp2[j]   
        
      }}
  } ## for n
  return(irf_Dimo)
}




### function for rotation and checking admissibility of proposed models: (needed for SR_regime)
signcheck <- function (a, hfirst, hlast, sr, impulses){
  
  QR <- qr(a)
  R <- qr.R(QR)
  Q <- qr.Q(QR)
  if (any(diag(R) < 0)) {
    Q[,which(diag(R)<0)] <- (-1)*Q[,which(diag(R)<0)]
  }
  impdraw        <- unname(do.call(cbind,lapply(seq_len(hlast),function(j) impulses[,(((j-1)*K)+1):(j*K)]%*%Q)))
  impeval        <- impdraw >=0
  constreval     <- sr >=0
  constrevalinv  <- (-1*sr) >=0
  #indemp <- vector("numeric", n=hlast)
  indimp <- matrix(0, nrow=K, ncol=hlast)
  indimpneg <- matrix(0, nrow=K, ncol=hlast)
  for (k in 1:K){ # for shocks
    for (i in 1:hlast) { # for horizons
      temp <- impeval[,k+((i-1)*K)]
      indimp[k,i]     <- ifelse(identical(temp, constreval)==TRUE, 1, 0)
      indimpneg[k,i]  <- ifelse(identical(temp, constrevalinv)==TRUE,1,0)
    }
  }
  zeilensum <- rowSums(indimp)
  zeilensum2 <- rowSums(indimpneg)
  yesno      <- which(zeilensum==hlast) 
  yesnoneg   <- which(zeilensum2==hlast) 
  if (length(c(yesno, yesnoneg))==1) {
    if(length(yesnoneg)==1){ 
      Q[,yesnoneg] <- (-1)*Q[,yesnoneg]
      yesno <- yesnoneg
    }
    wrzreturn <- list(Q=Q, q=impdraw, smat=indimp, smat2=indimpneg, shockind=yesno, keep=1)
  } else {
    
    wrzreturn <- list(Q=Q, impulse=impeval, keep=0)
  }
  #print(zeilensum,length(c(yesno, yesnoneg)))
  #print(length(c(yesno, yesnoneg)))
  return(wrzreturn)
}    


SR_regime <- function(parainp, nkeep=500, ndraws=1000, nsubdraws=100, hf=1, hl=1,  B1, B2, covar1, covar2, 
                      restrictions, names, argnorm, Qsingle=FALSE, impnames, transcode, stdX, 
                      compound = FALSE, f_low, f_high){
  
  ### for normalization, use argnorm = numeric value, if SR are negative, ie. a negative shock is considered, argnorm should be a negative 
  # number, because the entry in B corresponding to the specific shock 
  #(e.g. impact of MP shock on ffr) will be negative and neg/neg=pos :)
  
  #environment(irf.favar) <- environment()
  
  ### set global parameters
  K <- parainp$K  #dim(cholomega1)[1]  #or sth else
  p <- parainp$p                   # as chosen for baseline STFAVAR
  H <- parainp$H                  # immer h_max + 1 wählen, wegen h=0
  #nkeep <- 500
  #nsubdraws <- 1000
  #ndraws <- 1000
  N <- length(impnames)
  inparmat1 <- B1    ## in this spec, the MC simulations (draws) are in the rows!
  inparmat2 <- B2
  inpsigma1 <- t(covar1)
  inpsigma2 <- t(covar2)
  samplerowind <- sample(nrow(inparmat1), size=ndraws,replace=T)
  inparmat1 <- inparmat1[samplerowind,]     
  inparmat2 <- inparmat2[samplerowind,]
  inpsigma1 <- inpsigma1[samplerowind,]
  inpsigma2 <- inpsigma2[samplerowind,]
  #### storage
  MCMCA1 <- vector("list", ndraws)
  MCMCA2 <- vector("list", ndraws)
  MCMCsigmad1 <- vector("list", ndraws)
  MCMCsigmad2 <- vector("list", ndraws)
  Qcholaccept1 <- vector("list", nkeep)
  Qcholaccept2 <- vector("list", nkeep)
  acceptedsubdraws1 <- array(NA, dim=c(nkeep,H,N)) # nkeep = rows, H = columns, N = Number of matrices in array
  acceptedsubdraws2 <- array(NA, dim=c(nkeep,H,N))
  subacceptvec <- vector("numeric", ndraws)
  
  accept <- 0
  for (draws in 1:ndraws){
    
    if(compound==TRUE){
      
      ## ARmats for lower quantile and higher quantile of IU and therewith transition function
      MCMCA1[[draws]]  <- (1-f_high)*matrix(inparmat1[draws,], nrow=K, ncol=K*p) + f_high*matrix(inparmat2[draws,], nrow=K, ncol=K*p)
      MCMCA2[[draws]]  <- (1-f_low)*matrix(inparmat1[draws,], nrow=K, ncol=K*p) + f_low*matrix(inparmat2[draws,], nrow=K, ncol=K*p)
      ## Covariance matrices for lower and higher quantiles 
      MCMCsigmad1[[draws]] <- (1-f_high)*matrix(inpsigma1[draws,], ncol=K, nrow=K) + f_high*matrix(inpsigma2[draws,], ncol=K, nrow=K)
      MCMCsigmad2[[draws]]  <- (1-f_low)*matrix(inpsigma1[draws,], ncol=K, nrow=K) + f_low*matrix(inpsigma2[draws,], ncol=K, nrow=K)
      cholomega1 <- t(chol(MCMCsigmad1[[draws]]))
      cholomega2 <- t(chol(MCMCsigmad2[[draws]])) 
      
    } else {
    
    #####
    MCMCA1[[draws]]  <- matrix(inparmat1[draws,], nrow=K, ncol=K*p) # convert a row containing AR parameters into a matrix
    MCMCA2[[draws]]  <- matrix(inparmat2[draws,], nrow=K, ncol=K*p)
    MCMCsigmad1[[draws]] <- matrix(inpsigma1[draws,], ncol=K, nrow=K) # convert row containing RF covar paras into a matrix
    MCMCsigmad2[[draws]] <- matrix(inpsigma2[draws,], ncol=K, nrow=K)
    cholomega1 <- t(chol(MCMCsigmad1[[draws]]))
    cholomega2 <- t(chol(MCMCsigmad2[[draws]]))
    }
    
    ## compute irfs for both regimes:
    
    # high iu regime
    irf1 <- irf.favar(armat=MCMCA1[[draws]], B=cholomega1, loads=Loadings , neg=FALSE, 
                      hor=H, transcode = transcode, impnames = impnames, stdX = stdX)
    #irf1 <- irf.favar(armat=MCMCA1[[draws]], B=cholomega1, loads=Loadings , neg=FALSE, 
    #                  hor=H, normalize = argnorm)
    rownames(irf1) <- impnames
    # low iu regime 
    irf2 <- irf.favar(armat=MCMCA2[[draws]], B=cholomega2, loads=Loadings , neg=FALSE, hor=H,
                      transcode = transcode, impnames = impnames, stdX = stdX)
    #irf2 <- irf.favar(armat=MCMCA2[[draws]], B=cholomega2, loads=Loadings , neg=FALSE,
    #                 hor=H, normalize = argnorm)
    
    rownames(irf2) <- impnames
    ### select irfs for both regimes:
    #ind <- which(impnames %in% names)  
    irfselect1 <- irf1[names, ]  # had to be changes since which does not preserve the original order of variables
    irfselect2 <- irf2[names, ]
    
    i <- 1
    subaccept <- 0
    #imp1 <- irfselect1[,1:6]
    #imp2 <- irfselect2[,1:6]
    while (i<=nsubdraws & accept < nkeep){
      # while (accept < 1){  
      W1 <- matrix(rnorm(K^2, mean = 0, sd = 1), K, K)
      W2 <- matrix(rnorm(K^2, mean = 0, sd = 1), K, K)
      if (Qsingle==TRUE){
        W1 <- matrix(rnorm(K^2, mean = 0, sd = 1), K, K)
        W2 <- W1
      }
      
      draw1 <- signcheck(a=W1, hfirst = hf, hlast = hl, sr=restrictions, impulses = irfselect1)
      draw2 <- signcheck(a=W2, hfirst = hf, hlast = hl, sr=restrictions, impulses = irfselect2)
      
      #print(draw1$keep)
      #& draw2$keep==1
      if (draw1$keep==1 & draw2$keep==1){
        accept   <- accept + 1
        print(accept)
        subaccept <- subaccept + 1
        # if (accept >= nkeep) {
        #  break
        # } else {
        Bdraw1   <- cholomega1%*%draw1$Q
        Bdraw2   <- cholomega2%*%draw2$Q
        irfdraw1 <- irf.favar(armat=MCMCA1[[draws]], B=Bdraw1, loads=Loadings , neg=F,
                              hor=H, shockpos = draw1$shockind, normalize = argnorm, 
                              transcode = transcode, impnames = impnames, stdX = stdX) 
        irfdraw1 <- do.call(cbind, lapply(seq_len(H),function(j) irfdraw1[,(draw1$shockind+((j-1)*K))]))
        rownames(irfdraw1) <- impnames
        irfdraw2 <- irf.favar(armat=MCMCA2[[draws]], B=Bdraw2, loads=Loadings , neg=F, 
                              hor=H, shockpos = draw2$shockind, normalize = argnorm,
                              transcode = transcode, impnames = impnames, stdX = stdX) 
        irfdraw2 <- do.call(cbind, lapply(seq_len(H),function(j) irfdraw2[,(draw2$shockind+((j-1)*K))]))
        rownames(irfdraw2) <- impnames
        acceptedsubdraws1[accept,,] <-  unlist(lapply(seq_len(N), function(j) irfdraw1[j,] )) 
        ### -> dim(acceptedsubdraws) = (draws x h x N), so the IRFs of the shock on each variable is stored as 
        #      a matrix with horizon in the columns and draws in the rows, ie. for draws = 1000 and H=50, there are
        #      e.g., N=148 matrices of dimension (1000 x 50) in this array!
        acceptedsubdraws2[accept,,] <-  unlist(lapply(seq_len(N), function(j) irfdraw2[j,] ))
        
        
      } 
      
      i <- i+1 
    }
    
    subacceptvec[draws] <- subaccept
    #print(subaccept)
    if (accept >= nkeep) {
      break
      ldraw <- draws
      if (ldraw < ndraws) {
        MCMCA1      <- MCMCA1[1:ldraw]
        MCMCA1      <- MCMCA1[1:ldraw]
        MCMCsigmad1 <- MCMCsigmad1[1:ldraw]
        MCMCsigmad2 <- MCMCsigmad2[1:ldraw]
        
      } }
  } ##for 1. loop  
  #  if (accept < nkeep){
  
  if (accept == 0) {
    message("\n Not enough accepted draws to proceed!")
    return(NA)
    
  } else {
    
    #}
    out <- tryCatch( 
      {
    acceptedsubdraws1 <- acceptedsubdraws1[1:accept,,] # keep only those rows which belong to successful draws!
    dimnames(acceptedsubdraws1)[[3]] <- impnames      ## doesnt work if only one (1) model has been accepted!
    acceptedsubdraws2 <- acceptedsubdraws2[1:accept,,]
    dimnames(acceptedsubdraws2)[[3]] <- impnames
      }, error = function(e) NULL
    )
    
    if(is.null(out)==TRUE){
      
      return(out)
    }
    
    returnlist <- list(acceptedsubdraws1=acceptedsubdraws1,  acceptedsubdraws2=acceptedsubdraws2, MCMCA1=MCMCA1 , MCMCA2=MCMCA1 , 
                       MCMCsigmad1=MCMCsigmad1,  MCMCsigmad2=MCMCsigmad2, accept=accept, 
                       nsubdraws=nsubdraws, ndraws=ndraws, subacceptvec=subacceptvec, restrictions = names )
    
    #return(returnlist)
  } 
  
} ## for SR_regime function 


#--------------------------------------
#   plot regime dependent IRFs: 



plot_irf_sr <- function(inp1, inp2, realnames = FALSE, varselect, H=50, regime="both"){
  # inp1 = input in form of an array, as in sr_regime. Must have named dimension (A,B, Named)
  # varselect needs to be a vector of names included in the dim-names of inp1 and inp2
  if(realnames[1]==FALSE){
    nameind <- varselect
  } else {
    nameind <- realnames 
    } 
  
  plotirf1 <- data.frame(matrix(0, nrow=H*length(varselect), ncol=6))
  colnames(plotirf1) <- c("min","median", "max", "periode", "idvar", "idvar2")
  if (regime=="both" | regime=="1"){
    plotirf1[,"periode"] <- as.numeric(rep(0:(H-1), times=length(varselect)))
      plotirf1[,"idvar"]   <- rep(nameind, each=(H))
      plotirf1[,"idvar2"]  <- factor(plotirf1[,"idvar"], levels=nameind)
    irfquant1 <- apply(inp1, c(2,3), quantile, probs=c(0.16,0.5,0.84))
    for (i in 1:3){
      plotirf1[, i]  <- as.numeric(irfquant1[i,1:H,varselect])
    }
    x <- plotirf1
  } 
  if (regime=="both" | regime=="2"){
    plotirf2 <- data.frame(matrix(0, nrow=H*length(varselect), ncol=6))
    colnames(plotirf2) <- c("min","median", "max", "periode", "idvar", "idvar2")
    plotirf2[,"periode"] <- as.numeric(rep(0:(H-1), times=length(varselect)))
     plotirf2[,"idvar"]   <- rep(nameind, each=(H))
     plotirf2[,"idvar2"]  <- factor(plotirf2[,"idvar"], levels=nameind)
    irfquant2 <- apply(inp2, c(2,3), quantile, probs=c(0.16,0.5,0.84))
    for (i in 1:3){
      plotirf2[, i]  <- as.numeric(irfquant2[i,1:H,varselect])
    }
    x <- plotirf2
  }
  ### quantiles and median: 
  
  #irfquant1 <- apply(inp1, c(2,3), quantile, probs=c(0.16,0.5,0.84))
  #irfquant2 <- apply(inp2, c(2,3), quantile, probs=c(0.16,0.5,0.84))
  #for (i in 1:3){
  #  plotirf1[, i]  <- as.numeric(irfquant1[i,,varselect])
  #  plotirf2[, i]  <- as.numeric(irfquant2[i,,varselect])
  #  }
  # Quantile finden
  # Median als Punktschätzer?
  # Nur ein shock, d.h, beschriftung können variablennamen sein
  # 
  
  gg <- (ggplot(x, aes(x=periode, median, group=1)) +scale_x_continuous( breaks=c(seq(0,H, 10)), labels=c(seq(0,H, 10)), expand = c(0, 0), name="") 
         +scale_y_continuous(  name="") 
         #+theme(panel.background = element_rect(fill="white"))
         +theme_bw() 
         +theme(panel.background = element_rect(fill="white"),
                axis.title.x=element_blank()
                , panel.grid.major=element_line(size=0.3, color="#cccccc")
                , panel.grid.minor=element_line(size=0.5, color="#e6e6e6"))
         +theme(axis.text.y = element_text(colour="black", size=7, face="bold"), 
                axis.text.x = element_text(colour="black", size=8, face="bold"),
                strip.text=element_text(size=10, lineheight=0.5))
         +facet_wrap(~idvar2, scales="free"))
  #nrow=length(shocks)
  if (regime=="both"){
    for (i in nameind){
      # + geom_point(data=subset(dfp, idvar==d[i]))
      #
      gg <- gg + list(geom_line(data=subset(plotirf1, idvar==i), aes(y=median, group=1), linetype=1),      
                      geom_line(data=subset(plotirf2, idvar==i), aes(y=median, group=1), linetype=2, color="red")) 
      # red is low IU
      # blue is high IU 
      
    } 
    
    for (i in nameind){
      
      gg <- (gg +list(geom_ribbon(data=subset(plotirf1, idvar==i),
                                  aes(ymin=min, ymax=max, group=1), alpha=0.6,  fill="#74C2E1"),
                      geom_ribbon(data=subset(plotirf2, idvar==i),
                                  aes(ymin=min, ymax=max, group=1), alpha=0.6,  fill="grey") ))
    }
  }
  
  ##### for regime 1 or 2 only
  if (regime%in%c("1","2")){
    for (i in nameind){
      # + geom_point(data=subset(dfp, idvar==d[i]))
      #
      gg <- gg + list(geom_line(data=subset(x, idvar==i), aes(y=median, group=1), linetype=1))
      
      #gg <- (gg  + geom_line(data=subset(irfplot1, idvar==i, ...))  + geom_hline(yintercept=0, colour="blue"))
    } 
    
    for (i in nameind){
      
      gg <- (gg +list(geom_ribbon(data=subset(x, idvar==i),
                                  aes(ymin=min, ymax=max, group=1), alpha=0.6,  fill="#74C2E1") ))
    }
  }
  return(gg)
  
}  ## for plot_function_sr, main



plot_irf_diff <- function(inp1, inp2, realnames = FALSE, varselect, arraynames = list(NULL, NULL, all_names_sr) , H=50){
  # inp1 = input in form of an array, as in sr_regime. Must have named dimension (A,B, Named)
  # varselect needs to be a vector of names included in the dim-names of inp1 and inp2
  if(realnames[1]==FALSE){
    nameind <- varselect
  } else {
    nameind <- realnames 
  } 
  
  inpdiff <- array(dim = dim(inp1), dimnames = arraynames)
  
  for (i in 1:dim(inp1)[3]){
    
    inpdiff[,,i] <- inp2[,,i] - inp1[,,i]
    
  }
  
  plotirf1 <- data.frame(matrix(0, nrow=H*length(varselect), ncol=6))
  colnames(plotirf1) <- c("min","median", "max", "periode", "idvar", "idvar2")
  plotirf1[,"periode"] <- as.numeric(rep(0:(H-1), times=length(varselect)))
  plotirf1[,"idvar"]   <- rep(nameind, each=(H))
  plotirf1[,"idvar2"]  <- factor(plotirf1[,"idvar"], levels=nameind)
  irfquant1 <- apply(inpdiff, c(2,3), quantile, probs=c(0.16,0.5,0.84))
  for (i in 1:3){
    plotirf1[, i]  <- as.numeric(irfquant1[i,1:H,varselect])
  }
  x <- plotirf1
  
  
  
  gg <- (ggplot(x, aes(x=periode, median, group=1)) +scale_x_continuous( breaks=c(seq(0,H, 10)), labels=c(seq(0,H, 10)), expand = c(0, 0), name="") 
         +scale_y_continuous(  name="") 
         #+theme(panel.background = element_rect(fill="white"))
         +theme_bw() 
         +theme(panel.background = element_rect(fill="white"),
                axis.title.x=element_blank()
                , panel.grid.major=element_line(size=0.3, color="#cccccc")
                , panel.grid.minor=element_line(size=0.5, color="#e6e6e6"))
         +theme(axis.text.y = element_text(colour="black", size=7, face="bold"), 
                axis.text.x = element_text(colour="black", size=8, face="bold"),
                strip.text=element_text(size=10, lineheight=0.5))
         +facet_wrap(~idvar2, scales="free"))
  
  
  for (i in nameind){
    # + geom_point(data=subset(dfp, idvar==d[i]))
    #
    gg <- gg + list(geom_line(data=subset(x, idvar==i), aes(y=median, group=1), linetype=1))
    
    #gg <- (gg  + geom_line(data=subset(irfplot1, idvar==i, ...))  + geom_hline(yintercept=0, colour="blue"))
  } 
  
  for (i in nameind){
    
    gg <- (gg +list(geom_ribbon(data=subset(x, idvar==i),
                                aes(ymin=min, ymax=max, group=1), alpha=0.6,  fill="#74C2E1") ))
  }
  
  return(gg)
  
}  ## for plot_function_sr, main

########################################
# Lag selection (AIC + BIC):
p_opt <- function(inp, pmax, det_term="trend"){
  
  #inp: data as (K x T) 
  selecA <- matrix(nrow=pmax, ncol=2)
  selecS <- matrix(nrow=pmax, ncol=2)
  y <- inp[,(pmax+1):ncol(inp)]
  
  
  for (p in 1:pmax) {                                    #### main loop for calculating IC
    
    z <- inp[,pmax:(ncol(inp)-1)]
    
    i=2                                                    #### loop for adjusting reg-mat according to ar-order tested...
    while (i <= p) {
      
      z <- rbind(z, inp[, (pmax-i+1):(ncol(inp)-i)])           
      
      i=i+1
    }
    
    
    znone <- z
    
    zconst <-  rbind(as.vector(rep(1, ncol(z))), z)
    
    zlinear <- rbind(as.vector(rep(1, ncol(z))), seq(1,ncol(z),1), z)
    
    zquad <- rbind(as.vector(rep(1, ncol(z))), seq(1,ncol(z),1), seq(1, ncol(z),1)^2, z)
    
    detterms <- list(none=znone, constant=zconst, trend=zlinear, quad=zquad)
    
    z <- detterms[[det_term]]             #### choose det terms from function argument
    
    
    arhat<-(y%*%t(z))%*%solve(z%*%t(z))                    ## in armat K^2 + 1 column, due to constant, the first column contains constant estimates!!!
    uhat = y - arhat%*%z
    sighat = uhat%*%t(uhat)/ncol(uhat)                     ##### estimation of covar-mat
    
    AIC <- log(det(sighat))+2*ncol(arhat)*nrow(arhat)/ncol(uhat)
    
    SC <- log(det(sighat))+log(ncol(uhat))*ncol(arhat)*nrow(arhat)/ncol(uhat)
    
    
    selecA[p,1] <- AIC
    selecA[p,2] <- p
    
    selecS[p,1] <- SC
    selecS[p,2] <- p
  }
  
  selecS <- selecS[order(selecS[,1]),]
  selecA <- selecA[order(selecA[,1]),]
  
  selecG <- matrix(nrow=1, ncol=2)
  
  selecG[1,1] <- selecA[1,2]
  selecG[1,2] <- selecS[1,2]
  colnames(selecG) <- c("AIC", "SC")
  
  
  
  out <- list(selecG, selecS, selecA)
  
  names(out) <- c("seclecG", "selecS", "selecA")
  
  return(out)
  
}


# select det. terms by IC:

detterm.fun <- function(inp, p){
  
  selecA <- as.data.frame(matrix(nrow=4, ncol=2))
  selecS <- as.data.frame(matrix(nrow=4, ncol=2))
  #inp <- inp
  
  y <- inp[,(p+1):ncol(inp)]
  
  z <- inp[,p:(ncol(inp)-1)]
  
  i=2                                                    #### loop for adjusting reg-mat according to ar-order tested...
  while (i <= p) {
    
    z <- rbind(z, inp[, (p-i+1):(ncol(inp)-i)])           
    
    i=i+1
  }                                   #### main loop for calculating IC
  
  znone <- z
  
  zconst <-  rbind(as.vector(rep(1, ncol(z))), z)
  
  zlinear <- rbind(as.vector(rep(1, ncol(z))), seq(1,ncol(z),1), z)
  
  zquad <- rbind(as.vector(rep(1, ncol(z))), seq(1,ncol(z),1), seq(1, ncol(z),1)^2, z)
  
  detterms <- list(none=znone, constant=zconst, trend=zlinear, quad=zquad)
  
  for (i in 1:4) {
    
    z <- detterms[[i]]             #### choose det terms from function argument
    
    
    arhat<-(y%*%t(z))%*%solve(z%*%t(z))                    ## in armat K^2 + 1 column, due to constant, the first column contains constant estimates!!!
    uhat = y - arhat%*%z
    sighat = uhat%*%t(uhat)/ncol(uhat)                     ##### estimation of covar-mat
    
    AIC <- log(det(sighat))+2*ncol(arhat)*nrow(arhat)/ncol(uhat)
    
    BIC <- log(det(sighat))+log(ncol(uhat))*ncol(arhat)*nrow(arhat)/ncol(uhat)
    
    
    selecA[i,1] <- AIC
    
    selecS[i,1] <- BIC
    
  }
  
  
  selecA[,2] <- c("none", "constant", "trend", "quad")
  selecS[,2] <- c("none", "constant", "trend", "quad")
  
  selecS <- selecS[order(selecS[,1]),]
  selecA <- selecA[order(selecA[,1]),]
  
  selecG <- matrix(nrow=1, ncol=2)
  
  selecG[1,1] <- selecA[1,2]
  selecG[1,2] <- selecS[1,2]
  colnames(selecG) <- c("AIC", "BIC")
  
  
  
  out <- list(selecG, selecS, selecA)
  
  names(out) <- c("min", "BIC", "AIC")
  
  return(out)
  
}  


## test stability condition: 

stability.test <- function(inp, p,K){
  #inp <- beta0Est
  #p <- as.numeric(specs["p",2])
  Ap <- inp
  
  Kp <- K*p
  
  if(p==1){
    
    A <- Ap
    
  } else {
    
    A <- rbind(Ap, diag(Kp)[1:(Kp-K),])  
    
  }
  
  root <- Mod(eigen(A)$values)
  
  return(root)
}

###############################
# olssvd computes the OLS estimator (Jean Boivin (11/18/01))
olssvd <- function(y, z){ 
  
  # Input: Y and Z matrices
  # Output: Esimtators
  
  y <- t(y)
  z <- t(z)
  
  # for this the input must be data(KxT)
  result <- (y%*%t(z))%*%solve(z%*%t(z)) 
  return(t(result))
}


# function that computes a matrix containing all lags selected
lag_fun <- function(data, nlags){
  
  # Input:  data: T x k data matrix
  #         nlags: number of lags
  
  # Output: matrix containing all n lags in its columns including lag zero
  
  ny <- nlags + 1 # number of lags plus t = 0
  nv <- dim(data)[2] # number of variables
  nt <- dim(data)[1] # number of time observations
  
  lag_array <- array(0, dim = c((nt - nlags), nv, ny)) 
  
  for(j in 0:nlags){
    lag_mat <- matrix(0, ncol = nv, nrow = (nt - nlags))
    lag_mat <- data[(ny-j):(nt-j),]
    
    lag_array[,,j+1] <- lag_mat
  }
  all_lags <- matrix(lag_array, nrow = (nt - nlags), ncol = ny * nv)
  
  return(all_lags)
}


###############################################################
#### factor extraction: 

pcfac <- function(data, anzahlfac){
  
  ## PC factor extraction
  ## data must be (TxN)  
  nvars <- dim(data)[2]
  
  xx <- t(data)%*%data# X'X is proportional to covariance matrix
  eigen_x <- eigen(xx) 
  
  # applying the factor restricition
  lam <- sqrt(nvars)*eigen_x[[2]][,c(1:anzahlfac)] # lam'lam/n = I. loadings are, where n = number of time series in x
  factors <- data%*%lam/nvars     
  
  return(list(factors, lam))
}

fac.fun.bgh <- function(data = x_sd, beob, nfac = 5, niter=10, ctol = 10^(-8), standard = TRUE, yindex=c(1)) {
  # wenn nicht standardisiert gibt es probleme (negatives R^2) mit der Berechnung des marginalen R^2, da die Residuen nicht
  # mittelwert null haben... Denn in die OLS schätzung wurde ohne Konstante durchgeführt.
  # data has to be of the form: TxN (T = # of observations, N = # of variables)
  xvars <- dim(x_sd)[2]
  TL <- dim(x_sd)[1]
  olsfac <- function(y, x){ 
    
    # Input: Y and Z matrices
    # Output: Esimtators
    
    # for this the input must be data(KxT)
    result <- solve(crossprod(x),crossprod(x,y))
    return(result)
  }
  pcfac <- function(data, anzahlfac){
    ## data must be (TxN)  
    nvars <- dim(data)[2]
    
    xx <- t(data)%*%data# X'X is proportional to covariance matrix
    eigen_x <- eigen(xx) 
    
    # applying the factor restricition
    lam <- sqrt(xvars)*eigen_x[[2]][,c(1:anzahlfac)] # lam'lam/n = I. loadings are, where n = number of time series in x
    factors <- data%*%lam/xvars     
    
    return(list(factors, lam))
  }
  ifelse(is.vector(beob)==TRUE, r <- 1, r <- dim(beob)[2])
  #r <- 1 # number of observable factors.
  conv_tol <- ctol*dim(data)[1]*dim(data)[2]
  
  #******************************#
  
  temp <- pcfac(data=data, k_numb)
  factors <- temp[[1]]
  lam <- temp[[2]]
  
  # applying the factor restricition
  
  ifelse(standard==TRUE, ob <- scale(beob), ob <- beob)
  
  designmat <- cbind(ob,rep(1,TL), factors)
  
  loadtemp <- t(olsfac(data,designmat)) # dim = (NxK), K number of all factors
  #View(loadtemp)
  
  loadobservable <- loadtemp[,1:r] # (Nxr)
  
  eobs <- data - ob%*%t(loadobservable)
  ### first estimates needed for convergence evaluation:
  #e <- Xx - designmat%*%t(loadtemp)
  SSRold <- 0 # sum(diag(t(e)%*%e))
  factorsold <- factors              
  lamold <- lam   
  
  conv <- 100
  counter <- 1
  while (conv > conv_tol & counter<=niter){
   # while (counter < 11){  
    
    
    #xx <- t(eobs)%*%eobs# X'X is proportional to covariance matrix
    #eigen_x <- eigen(xx) # get the eigenvalues  --> which are already in decreasing order
    #eigen_x[[1]]
    # applying the factor restricition
    #lam <- sqrt(xvars)*eigen_x[[2]][,c(1:k_numb)] # lam'lam/n = I. loadings are, where n = number of time series in x
    #factors <- eobs%*%lam/xvars 
    temp <- pcfac(data=eobs, nfac)
    factors <- temp[[1]]
    
    #designmat <- cbind(ob, factors)  # dim = (TxK)
    
    designmat <- cbind(ob,rep(1,TL), factors) # dim = (Nxr)
    #loadtemp <- loadtemp[,-(r+1)]
    loadtemp <- t(olsfac(data,designmat)) # dim = (NxK), K number of all factors
    
    loadobservable <- loadtemp[,1:r]
    
    eobs <- data - ob%*%t(loadobservable)  #dim(eobs)=TxN
    e <- data - designmat%*%t(loadtemp)
    SSR <- sum(diag(t(e)%*%e))
    
    #ifelse(abs(SSR-SSRold) <  conv_tol | , conv <- 1, conv <- 0)  # new conv criterion: |SSRnew - SSRold|
    conv <- abs(SSRold - SSR)
    print(conv)
    SSRold <- SSR
    #conv <- abs(sum(colSums(lam-lamold)))       # absolute difference of SSRnew and SSRold
    #conv <- abs(sum(colSums(lam-lamold))/sum(colSums(lamold))) # old relative convergence criterion
    print(cor(designmat[,-(r+1)])[,1:2])
    if (conv < conv_tol){
      #lamlast <- lam 
    } else {
      #lamold <- lam 
    }
    #factorsold <- factors
    counter <- counter+1
  }
  
  fac_est <- designmat[,-(r+1)] 
  #load_est <- loadtemp
  load_est <- t(olsfac(data,fac_est)) # new according to BGM
  Loadings <- rbind(diag(var_numb)[1:r,], load_est)
  
  ####### calculate R^2 for each series + overall: 
  #dim(fac_est)
  #dim(load_est)
  residtot <- data - fac_est%*%t(load_est) # (TxK) (KxN)
  ssrvec <- diag(t(residtot)%*%residtot)
  sstvec <- diag(t(data)%*%data)
  
  #var(Xx[,1:10])
  
  R2vec  <- 1-ssrvec/sstvec  # communality for each series
  names(R2vec) <- all_names[-yindex]
  R2vecord <- R2vec[order(R2vec, decreasing = TRUE)]
  #View(R2vecord)
  names(R2vec) <- all_names[-yindex]
  R2tot  <- 1-sum(ssrvec)/sum(sstvec) # R2 total
  
  ###### compute R^2 for each factor for all series and for each series: 
  
  sfacR2 <- matrix(NA, nrow=dim(fac_est)[2], ncol=2)
  sfacR2[,1] <- c(paste(all_names[yindex],sep=""), paste("factor"," ",(r+1):dim(fac_est)[2], sep=""))
  for (i in 1:ncol(fac_est)){
    ee     <- data - fac_est[,i]%*%t(load_est[,i])
    ssrone <- sum(diag(t(ee)%*%ee))
    sfacR2[i,2] <- 1 - ssrone/sum(sstvec)
  }
  
  return(list(iteration=counter, factors = fac_est, loadings = load_est, Loadings = Loadings,
              R2tot = R2tot, R2vec = R2vec, R2vecord = R2vecord, tablemR2 = sfacR2  ))
  
} # for function

### two step factor estimation, following Hwang_2009 (economics letters)

twostepfac.fun <- function(data = x_sd, beob, numfac, standard = TRUE, yindex=c(1)) {
  
  # data has to be of the form: TxN (T = # of observations, N = # of variables)
  ifelse(is.vector(beob)==TRUE, r <- 1, r <- dim(beob)[2])
  ifelse(standard==TRUE, ob <- scale(beob), ob <- beob)
  B <- t(olssvd(data, ob)) # dim = (NxK), K number of all factors
  
  U <- data - ob%*%t(B)  # U is TxN
  
  #dim(U)
  
  temp <- pcfac(data=U, numfac)
  factors <- temp[[1]]
  load <- temp[[2]]
  
  fac_est <- cbind(ob, factors)
  load_est <- cbind(B,load)                        # dim F (TxK), loadings = (NxK)
  Loadings <- rbind(diag(ncol(fac_est))[1:r,], load_est)  
  corfac <- round(cor(fac_est), digits = 4)
  
  
  
  
  return(list(factors = fac_est, loadings = load_est, Loadings = Loadings, corr = corfac, B=B ))
  
} # for function


# Factor rotation
# Bernanke, Boivin and Eliasz (2005) differenciate between slow and fast moving
# variables and then extract the factors using this information.
facrot <- function(Fest, Ffast, Fslow){ 
  
  # Input: Fest: unrestricted PC estimates
  #        Ffast: fast moving factors, monetary policy
  #        Fslow = slow moving factors
  
  # Output: Rotated factor estimates
  
  k1 <- dim(as.matrix(Ffast))[2]
  b <- olssvd(Fest, cbind(rep(1, nobs), Ffast, Fslow))
  Fr <- Fest - as.matrix(Ffast)%*%b[2:(k1+1),] # new estimate for factors
  return(Fr)
}


# Function returning the variables with the highest factor loadings
lmr2 <- function(fac_num){
  
  # Input: number of the factor which should be analyzed
  # Output: names and loadings for which this factor loads highest
  
  orderings <- order(abs(loading_int[,fac_num]), decreasing = T)
  
  r2 <- rep(NA, 10)
  for(i in 1:10){
    r2[i] <- summary(lm(x_sd[,orderings[i]] ~ rot_fac[,fac_num]))$r.squared
  }
  
  r2names <- x_names[orderings][1:10]
  return(cbind(r2names, r2))
}


# Automatization of factor extraction
fa_estimate <- function(fac_numb){
  
  # Input: Number of Factors
  # Output: Factor estimates and error variance
  
  xx <- t(x_sd)%*%x_sd 
  eigen_x <- eigen(xx) # eigenvalues
  
  # lambda and factors
  lam <- sqrt(xvars)*eigen_x[[2]][,c(1:fac_numb)] 
  factors <- x_sd%*%lam/xvars 
  
  # extract principal components from slow moving factors
  slowindex <- which(adj_speed==1)   
  x_mean <- (x_data -  mean_mat)     
  slow_fac <- x_mean[,c(slowindex)]
  xsvars <- dim(slow_fac)[2]
  
  xx_s <- t(slow_fac) %*% slow_fac 
  eigen_xs <- eigen(xx_s)
  lam_s <- sqrt(xsvars) * eigen_xs[[2]][,c(1:fac_numb)]  
  factors_s <- slow_fac%*%lam_s/xsvars
  
  # rotate factors. for factor rotation ffr has to be standardized
  std_ffr <- (ffr- mean(ffr))/sd(ffr)
  rot_fac <- facrot(factors, std_ffr, factors_s) 
  
  # calculate the loadings matrix
  XY <- cbind(x_sd, obser)
  FY <- cbind(rot_fac, obser) 
  loading = t(olssvd(XY,FY))
  
  e <- XY - FY %*% t(loading)
  SIGMA1 <- t(e) %*% e / nobs
  SIGMA <- diag(diag(SIGMA1))
  SIGMA_BAI1 <- sum(SIGMA)/(nobs*xvars)
  return(list(FY, SIGMA_BAI1))
}



#####################################################################
#####################################################################


# Information Criterion Bai & Ng
IC <- function(fac_numb){
  
  # Input: Number of factors to be included in factor extraction
  # Output: Value of the Bai & Ng Criterion
  
  errors <- fa_estimate(fac_numb)[[2]]
  IC_value <- log(errors) + fac_numb*((xvars+nobs)/(xvars*nobs))*(log(min(xvars, nobs)))
  return(IC_value)
}


# Vectorize Covariance matrices
vech <- function(covm){
  
  # puts covariance matrix into vector
  
  vect <- list()
  
  for(i in 1:var_numb){
    vect[[i]] <- covm[c(i:var_numb), i]
    vectors <- unlist(vect)
  }
  return(vectors)
}


# Unvectorize Estimator vector
unvech <- function(raws){
  
  # Input: raws - vector of guesses containing lower triangular entries
  # of both omegas and the guess for gamma.
  
  # Output: list containing [[1]] the guess for gamma, [[2]] omega_l
  # matric and [[3]] omega_h matrix
  
  nr <- raws[c(1:(length(raws)-1))]
  raws1 <- nr[c(1:(length(nr)/2))]
  raws2 <- nr[c(((length(nr)/2)+1):length(nr))]
  
  covm1 <- matrix(0, ncol = var_numb, nrow = var_numb)
  covm1[lower.tri(covm1, diag = T)] <- raws1
  
  covm2 <- matrix(0, ncol = var_numb, nrow = var_numb)
  covm2[lower.tri(covm2, diag = T)] <- raws2
  
  omega0 <- covm1 %*% t(covm1)
  omega1 <- covm2 %*% t(covm2)
  
  return(list(raws[length(raws)], omega0, omega1))
}


# function used to unvectorize the parameter guess vector and thereby selecting
# one of the regime dependend covariance matrix and storing them in a matrix
unvech1 <- function(guess, l){
  unvech(guess)[[l]]  
}


# Put vector of the pi's into matrix
unvec <- function(pi_vec, n_est, n_var){
  
  # Input: pi_vec - verctor containing the entries of both pi's
  #        n_est: total number of estimates
  #        n_var: number of variables
  
  # Output: matrix with pi estimators
  
  v <- matrix(NA, n_est, n_var)
  for (i in 1:n_var){
    v[,i] <- pi_vec[((i-1)*n_est+1):(i*n_est)]
  }
  return(v)
}


# Maximum lieklihood for starting values
first_log <- function(param){
  
  # returns log-likelihood values for each observation. This function will be 
  # maximized using maxBHHH
  # Input are the guesses for theta and the guesses for the lower triangular
  # matrix of the cholesky covariance matrices
  
  logL <- vector("numeric")
  
  theta <- param[length(param)]
  om1 <- unvech(param)[[2]]
  om2 <- unvech(param)[[3]]
  N <- length(z_lag)
  
  if(theta < 0){
    
    if(min(eigen(om1)$values) > 0 && min(eigen(om2)$values) >0){
      
      # loglikelihood
      M_Z_t <- exp(theta*z_lag)/(1 + exp(theta*z_lag))
      
      for(i in 1:N){
        omt <- om1 * (1 - M_Z_t[i]) + om2 * M_Z_t[i]
        e_t <- app_res[i,]
        logL[i] <- - 0.5 *((log(det(omt))) + (e_t %*% solve(omt) %*% e_t))
      }
      
    }else{
      logL <- rep(-1e12, N)
    }
  }else{
    logL <- rep(-1e12,N)
  }
 return( sum((-1)*logL))
}    ### --> negative log-lik now (and sumed...)


# Return negative log likelihood values. used to minimize the loglikelihood
# function in the mcmc algorithm
second_log <- function(param){
  
  
  # Input: Parameter guesses for theta and cholesky omegas
  # Output: log likelihood value
  # This time the negative log likelihood should be minimized!
  
  guesses <- unvech(param)
  
  theta <- guesses[[1]]
  om1 <- guesses[[2]]
  om2 <- guesses[[3]]
  
  if(theta<0){
    if(min(eigen(om1)$values) > 0 && min(eigen(om2)$values) >0){
      # curvature of the transition function should be the same for 
      # contemporaneous response and dynamic response
      gamma <- theta
      
      # transition equation
      M_Z_t <- exp(theta*Z0)/(1 + exp(theta*Z0))
      F_Z_t <- exp(gamma*Z0)/(1 + exp(gamma*Z0))
      
      low <- matrix(rep((1 - F_Z_t), var_numb*nlags), ncol = var_numb*nlags) ## wurde vertauscht: wenn z_t -> infinity --> F_zt -> 0
      high <- matrix(rep(F_Z_t, var_numb*nlags), ncol = var_numb*nlags)    ## daher ist 1-F_z_t HIGH 
      
      XM <- cbind(X*low, X*high, rep(1, n_obs), ex)
      
      wgtM <- list()
      dims <- c(var_numb*dim(XM)[2], var_numb*dim(XM)[2], N)
      krons <- matrix(0, dims[1], dims[2])
      kronsum <- krons
      
      vecs <- rep(0, dims[1])
      vecsum <- vecs
      
      for(i in 1:N){ 
        wgtM <- solve(om1 * (1-M_Z_t[i]) + om2 * M_Z_t[i])
        krons <- kronecker(wgtM, (XM[i,])%*% t(XM[i,]))
        kronsum <- kronsum + krons
        vecs <- as.vector(XM[i,] %*% t(Y0[i,]) %*% wgtM)
        vecsum <- vecsum + vecs
      }
      
      pis <- solve(kronsum)%*%vecsum
      pi_mat <- unvec(pis, dim(X)[2]*2+1+Nforecs, var_numb)
      
      pi_l <- t(pi_mat[c(1:dim(X)[2]),])
      pi_h <- t(pi_mat[c((dim(X)[2] + 1):(2*dim(X)[2])),])
      const <- t(pi_mat[(2*dim(X)[2]+1),])
      others <- t(pi_mat[c((2*dim(X)[2]+2):(dim(pi_mat)[1])),])
      
      residsM <- Y0 - XM%*%pi_mat
      
      # loglikelihoodvalue
      log_value <- 0
      
      for (i in 1: N){
        wgtM <- om1 * (1-M_Z_t[i]) + om2 * M_Z_t[i]
        log_value <- log_value-0.5*((log(det(wgtM))) + 
                                      (residsM[i,] %*%solve(wgtM) %*% residsM[i,]))
      }
      
    }else{
      log_value <- -1e12
    }}else{
      log_value <- -1e12
    }
  
  return(-log_value)
}


# Returns parameters like log likelihood value, pi values and residuals
# for a given guess for the omegas
return_paras <- function(param){
  
  # Input: Parameter guesses for theta and cholesky omegas
  
  # Output: log likelihood value
  # With this function, calculated values for the given guesses for theta
  # and the cholesky omegas will be returned
  
  guesses <- unvech(param)
  
  theta <- guesses[[1]] 
  om1 <- guesses[[2]]
  om2 <- guesses[[3]]
  
  gamma <- theta 
  
  M_Z_t <- exp(theta*Z0)/(1 + exp(theta*Z0)) 
  F_Z_t <- exp(gamma*Z0)/(1 + exp(gamma*Z0)) 
  
  low <- matrix(rep((1 - F_Z_t), var_numb*nlags), ncol = var_numb*nlags)
  high <- matrix(rep(F_Z_t, var_numb*nlags), ncol = var_numb*nlags)
  
  XM <- cbind(X*low, X*high, rep(1, n_obs), ex)
  
  wgtM <- list()
  dims <- c(var_numb*dim(XM)[2], var_numb*dim(XM)[2], N)
  krons <- matrix(0, dims[1], dims[2])
  kronsum <- krons
  
  vecs <- rep(0, dims[1])
  vecsum <- vecs
  
  for(i in 1:N){ 
    wgtM <- solve(om1 * (1-M_Z_t[i]) + om2 * M_Z_t[i])
    krons <- kronecker(wgtM, (XM[i,])%*% t(XM[i,]))
    kronsum <- kronsum + krons
    vecs <- as.vector(XM[i,] %*% t(Y0[i,]) %*% wgtM)
    vecsum <- vecsum + vecs
  }
  
  pis <- solve(kronsum)%*%vecsum
  pi_mat <- unvec(pis, dim(X)[2]*2+1+Nforecs, var_numb)
  
  pi_l <- t(pi_mat[c(1:dim(X)[2]),])
  pi_h <- t(pi_mat[c((dim(X)[2] + 1):(2*dim(X)[2])),])
  const <- t(pi_mat[(2*dim(X)[2]+1),])
  others <- t(pi_mat[c((2*dim(X)[2]+2):(dim(pi_mat)[1])),])
  
  residsM <- Y0 - XM%*%pi_mat
  
  # loglikelihoodvalue
  log_value <- 0
  
  for (i in 1: N){
    wgtM <- om1 * (1-M_Z_t[i]) + om2 * M_Z_t[i]
    log_value <- log_value-0.5*((log(det(wgtM))) + 
                                  (residsM[i,] %*%solve(wgtM) %*% residsM[i,]))
  }
  output <-  list(pi_l, pi_h, om1, om2, theta, gamma, residsM, const, others, -log_value)
  names(output) <- c("Pi_low", "Pi_high", "Omega_low", "Omega_high", "Theta",
                     "Gamma", "Residuals", "Constants", "Other_Slopes", "Log_Likelihood_Value")
  return(output)
}


# Function to plot the values of the draws -> Trace Plots
plots <- function(mat, columns){
  par(mfrow = c(2,2))
  for(j in 1:4){
    plot(c(1:iters), mat[,columns[j]], type = "l")
    
  }
}


# Impulse response function and Confidence Invervals
irf_fa_CIsmall <- function(betamat, betaNLE, Omegamat, OmegaNL, loads, regime, shockpos, neg = FALSE){
  
  # function to calculate the impulse responses considering factor 
  # augmentation:
  
  # Input: betamat: Array of all the accepted beta draws in the MCMC
  #        betaNLE: Estimated VAR coefficients in matrix form
  #        Omegamat: Array of all the excepted Omega-guesses from the MCMC
  #        OmegaNL: Covariance matrix of estimated error terms
  #        loads: loadings of the factor augmentation
  #        regime: Set regime for inflation uncertainty: 1 = low iu, 2 = high iu
  #        shockpos: The position of the shockvariable in the VAR equation
  
  # Output: irf_imp: impulse response impact matrix
  #         CIup: upper confidence interval
  #         CIlow: lower confidence interval
  #         IRFse: standard error for confidence interval
  #         IRF_CImean: mean value of irf across draws
  
  
  # number of draws to compute confidence bands and standard error for the irfs
  CIiters <- 500
  
  # confidence intervals
  up_percentile <- 95
  low_percentile <- 5
  
  reg_cov <- regime +1
  
  # companion form
  betaNL <- rbind(betaNLE, 
                  cbind(diag(var_numb * (nlags-1)), 
                        matrix(0, nrow = (var_numb * (nlags-1)), 
                               ncol = (dim(betaNLE)[2]-(var_numb * (nlags-1)))))) 
  
  
  CholOmegaNL <- t(chol(OmegaNL))
  shockvec = CholOmegaNL[,shockpos]/CholOmegaNL[shockpos, shockpos]
  
  if(neg == TRUE){
    shockvec <- -shockvec
  }
  
  # irf with a shock
  irf_lin <- c(shockvec, rep(0, var_numb*(nlags-1)))
  irf_imp <- irf_lin
  
  for(i in 2:irf_hor){
    irf_lin <- betaNL %*% irf_lin
    irf_imp <- cbind(irf_imp, irf_lin)
  }
  
  # point estimate of the ir
  irf_Dimo <- loads %*% irf_imp[1:var_numb,]
  
  
  # confidence intervals for the irf using MCMC draws
  irf_Dimp0mat <- array(NA, dim =c(CIiters, dim(irf_imp)))
  
  beta_NL_vec <- as.vector(betaNLE)
  set.seed(234)
  
  OmegaVARLIN <- OmegaNL
  Dn <- duplication.matrix(var_numb)
  Dplus <- (solve(t(Dn)%*% Dn)) %*% t(Dn)
  OmegaVARLIN_var <- (2*(Dplus %*% (kronecker(OmegaNL, OmegaNL)) %*% t(Dplus)))/N
  
  for(i in 1:CIiters){
    
    # sample each step a random position from the available draws
    position <- max(1, min(round(runif(1)*dim(betamat)[1]), dim(betamat)[1])) # replace beta1mat with function argument
    
    # dynamic responses draw
    beta_draw0 <- betamat[position,]
    beta_draw <- unvec(beta_draw0, var_numb, var_numb*nlags)
    
    # covariance matrix draw
    OmegaNLdraw2 <- vech(OmegaVARLIN) + t(chol(OmegaVARLIN_var)) %*% rnorm(length(vech(OmegaVARLIN)))
    unvech_om3 <- matrix(0, ncol = var_numb, nrow = var_numb)
    unvech_om3[lower.tri(unvech_om3, diag = T)] <- OmegaNLdraw2
    
    for(k in 1:var_numb){
      for(h in k:var_numb){
        unvech_om3[k,h] <- unvech_om3[h,k]
      }
    }
    OmegaNLdraw3 <- unvech_om3
    
    CholOmegaNL0 <- t(chol(OmegaNLdraw3))
    shockvec0 <- CholOmegaNL0[,shockpos]/CholOmegaNL0[shockpos,shockpos]
    
    if(neg == TRUE){
      shockvec0 <- -shockvec0
    }
    
    betaNLE <- beta_draw
    betaNL <- rbind(betaNLE,  # companion form 
                    cbind(diag(var_numb * (nlags-1)), 
                          matrix(0, nrow = (var_numb * (nlags-1)), 
                                 ncol = (dim(betaNLE)[2]-(var_numb * (nlags-1)))))) 
    
    irf_lin_draw <- c(shockvec0, rep(0, var_numb*(nlags-1)))
    irf_imp_draw <- irf_lin_draw
    
    for(j in 2:irf_hor){
      irf_lin_draw <- betaNL %*% irf_lin_draw
      irf_imp_draw <- cbind(irf_imp_draw, irf_lin_draw)
    }
    
    irf_Dimp0mat[i,,] <- irf_imp_draw
    
  }
  
  irf_Dimp0mat1 <- array(NA, c(CIiters, var_numb_FA, irf_hor))
  
  for(k in 1: CIiters){
    irf_Dimp0mat1[k,,] <- loads %*% irf_Dimp0mat[k,c(1:var_numb),] 
  }
  
  
  # compute confidence intervals
  
  CIup <- matrix(NA, nrow = var_numb_FA, ncol = irf_hor)
  CIlow <- matrix(NA, nrow = var_numb_FA, ncol = irf_hor)
  IRFse <- matrix(NA, nrow = var_numb_FA, ncol = irf_hor)
  IRF_CImean <- matrix(NA, nrow = var_numb_FA, ncol = irf_hor)
  
  
  for(i in 1:var_numb_FA){
    CIup[i,] <- apply(irf_Dimp0mat1[,i,], 2, quantile, probs = up_percentile/100)
    CIlow[i,] <- apply(irf_Dimp0mat1[,i,], 2, quantile, probs = low_percentile/100)
    IRFse[i,] <- apply(irf_Dimp0mat1[,i,], 2, sd)
    IRF_CImean[i,] <- apply(irf_Dimp0mat1[,i,], 2, median)
  }
  
  results <- list(irf_Dimo, CIup, CIlow, IRFse, IRF_CImean)
  names(results) <- c("irf_imp", "CIup", "CIlow", "IRFse", "IRF_CImean")
  return(results)
  
}

# Plot impulse responses
irf_plots3_sub <- function(impulse, name, impulse2){
  
  # function to plot the impulse responses.
  # Input: impulse response of all variables, name refers to the id of the 
  #        respective data series. Names can be chosen from all_names above. A 
  #        different title can be added instead of using the id. Second impulse
  #        respond can be added (e.g. Other regime) for comparison
  
  ind <- which(all_names == name)
  
  imp1 <- impulse$irf_imp
  imp2 <- impulse$CIup
  imp3 <- impulse$CIlow
  
  imp12 <- impulse2$irf_imp
  imp22 <- impulse2$CIup
  imp32 <- impulse2$CIlow
  
  over_min <- min(imp1[ind,], imp2[ind,], imp3[ind,], 
                  imp12[ind,], imp22[ind,], imp32[ind,])
  over_max <- max(imp1[ind,], imp2[ind,], imp3[ind,],
                  imp12[ind,], imp22[ind,], imp32[ind,])
  
  x_axis <- c(1:irf_hor, irf_hor:1)
  plot(imp1[ind,], type = "l", ylim = c(over_min, over_max),
       #xlim = c(1.65, irf_hor-1.75), xlab = NA, ylab = NA)#, col = col1)
       xlim = c(1.75, irf_hor-1.75), xlab = NA, ylab = NA)#, col = col1)
  title(sub = descriptions[ind], mgp = c(1,1,0))#, axes = F)
  
  polygon(x_axis, c(imp22[ind,], rev(imp32[ind,])), 
          col = "grey", border = NA)
  polygon(x_axis, c(imp2[ind,], rev(imp3[ind,])), 
          col = "grey", border = NA)
  
  lines(imp12[ind,], type = "l", col = "blue")#, lwd = 2)
  lines(imp1[ind,], type = "l", col = "red")#, lwd = 2)
  
  
  abline(h = 0, lty = "dotted")
  
}
