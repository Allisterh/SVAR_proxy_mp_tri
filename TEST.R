# Load functions ----------------------------------------------------------
rm(list = ls())
Functions <- devtools::as.package("../Functions")
devtools::load_all("../Functions")
library(ggfortify)
library(gridExtra)
library(magrittr)
library(svars)
library(readr)

#devtools::install_github("ropenscilabs/ochRe")
library(ochRe)

# Import data -------------------------------------------------------------
USA_Tri_raw   <- read.csv("../Applications/Data/USA_Tri.csv")
USA_Tri       <- ts(USA_Tri_raw[,-1], frequency = 4, start = c(floor(USA_Tri_raw$Year[1]), (USA_Tri_raw$Year[1] %% 1) / 0.25 + 1 ))
#autoplot(USA_Tri, ncol = 1) + theme_bw()

# VAR Spesicifation -------------------------------------------------------
## lag length
VARselect(USA_Tri, lag.max = 8)
var4 <- VAR(USA_Tri, p = 4, type = "const")
#var4 %>% summary
#var4 %>% resid %>% t() %>% write.table(file = "NonFundamentalness/resid.txt", col.names = F, row.names = F)
#apply(USA_Tri, 1, function(x){x-colMeans(USA_Tri)}) %>% t() %>% write.table(file = "NonFundamentalness/dat.txt", col.names = F, row.names = F)
## LM test for serial correlation
#serial.test(var4, type = "BG")
#serial.test(var4, type = "ES")

## detecting structural break
#set.seed(1234)
#chow.test(var4, SB = c(1984, 1), frequency = 4, nboot = 2000) %>% summary
#set.seed(1234)
#chow.test(var4, SB = c(1979, 3), frequency = 4, nboot = 2000) %>% summary

## Gaussianity test
#ICtest::FOBIboot(var4 %>% resid, n.boot = 2000, 1) # test against the last 2 components are gaussian: reject
#ICtest::FOBIboot(var4 %>% resid, n.boot = 2000, 2) # test against the last 1 component is gaussian: reject

#tseries::jarque.bera.test(resid(var4)[,1])
#tseries::jarque.bera.test(resid(var4)[,2])
#tseries::jarque.bera.test(resid(var4)[,3])


# Identification ----------------------------------------------------------
## Sign restrictions
set.seed(1234)
Signmat     <- matrix(c(1,1,1,-1,1,1,NA,-1,1), nrow = 3, ncol = 3)
SR          <- get.id.Sign(var4, Iter = 1000, Itermax = Inf, Signmat = Signmat, RHorizont = 0, Normalize = F, Step = 15, MT = T, Ci = 0.68, Plot = T, Hist = F, Epsname =  c("d", "s", "mp"))

## Identification by change in volatility
set.seed(1234)
CV          <- id.cv(var4, SB = c(1984, 1), frequency = 4)
CV          %>% summary


## Proxy SVAR
# Smets-Wouters (2007)
SW.raw      <- read.csv("../Applications/Data/Instruments/SW.csv")
SW          <- ts(SW.raw[,-1], frequency = 4, start = c(floor(SW.raw$Year[1]), (SW.raw$Year[1] %% 1) / 0.25 + 1 ))
IV.SW       <- get.id.iv(var4, instruments = SW, Synchro = T)

Tn = 1966
T1 = 1979.5
T2 = 1982.5

x <- CV
series = 2
Partial = 3
Epsname = c("d", "s", "mp")
Start = Tn
End = T2
Freq = 4



temp_plot <- y_all %>% select(time, y_hat, y_true) %>% reshape2::melt(id = "time")
temp_plot %>% ggplot(aes(x = time, y =  value, group = variable)) +
  geom_line(aes(linetype = variable))

hd.temp <- hd(CV, series = 2)
temp_plot2 <- hd.temp$hidec %>% as.data.frame() %>% select(1,2) %>% tibble::add_column("time" = y_all$time)%>% reshape2::melt(id = "time")
temp_plot2 %>% ggplot(aes(x = time, y =  value, group = variable)) +
  geom_line(aes(linetype = variable))


compare <- data.frame(y_all %>% select(time, y_hat, y_true), "true" =  hd.temp$hidec[,1], "construct" =  hd.temp$hidec[,2])

impulse <- irf.svars(CV, n.ahead = CV$VAR$obs)
impulse %>% plot
impulse$irf[1,-1] %>% matrix(nrow = 3, ncol = 3, byrow = T)
CV$B


# historic decomp ---------------------------------------------------------
# Function to calculate matrix potence
"%^%" <- function(A, n){
  if(n == 1){
    A
  }else{
    A %*% (A %^% (n-1))
  }
}


# function to calculate impulse response
IrF <- function(A_hat, B_hat, horizon){
  k <- nrow(A_hat)
  p <- ncol(A_hat)/k
  if(p == 1){
    irfa <- array(0, c(k, k, horizon))
    irfa[,,1] <- B_hat
    for(i in 1:horizon){
      irfa[,,i] <- (A_hat%^%i)%*%B_hat
    }
    return(irfa)
  }else{
    irfa <- array(0, c(k, k, horizon))
    irfa[,,1] <- B_hat
    Mm <- matrix(0, nrow = k*p, ncol = k*p)
    Mm[1:k, 1:(k*p)] <- A_hat
    Mm[(k+1):(k*p), 1 : ((p-1)*k)] <- diag(k*(p-1))
    Mm1 <- diag(k*p)
    for(i in 1:(horizon-1)){
      Mm1 <- Mm1%*%Mm
      irfa[,,(i+1)] <- Mm1[1:k, 1:k]%*%B_hat
    }
    return(irfa)
  }
}

## proxy SVAR
x <- IV.SW

Varname <- x$Varname
k <- length(x$Varname)
A_hat <- t(x$Beta)[,  1 : (x$p * k)]

B_hat <- matrix(x$B, nrow = k, ncol = k)
horizon <- length(x$epsilon)

# MA coeff
IR <- IrF(A_hat, B_hat, horizon)
if(is.null(Epsname)){Epsname <- Varname}
impulse <- IR[series,,]

## struct shocks
p <- x$p
s.time <- x$eps.ts$time
y <-  ts2df(x$dat, Varname) %>% 
  dplyr::filter(time >= s.time[1] & time <= s.time[length(s.time)]) %>% 
  dplyr::select(-1)

s.errors <- x$epsilon
y_hat <- rep(NA, horizon)
for (i in 1:horizon) {
  y_hat[i] <- t(impulse[Partial,1:i]) %*% s.errors[i:1]
}
yhat <- data.frame(s.time, y_hat)
colnames(yhat) <- c("t",
                    paste("Cumulative effect of ", Epsname[Partial], "shock on ", Varname[series]))




## sign restriction
x <- SR
k <- length(x$Varname)
A_hat <- t(x$Beta)[,  1 : (x$p * k)]

B_hat <- x$IRF.MT[,,1] #median target
horizon <- x$Tob

IR <- IrF(A_hat, B_hat, horizon)
if(is.null(Epsname)){Epsname <- x$Varname}
impulse <- IR[series,,]

Epsname = c("d", "s", "mp")



p <- x$p
s.time <- time(x$dat)[-(1:p)]
y <- x$dat[-c(1:p), series]
if (is.null(Partial)){
  s.errors <- x$epsilon.MT
  y_hat <- matrix(NA, nrow = horizon, ncol = k)
  for (i in 1:horizon) {
    for (j in 1:k) {
      y_hat[i,j] <- t(impulse[j,1:i]) %*% s.errors[j,i:1]
    }
  }

  y_hat_a <- rowSums(y_hat)
  
  yhat <- data.frame(s.time, (y - mean(y)), y_hat_a, y_hat)
  colnames(yhat)[1:3]<- c("t", 
                      paste("Demeaned series ", x$Varname[series]),
                      paste("Constructed series ", x$Varname[series]))
  for(i in 4:ncol(yhat)){
    colnames(yhat)[i] <- paste("Cumulative effect of ", Epsname[i-3], "shock on ", x$Varname[series])
  }
  

  
} else {
  s.errors <- x$epsilon.MT[Partial,]
  y_hat <- rep(NA, horizon)
  for (i in 1:horizon) {
    y_hat[i] <- t(impulse[Partial,1:i]) %*% s.errors[i:1]
  }
  yhat <- data.frame(s.time, y_hat)
  colnames(yhat) <- c("t", 
                      paste("Cumulative effect of ", Epsname[Partial], "shock on ", Varname[series]))
  
}



# CV ----------------------------------------------------------------------
x <- CV

Varname <- colnames(x$y)
k <- x$K
p <- x$p
horizon <- x$n
B_hat <- x$B

if(x$type == "const"){
  A_hat <- x$A_hat[,-1]
}else if(x$type == "trend"){
  A_hat <- x$A_hat[,-1]
}else if(x$type == "both"){
  A_hat <- x$A_hat[,-c(1,2)]
}else{
  A_hat <- x$A_hat
}

IR <- IrF(A_hat, B_hat, horizon)

if(is.null(Epsname)){Epsname <- Varname}
impulse <- IR[series,,]

## struct shocks
s.time <- time(x$y)[-(1:p)]
p <- x$p
y <- x$y[-c(1:p), series]

u <- t(resid(x$VAR))

s.errors <- solve(B_hat)%*%u

if (is.null(Partial)){
  y_hat <- matrix(NA, nrow = horizon, ncol = k)
  for (i in 1:horizon) {
    for (j in 1:k) {
      y_hat[i,j] <- t(impulse[j,1:i]) %*% s.errors[j,i:1]
    }
  }
  
  

  y_hat_a <- rowSums(y_hat)
  
  yhat <- data.frame(s.time, (y - mean(y)), y_hat_a, y_hat)
  colnames(yhat)[1:3]<- c("t",
                          paste("Demeaned series ", Varname[series]),
                          paste("Constructed series ", Varname[series]))
  for(i in 4:ncol(yhat)){
    colnames(yhat)[i] <- paste("Cumulative effect of ", Epsname[i-3], "shock on ", Varname[series])
  }
  
}else {
  s.errors <- s.errors[Partial,]
  y_hat <- rep(NA, horizon)
  for (i in 1:horizon) {
    y_hat[i] <- t(impulse[Partial,1:i]) %*% s.errors[i:1]
  }
  yhat <- data.frame(s.time, y_hat)
  colnames(yhat) <- c("t",
                      paste("Cumulative effect of ", Epsname[Partial], "shock on ", Varname[series]))
}

plot.df <- reshape2::melt(yhat, id = "t")
ggplot(plot.df, aes(x = t, y = value)) + geom_line() +
  facet_wrap(~variable, ncol = 1) +
  xlab("Time") + theme_bw()

