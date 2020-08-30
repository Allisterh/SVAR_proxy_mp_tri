# Load functions ----------------------------------------------------------
rm(list = ls())
Functions <- devtools::as.package("../Functions")
devtools::load_all("../Functions")


# Solving DSGE to VAR -----------------------------------------------------
## DSGE Calibration
DSGE.Par <- DSGE2VAR.DGP1(
  alpha   = 0.5,
  beta    = 0.99,
  k       = 0.05,
  gamma   = 0.5,
  delta_x = 0.1,
  tau_x   = 0.5,
  tau_pi  = 1.8,
  tau_r   = 0.6,
  rho_x   = 0.5,
  rho_pi  = 0.5,
  rho_r   = 0.5)


# True IRFs and AoA -------------------------------------------------------
## IRF
IRF.TRUE <- get.my.irf(Acoef = array(c(DSGE.Par$A_1, DSGE.Par$A_2), c(3,3,2)), Step = 15,
                       Normalize = TRUE, Bmat = DSGE.Par$B, Varname = c("x", "pi", "r"))
## AoA
#set.seed(1234)
#AOA <- AOA.DGP1(Simu = 10000, Step = 15, Ci = 0.9)
AOA <- readRDS("Replica/DGP1.AOA.rds")
## plot
# plot.AOA(IRF.TRUE, AOA)


# Monte Carlo configuration -----------------------------------------------
## Note: Number of experiments (Size) and observations (Tob) are reduced for show case !!
Size <- 200  # Number of Monte Carlo experiments
Tob  <- 200  # Length of time series for each replicate

## for questions: run "?DGP1"
Simu.par.DGP <- list(
  "Par"       = DSGE.Par,
  "Size"      = Size,
  "Tob"       = Tob,
  "Break"     = c(Tob/3, Tob*2/3),
  "Volaregim" = c(c(1.2, 1.1, 1),c(3.6, 1.1, 0.5)),
  "VolaRand"  = NULL,
  "Edist"     = "normdist",
  "v"         = 5,
  "Combi"     = matrix(c(0.4,  0.3, 0.3,
                         0.25, 0.5, 0.25,
                         0.2,  0.2, 0.6), nrow = 3, ncol = 3, byrow = T),
  "ME"        = NULL,
  "MERand"    = NULL)


# Data generating ---------------------------------------------------------
Y <- simu.datas(Config = Simu.par.DGP, DGP = 1)


# Identification and evaluation -------------------------------------------
## Sign restrictions
Signmat   <- matrix(c(1,1,1,-1,1,1,NA,-1,1), nrow = 3, ncol = 3) # agnostic sign restrictions
TrueSign  <- matrix(c(1,1,1,-1,1,1,-1,-1,1), nrow = 3, ncol = 3) # true signs
SR        <- simulate.sign(Data = Y, AOA = AOA, Iter = 50, Itermax = Inf, Signmat = Signmat, TrueSign = TrueSign, Bmat = DSGE.Par$B, RHorizont = 0, MT = T, Ci = 0.68)

## Identification by change in volatility
CV        <- simulate.cv(Data = Y, AOA = AOA, Break = Simu.par.DGP$Break, Bmat = DSGE.Par$B, Guess = TRUE)

## Identification by minimizing distance covariance
dCov      <- simulate.dcov(Data = Y, AOA = AOA, Bmat = DSGE.Par$B)

## Identification by minimizing Cramer-von Mises distance (VERY SLOW! Trust me, it works.)
#CvM       <- simulate.cvm(Data = Y, AOA = AOA, Iter = 1000, Steptol = 200, Stage2 = 500)

## Identification by valid instruments
IV        <- simulate.iv(Data = Y, AOA = AOA, Bmat = DSGE.Par$B)

## Identification by weak instruments
IV.W      <- simulate.iv(Data = Y, AOA = AOA, Weak = T, Bmat = DSGE.Par$B)

## Identification by invalid instruments
IV.P      <- simulate.iv(Data = Y, AOA = AOA, Combi = T, Bmat = DSGE.Par$B)


# Tabulate ----------------------------------------------------------------

tabulate.simu.aoa(list(SR$M, CV$aoa, dCov$aoa, IV$aoa, IV.W$aoa, IV.P$aoa),
                  Label = c("SR", "CV", "dCov", "IV", "IV-W", "IV-P"))

tabulate.simu.ump(list(SR, CV, dCov, IV, IV.W, IV.P),
                  Label = c("SR", "CV", "dCov", "IV", "IV-W", "IV-P"))

# Plots -------------------------------------------------------------------
plot.simu.multi(list(SR$M,  CV$aoa, dCov$aoa,
                     # CvM,
                     IV$aoa, IV.W$aoa, IV.P$aoa),
                Label = c("SR", "CV", "dCov",
                          # "CvM",
                          "IV", "IV-W", "IV-P"), shape_manual = c(16:17, 5, 1:4))


# Other distribution scenarios? -------------------------------------------
Simu.par.DGP$Edist <- "tdist"   # Run code line 56-91
Simu.par.DGP$Edist <- "chidist" # Run code line 56-91
Simu.par.DGP$Edist <- "mixed_a" # Run code line 56-91
Simu.par.DGP$Edist <- "mixed_b" # Run code line 56-91


# Subsample performance ---------------------------------------------------
## sign restrictions
SR_1 <- simulate.sign(Data = Y, AOA = AOA, Iter = 200, Itermax = Inf,  Subsample = c(1,79),
                             Signmat, TrueSign, RHorizont = 0, MT = T, Ci = 0.68)
SR_2 <- simulate.sign(Data = Y, AOA = AOA, Iter = 200, Itermax = Inf,  Subsample = c(80,159),
                             Signmat, TrueSign, RHorizont = 0, MT = T, Ci = 0.68)
SR_3 <- simulate.sign(Data = Y, AOA = AOA, Iter = 200, Itermax = Inf,  Subsample = c(160,240),
                             Signmat, TrueSign, RHorizont = 0, MT = T, Ci = 0.68)
## plots
plot.simu.multi(list(SR_1$MT, SR_2$MT, SR_3$MT), Label = c("I", "II", "III"), shape_manual = c(0,7,15)) +
  labs(color = "Subsample",shape = "Subsample", linetype = "Subsample")

## Identification by invalid instruments
IV.P_1 <- simulate.iv(Data = Y, AOA = AOA, Combi = T, Subsample = c(1,79))
IV.P_2 <- simulate.iv(Data = Y, AOA = AOA, Combi = T, Subsample = c(80,159))
IV.P_3 <- simulate.iv(Data = Y, AOA = AOA, Combi = T, Subsample = c(160,240))
## plots
plot.simu.multi(list(IV.P_1, IV.P_2, IV.P_3), Label = c("I", "II", "III"), shape_manual = c(0,7,15)) +
  labs(color = "Subsample",shape = "Subsample", linetype = "Subsample")
