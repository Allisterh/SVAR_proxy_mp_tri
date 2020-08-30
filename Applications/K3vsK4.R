# Load functions ----------------------------------------------------------
rm(list = ls())
Functions <- devtools::as.package("../Functions")
devtools::load_all("../Functions")
library(ggfortify)
library(gridExtra)
library(magrittr)
library(svars)

# K = 3 -------------------------------------------------------------------
USA_Tri_raw   <- read.csv("data/USA_Tri.csv")
USA_Tri       <- ts(USA_Tri_raw[,-1], frequency = 4, start = c(floor(USA_Tri_raw$Year[1]), (USA_Tri_raw$Year[1] %% 1) / 0.25 + 1 ))
var4 <- VAR(USA_Tri, p = 4, type = "const")

## Sign restrictions
set.seed(1234)
Signmat3     <- matrix(c(1,1,1,-1,1,1,NA,-1,1), nrow = 3, ncol = 3)
SR3          <- get.id.Sign(var4, Iter = 1000, Itermax = Inf, Signmat = Signmat3, RHorizont = 0, Normalize = F, Step = 15, MT = T, Ci = 0.68, Plot = T, Hist = F, Epsname =  c("d", "s", "mp"))

## Identification by change in volatility
set.seed(1234)
CV3          <- id.cv(var4, SB = c(1984, 1), frequency = 4)

## Identification by minimizing distance covariance
set.seed(1234)
dCov3        <- id.dc(var4)
dCov3$B      <- dCov3$B[,c(3,2,1)]
dCov3$B[,3]  <- dCov3$B[,3]*-1


## Proxy SVAR
# Smets-Wouters (2007)
SW.raw      <- read.csv("Data/Instruments/SW.csv")
SW          <- ts(SW.raw[,-1], frequency = 4, start = c(floor(SW.raw$Year[1]), (SW.raw$Year[1] %% 1) / 0.25 + 1 ))
IV.SW3       <- get.id.iv(var4, instruments = SW, Synchro = T)

# Romer Romer (2004)
RR.raw      <- read.csv("Data/Instruments/RR.csv")
RR          <- ts(RR.raw[,-1], frequency = 4, start = c(floor(RR.raw$Year[1]), (RR.raw$Year[1] %% 1) / 0.25 + 1 ))
RR          <- Remove_AC(RR, lag = 1) 
IV.RR3       <- get.id.iv(var4, instruments = RR, Synchro = T)

# Sims Zha (2006)
SZ.raw      <- read.csv("Data/Instruments/SZ.csv")
SZ          <- ts(SZ.raw[,-1], frequency = 4, start = c(floor(SZ.raw$Year[1]), (SZ.raw$Year[1] %% 1) / 0.25 + 1 ))
SZ          <- Remove_AC(SZ, 2)
IV.SZ3       <- get.id.iv(var4, instruments = SZ, Synchro = T)


# quasi HD ----------------------------------------------------------------
Start = 1979.5
End1 = 1982.5
End2 = 1983

df3 <- data.frame("period" = c("1979Q3 - 1982Q3", "1979Q3 - 1983Q1"),
                  "K" = "K = 3",
                  "SR" = c(quasi.hd(SR3, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End1, Freq = 4),
                           quasi.hd(SR3, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End2, Freq = 4)),
                  "dCov" = c(quasi.hd(dCov3, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End1, Freq = 4),
                             quasi.hd(dCov3, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End2, Freq = 4)),
                  "CV" = c(quasi.hd(CV3, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End1, Freq = 4),
                        quasi.hd(CV3, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End2, Freq = 4)),
                  "SW" = c(quasi.hd(IV.SW3, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End1, Freq = 4),
                           quasi.hd(IV.SW3, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End2, Freq = 4)),
                  "SZ" = c(quasi.hd(IV.SZ3, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End1, Freq = 4),
                           quasi.hd(IV.SZ3, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End2, Freq = 4)),
                  "RR" = c(quasi.hd(IV.RR3, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End1, Freq = 4),
                           quasi.hd(IV.RR3, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End2, Freq = 4)))




# K = 4 -------------------------------------------------------------------

# Import data -------------------------------------------------------------
USA_Qua.raw <- read.csv("Data/USA_Qua.csv")
USA_Qua <- ts(USA_Qua.raw[,-1], frequency = 4, start = c(floor(USA_Qua.raw$Year[1]), (USA_Qua.raw$Year[1] %% 1) / 0.25 + 1 ))
var3 <- VAR(USA_Qua, p = 3, type = "const")

## Sign restrictions
set.seed(1234)
Signmat4     <- matrix(c(1,1,1,NA,-1,1,1,NA,NA,-1,1,rep(NA,5)), nrow = 4, ncol = 4)
SR4          <- get.id.Sign(var3, Iter = 1000, Itermax = Inf, Signmat = Signmat4, RHorizont = 0, Normalize = F, Step = 15, MT = T, Ci = 0.68, Plot = T, Hist = F, Epsname =  c("d", "s", "mp", "f"))

## Identification by change in volatility
set.seed(1234)
CV4          <- id.cv(var3, SB = c(1984, 1), frequency = 4)

## Identification by minimizing distance covariance
set.seed(1234)
dCov4        <- id.dc(var3)

## Proxy SVAR
# Smets-Wouters (2007)
SW.raw      <- read.csv("Data/Instruments/SW.csv")
SW          <- ts(SW.raw[,-1], frequency = 4, start = c(floor(SW.raw$Year[1]), (SW.raw$Year[1] %% 1) / 0.25 + 1 ))
IV.SW4       <- get.id.iv(var3, instruments = SW, Synchro = T)

# Romer Romer (2004)
RR.raw      <- read.csv("Data/Instruments/RR.csv")
RR          <- ts(RR.raw[,-1], frequency = 4, start = c(floor(RR.raw$Year[1]), (RR.raw$Year[1] %% 1) / 0.25 + 1 ))
RR          <- Remove_AC(RR, lag = 1) 
IV.RR4       <- get.id.iv(var3, instruments = RR, Synchro = T)

# Sims Zha (2006)
SZ.raw      <- read.csv("Data/Instruments/SZ.csv")
SZ          <- ts(SZ.raw[,-1], frequency = 4, start = c(floor(SZ.raw$Year[1]), (SZ.raw$Year[1] %% 1) / 0.25 + 1 ))
SZ          <- Remove_AC(SZ, 2)
IV.SZ4       <- get.id.iv(var3, instruments = SZ, Synchro = T)


# quasi HD ----------------------------------------------------------------
Start = 1979.5
End1 = 1982.5
End2 = 1983

df4 <- data.frame("period" = c("1979Q3 - 1982Q3", "1979Q3 - 1983Q1"),
                  "K" = "K = 4",
                  "SR" = c(quasi.hd(SR4, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End1, Freq = 4),
                           quasi.hd(SR4, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End2, Freq = 4)),
                  "dCov" = c(quasi.hd(dCov4, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End1, Freq = 4),
                             quasi.hd(dCov4, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End2, Freq = 4)),
                  "CV" = c(quasi.hd(CV4, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End1, Freq = 4),
                           quasi.hd(CV4, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End2, Freq = 4)),
                  "SW" = c(quasi.hd(IV.SW4, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End1, Freq = 4),
                           quasi.hd(IV.SW4, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End2, Freq = 4)),
                  "SZ" = c(quasi.hd(IV.SZ4, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End1, Freq = 4),
                           quasi.hd(IV.SZ4, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End2, Freq = 4)),
                  "RR" = c(quasi.hd(IV.RR4, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End1, Freq = 4),
                           quasi.hd(IV.RR4, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End2, Freq = 4)))


# plot --------------------------------------------------------------------
df.all <- rbind(df3, df4)

df.plot <- df.all %>% melt(id = c("period", "K"))

df.plot <- df.plot %>% mutate(value = value / 4)

df.plot %>% ggplot(aes(x = variable, y = value, group = K)) + 
  facet_grid(rows = vars(K), cols = vars(period)) +
  geom_hline(yintercept = 0, color = "red", alpha = 0.8) +
  geom_bar(stat="identity", alpha = .7, width=0.5, aes(fill = variable)) + 
  xlab("") + ylab(" ") + theme_bw() + theme(legend.position = "none")

ggsave("../Plots/inflation_fighter.pdf", plot = last_plot(), width = 12, height = 6)


delta_pi <- c(window(USA_Tri, End1, End1)[,"pi"] %>% as.numeric() - window(USA_Tri, Start, Start)[,"pi"] %>% as.numeric() ,
              window(USA_Tri, End2, End2)[,"pi"] %>% as.numeric() - window(USA_Tri, Start, Start)[,"pi"] %>% as.numeric() ,
              window(USA_Qua, End1, End1)[,"pi"] %>% as.numeric() - window(USA_Qua, Start, Start)[,"pi"] %>% as.numeric() ,
              window(USA_Qua, End2, End2)[,"pi"] %>% as.numeric() - window(USA_Qua, Start, Start)[,"pi"] %>% as.numeric() )
cbind(df.all, delta_pi) %>% kableExtra::kable(format = "latex", booktabs = T, digits = 3, linesep = "")




# Ivertibility ------------------------------------------------------------

GC_inv <- function(Datas, IV, K){
  mod1 <- suppressWarnings( window(Datas, start = time(IV)[1], end = time(SW)[length(IV)], frequency = 4))
  mod1 <- cbind(mod1, IV) %>% as.data.frame() %>% tidyr::drop_na()
  mod1 %>% VARselect(lag.max = 8)
  temp <- mod1 %>% VAR(lag.max = 8, ic = "AIC", type = "const")
  temp1 <- temp %>% summary()
  for (i in 1:K) {
    print(temp1$varresult[[i]]$coeff[(K+1),4])
  }
}

GC_inv(USA_Tri, SW, 3)
GC_inv(USA_Tri, SZ, 3)
GC_inv(USA_Tri, RR, 3)

GC_inv(USA_Qua, SW, 4)
GC_inv(USA_Qua, SZ, 4)
GC_inv(USA_Qua, RR, 4)


