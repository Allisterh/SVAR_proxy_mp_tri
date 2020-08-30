# Load functions ----------------------------------------------------------
rm(list = ls())
Functions <- devtools::as.package("../Functions")
devtools::load_all("../Functions")
library(ggfortify)
library(gridExtra)
library(magrittr)
library(svars)
library(tidyr)


# Import data -------------------------------------------------------------
USA_Qua.raw <- read.csv("Data/USA_Qua.csv")
USA_Qua <- ts(USA_Qua.raw[,-1], frequency = 4, start = c(floor(USA_Qua.raw$Year[1]), (USA_Qua.raw$Year[1] %% 1) / 0.25 + 1 ))
autoplot(USA_Qua, ncol = 1) + theme_bw()

#urca::ur.df(USA_Qua[,1], type = "drift", selectlags = "AIC")
#urca::ur.df(USA_Qua[,2], type = "drift", selectlags = "AIC")
#urca::ur.df(USA_Qua[,3], type = "drift", selectlags = "AIC")
#urca::ur.df(USA_Qua[,4], type = "drift", selectlags = "AIC")


# VAR Spesicifation -------------------------------------------------------
## lag length
VARselect(USA_Qua, lag.max = 8)
var3 <- VAR(USA_Qua, p = 3, type = "const")
#var3 %>% summary

## LM test for serial correlation
#serial.test(var3, type = "BG")
#serial.test(var3, type = "ES")

## detecting structural break
#chow.test(var3, SB = c(1984, 1), frequency = 4, nboot = 2000) %>% summary 
#chow.test(var3, SB = c(1979, 3), frequency = 4, nboot = 2000) %>% summary 

## Gaussianity test
#ICtest::FOBIboot(var3 %>% resid, n.boot = 2000, 1) # test against the last 3 components are gaussian: reject
#ICtest::FOBIboot(var3 %>% resid, n.boot = 2000, 2) # test against the last 2 components are gaussian: reject
#ICtest::FOBIboot(var3 %>% resid, n.boot = 2000, 3) # test against the last 1 component is gaussian components: reject

#tseries::jarque.bera.test(resid(var3)[,1])
#tseries::jarque.bera.test(resid(var3)[,2])
#tseries::jarque.bera.test(resid(var3)[,3])
#tseries::jarque.bera.test(resid(var3)[,4])


# Identification ----------------------------------------------------------
## Sign restrictions
set.seed(1234)
Signmat     <- matrix(c(1,1,1,NA,-1,1,1,NA,NA,-1,1,rep(NA,5)), nrow = 4, ncol = 4)
SR          <- get.id.Sign(var3, Iter = 1000, Itermax = Inf, Signmat = Signmat, RHorizont = 0, Normalize = F, Step = 15,
                           MT = T, Ci = 0.68, Plot = T, Hist = F, Epsname =  c("d", "s", "mp", "f"))

## Identification by change in volatility
set.seed(1234)
CV          <- id.cv(var3, SB = c(1984, 1), frequency = 4)
CV          %>% summary

## Identification by minimizing distance covariance
set.seed(1234)
dCov        <- id.dc(var3)
dCov        %>% summary

## Identification by minimizing CvM distance
#set.seed(1234)
#CvM        <- id.cvm(var3, itermax = 2000, steptol = 200, iter2 = 200)
#CvM$B[,c(3)] <- CvM$B[,c(3)]  *-1
#CvM        %>% summary
#CvM %>% irf() %>% plot
#ggsave("../Plots/Qua_CvM.pdf", plot = last_plot(), width = 12, height = 8)



## Proxy SVAR
# Smets-Wouters (2007)
SW.raw      <- read.csv("Data/Instruments/SW.csv")
SW          <- ts(SW.raw[,-1], frequency = 4, start = c(floor(SW.raw$Year[1]), (SW.raw$Year[1] %% 1) / 0.25 + 1 ))
IV.SW       <- get.id.iv(var3, instruments = SW, Synchro = T)
IV.SW$B
IV.SW$F_test

# Guerkaynak, Sack and Swanson (2005)
GSS.raw     <- read.csv("Data/Instruments/GSS.csv")
GSS         <- ts(GSS.raw[,-1], frequency = 4, start = c(floor(GSS.raw$Year[1]), (GSS.raw$Year[1] %% 1) / 0.25 + 1 ))
IV.GSS      <- get.id.iv(var3, instruments = GSS, Synchro = T)
IV.GSS$B
IV.GSS$F_test

# Romer Romer (2004)
RR.raw      <- read.csv("Data/Instruments/RR.csv")
RR          <- ts(RR.raw[,-1], frequency = 4, start = c(floor(RR.raw$Year[1]), (RR.raw$Year[1] %% 1) / 0.25 + 1 ))
RR          <- Remove_AC(RR, lag = 1) 
IV.RR       <- get.id.iv(var3, instruments = RR, Synchro = T)
IV.RR$B
IV.RR$F_test

# Sims Zha (2006)
SZ.raw      <- read.csv("Data/Instruments/SZ.csv")
SZ          <- ts(SZ.raw[,-1], frequency = 4, start = c(floor(SZ.raw$Year[1]), (SZ.raw$Year[1] %% 1) / 0.25 + 1 ))
SZ          <- Remove_AC(SZ, 2)
IV.SZ       <- get.id.iv(var3, instruments = SZ, Synchro = T)
IV.SZ$B
IV.SZ$F_test


# PCA ---------------------------------------------------------------------
iv2ts <- function(x, varname){
  erg <- data.frame(x %>% time %>% as.numeric(),  x)
  colnames(erg) <- c("time", varname)
  return(erg)}
SW.ts <- iv2ts(SW, "SW")
GSS.ts <- iv2ts(GSS, "GSS")
RR.ts <- iv2ts(RR, "RR")
SZ.ts <- iv2ts(SZ, "SZ")


# Inference ---------------------------------------------------------------
#set.seed(1234)
#CV.boot     <- svars::wild.boot(CV, n.ahead = 20, nboot = 2000)
#saveRDS(CV.boot, file = "Output/CV_boot4.rds")
#set.seed(1234)
#dCov.boot   <- svars::wild.boot(dCov, n.ahead = 20, nboot = 2000)
#saveRDS(dCov.boot, file = "Output/dCov_boot4.rds")

CV.boot <- readRDS("Output/CV_boot4.rds")
dCov.boot <- readRDS("Output/dCov_boot4.rds")

CV.boot %>% signcheck.sboot()
dCov.boot %>% signcheck.sboot()



set.seed(1234)
IV.SW.boot  <- get.MBB.fixed(IV.SW, design = "fixed", Step = 15, nboot = 2000, Ci = 0.68)
#IV.SW.boot  <- get.wild.fixed.Boots(IV.SW, Step = 15, nboot = 2000, Ci = 0.9)
#set.seed(1234)
#IV.GSS.boot <- get.MBB.fixed(IV.GSS, design = "fixed", Step = 15, nboot = 2000, Ci = 0.68)
#IV.GSS.boot <- get.wild.fixed.Boots(IV.GSS, Step = 15, nboot = 2000, Ci = 0.9)
set.seed(1234)
IV.RR.boot <- get.MBB.fixed(IV.RR, design = "fixed", Step = 15, nboot = 2000, Ci = 0.68)
#IV.RR.boot  <- get.wild.fixed.Boots(IV.RR, Step = 15, nboot = 2000, Ci = 0.9)
set.seed(1234)
IV.SZ.boot  <- get.MBB.fixed(IV.SZ, design = "fixed", Step = 15, nboot = 1000, Ci = 0.68)
#IV.SZ.boot  <- get.wild.fixed.Boots(IV.SZ, Step = 15, nboot = 1000, Ci = 0.9)

# Tab 
IV_id.tab <- function(ivid, ivboots){
  temp <- cbind(ivboots$Bmat$sd[,3], ivid$B)
  return(rbind(temp, rep(ivid$F_test,2))) 
}
tab_iv_id <- cbind(IV_id.tab(IV.SW, IV.SW.boot),
                   IV_id.tab(IV.RR, IV.RR.boot),
                   IV_id.tab(IV.SZ, IV.SZ.boot))
                   
rownames(tab_iv_id) <- NULL
kableExtra::kable(tab_iv_id, format = "latex", booktabs = T, digits = 3, linesep = "")



# Plot --------------------------------------------------------------------
# standardize IRFs plots
st_fct <- function(sboot, h){
  sboot$true$irf <- sboot$true$irf[1:(h+1),]
  sboot$true$irf[,1] <- 0:h
  for (i in 1:length(sboot[["bootstrap"]])) {
    sboot[["bootstrap"]][[i]]$irf <- sboot[["bootstrap"]][[i]]$irf[1:(h+1),]
    sboot[["bootstrap"]][[i]]$irf[,1] <- 0:h
  }
  return(sboot)}
CV.boot    %<>% st_fct(h = 15)
dCov.boot  %<>% st_fct(h = 15)


#change shock names
colnames(CV.boot$true$irf)[4] <- "epsilon[ mp ] %->% x"
colnames(CV.boot$true$irf)[8] <- "epsilon[ mp ] %->% pi"
colnames(CV.boot$true$irf)[12] <- "epsilon[ mp ] %->% GBR1"
colnames(CV.boot$true$irf)[16] <- "epsilon[ mp ] %->% sp"

colnames(dCov.boot$true$irf)[4] <- "epsilon[ mp ] %->% x"
colnames(dCov.boot$true$irf)[8] <- "epsilon[ mp ] %->% pi"
colnames(dCov.boot$true$irf)[12] <- "epsilon[ mp ] %->% GBR1"
colnames(dCov.boot$true$irf)[16] <- "epsilon[ mp ] %->% sp"


# Normalization 1MP->25bpsGBR1 --------------------------------------------
SR$IRF.plot %<>% mutate(value = ifelse(Label == "median", 0.25*value/SR$IRF.M[3,3,1], ifelse(Label == "median target (MT)", 0.25*value/SR$IRF.MT[3,3,1], NA)))
SR$IRF.plot %<>% mutate(L =  0.25*L/SR$IRF.M[3,3,1])
SR$IRF.plot %<>% mutate(U =  0.25*U/SR$IRF.M[3,3,1])

SR.MP     <- plot.my.irf(My.id = SR, Partial = 3, Hide.Legend = T, Color.manual = c("black", "blue"), Lty.manual = c("solid", "longdash"), Epsname =  c("d", "s", "mp", "f")) + xlab("") + ylab("")
CV.MP     <- CV.boot %>% plot.sboot.mod(lowerq = 0.16, upperq = 0.84, Partial = 3, Normalize = CV$B[3,3]/0.25) + xlab("") + ylab("")
dCov.MP   <- dCov.boot %>% plot.sboot.mod(lowerq = 0.16, upperq = 0.84, Partial = 3, Normalize = dCov$B[3,3]/0.25) + xlab("") + ylab("")
IRF_nonIV <- grid.arrange(SR.MP, CV.MP, dCov.MP, ncol = 3)
ggsave("../Plots/Qua_A.pdf", plot = IRF_nonIV, width = 12, height = 8)

SW.MP     <- IV.SW.boot %>% plot.my.irf(Partial = 3, Epsname =  c("d", "s", "mp", "f"), Normalize = IV.SW$B[3]/0.25) + xlab("") + ylab("")
RR.MP     <- IV.RR.boot %>% plot.my.irf(Partial = 3, Epsname =  c("d", "s", "mp", "f"), Normalize = IV.RR$B[3]/0.25) + xlab("") + ylab("")
SZ.MP     <- IV.SZ.boot %>% plot.my.irf(Partial = 3, Epsname =  c("d", "s", "mp", "f"), Normalize = IV.SZ$B[3]/0.25) + xlab("") + ylab("")
IRF_IV <- grid.arrange(SW.MP, RR.MP, SZ.MP, ncol = 3)
ggsave("../Plots/Qua_B_MBB_cum.pdf", plot = IRF_IV, width = 12, height = 8)


# SP2MP plots -------------------------------------------------------------
SP2MP.SR <- SR.MP$data %>% filter(variable == "epsilon[ mp ] %->% sp")
# delet MT
SP2MP.SR <- SP2MP.SR %>% filter(Label == "median") %>% select(-Label)

SP2MP.SR$variable <- "SR"

SP2MP.CV     <- CV.MP$data %>% filter(variable == "epsilon[ mp ] %->% sp") %>% select(-probs)
SP2MP.CV <- SP2MP.CV[,c(1,2,5,3,4)]
colnames(SP2MP.CV) <- colnames(SP2MP.SR)
SP2MP.CV$variable <- "CV"

SP2MP.dCov   <- dCov.MP$data %>% filter(variable == "epsilon[ mp ] %->% sp") %>% select(-probs)
SP2MP.dCov <- SP2MP.dCov[,c(1,2,5,3,4)]
colnames(SP2MP.dCov) <- colnames(SP2MP.SR)
SP2MP.dCov$variable <- "dCov"

SP2MP.nonIV <- rbind(SP2MP.SR, SP2MP.CV, SP2MP.dCov)
SP2MP.nonIV$variable <- factor(SP2MP.nonIV$variable, levels = c("SR", "CV", "dCov"))


SP2MP.SW     <- SW.MP$data %>% filter(variable == "epsilon[ mp ] %->% sp")
SP2MP.RR     <- RR.MP$data %>% filter(variable == "epsilon[ mp ] %->% sp")
SP2MP.SZ     <- SZ.MP$data %>% filter(variable == "epsilon[ mp ] %->% sp")

SP2MP.SW$variable <- "SW"
SP2MP.RR$variable <- "RR"
SP2MP.SZ$variable <- "SZ"

SP2MP.IV <- rbind(SP2MP.SW, SP2MP.RR, SP2MP.SZ)
SP2MP.IV$variable <- factor(SP2MP.IV$variable, levels = c("SW", "RR", "SZ"))



SP2MP_all <- rbind(SP2MP.nonIV, SP2MP.IV)


#saveRDS(SP2MP_all, "SP2MP_all.rds")


# cummulation -------------------------------------------------------------
cumfuc <- function(x){exp(cumsum(x) - 1)*100}
#cumfuc <- function(x){cumsum(x)}
cumfuc <- function(x){exp(cumsum(x)/100)-1}

cumfuc <- function(x){exp(cumsum(x/100))-1}

SP2MP_all %>% filter(variable == "SR") %>% select(value) %>% cumfuc



 
plot1.1 <- SP2MP_all %>% filter(variable == "SR") %>%
  mutate(value = cumfuc(value), L = cumfuc(L), U = cumfuc(U)) %>% 
  ggplot(aes(x = h, y = value)) + geom_line() + 
  geom_hline(yintercept = 0, color = 'red') + geom_ribbon(aes(ymin = L, ymax = U), alpha=0.3) +
  xlab("") + ylab("") + theme_bw() + facet_grid(. ~ variable)

plot1.2 <- SP2MP_all %>% filter(variable == "CV") %>% 
  mutate(value = cumfuc(value), L = cumfuc(L), U = cumfuc(U)) %>% 
  ggplot(aes(x = h, y = value)) + geom_line() + 
  geom_hline(yintercept = 0, color = 'red') + geom_ribbon(aes(ymin = L, ymax = U), alpha=0.3) +
  xlab("") + ylab("") + theme_bw() + facet_grid(. ~ variable) #+ ylim(yl, yu)

plot1.3 <- SP2MP_all %>% filter(variable == "dCov") %>% 
  mutate(value = cumfuc(value), L = cumfuc(L), U = cumfuc(U)) %>% 
  ggplot(aes(x = h, y = value)) + geom_line() + 
  geom_hline(yintercept = 0, color = 'red') + geom_ribbon(aes(ymin = L, ymax = U), alpha=0.3) +
  xlab("") + ylab("") + theme_bw() + facet_grid(. ~ variable) #+ ylim(yl, yu)

plot2.1 <- SP2MP_all %>% filter(variable == "SW") %>%   
  mutate(value = cumfuc(value), L = cumfuc(L), U = cumfuc(U)) %>% 
  ggplot(aes(x = h, y = value)) + geom_line() + 
  geom_hline(yintercept = 0, color = 'red') + geom_ribbon(aes(ymin = L, ymax = U), alpha=0.3) +
  xlab("") + ylab("") + theme_bw() + facet_grid(. ~ variable) #+ ylim(yl, yu)

plot2.2 <- SP2MP_all %>% filter(variable == "RR") %>% 
  mutate(value = cumfuc(value), L = cumfuc(L), U = cumfuc(U)) %>% 
  ggplot(aes(x = h, y = value)) + geom_line() + 
  geom_hline(yintercept = 0, color = 'red') + geom_ribbon(aes(ymin = L, ymax = U), alpha=0.3) +
  xlab("") + ylab("") + theme_bw() + facet_grid(. ~ variable) #+ ylim(yl, yu)

plot2.3 <- SP2MP_all %>% filter(variable == "SZ") %>% 
  mutate(value = cumfuc(value), L = cumfuc(L), U = cumfuc(U)) %>% 
  ggplot(aes(x = h, y = value)) + geom_line() + 
  geom_hline(yintercept = 0, color = 'red') + geom_ribbon(aes(ymin = L, ymax = U), alpha=0.3) +
  xlab("") + ylab("") + theme_bw() + facet_grid(. ~ variable) #+ ylim(yl, yu)

plot_SP2MP_all <- grid.arrange(plot1.1, plot1.2, plot1.3, plot2.1, plot2.2, plot2.3, ncol = 3)


ggsave("../Plots/SP2MP.pdf", plot = plot_SP2MP_all, width = 12, height = 6)



# Correlation of MP shocks ------------------------------------------------
mp_sr     <- solve(SR$IRF.M[,,1]) %*% t(resid(var3)) %>% t() %>% as.data.frame() %>% select(3)
mp_cv     <- solve(CV$B)          %*% t(resid(var3)) %>% t() %>% as.data.frame() %>% select(3)
mp_dcov   <- solve(dCov$B)        %*% t(resid(var3)) %>% t() %>% as.data.frame() %>% select(3)
mp_noniv  <- cbind(mp_sr, mp_cv, mp_dcov); colnames(mp_noniv) <- c("SR", "CV", "dCov")
mp_iv     <- left_join(IV.SW$eps.ts, IV.SZ$eps.ts, by = "time")
mp_iv     <- left_join(mp_iv, IV.RR$eps.ts, by = "time") %>% tidyr::drop_na()
colnames(mp_iv) <- c("time", "SW", "SZ", "RR")
mp_noniv  <- cbind(time(var3$y)[-c(1:3)], mp_noniv); colnames(mp_noniv)[1] <- "time"
mp_all    <- left_join(mp_noniv, mp_iv, by = "time") %>% tidyr::drop_na()

cor(mp_all %>% select(-1)) %>% kableExtra::kable(format = "latex", booktabs = T, digits = 2, linesep = "")
cor(mp_all %>% filter(time < 1987)%>% select(-1))
cor(mp_all %>% filter(time > 1987)%>% select(-1))

mp_sds <- mp_all %>% select(-1) %>% apply( 2, sd)

mp_all_stand <- mp_all

for (j in 1:dim(mp_all)[1]) {
  for (i in 2:7) {
    mp_all_stand[j,i] <- mp_all[j,i] / mp_sds[i-1]
  } 
}

mp_all_stand %>% select(-1) %>% apply( 2, sd)


mp_all_plot4 <- reshape2::melt( mp_all_stand, id.var = "time")

library(ggplot2)

MP_shock_plot4 <- ggplot(mp_all_plot4, aes(x = time, y = value, group = variable)) +
  geom_line(aes(linetype = variable, color = variable)) + 
  scale_x_continuous(breaks = seq(mp_all$time[1] %>% ceiling(),mp_all$time[nrow(mp_all)], 2)) +
  xlab("") + ylab(" ")  + theme_bw() + theme(legend.title = element_blank())
MP_shock_plot <- grid.arrange(MP_shock_plot3, MP_shock_plot4, ncol = 2)

# K = 3 vs K = 4
(mp_all_plot4) %>% dim
(mp_all_plot3) %>% dim

mp_all_plot3 <- mp_all_plot3 %>% add_column("K" = "K = 3")
mp_all_plot4 <- mp_all_plot4 %>% add_column("K" = "K = 4")

mp_all_plot_all <- rbind(mp_all_plot3, mp_all_plot4)

# with SR (delet %>% select(-SR) before!)
#mp_all_plot_all$variable <- factor(  mp_all_plot_all$variable, levels = c(   "CV",   "dCov", "SW",   "SZ",   "RR","SR" ))
MP_shock_plot_all <- ggplot(mp_all_plot_all, aes(x = time, y = value, group = variable)) +
  geom_line(aes(linetype = variable, color = variable)) + 
  facet_wrap(~K, scales = "free_y", labeller = label_value, ncol = 1) +
  scale_x_continuous(breaks = seq(mp_all$time[1] %>% ceiling(),mp_all$time[nrow(mp_all)], 2)) +
  xlab("") + ylab(" ")  + theme_bw() + theme(legend.title = element_blank()) + geom_hline(yintercept = c(-1, 1), linetype = "dashed", alpha = 0.3) 
ggsave("../Plots/MP_shocks.pdf", plot = MP_shock_plot_all, width = 12, height = 6)



# Event analysis ----------------------------------------------------------
# April 1974, October 1979, December 1988 and February 1994 - were contractionary
# December 1990, October 1998, April 2001, and November 2002  -were expansionary 
mp_all_stand_tab <- mp_all_stand %>% filter(time == 1974.25 | time == 1979.75 | time == 1988.75 | time == 1994 | time == 1990.75 | time == 1998.75 | time == 2001.25| time == 2002.75)
#mp_all_stand_tab <- mp_all_stand_tab %>% select(-SR)
mp_all_stand_tab <- mp_all_stand_tab %>% add_column(SUM = rowSums(mp_all_stand_tab[,-1]))

mp_all_stand_tab[c(4,5),] <- mp_all_stand_tab[c(5,4),]
mp_all_stand_tab %>% t() %>% kableExtra::kable(format = "latex", booktabs = T, digits = 3, linesep = "")



# ample -------------------------------------------------------------------
PC1.ts <- iv2ts(PC1, "PC1")
Min.ts <- iv2ts(proxy.min, "Min")
Max.ts <- iv2ts(proxy.max, "Max")
Med.ts <- iv2ts(proxy.med, "Med")

Ample <- function(NonMP, iv){
  K <- dim(NonMP)[2]
  temp <- left_join(NonMP, iv, by = "time") %>% tidyr::drop_na() %>% select(-1)
  temp2 <- lm(temp[,K] ~ temp[,-K] %>% as.matrix()-1) %>% summary()
  return(c(cor(temp)[K, -K], temp2$r.squared,temp2$fstatistic[1])) 
}

Conf_bound <- function(NonMP, iv){
  K <- dim(NonMP)[2]
  temp <- left_join(NonMP, iv, by = "time") %>% tidyr::drop_na() %>% select(-1)
  return(2/sqrt(nrow(temp)))
}


demand_cv   <- solve(CV$B)        %*% t(resid(var3)) %>% t() %>% as.data.frame() %>% select(1)
demand_cv <- cbind(time(var3$y)[-c(1:3)], demand_cv); colnames(demand_cv) <- c("time", "demand_cv")
supply_cv   <- solve(CV$B)        %*% t(resid(var3)) %>% t() %>% as.data.frame() %>% select(2)
supply_cv <- cbind(time(var3$y)[-c(1:3)], supply_cv); colnames(supply_cv) <- c("time", "supply_cv")
finan_cv   <- solve(CV$B)        %*% t(resid(var3)) %>% t() %>% as.data.frame() %>% select(4)
finan_cv <- cbind(time(var3$y)[-c(1:3)], finan_cv); colnames(finan_cv) <- c("time", "finan_cv")

non_mp_cv <- left_join(demand_cv, supply_cv, by = "time") %>% left_join(finan_cv,  by = "time")


ample_tab_cv <- rbind(Ample(non_mp_cv, SW.ts),
                   Ample(non_mp_cv, RR.ts),
                   Ample(non_mp_cv, SZ.ts))
row.names(ample_tab_cv) <- c("SW",  "RR", "SZ")
colnames(ample_tab_cv) <- NULL


ample_tab_cv <- rbind(Ample(non_mp_cv, IV.SW$eps.ts),
                      Ample(non_mp_cv, IV.RR$eps.ts),
                      Ample(non_mp_cv, IV.SZ$eps.ts))
row.names(ample_tab_cv) <- c("SW",  "RR", "SZ")
colnames(ample_tab_cv) <- NULL

kableExtra::kable(t(ample_tab_cv), format = "latex", booktabs = T, digits = 3, linesep = "")


# ample kurz --------------------------------------------------------------

non_mp_cv_s <- left_join(non_mp_cv, IV.SW$eps.ts, by = "time") %>% left_join(IV.RR$eps.ts, by = "time") %>% left_join(IV.SZ$eps.ts, by = "time")  %>% tidyr::drop_na() %>% select(-1)

non_mp_cv_s <- left_join(non_mp_cv, SW.ts, by = "time") %>% 
  left_join(RR.ts, by = "time") %>% left_join(SZ.ts, by = "time") %>% tidyr::drop_na() %>% select(-1)

non_mp_cv_s_tab <- data.frame(cor(non_mp_cv_s[,4], non_mp_cv_s$demand_cv), 
                              cor(non_mp_cv_s[,4], non_mp_cv_s$supply_cv), 
                              cor(non_mp_cv_s[,4], non_mp_cv_s$finan_cv), 
                              summary(lm(non_mp_cv_s[,4] ~ non_mp_cv_s[,1:3] %>% as.matrix()-1))$r.squared,
                              summary(lm(non_mp_cv_s[,4] ~ non_mp_cv_s[,1:3] %>% as.matrix()-1))$fstatistic[1])
colnames(non_mp_cv_s_tab) <- 1:5
for (i in 5:dim(non_mp_cv_s)[2]) {
  temp <- data.frame(cor(non_mp_cv_s[,i], non_mp_cv_s$demand_cv), 
                     cor(non_mp_cv_s[,i], non_mp_cv_s$supply_cv), 
                     cor(non_mp_cv_s[,i], non_mp_cv_s$finan_cv), 
                     summary(lm(non_mp_cv_s[,i] ~ non_mp_cv_s[,1:3] %>% as.matrix()-1))$r.squared,
                     summary(lm(non_mp_cv_s[,i] ~ non_mp_cv_s[,1:3] %>% as.matrix()-1))$fstatistic[1])
  colnames(temp) <- 1:5
  non_mp_cv_s_tab <- rbind(non_mp_cv_s_tab, temp)
}


rownames(non_mp_cv_s_tab) <- NULL
kableExtra::kable(t(non_mp_cv_s_tab), format = "latex", booktabs = T, digits = 3, linesep = "")




# Endogeneity test --------------------------------------------------------
# SW
USA_Qua_SW <- data.frame("time" = USA_Qua %>% time %>% as.numeric, USA_Qua) %>% left_join(iv2ts(SW, "SW"), by = "time") %>% drop_na()
colnames(USA_Qua_SW)[6] <- "SW"
USA_Qua_SW <- ts(USA_Qua_SW %>% select(-1), frequency = 4, start = c(floor(USA_Qua_SW$time[1]), (USA_Qua_SW$time[1] %% 1) / 0.25 + 1 ))
VARselect(USA_Qua_SW, lag.max = 8)
var3_SW <- VAR(USA_Qua_SW, p = 3, type = "const")
CV_SW <- id.cv(var3_SW, SB = c(1984, 1), frequency = 4)
CV_SW
CV
rest_mat <- matrix(c(rep(NA, 20), 0, 0, NA, 0, NA), nrow = 5, ncol = 5, byrow = T)
rest_mat
CV_SW_rest <- id.cv(var3_SW, SB = c(1984, 1), frequency = 4, restriction_matrix = rest_mat)
CV_SW_rest %>% summary
CV_SW_rest$lRatioTest

# RR
USA_Qua_RR <- data.frame("time" = USA_Qua %>% time %>% as.numeric, USA_Qua) %>% left_join(iv2ts(RR, "RR"), by = "time") %>% drop_na()
colnames(USA_Qua_RR)[6] <- "RR"
USA_Qua_RR <- ts(USA_Qua_RR %>% select(-1), frequency = 4, start = c(floor(USA_Qua_RR$time[1]), (USA_Qua_RR$time[1] %% 1) / 0.25 + 1 ))
VARselect(USA_Qua_RR, lag.max = 8)
var3_RR <- VAR(USA_Qua_RR, p = 3, type = "const")
CV_RR <- id.cv(var3_RR, SB = c(1984, 1), frequency = 4)
CV_RR
CV
rest_mat <- matrix(c(rep(NA, 20), 0, 0, NA, 0, NA), nrow = 5, ncol = 5, byrow = T)
rest_mat
CV_RR_rest <- id.cv(var3_RR, SB = c(1984, 1), frequency = 4, restriction_matrix = rest_mat)
CV_RR_rest %>% summary
CV_RR_rest$lRatioTest

# SZ
USA_Qua_SZ <- data.frame("time" = USA_Qua %>% time %>% as.numeric, USA_Qua) %>% left_join(iv2ts(SZ, "SZ"), by = "time") %>% drop_na()
colnames(USA_Qua_SZ)[6] <- "SZ"
USA_Qua_SZ <- ts(USA_Qua_SZ %>% select(-1), frequency = 4, start = c(floor(USA_Qua_SZ$time[1]), (USA_Qua_SZ$time[1] %% 1) / 0.25 + 1 ))
VARselect(USA_Qua_SZ, lag.max = 8)
var3_SZ <- VAR(USA_Qua_SZ, p = 3, type = "const")
CV_SZ <- id.cv(var3_SZ, SB = c(1984, 1), frequency = 4)
CV_SZ
CV
rest_mat <- matrix(c(rep(NA, 20), 0, 0, NA, 0, NA), nrow = 5, ncol = 5, byrow = T)
rest_mat
CV_SZ_rest <- id.cv(var3_SZ, SB = c(1984, 1), frequency = 4, restriction_matrix = rest_mat)
CV_SZ_rest %>% summary
CV_SZ_rest$lRatioTest

Endogen_test <- rbind(CV_SW_rest$lRatioTest, 
                      CV_RR_rest$lRatioTest,
                      CV_SZ_rest$lRatioTest,
                      CV_PC_rest$lRatioTest,
                      CV_min_rest$lRatioTest,
                      CV_max_rest$lRatioTest,
                      CV_med_rest$lRatioTest) %>% t()
colnames(Endogen_test) <- c("SW", "RR", "SZ", "PC1", "Min", "Max", "Median")
kableExtra::kable(Endogen_test, format = "latex", booktabs = T, digits = 3, linesep = "")


# Hist Decomp -------------------------------------------------------------
# pi
hd_pi <- data.frame("SR" = SR %>% cf.mod(series = 2, Partial = 3, Epsname =  c("d", "s", "mp", "f")) %>% select(-y_cf),
                    "dCov" = hd(dCov, series = 2)$hidec[,"Cumulative effect of flow  GBR1 shock on  pi"], 
                    "CV" =   hd(CV, series = 2)$hidec[,"Cumulative effect of flow  GBR1 shock on  pi"])
hd_pi[,c(2,3)] <- hd_pi[,c(3,2)]
colnames(hd_pi)[1:3] <- c("time", "Actual", "SR")

hd_pi <- hd_pi %>% left_join(cf.mod(IV.SW, series = 2, Epsname =  c("d", "s", "mp", "f")) %>% select(-c(y_true, y_cf)), by = "time") %>% 
  left_join(cf.mod(IV.SZ, series = 2, Epsname =  c("d", "s", "mp", "f")) %>% select(-c(y_true, y_cf)) , by = "time") %>% 
  left_join(cf.mod(IV.RR, series = 2, Epsname =  c("d", "s", "mp", "f")) %>% select(-c(y_true, y_cf)), by = "time") %>% tidyr::drop_na()

colnames(hd_pi)[6:8] <- c("SW", "SZ", "RR")

hd_pi <- hd_pi %>% filter(time >= 1979.5 & time <= 1982.5)

hd_pi_plot <- melt(hd_pi[,-2], id = "time")
hd_pi_plot <- hd_pi_plot %>% add_column("actual" = rep( hd_pi[,2], 6))

hd_pi_plot$variable <- factor(hd_pi_plot$variable, levels = c("CV", "dCov", "SW", "SZ", "RR", "SR"))

hd_pi_plot %>% ggplot(aes(x = time, y = value, group = variable)) +
  geom_line(linetype = "longdash", alpha = 0.6) + 
  geom_line(aes(x = time, y = actual, linetype = "solid")) +
  facet_wrap(~variable, scales = "free_x", ncol = 3) + geom_hline(yintercept = 0, color = "red", alpha = 0.7) +
  xlab("") + ylab(" ") + theme_bw() + theme(legend.position = "none")

ggsave("../Plots/HD_pi_4sub_volcker.pdf", plot = last_plot(), width = 12, height = 8)


# bar plot
hd_pi_bar <- hd_pi %>% filter(time >= 1979.5 & time <= 1982.5)

#(hd_pi_bar <- apply(hd_pi_bar, 2, diff) %>% colSums)

(hd_pi_bar <- hd_pi_bar[nrow(hd_pi_bar),]- hd_pi_bar[1,])

SUM <- hd_pi_bar %>% select(-c(time, Actual)) %>% sum

hd_pi_bar <- hd_pi_bar %>% mutate(SR = SR/SUM, 
                                  dCov = dCov/SUM,
                                  CV = CV/SUM,
                                  SW = SW/SUM,
                                  SZ = SZ/SUM,
                                  RR = RR/SUM) %>% select(-c(time, Actual))

melt(hd_pi_bar, id = NULL) %>% ggplot(aes(x = variable, y = value, fill = variable)) + 
  geom_bar(stat="identity", alpha = .8) + scale_fill_brewer(palette="Accent") +
  xlab("") + ylab(" ") + theme_bw() + theme(legend.position = "none")

ggsave("../Plots/HD_pi_4bar_volcker.pdf", plot = last_plot(), width = 12, height = 6)


# Other instruments -------------------------------------------------------


# RR83b
RR83b.raw      <- read.csv("Data/Instruments/RR83b.csv")
RR83b          <- ts(RR83b.raw[,-1], frequency = 4, start = c(floor(RR83b.raw$Year[1]), (RR83b.raw$Year[1] %% 1) / 0.25 + 1 ))
IV.RR83b       <- get.id.iv(var3, instruments = RR83b, Synchro = T)
IV.RR83b$B
IV.RR83b$F_test

# FF1vr
FF1vr.raw      <- read.csv("Data/Instruments/FF1vr.csv")
FF1vr          <- ts(FF1vr.raw[,-1], frequency = 4, start = c(floor(FF1vr.raw$Year[1]), (FF1vr.raw$Year[1] %% 1) / 0.25 + 1 ))
IV.FF1vr       <- get.id.iv(var3, instruments = FF1vr, Synchro = T)
IV.FF1vr$B
IV.FF1vr$F_test

# FF4vr
FF4vr.raw      <- read.csv("Data/Instruments/FF4vr.csv")
FF4vr          <- ts(FF4vr.raw[,-1], frequency = 4, start = c(floor(FF4vr.raw$Year[1]), (FF4vr.raw$Year[1] %% 1) / 0.25 + 1 ))
IV.FF4vr       <- get.id.iv(var3, instruments = FF4vr, Synchro = T)
IV.FF4vr$B
IV.FF4vr$F_test

# FF1gb
FF1gb.raw      <- read.csv("Data/Instruments/FF1gb.csv")
FF1gb          <- ts(FF1gb.raw[,-1], frequency = 4, start = c(floor(FF1gb.raw$Year[1]), (FF1gb.raw$Year[1] %% 1) / 0.25 + 1 ))
IV.FF1gb       <- get.id.iv(var3, instruments = FF1gb, Synchro = T)
IV.FF1gb$B
IV.FF1gb$F_test


# FF4gb
FF4gb.raw      <- read.csv("Data/Instruments/FF4gb.csv")
FF4gb          <- ts(FF4gb.raw[,-1], frequency = 4, start = c(floor(FF4gb.raw$Year[1]), (FF4gb.raw$Year[1] %% 1) / 0.25 + 1 ))
IV.FF4gb       <- get.id.iv(var3, instruments = FF4gb, Synchro = T)
IV.FF4gb$B
IV.FF4gb$F_test


# ED2vr
ED2vr.raw      <- read.csv("Data/Instruments/ED2vr.csv")
ED2vr          <- ts(ED2vr.raw[,-1], frequency = 4, start = c(floor(ED2vr.raw$Year[1]), (ED2vr.raw$Year[1] %% 1) / 0.25 + 1 ))
IV.ED2vr       <- get.id.iv(var3, instruments = ED2vr, Synchro = T)
IV.ED2vr$B
IV.ED2vr$F_test


# ED2gb
ED2gb.raw      <- read.csv("Data/Instruments/ED2gb.csv")
ED2gb          <- ts(ED2gb.raw[,-1], frequency = 4, start = c(floor(ED2gb.raw$Year[1]), (ED2gb.raw$Year[1] %% 1) / 0.25 + 1 ))
IV.ED2gb       <- get.id.iv(var3, instruments = ED2gb, Synchro = T)
IV.ED2gb$B
IV.ED2gb$F_test


# MA
MA.raw      <- read.csv("Data/Instruments/MA.csv")
MA          <- ts(MA.raw[,-1], frequency = 4, start = c(floor(MA.raw$Year[1]), (MA.raw$Year[1] %% 1) / 0.25 + 1 ))
IV.MA       <- get.id.iv(var3, instruments = MA, Synchro = T)
IV.MA$B
IV.MA$F_test


# data driven boot sign ---------------------------------------------------


extract_B0_boot <- function(x){
  irfs <- lapply(x$bootstrap, '[[', 'irf')
  B0 <- irfs[[1]] %>% filter(V1 == 1)
  K <- sqrt(ncol(B0)-1)
  for (i in 2:length(irfs)) {
    B0 <- rbind(B0, irfs[[i]] %>% filter(V1 == 1))
  }
  return(B0)
}


CV.B0.boot <- CV.boot %>% extract_B0_boot
CV.B0.boot %>% colnames() 
colnames(CV.B0.boot)[c(8, 12)] <- c("p", "r")
CV.B0.boot <- CV.B0.boot %>% mutate(labling = ifelse(p*r<0, 1, 0))
CV.B0.boot$labling %>% mean


dCov.B0.boot <- dCov.boot %>% extract_B0_boot
dCov.B0.boot %>% colnames() 
colnames(dCov.B0.boot)[c(8, 12)] <- c("p", "r")
dCov.B0.boot <- dCov.B0.boot %>% mutate(labling = ifelse(p*r<0, 1, 0))
dCov.B0.boot$labling %>% mean



