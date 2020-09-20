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
USA_Tri_raw   <- read.csv("Data/USA_Tri.csv")
USA_Tri       <- ts(USA_Tri_raw[,-1], frequency = 4, start = c(floor(USA_Tri_raw$Year[1]), (USA_Tri_raw$Year[1] %% 1) / 0.25 + 1 ))
#autoplot(USA_Tri, ncol = 1) + theme_bw()

# VAR Spesicifation -------------------------------------------------------
## lag length
#VARselect(USA_Tri, lag.max = 8)
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

## Identification by minimizing distance covariance
set.seed(1234)
dCov        <- id.dc(var4)
dCov$B      <- dCov$B[,c(3,2,1)]
dCov$B[,3]  <- dCov$B[,3]*-1
dCov        %>% summary

## Identification by minimizing CvM distance
#set.seed(1234)
#CvM        <- id.cvm(var4, itermax = 2000, steptol = 200, iter2 = 200)
#CvM$B <- CvM$B[,c(3,2,1)]
#CvM$B[,1] <- CvM$B[,1] *-1
#CvM        %>% summary
#CvM %>% irf() %>% plot
#ggsave("../Plots/Tri_CvM.pdf", plot = last_plot(), width = 12, height = 8)



## Proxy SVAR
# Smets-Wouters (2007)
SW.raw      <- read.csv("Data/Instruments/SW.csv")
SW          <- ts(SW.raw[,-1], frequency = 4, start = c(floor(SW.raw$Year[1]), (SW.raw$Year[1] %% 1) / 0.25 + 1 ))
IV.SW       <- get.id.iv(var4, instruments = SW, Synchro = T)
IV.SW$B
IV.SW$F_test
#Box.test(SW, lag = 1, type = "Ljung-Box")
#Box.test(SW, lag = 2, type = "Ljung-Box")
#Box.test(SW, lag = 3, type = "Ljung-Box")

# Guerkaynak, Sack and Swanson (2005)
GSS.raw     <- read.csv("Data/Instruments/GSS.csv")
GSS         <- ts(GSS.raw[,-1], frequency = 4, start = c(floor(GSS.raw$Year[1]), (GSS.raw$Year[1] %% 1) / 0.25 + 1 ))
IV.GSS      <- get.id.iv(var4, instruments = GSS, Synchro = T)
IV.GSS$B
IV.GSS$F_test

# Romer Romer (2004)
RR.raw      <- read.csv("Data/Instruments/RR.csv")
RR          <- ts(RR.raw[,-1], frequency = 4, start = c(floor(RR.raw$Year[1]), (RR.raw$Year[1] %% 1) / 0.25 + 1 ))
RR          <- Remove_AC(RR, lag = 1) 
IV.RR       <- get.id.iv(var4, instruments = RR, Synchro = T)
IV.RR$B
IV.RR$F_test
#Box.test(RR, lag = 1, type = "Ljung-Box")
#Box.test(RR, lag = 2, type = "Ljung-Box")
#Box.test(RR, lag = 3, type = "Ljung-Box")

# Sims Zha (2006)
SZ.raw      <- read.csv("Data/Instruments/SZ.csv")
SZ          <- ts(SZ.raw[,-1], frequency = 4, start = c(floor(SZ.raw$Year[1]), (SZ.raw$Year[1] %% 1) / 0.25 + 1 ))
SZ          <- Remove_AC(SZ, 2)
IV.SZ       <- get.id.iv(var4, instruments = SZ, Synchro = T)
IV.SZ$B
IV.SZ$F_test
#Box.test(SZ, lag = 1, type = "Ljung-Box")
#Box.test(SZ, lag = 2, type = "Ljung-Box")
#Box.test(SZ, lag = 3, type = "Ljung-Box")
#Box.test(SZ, lag = 4, type = "Ljung-Box")


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
sisgnrest <- list(demand = c(1,1,1), supply = c(-1,1,1), money = c(-1,-1,1))
set.seed(1234)
CV.boot     <- svars::wild.boot(CV, n.ahead = 20, nboot = 2000, signrest = signrest)
#saveRDS(CV.boot, file = "Output/CV_boot3.rds")
set.seed(1234)
dCov.boot   <- svars::wild.boot(dCov, n.ahead = 20, nboot = 2000, signrest = signrest)
#saveRDS(dCov.boot, file = "Output/dCov_boot3.rds")

CV.boot <- readRDS("Output/CV_boot3.rds")
dCov.boot <- readRDS("Output/dCov_boot3.rds")

CV.boot %>% summary
dCov.boot %>% summary

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
colnames(CV.boot$true$irf)[7] <- "epsilon[ mp ] %->% pi"
colnames(CV.boot$true$irf)[10] <- "epsilon[ mp ] %->% GBR1"

colnames(dCov.boot$true$irf)[4] <- "epsilon[ mp ] %->% x"
colnames(dCov.boot$true$irf)[7] <- "epsilon[ mp ] %->% pi"
colnames(dCov.boot$true$irf)[10] <- "epsilon[ mp ] %->% GBR1"


# Normalization 1MP->25bpsGBR1 --------------------------------------------
SR$IRF.plot %<>% mutate(value = ifelse(Label == "median", 0.25*value/SR$IRF.M[3,3,1], ifelse(Label == "median target (MT)", 0.25*value/SR$IRF.MT[3,3,1], NA)))
SR$IRF.plot %<>% mutate(L =  0.25*L/SR$IRF.M[3,3,1])
SR$IRF.plot %<>% mutate(U =  0.25*U/SR$IRF.M[3,3,1])


SR.MP     <- plot.my.irf(My.id = SR, Partial = 3, Hide.Legend = T, Color.manual = c("black", "blue"), Lty.manual = c("solid", "longdash"), 
                         Epsname =  c("d", "s", "mp")) + xlab("") + ylab("")
CV.MP     <- CV.boot %>% plot.sboot.mod(lowerq = 0.16, upperq = 0.84, Partial = 3, Normalize = CV$B[3,3]/0.25) + xlab("") + ylab("") 
dCov.MP   <- dCov.boot %>% plot.sboot.mod(lowerq = 0.16, upperq = 0.84, Partial = 3, Normalize = dCov$B[3,3]/0.25) + xlab("") + ylab("")
IRF_nonIV <- grid.arrange(SR.MP, CV.MP, dCov.MP, ncol = 3)
ggsave("../Plots/Tri_A.pdf", plot = IRF_nonIV, width = 12, height = 8)

SW.MP     <- IV.SW.boot %>% plot.my.irf(Partial = 3, Epsname = c("d", "s", "mp"), Normalize = IV.SW$B[3]/0.25) + xlab("") + ylab("")
RR.MP     <- IV.RR.boot %>% plot.my.irf(Partial = 3, Epsname = c("d", "s", "mp"), Normalize = IV.RR$B[3]/0.25) + xlab("") + ylab("")
SZ.MP     <- IV.SZ.boot %>% plot.my.irf(Partial = 3, Epsname = c("d", "s", "mp"), Normalize = IV.SZ$B[3]/0.25) + xlab("") + ylab("")
IRF_IV <- grid.arrange(SW.MP, RR.MP, SZ.MP, ncol = 3)
ggsave("../Plots/Tri_B_MBB.pdf", plot = IRF_IV, width = 12, height = 8)

IV.SW_pro.boot %>% plot.my.irf(Partial = 3) + xlab("") + ylab("")

# Correlation of MP shocks ------------------------------------------------
mp_sr     <- solve(SR$IRF.MT[,,1]) %*% t(resid(var4)) %>% t() %>% as.data.frame() %>% select(3)
mp_cv     <- solve(CV$B)          %*% t(resid(var4)) %>% t() %>% as.data.frame() %>% select(3)
mp_dcov   <- solve(dCov$B)        %*% t(resid(var4)) %>% t() %>% as.data.frame() %>% select(3)
mp_noniv  <- cbind(mp_sr, mp_cv, mp_dcov); colnames(mp_noniv) <- c("SR", "CV", "dCov")
mp_iv     <- left_join(IV.SW$eps.ts, IV.SZ$eps.ts, by = "time")
mp_iv     <- left_join(mp_iv, IV.RR$eps.ts, by = "time") %>% tidyr::drop_na()
colnames(mp_iv) <- c("time", "SW", "SZ", "RR")
mp_noniv  <- cbind(time(var4$y)[-c(1:4)], mp_noniv); colnames(mp_noniv)[1] <- "time"
mp_all    <- left_join(mp_noniv, mp_iv, by = "time") %>% tidyr::drop_na() 
#mp_all <- mp_all %>% select(-1)

cor(mp_all) %>% kableExtra::kable(format = "latex", booktabs = T, digits = 2, linesep = "")

mp_sds <- mp_all %>% select(-1) %>% apply( 2, sd)

mp_all_stand <- mp_all

for (j in 1:dim(mp_all)[1]) {
  for (i in 2:7) {
    mp_all_stand[j,i] <- mp_all[j,i] / mp_sds[i-1]
  } 
}

mp_all_stand %>% select(-1) %>% apply( 2, sd)


mp_all_plot3 <- reshape2::melt( mp_all_stand, id.var = "time")

library(ggplot2)
mp_all_plot3$variable <- factor(  mp_all_plot3$variable, levels = c(   "CV",   "dCov", "SW",   "SZ",   "RR","SR" ))
MP_shock_plot3 <- ggplot(mp_all_plot3, aes(x = time, y = value, group = variable)) +
  geom_line(aes(linetype = variable, color = variable), size = 0.7) + 
  scale_x_continuous(breaks = seq(mp_all_stand$time[1] %>% ceiling(),mp_all_stand$time[nrow(mp_all_stand)], 2)) + 
  xlab("") + ylab(" ")  + theme_bw() + theme(legend.title = element_blank()) + 
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", alpha = 0.3) + scale_color_ochre(palette="lorikeet") 
  #annotate("rect", xmin = Start, xmax = End2, ymin = -Inf, ymax = Inf, alpha = 0.3) + ylim(-4.5, 4.5)
MP_shock_plot3
ggsave("../Plots/MP_shocks3.pdf", plot = MP_shock_plot3, width = 12, height = 4)


MP_shocks_plotSW_Fa_df <- IV.SW$eps.ts %>% left_join(FaIV.SW$eps.ts, by = "time") %>% filter(time >= mp_all_stand$time[1] & time <= mp_all_stand$time[nrow(mp_all_stand)])
MP_shocks_plotSW_Fa_df[,-1] <- scale(MP_shocks_plotSW_Fa_df[,-1], center = F)
colnames(MP_shocks_plotSW_Fa_df)[-1] <- c("VAR", "FAVAR")
MP_shocks_plotSW_Fa_df_long <-  reshape2::melt(MP_shocks_plotSW_Fa_df, id.var = "time")

MP_shocks_plotSW_Fa <- ggplot(MP_shocks_plotSW_Fa_df_long, aes(x = time, y = value, group = variable)) +
  geom_line(aes(linetype = variable)) + scale_linetype_manual(values = c("solid","longdash")) +
  scale_x_continuous(breaks = seq(mp_all_stand$time[1] %>% ceiling(),mp_all_stand$time[nrow(mp_all_stand)], 2)) +
  xlab("") + ylab(" ")  + theme_bw()  + geom_hline(yintercept = c(-1, 1), linetype = "dashed", alpha = 0.3) + 
  #annotate("rect", xmin = Start, xmax = End2, ymin = -Inf, ymax = Inf, alpha = 0.3) + 
  theme(legend.title = element_blank()) + ylim(-4.5, 4.5)

new_shock_plots <- grid.arrange(MP_shock_plot3, MP_shocks_plotSW_Fa, ncol = 1)

ggsave("../Plots/MP_shocks3_Fa3_new.pdf", plot = new_shock_plots, width = 12, height = 6)


# Event analysis ----------------------------------------------------------
# April 1974, October 1979, December 1988 and February 1994 - were contractionary
# December 1990, October 1998, April 2001, and November 2002  -were expansionary 
mp_all_tab <- mp_all %>% filter(time == 1974.25 | time == 1979.75 | time == 1988.75 | time == 1994 | time == 1990.75 | time == 1998.75 | time == 2001.25| time == 2002.75)
#mp_all_stand_tab <- mp_all_stand_tab %>% select(-SR)
mp_all_tab <- mp_all_tab %>% add_column(SUM = rowSums(mp_all_tab[,-1]))
  
mp_all_tab[c(4,5),] <- mp_all_tab[c(5,4),]
mp_all_tab %>% t() %>% kableExtra::kable(format = "latex", booktabs = T, digits = 3, linesep = "")


# ample -------------------------------------------------------------------

demand_cv   <- solve(CV$B)        %*% t(resid(var4)) %>% t() %>% as.data.frame() %>% select(1)
demand_cv <- cbind(time(var4$y)[-c(1:4)], demand_cv); colnames(demand_cv) <- c("time", "demand_cv")
supply_cv     <- solve(CV$B)          %*% t(resid(var4)) %>% t() %>% as.data.frame() %>% select(2)
supply_cv <- cbind(time(var4$y)[-c(1:4)], supply_cv); colnames(supply_cv) <- c("time", "supply_cv")


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


non_mp_cv <- left_join(demand_cv, supply_cv, by = "time") 
#PC1.ts <- iv2ts(PC1, "PC1")

ample_tab_cv <- rbind(Ample(non_mp_cv, IV.SW$eps.ts),
                      Ample(non_mp_cv, IV.RR$eps.ts),
                      Ample(non_mp_cv, IV.SZ$eps.ts))
Conf_tab <- rbind(Conf_bound(non_mp_cv, IV.SW$eps.ts),
                  Conf_bound(non_mp_cv, IV.RR$eps.ts),
                  Conf_bound(non_mp_cv, IV.SZ$eps.ts))



row.names(ample_tab_cv) <- c("SW", "RR", "SZ")
colnames(ample_tab_cv) <- NULL

ample_tab <- cbind(Conf_tab, ample_tab_cv)
kableExtra::kable(t(ample_tab), format = "latex", booktabs = T, digits = 3, linesep = "")




# ample kurz --------------------------------------------------------------

#non_mp_cv_s <- left_join(non_mp_cv, SW.ts, by = "time") %>% left_join(RR.ts, by = "time") %>% left_join(SZ.ts, by = "time")  %>% tidyr::drop_na() %>% select(-1)
non_mp_cv_s <- left_join(non_mp_cv, IV.SW$eps.ts, by = "time") %>% left_join(IV.RR$eps.ts, by = "time") %>% left_join(IV.SZ$eps.ts, by = "time")  %>% tidyr::drop_na() %>% select(-1)


non_mp_cv_s_tab <- data.frame(cor(non_mp_cv_s[,3], non_mp_cv_s$demand_cv), 
                              cor(non_mp_cv_s[,3], non_mp_cv_s$supply_cv), 
                              summary(lm(non_mp_cv_s[,3] ~ non_mp_cv_s[,1:2] %>% as.matrix()-1))$r.squared,
                              summary(lm(non_mp_cv_s[,3] ~ non_mp_cv_s[,1:2] %>% as.matrix()-1))$fstatistic[1])
colnames(non_mp_cv_s_tab) <- 1:4
for (i in 4:dim(non_mp_cv_s)[2]) {
  temp <- data.frame(cor(non_mp_cv_s[,i], non_mp_cv_s$demand_cv), 
                     cor(non_mp_cv_s[,i], non_mp_cv_s$supply_cv), 
                     summary(lm(non_mp_cv_s[,i] ~ non_mp_cv_s[,1:2] %>% as.matrix()-1))$r.squared,
                     summary(lm(non_mp_cv_s[,i] ~ non_mp_cv_s[,1:2] %>% as.matrix()-1))$fstatistic[1])
  colnames(temp) <- 1:4
  non_mp_cv_s_tab <- rbind(non_mp_cv_s_tab, temp)
}


rownames(non_mp_cv_s_tab) <- NULL
kableExtra::kable(t(non_mp_cv_s_tab), format = "latex", booktabs = T, digits = 3, linesep = "")


# Waschmaschine -----------------------------------------------------------
demand_cv   <- solve(CV$B)        %*% t(resid(var4)) %>% t() %>% as.data.frame() %>% select(1)
demand_cv <- cbind(time(var4$y)[-c(1:4)], demand_cv); colnames(demand_cv) <- c("time", "demand_cv")
supply_cv     <- solve(CV$B)          %*% t(resid(var4)) %>% t() %>% as.data.frame() %>% select(2)
supply_cv <- cbind(time(var4$y)[-c(1:4)], supply_cv); colnames(supply_cv) <- c("time", "supply_cv")

non_mp_cv <- left_join(demand_cv, supply_cv, by = "time")

## correlation btw. SR shocks and CV non_MP shocks
#temp = cbind(mp_sr, non_mp_cv %>% select(-1))
#colnames(temp)[1] = "SR"
#cor(temp)
#lm(SR~demand_cv + supply_cv - 1, data = temp) %>% summary


#temp <- left_join(non_mp_cv, SW.ts, by = "time") %>% tidyr::drop_na()
#SW_pro <- data.frame("time" = temp$time, "SW_pro" =  lm(temp$SW ~ temp[,2:3] %>% as.matrix()-1) %>% resid())
#SW_pro <- ts(SW_pro$SW_pro, frequency = 4, start = c(floor(SW_pro$time[1]), (SW_pro$time[1] %% 1) / 0.25 + 1 ))

#IV.SW_pro       <- get.id.iv(var4, instruments = SW_pro, Synchro = T)
#IV.SW_pro$B
#IV.SW_pro$F_test

# RR_pro
temp <- left_join(non_mp_cv, RR.ts, by = "time") %>% tidyr::drop_na()
lm1 <-  lm(temp$RR ~ temp[,2:3] %>% as.matrix()-1)

RR_pro <- data.frame("time" = temp$time, "RR_pro" =  lm1 %>% resid())
RR_pro <- ts(RR_pro$RR_pro, frequency = 4, start = c(floor(RR_pro$time[1]), (RR_pro$time[1] %% 1) / 0.25 + 1 ))



IV.RR_pro <- get.id.iv(var4, instruments = RR_pro, Synchro = T)
IV.RR_pro$B
IV.RR_pro$F_test

#SZ_pro
temp <- left_join(non_mp_cv, SZ.ts, by = "time") %>% tidyr::drop_na()
lm2 <- lm(temp$SZ ~ temp[,2:3] %>% as.matrix()-1)

SZ_pro <- data.frame("time" = temp$time, "SZ_pro" =  lm2 %>% resid())
SZ_pro <- ts(SZ_pro$SZ_pro, frequency = 4, start = c(floor(SZ_pro$time[1]), (SZ_pro$time[1] %% 1) / 0.25 + 1 ))

IV.SZ_pro <- get.id.iv(var4, instruments = SZ_pro, Synchro = T)
IV.SZ_pro$B
IV.SZ_pro$F_test

#xtable(summary(lm1), type="latex")
#xtable(summary(lm2), type="latex")


#sjPlot::tab_model(lm1,lm2, file="output.tex")

# Cormat
new_cormat <- mp_all %>% left_join(IV.RR_pro$eps.ts, by = "time") %>% left_join(IV.SZ_pro$eps.ts, by = "time") 
colnames(new_cormat)[c(8,9)] <- c("RR_pro", "SZ_pro")
new_cormat %>% select(-1) %>% cor %>% round(3) %>% kableExtra::kable(format = "latex", booktabs = T, digits = 3, linesep = "")





# shock plot
mp_all_pro <- mp_all %>% left_join(IV.RR_pro$eps.ts, by = "time") %>% left_join(IV.SZ_pro$eps.ts, by = "time")
colnames(mp_all_pro)[c(8,9)] <- c("RR_pro", "SZ_pro")

# including cv demand and supply shocks
mp_all_pro %<>% left_join(non_mp_cv, by = "time") 

#temp <- mp_all_pro %>% select(-c(SR, dCov, CV)) %>%  filter(time >= Start & time <= End2) #%>% kableExtra::kable(format = "latex", booktabs = T, digits = 3, linesep = "")
#temp %<>% abs  
#temp_logi <- temp>qnorm(0.95) 
#temp_logi %>% colSums()

mp_sds_pro <- mp_all_pro %>% select(-1) %>% apply( 2, sd)
mp_all_pro_stand <- mp_all_pro
for (j in 1:dim(mp_all_pro)[1]) {
  for (i in 2:ncol(mp_all_pro)) {
    mp_all_pro_stand[j,i] <- mp_all_pro[j,i] / mp_sds_pro[i-1]
  } 
}

mp_all_pro_stand %>% select(-1) %>% apply( 2, sd)



# Event analysis ----------------------------------------------------------
# April 1974, October 1979, December 1988 and February 1994 - were contractionary
# December 1990, October 1998, April 2001, and November 2002  -were expansionary 
mp_pros <- IV.RR_pro$eps.ts %>% left_join(IV.SZ_pro$eps.ts, by = "time")
colnames(mp_pros)[-1] <- c("RR", "SZ")
mp_pro_tab <- mp_pros %>% filter(time == 1974.25 | time == 1979.75 | time == 1988.75 | time == 1994 | time == 1990.75 | time == 1998.75 | time == 2001.25| time == 2002.75)
#mp_all_stand_tab <- mp_all_stand_tab %>% select(-SR)
#mp_all_pro_stand <- mp_all_pro_stand %>% add_column(SUM = rowSums(mp_all_pro_stand[,-1]))

mp_pro_tab[c(4,5),] <- mp_pro_tab[c(5,4),]
mp_pro_tab %>% t() %>% kableExtra::kable(format = "latex", booktabs = T, digits = 3, linesep = "")



mps_RR_pro <- mp_all_pro_stand %>% select(time, old=RR, new=RR_pro) %>% melt(id.var = "time")
mps_SZ_pro <- mp_all_pro_stand %>% select(time, old=SZ, new=SZ_pro) %>% melt(id.var = "time")


mps_plot.SW <- ggplot(mp_all_pro_stand, aes(x = time, y = SW)) + geom_line() + ylim(-4.5,  4.5)+
  scale_x_continuous(breaks = seq(mp_all_pro_stand$time[1] %>% ceiling(),mp_all_pro_stand$time[nrow(mp_all)], 2)) +
  xlab("") + ylab(" ")  + theme_bw()  + geom_hline(yintercept = c(-1, 1), linetype = "dashed", alpha = 0.3) + 
  annotate("rect", xmin = 1979.5, xmax = 1983, ymin = -Inf, ymax = Inf, alpha = 0.3)

mps_plot.RR <- ggplot(mps_RR_pro, aes(x = time, y = value, group = variable)) + ylim(-4.5,  4.5)+
  geom_line(aes(linetype = variable)) + scale_linetype_manual(values = c("longdash", "solid")) +
  scale_x_continuous(breaks = seq(mp_all_pro_stand$time[1] %>% ceiling(),mp_all_pro_stand$time[nrow(mp_all)], 2)) +
  xlab("") + ylab(" ")  + theme_bw()  + geom_hline(yintercept = c(-1, 1), linetype = "dashed", alpha = 0.3) + 
  theme(legend.position = "none") + annotate("rect", xmin = 1979.5, xmax = 1983, ymin = -Inf, ymax = Inf, alpha = 0.3)

mps_plot.SZ <- ggplot(mps_SZ_pro, aes(x = time, y = value, group = variable)) + ylim(-4.5,  4.5)+
  geom_line(aes(linetype = variable)) + scale_linetype_manual(values = c("longdash", "solid")) +
  scale_x_continuous(breaks = seq(mp_all_pro_stand$time[1] %>% ceiling(),mp_all_pro_stand$time[nrow(mp_all)], 2)) +
  xlab("") + ylab(" ")  + theme_bw()  + geom_hline(yintercept = c(-1, 1), linetype = "dashed", alpha = 0.3) + 
  theme(legend.position = "none") + annotate("rect", xmin = 1979.5, xmax = 1983, ymin = -Inf, ymax = Inf, alpha = 0.3)

mps_plot_pros <- grid.arrange(mps_plot.SW, mps_plot.RR, mps_plot.SZ)
ggsave("../Plots/MP_shocks_pro3.pdf", plot = mps_plot_pros, width = 12, height = 6)


non_mp_cv_stand <- data.frame("time" =  non_mp_cv$time, scale(non_mp_cv[,-1]))
all_shocks.CV <- mp_all_pro_stand %>% left_join(non_mp_cv_stand, by = "time")

all_shocks_plot.CV <- ggplot(all_shocks.CV %>% select(time, AD = demand_cv, AS = supply_cv) %>% melt(id.var = "time"), aes(x = time, y = value, group = variable)) + ylim(-4.5,  4.5)+
  geom_line(aes(linetype = variable)) + scale_linetype_manual(values = c("longdash", "dotted")) +
  scale_x_continuous(breaks = seq(all_shocks.CV$time[1] %>% ceiling(),mp_all_pro_stand$time[nrow(mp_all)], 2)) +
  xlab("") + ylab(" ")  + theme_bw()  + geom_hline(yintercept = c(-1, 1), linetype = "dashed", alpha = 0.3) + 
  theme(legend.position = "none") + annotate("rect", xmin = 1979.5, xmax = 1983, ymin = -Inf, ymax = Inf, alpha = 0.3)
all_plot_pros <- grid.arrange(mps_plot.SW, mps_plot.RR, mps_plot.SZ, all_shocks_plot.CV, ncol = 1)
ggsave("../Plots/MP_shocks_pro3_mit_CVnonMP.pdf", plot = all_plot_pros, width = 12, height = 8)

# Volcker -----------------------------------------------------------------
Start = 1979.75
End1 = 1982.5
End2 = 1983.5

non_mp_cv %>% apply(2, sd)

all_shocks.CV %>% filter(time >= Start & time <= End2) %>% select(time, SW, RR, RR_pro, SZ, SZ_pro, demand_cv, supply_cv) %>% kableExtra::kable(format = "latex", booktabs = T, digits = 3, linesep = "")


df_pro <- data.frame("period" = c("1979Q4 - 1982Q3", "1979Q4 - 1983Q3"),
                  "SR" = c(quasi.hd(SR, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End1, Freq = 4),
                           quasi.hd(SR, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End2, Freq = 4)),
                  "dCov" = c(quasi.hd(dCov, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End1, Freq = 4),
                             quasi.hd(dCov, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End2, Freq = 4)),
                  "CV" = c(quasi.hd(CV, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End1, Freq = 4),
                           quasi.hd(CV, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End2, Freq = 4)),
                  "SW" = c(quasi.hd(IV.SW, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End1, Freq = 4),
                           quasi.hd(IV.SW, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End2, Freq = 4)),
                  "SZ" = c(quasi.hd(IV.SZ, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End1, Freq = 4),
                           quasi.hd(IV.SZ, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End2, Freq = 4)),
                  "RR" = c(quasi.hd(IV.RR, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End1, Freq = 4),
                           quasi.hd(IV.RR, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End2, Freq = 4)),
                  "SZ_pro" = c(quasi.hd(IV.SZ_pro, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End1, Freq = 4),
                           quasi.hd(IV.SZ_pro, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End2, Freq = 4)),
                  "RR_pro" = c(quasi.hd(IV.RR_pro, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End1, Freq = 4),
                           quasi.hd(IV.RR_pro, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End2, Freq = 4))) %>% melt(id = c("period"))
df_pro %<>% add_column(lab = "old") %>% mutate(lab = ifelse(variable == "SZ_pro" | variable == "RR_pro", "new", "old"))
df_pro[13:14,]$variable <- "SZ"
df_pro[15:16,]$variable <- "RR"

df_pro %<>% mutate(value = value/4)

df_pro %>% ggplot(aes(x = variable, y = value)) + 
  geom_hline(yintercept = 0, color = "red", alpha = 0.8) +
  facet_wrap(~period) +
  geom_bar(data=subset(df_pro, lab == "old"), aes(fill = variable),
           stat="identity", alpha = .7, width=0.5,  color="#999999", linetype = "solid", size=0.8) + 
  geom_bar(data=subset(df_pro, lab == "new"), aes(fill = variable),
           stat="identity", alpha = .4, width=0.5,  color="#999999", linetype = "dashed", size=0.8) + 
  xlab("") + ylab(" ") + theme_bw() + theme(legend.position = "none") + scale_fill_ochre(palette="parliament")

ggsave("../Plots/inflation_fighter_pro3.pdf", plot = last_plot(), width = 12, height = 4)






# exogen
ample_tab_pro <- rbind(Ample(non_mp_cv, IV.RR_pro$eps.ts),
                       Ample(non_mp_cv, IV.SZ_pro$eps.ts))



# MBB
set.seed(1234)
MP_RR_pro.boot <-  get.MBB.fixed(IV.RR_pro, design = "fixed", Step = 15, nboot = 2000, Ci = 0.68)

set.seed(1234)
MP_SZ_pro.boot <-  get.MBB.fixed(IV.SZ_pro, design = "fixed", Step = 15, nboot = 2000, Ci = 0.68) 


# IRFs
MP_RR_pro <-  MP_RR_pro.boot %>% plot.my.irf(Partial = 3, Epsname = c("d", "s", "mp"), Normalize = IV.RR_pro$B[3]/0.25) + xlab("") + ylab("")
MP_SZ_pro <-  MP_SZ_pro.boot %>% plot.my.irf(Partial = 3, Epsname = c("d", "s", "mp"), Normalize = IV.RR_pro$B[3]/0.25) + xlab("") + ylab("")


pros <- grid.arrange(CV.MP, SW.MP, MP_RR_pro, MP_SZ_pro, ncol = 4)
ggsave("../Plots/pros.pdf", plot = pros, width = 12, height = 6)

# IV results
tab_iv_id <- cbind(IV_id.tab(IV.SW, IV.SW.boot),
                   IV_id.tab(IV.RR, IV.RR.boot),
                   IV_id.tab(IV.SZ, IV.SZ.boot),
                   IV_id.tab(IV.RR_pro, MP_RR_pro.boot),
                   IV_id.tab(IV.SZ_pro, MP_SZ_pro.boot))

rownames(tab_iv_id) <- NULL
kableExtra::kable(tab_iv_id, format = "latex", booktabs = T, digits = 3, linesep = "")






# SW_pro
USA_Tri_SW_pro <- data.frame("time" = USA_Tri %>% time %>% as.numeric, USA_Tri) %>% left_join(iv2ts(SW_pro, "SW_pro"), by = "time") %>% drop_na()
colnames(USA_Tri_SW_pro)[5] <- "SW_pro"
USA_Tri_SW_pro <- ts(USA_Tri_SW_pro %>% select(-1), frequency = 4, start = c(floor(USA_Tri_SW_pro$time[1]), (USA_Tri_SW_pro$time[1] %% 1) / 0.25 + 1 ))

VARselect(USA_Tri_SW_pro, lag.max = 8)

var_SW_pro <- VAR(USA_Tri_SW_pro, p = 4, type = "const")
#chow.test(var_SW_pro, SB = c(1984, 1), frequency = 4, nboot = 2000) %>% summary 

CV_SW_pro <- id.cv(var_SW_pro, SB = c(1984, 1), frequency = 4)
CV_SW_pro

rest_mat <- matrix(c(rep(NA, 12), 0, 0, NA, NA), nrow = 4, ncol = 4, byrow = T)
rest_mat
CV_SW_pro_rest <- id.cv(var_SW_pro, SB = c(1984, 1), frequency = 4, restriction_matrix = rest_mat)
#CV_SW_pro_rest %>% summary
CV_SW_pro_rest$lRatioTest


var_SW_pro_2 <- VAR(USA_Tri_SW_pro[,c(4,1,2,3)], p = 4, type = "const")
Chol_SW_pro <- id.chol(var_SW_pro_2)
Chol_SW_pro$B

set.seed(1234)
MP_IV_SW_pro <-  get.MBB.fixed(IV.SW_pro, design = "fixed", Step = 15, nboot = 2000, Ci = 0.68) %>% plot.my.irf(Partial = 3) + xlab("") + ylab("")
ggsave("../Plots/Tri_IV_SW_pro.pdf", plot = MP_IV_SW_pro, width = 12, height = 8)

set.seed(1234)
MP_Chol_SW_pro <- wild.boot(Chol_SW_pro, design = "fixed", distr = "rademacher", n.ahead = 20, nboot = 2000)  %>% plot(lowerq = 0.16, upperq = 0.84)
ggsave("../Plots/Tri_Chol_SW_pro.pdf", plot = MP_Chol_SW_pro, width = 12, height = 8)

set.seed(1234)
MP_CV_SW_pro <- wild.boot(CV_SW_pro, n.ahead = 20, nboot = 2000) %>% plot
ggsave("../Plots/Tri_CV_SW_pro.pdf", plot = MP_CV_SW_pro, width = 12, height = 8)

set.seed(1234)
MP_CV_SW_pro_rest <- wild.boot(CV_SW_pro_rest, n.ahead = 20, nboot = 2000) %>% plot
ggsave("../Plots/Tri_CV_SW_pro_rest.pdf", plot = MP_CV_SW_pro_rest, width = 12, height = 8)



# Decompose Xi ------------------------------------------------------------
demand_cv   <- solve(CV$B)        %*% t(resid(var4)) %>% t() %>% as.data.frame() %>% select(1)
demand_cv <- cbind(time(var4$y)[-c(1:4)], demand_cv); colnames(demand_cv) <- c("time", "demand_cv")
supply_cv     <- solve(CV$B)          %*% t(resid(var4)) %>% t() %>% as.data.frame() %>% select(2)
supply_cv <- cbind(time(var4$y)[-c(1:4)], supply_cv); colnames(supply_cv) <- c("time", "supply_cv")

non_mp_cv <- left_join(demand_cv, supply_cv, by = "time")
non_mp_cv_stand <- data.frame("time" = non_mp_cv$time, scale(non_mp_cv[,-1]))
non_mp_cv_stand %>% apply(2,sd)

non_mp_cv_stand %>% filter(time >= Start & time <= End2) %>% apply(2,function(x){ifelse(abs(x)>2,1,0)}) %>% colSums()

  
non_mp_cv %>% filter(time >= Start & time <= End2)# %>% colSums()

tempdat <- cbind(non_mp_cv, mp_cv) %>% left_join(IV.SW$eps.ts, by = "time") %>% left_join(IV.RR$eps.ts, by = "time") %>% left_join(IV.SZ$eps.ts, by = "time") 
colnames(tempdat)[c(4,5,6,7)] <- c("mp_cv", "mp_sw", "mp_rr", "mp_sz")

##### SW, bessere Idee!
Start = 1979.75
End1 = 1982.5
End2 = 1983.5


# RR
topmod1 <-  lm(mp_sw ~ mp_rr-1, data = tempdat %>% tidyr::drop_na())
topmod2 <-  lm(mp_rr ~ demand_cv + supply_cv-1, data = tempdat %>% tidyr::drop_na())

summary(topmod1)$coeff %>% round(3)
summary(topmod2)$coeff %>% round(3)

Comp1 <- data.frame("time" = tempdat %>% tidyr::drop_na() %>% select(1), "Comp1" = topmod1 %>% resid)
Comp2 <- tempdat %>% select(time, demand_cv) %>% mutate(demand_cv = demand_cv*topmod1$coefficients*topmod2$coefficients[1])
Comp3 <- tempdat %>% select(time, supply_cv) %>% mutate(supply_cv = supply_cv*topmod1$coefficients*topmod2$coefficients[2])
Comp4 <- data.frame("time" = tempdat %>% tidyr::drop_na() %>% select(1), "Comp4" = topmod2 %>% resid) %>% mutate(Comp4 = Comp4*topmod1$coefficients)

# TRUE value
0.25 *c(quasi.hd(IV.SW, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End1, Freq = 4),
        quasi.hd(IV.SW, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End2, Freq = 4))

decomp1 <- 0.25 *c(decompose_Xi(IV.SW, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End1, Freq = 4, wrong_shocks = Comp1),
        decompose_Xi(IV.SW, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End1, Freq = 4, wrong_shocks = Comp2),
        decompose_Xi(IV.SW, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End1, Freq = 4, wrong_shocks = Comp3),
        decompose_Xi(IV.SW, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End1, Freq = 4, wrong_shocks = Comp4))  %>% t %>% t

decomp2 <- 0.25 *c(decompose_Xi(IV.SW, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End2, Freq = 4, wrong_shocks = Comp1),
        decompose_Xi(IV.SW, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End2, Freq = 4, wrong_shocks = Comp2),
        decompose_Xi(IV.SW, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End2, Freq = 4, wrong_shocks = Comp3),
        decompose_Xi(IV.SW, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End2, Freq = 4, wrong_shocks = Comp4))   %>% t %>% t

# SZ
topmod3 <-  lm(mp_sw ~ mp_sz-1, data = tempdat %>% tidyr::drop_na())
topmod4 <-  lm(mp_sz ~ demand_cv + supply_cv-1, data = tempdat %>% tidyr::drop_na())

summary(topmod3)$coeff %>% round(3)
summary(topmod4)$coeff %>% round(3)

Comp5 <- data.frame("time" = tempdat %>% tidyr::drop_na() %>% select(1), "Comp5" = topmod3 %>% resid)
Comp6 <- tempdat %>% select(time, demand_cv) %>% mutate(demand_cv = demand_cv*topmod3$coefficients*topmod4$coefficients[1])
Comp7 <- tempdat %>% select(time, supply_cv) %>% mutate(supply_cv = supply_cv*topmod3$coefficients*topmod4$coefficients[2])
Comp8 <- data.frame("time" = tempdat %>% tidyr::drop_na() %>% select(1), "Comp8" = topmod4 %>% resid) %>% mutate(Comp8 = Comp8*topmod3$coefficients)

decomp3 <- 0.25 *c(decompose_Xi(IV.SW, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End1, Freq = 4, wrong_shocks = Comp5),
                   decompose_Xi(IV.SW, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End1, Freq = 4, wrong_shocks = Comp6),
                   decompose_Xi(IV.SW, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End1, Freq = 4, wrong_shocks = Comp7),
                   decompose_Xi(IV.SW, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End1, Freq = 4, wrong_shocks = Comp8))  %>% t %>% t

decomp4 <- 0.25 *c(decompose_Xi(IV.SW, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End2, Freq = 4, wrong_shocks = Comp5),
                   decompose_Xi(IV.SW, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End2, Freq = 4, wrong_shocks = Comp6),
                   decompose_Xi(IV.SW, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End2, Freq = 4, wrong_shocks = Comp7),
                   decompose_Xi(IV.SW, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End2, Freq = 4, wrong_shocks = Comp8))   %>% t %>% t


rbind(cbind(decomp1, decomp3, decomp2, decomp4), c(sum(decomp1), sum(decomp3), sum(decomp2), sum(decomp4))) %>% round(3) %>% kableExtra::kable(format = "latex", booktabs = T, digits = 3, linesep = "")


# FAVAR -------------------------------------------------------------------
r <- 3
SDFM <- FALSE
source("Fa.R")
colnames(USA_FAVAR) <- c(colnames(USA_Tri), 1:r)
favar4 <- VAR(USA_FAVAR, p = 4, type = "const")

FaIV.SW      <- get.id.iv(favar4, instruments = SW, Synchro = T)
#FaIV.SW$B %>% round(3)
#FaIV.SW$F_test

FaIV.RR       <- get.id.iv(favar4, instruments = RR, Synchro = T)
#FaIV.RR$B %>% round(3)
#FaIV.RR$F_test

FaIV.SZ       <- get.id.iv(favar4, instruments = SZ, Synchro = T)
#FaIV.SZ$B %>% round(3)
#FaIV.SZ$F_test

mp_iv_fa     <- left_join(FaIV.SW$eps.ts, FaIV.SZ$eps.ts, by = "time")
mp_iv_fa     <- left_join(mp_iv_fa, FaIV.RR$eps.ts, by = "time") %>% tidyr::drop_na()
colnames(mp_iv_fa) <- c("time", "SW", "SZ", "RR")

#mp_iv_fa_stand <- mp_iv_fa
#mp_iv_fa_stand[-1] <- scale(mp_iv_fa, center = F)

mp_iv_fa_tab <-  mp_iv_fa %>% filter(time == 1974.25 | time == 1979.75 | time == 1988.75 | time == 1994 | time == 1990.75 | time == 1998.75 | time == 2001.25| time == 2002.75)
mp_iv_fa_tab[c(4,5),] <- mp_iv_fa_tab[c(5,4),]
mp_iv_fa_tab %>% t() %>% kableExtra::kable(format = "latex", booktabs = T, digits = 3, linesep = "")

#new_cormat2 <- new_cormat %>% left_join(FaIV.SW$eps.ts, by = "time") %>% left_join(FaIV.RR$eps.ts, by = "time") %>% left_join(FaIV.SZ$eps.ts, by = "time") 
#colnames(new_cormat2)[c(10,11,12)] <- c("SW_fa", "RR_fa", "SZ_fa")
#new_cormat2 %>% select(-1) %>% cor %>% round(3) %>% kableExtra::kable(format = "latex", booktabs = T, digits = 3, linesep = "")

# Inference
set.seed(1234)
MP_fa.SW  <- get.MBB.fixed(FaIV.SW, design = "fixed", Step = 15, nboot = 2000, Ci = 0.68) %>% plot.my.irf(Partial = 3, Epsname = c("d", "s", "mp", "X1", "X2", "X3"), Normalize = FaIV.SW$B[3]/0.25) + xlab("") + ylab("")



MP_fa.SW1 <- MP_fa.SW$data %>% filter(!stringr::str_detect(variable, "X")) %>% ggplot( aes(x = h, y = value)) + geom_line() + geom_hline(yintercept = 0, color = 'red') +
  facet_wrap(~variable, scales = "free_y", labeller = label_parsed, ncol = 1) +
  geom_ribbon(aes(ymin = L, ymax = U), alpha=0.3) +
  xlab("Horizon") + ylab("Response") +
  theme_bw()

MP_fa.SW2 <- MP_fa.SW$data %>% filter(!stringr::str_detect(variable, "X")) %>% ggplot( aes(x = h, y = value)) + geom_line() + geom_hline(yintercept = 0, color = 'red') +
  facet_wrap(~variable, scales = "free_y", labeller = label_parsed, ncol = 1) +
  geom_ribbon(aes(ymin = L, ymax = U), alpha=0.3) +
  xlab("Horizon") + ylab("Response") +
  theme_bw()

MP_fa.SW3 <- MP_fa.SW$data %>% filter(!stringr::str_detect(variable, "X")) %>% ggplot( aes(x = h, y = value)) + geom_line() + geom_hline(yintercept = 0, color = 'red') +
  facet_wrap(~variable, scales = "free_y", labeller = label_parsed, ncol = 1) +
  geom_ribbon(aes(ymin = L, ymax = U), alpha=0.3) +
  xlab(" ") + ylab(" ") +
  theme_bw() + 

ggsave("../Plots/MP_fa_SW3.pdf", plot = MP_fa.SW3, width = 12, height = 3)

MP_fa_SW <- grid.arrange(MP_fa.SW1, MP_fa.SW2, MP_fa.SW3, ncol = 3)
ggsave("../Plots/MP_fa_SW123.pdf", plot = MP_fa_SW, width = 12, height = 8)

# Structrual dynamic factor model -----------------------------------------
r <- 8
SDFM <- TRUE
source("Fa.R")

DFM <- VAR(Fact, p = 2, type = "none")

DFM_IV.SW      <- get.id.iv(DFM, instruments = SW, Synchro = T)
DFM_IV.SW$B %>% round(3)
DFM_IV.SW$F_test

DFM_IV.RR       <- get.id.iv(DFM, instruments = RR, Synchro = T)
DFM_IV.RR$B %>% round(3)
DFM_IV.RR$F_test

DFM_IV.SZ       <- get.id.iv(DFM, instruments = SZ, Synchro = T)
DFM_IV.SZ$B %>% round(3)
DFM_IV.SZ$F_test

new_cormat3 <- new_cormat2 %>% left_join(DFM_IV.SW$eps.ts, by = "time") %>% left_join(DFM_IV.RR$eps.ts, by = "time") %>% left_join(DFM_IV.SZ$eps.ts, by = "time") 
colnames(new_cormat3)[c(13, 14, 15)] <- c("SW_DFM", "RR_DFM", "SZ_DFM")
new_cormat3 %>% select(-1) %>% cor %>% round(3) %>% kableExtra::kable(format = "latex", booktabs = T, digits = 3, linesep = "")


# Exogeneity test ---------------------------------------------------------
# SW
USA_Tri_SW <- data.frame("time" = USA_Tri %>% time %>% as.numeric, USA_Tri) %>% left_join(iv2ts(SW, "SW"), by = "time") %>% tidyr::drop_na()
USA_Tri_SW <- ts(USA_Tri_SW %>% select(-1), frequency = 4, start = c(floor(USA_Tri_SW$time[1]), (USA_Tri_SW$time[1] %% 1) / 0.25 + 1 ))
VARselect(USA_Tri_SW, lag.max = 8)
var3_SW <- VAR(USA_Tri_SW, type = "const", p = 3)
CV_SW <- id.cv(var3_SW, SB = c(1984, 1), frequency = 4)
CV_SW
CV
rest_mat <- matrix(c(rep(NA, 12), 0, 0, NA, NA), nrow = 4, ncol = 4, byrow = T)
set.seed(1234)
CV_SW_rest <- id.cv(var3_SW, SB = c(1984, 1), frequency = 4, restriction_matrix = rest_mat)
CV_SW_rest$lRatioTest

# RR
USA_Tri_RR <- data.frame("time" = USA_Tri %>% time %>% as.numeric, USA_Tri) %>% left_join(iv2ts(RR, "RR"), by = "time") %>% tidyr::drop_na()
USA_Tri_RR <- ts(USA_Tri_RR %>% select(-1), frequency = 4, start = c(floor(USA_Tri_RR$time[1]), (USA_Tri_RR$time[1] %% 1) / 0.25 + 1 ))
VARselect(USA_Tri_RR, lag.max = 8)
var3_RR <- VAR(USA_Tri_RR, p = 3, type = "const")
CV_RR <- id.cv(var3_RR, SB = c(1984, 1), frequency = 4)
CV_RR
CV
rest_mat <- matrix(c(rep(NA, 12), 0, 0, NA, NA), nrow = 4, ncol = 4, byrow = T)
set.seed(1234)
CV_RR_rest <- id.cv(var3_RR, SB = c(1984, 1), frequency = 4, restriction_matrix = rest_mat)
CV_RR_rest$lRatioTest


# SZ
USA_Tri_SZ <- data.frame("time" = USA_Tri %>% time %>% as.numeric, USA_Tri) %>% left_join(iv2ts(SZ, "SZ"), by = "time") %>% tidyr::drop_na()
USA_Tri_SZ <- ts(USA_Tri_SZ %>% select(-1), frequency = 4, start = c(floor(USA_Tri_SZ$time[1]), (USA_Tri_SZ$time[1] %% 1) / 0.25 + 1 ))
VARselect(USA_Tri_SZ, lag.max = 8)
var3_SZ <- VAR(USA_Tri_SZ, p = 3, type = "const")
CV_SZ <- id.cv(var3_SZ, SB = c(1984, 1), frequency = 4)
CV_SZ
CV
rest_mat <- matrix(c(rep(NA, 12), 0, 0, NA, NA), nrow = 4, ncol = 4, byrow = T)
set.seed(1234)
CV_SZ_rest <- id.cv(var3_SZ, SB = c(1984, 1), frequency = 4, restriction_matrix = rest_mat)
CV_SZ_rest$lRatioTest


# Counterfactural ---------------------------------------------------------
cf_pi <- data.frame("SR" = SR %>% cf.mod(series = 2, Partial = 3, Epsname = c("d", "s", "mp")) %>% select(-y_hat),
                    "dCov" = cf(dCov, series = 2)$counter$`pi with and without cumulative effect of flow GBR1 shock`, 
                    "CV" =   cf(CV, series = 2)$counter$`pi with and without cumulative effect of flow GBR1 shock`)
colnames(cf_pi)[1:3] <- c("time", "Actual", "SR")

cf_pi <- cf_pi %>% left_join(cf.mod(IV.SW, series = 2, Epsname = c("d", "s", "mp")) %>% select(-c(y_true, y_hat)), by = "time") %>% 
  left_join(cf.mod(IV.SZ, series = 2, Epsname = c("d", "s", "mp")) %>% select(-c(y_true, y_hat)) , by = "time") %>% 
  left_join(cf.mod(IV.RR, series = 2, Epsname = c("d", "s", "mp")) %>% select(-c(y_true, y_hat)), by = "time") %>% tidyr::drop_na()

colnames(cf_pi)[6:8] <- c("SW", "SZ", "RR")

#cf_pi <- cf_pi %>% filter(time >= 1978.5 & time <= 1983.25)
cf_pi_plot <- melt(cf_pi[,-2], id = "time")
cf_pi_plot <- cf_pi_plot %>% add_column("actual" = rep( cf_pi[,2], 6))

cf_pi_plot$variable <- factor(cf_pi_plot$variable, levels = c("CV", "dCov", "SW", "SZ", "RR", "SR"))

library(RColorBrewer)

cf_pi_plot %>% ggplot(aes(x = time, y = value, group = variable)) +
  geom_line(linetype = "longdash", alpha = 0.6) + 
  geom_line(aes(x = time, y = actual, linetype = "solid")) +
  scale_x_continuous(breaks = seq(cf_pi$time[1] %>% ceiling(),cf_pi$time[nrow(cf_pi)], 4)) +
  facet_wrap(~variable, scales = "free_x", ncol = 3) +
  xlab("") + ylab(" ") + theme_bw() + theme(legend.position = "none")

ggsave("../Plots/CF_pi_3sub.pdf", plot = last_plot(), width = 12, height = 8)




# Hist Decomp -------------------------------------------------------------
# pi
hd_pi <- data.frame("SR" = SR %>% cf.mod(series = 2, Partial = 3, Epsname = c("d", "s", "mp")) %>% select(-y_cf),
                    "dCov" = hd(dCov, series = 2)$hidec[,"Cumulative effect of flow  GBR1 shock on  pi"], 
                    "CV" =   hd(CV, series = 2)$hidec[,"Cumulative effect of flow  GBR1 shock on  pi"])
hd_pi[,c(2,3)] <- hd_pi[,c(3,2)]
colnames(hd_pi)[1:3] <- c("time", "Actual", "SR")

hd_pi <- hd_pi %>% left_join(cf.mod(IV.SW, series = 2, Epsname = c("d", "s", "mp")) %>% select(-c(y_true, y_cf)), by = "time") %>% 
  left_join(cf.mod(IV.SZ, series = 2, Epsname = c("d", "s", "mp")) %>% select(-c(y_true, y_cf)) , by = "time") %>% 
  left_join(cf.mod(IV.RR, series = 2, Epsname = c("d", "s", "mp")) %>% select(-c(y_true, y_cf)), by = "time") %>% tidyr::drop_na()

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

ggsave("../Plots/HD_pi_3sub_volcker.pdf", plot = last_plot(), width = 12, height = 8)

# GBR1
hd_GBR1 <- data.frame("SR" = SR %>% cf.mod(series = 3, Partial = 3, Epsname = c("d", "s", "mp")) %>% select(-y_cf),
                      "dCov" = hd(dCov, series = 3)$hidec[,"Cumulative effect of flow  GBR1 shock on  GBR1"], 
                      "CV" =   hd(CV, series = 3)$hidec[,"Cumulative effect of flow  GBR1 shock on  GBR1"])
hd_GBR1[,c(2,3)] <- hd_GBR1[,c(3,2)]
colnames(hd_GBR1)[1:3] <- c("time", "Actual", "SR")

hd_GBR1 <- hd_GBR1 %>% left_join(cf.mod(IV.SW, series = 3, Epsname = c("d", "s", "mp")) %>% select(-c(y_true, y_cf)), by = "time") %>% 
  left_join(cf.mod(IV.SZ, series = 3, Epsname = c("d", "s", "mp")) %>% select(-c(y_true, y_cf)) , by = "time") %>% 
  left_join(cf.mod(IV.RR, series = 3, Epsname = c("d", "s", "mp")) %>% select(-c(y_true, y_cf)), by = "time") %>% tidyr::drop_na()

colnames(hd_GBR1)[6:8] <- c("SW", "SZ", "RR")

hd_GBR1_plot <- melt(hd_GBR1[,-2], id = "time")
hd_GBR1_plot <- hd_GBR1_plot %>% add_column("actual" = rep( hd_GBR1[,2], 6))

hd_GBR1_plot$variable <- factor(hd_GBR1_plot$variable, levels = c("CV", "dCov", "SW", "SZ", "RR", "SR"))

hd_GBR1_plot %>% ggplot(aes(x = time, y = value, group = variable)) +
  geom_line(linetype = "longdash", alpha = 0.6) + 
  
  geom_line(aes(x = time, y = actual, linetype = "solid")) +
  scale_x_continuous(breaks = seq(hd_GBR1$time[1] %>% ceiling(),hd_GBR1$time[nrow(hd_GBR1)], 4)) +
  facet_wrap(~variable, scales = "free_x", ncol = 3) +
  xlab("") + ylab(" ") + theme_bw() + theme(legend.position = "none")

ggsave("../Plots/HD_GBR1_3sub.pdf", plot = last_plot(), width = 12, height = 8)


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

ggsave("../Plots/HD_pi_3bar_volcker.pdf", plot = last_plot(), width = 12, height = 6)



# Other Instruments -------------------------------------------------------
# RR83b
RR83b.raw      <- read.csv("Data/Instruments/RR83b.csv")
RR83b          <- ts(RR83b.raw[,-1], frequency = 4, start = c(floor(RR83b.raw$Year[1]), (RR83b.raw$Year[1] %% 1) / 0.25 + 1 ))
RR83b          <- Remove_AC(RR83b, lag = 1) 
IV.RR83b       <- get.id.iv(var4, instruments = RR83b, Synchro = T)
IV.RR83b$B
IV.RR83b$F_test
Box.test(RR83b, lag = 1, type = "Ljung-Box")
Box.test(RR83b, lag = 2, type = "Ljung-Box")
Box.test(RR83b, lag = 3, type = "Ljung-Box")
Box.test(RR83b, lag = 4, type = "Ljung-Box")



# FF1vr
FF1vr.raw      <- read.csv("Data/Instruments/FF1vr.csv")
FF1vr          <- ts(FF1vr.raw[,-1], frequency = 4, start = c(floor(FF1vr.raw$Year[1]), (FF1vr.raw$Year[1] %% 1) / 0.25 + 1 ))
FF1vr          <- Remove_AC(FF1vr, lag = 1) 
IV.FF1vr       <- get.id.iv(var4, instruments = FF1vr, Synchro = T)
IV.FF1vr$B
IV.FF1vr$F_test
Box.test(FF1vr, lag = 1, type = "Ljung-Box")
Box.test(FF1vr, lag = 2, type = "Ljung-Box")
Box.test(FF1vr, lag = 3, type = "Ljung-Box")

# FF4vr
FF4vr.raw      <- read.csv("Data/Instruments/FF4vr.csv")
FF4vr          <- ts(FF4vr.raw[,-1], frequency = 4, start = c(floor(FF4vr.raw$Year[1]), (FF4vr.raw$Year[1] %% 1) / 0.25 + 1 ))
FF4vr          <- Remove_AC(FF4vr, lag = 1) 
IV.FF4vr       <- get.id.iv(var4, instruments = FF4vr, Synchro = T)
IV.FF4vr$B
IV.FF4vr$F_test
Box.test(FF4vr, lag = 1, type = "Ljung-Box")
Box.test(FF4vr, lag = 2, type = "Ljung-Box")
Box.test(FF4vr, lag = 3, type = "Ljung-Box")


# FF1gb
FF1gb.raw      <- read.csv("Data/Instruments/FF1gb.csv")
FF1gb          <- ts(FF1gb.raw[,-1], frequency = 4, start = c(floor(FF1gb.raw$Year[1]), (FF1gb.raw$Year[1] %% 1) / 0.25 + 1 ))
FF1gb          <- Remove_AC(FF1gb, lag = 1) 
IV.FF1gb       <- get.id.iv(var4, instruments = FF1gb, Synchro = T)
IV.FF1gb$B
IV.FF1gb$F_test
Box.test(FF1gb, lag = 1, type = "Ljung-Box")
Box.test(FF1gb, lag = 2, type = "Ljung-Box")
Box.test(FF1gb, lag = 3, type = "Ljung-Box")

# FF4gb
FF4gb.raw      <- read.csv("Data/Instruments/FF4gb.csv")
FF4gb          <- ts(FF4gb.raw[,-1], frequency = 4, start = c(floor(FF4gb.raw$Year[1]), (FF4gb.raw$Year[1] %% 1) / 0.25 + 1 ))
FF4gb          <- Remove_AC(FF4gb, lag = 1) 
IV.FF4gb       <- get.id.iv(var4, instruments = FF4gb, Synchro = T)
IV.FF4gb$B
IV.FF4gb$F_test
Box.test(FF4gb, lag = 1, type = "Ljung-Box")


# ED2vr
ED2vr.raw      <- read.csv("Data/Instruments/ED2vr.csv")
ED2vr          <- ts(ED2vr.raw[,-1], frequency = 4, start = c(floor(ED2vr.raw$Year[1]), (ED2vr.raw$Year[1] %% 1) / 0.25 + 1 ))
IV.ED2vr       <- get.id.iv(var4, instruments = ED2vr, Synchro = T)
IV.ED2vr$B
IV.ED2vr$F_test


# ED2gb
ED2gb.raw      <- read.csv("Data/Instruments/ED2gb.csv")
ED2gb          <- ts(ED2gb.raw[,-1], frequency = 4, start = c(floor(ED2gb.raw$Year[1]), (ED2gb.raw$Year[1] %% 1) / 0.25 + 1 ))
IV.ED2gb       <- get.id.iv(var4, instruments = ED2gb, Synchro = T)
IV.ED2gb$B
IV.ED2gb$F_test


# MA
MA.raw      <- read.csv("Data/Instruments/MA.csv")
MA          <- ts(MA.raw[,-1], frequency = 4, start = c(floor(MA.raw$Year[1]), (MA.raw$Year[1] %% 1) / 0.25 + 1 ))
IV.MA       <- get.id.iv(var4, instruments = MA, Synchro = T)
IV.MA$B
IV.MA$F_test
Box.test(MA, lag = 1, type = "Ljung-Box")
Box.test(MA, lag = 2, type = "Ljung-Box")
Box.test(MA, lag = 3, type = "Ljung-Box")



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
colnames(CV.B0.boot)[c(7, 10)] <- c("p", "r")
CV.B0.boot <- CV.B0.boot %>% mutate(labling = ifelse(p*r<0, 1, 0))
CV.B0.boot$labling %>% mean


dCov.B0.boot <- dCov.boot %>% extract_B0_boot
dCov.B0.boot %>% colnames() 
colnames(dCov.B0.boot)[c(7, 10)] <- c("p", "r")
dCov.B0.boot <- dCov.B0.boot %>% mutate(labling = ifelse(p*r<0, 1, 0))
dCov.B0.boot$labling %>% mean


# not good anymore --------------------------------------------------------


#tempmod1 <- lm(mp_rr ~ mp_sw + demand_cv + supply_cv -1, data = tempdat %>% tidyr::drop_na())
#tempmod2 <- lm(mp_sz ~ mp_sw + demand_cv + supply_cv -1, data = tempdat %>% tidyr::drop_na())#

#tempmod3 <- lm(mp_rr - mp_sw ~ + demand_cv + supply_cv -1, data = tempdat %>% tidyr::drop_na())
#tempmod4 <- lm(mp_sz - mp_sw ~ + demand_cv + supply_cv -1, data = tempdat %>% tidyr::drop_na())

#tempmod3 %>% summary
#tempmod4 %>% summary

#tempdat2 <- cbind(non_mp_cv, mp_cv) %>% left_join(mp_iv, by ="time") %>% tidyr::drop_na() %>% select(-1)
#colnames(tempdat2)[3] <- "mp_cv"

#tempmod1 <- lm(RR ~ SW + demand_cv + supply_cv -1, data = tempdat2)
#tempmod2 <- lm(SZ ~ SW + demand_cv + supply_cv -1, data = tempdat2)

#summary(tempmod1)$coeff %>% round(3)
#summary(tempmod2)$coeff %>% round(3)

#tempmod2 %>% summary
#1- crossprod(resid(tempmod1)) /var(tempdat2$RR)
#1-sigma(tempmod2)^2/var(tempdat2$SZ)

#lm(mp_rr ~ mp_cv + demand_cv + supply_cv -1, data = tempdat) %>% summary
#lm(mp_sz ~ mp_cv + demand_cv + supply_cv -1, data = tempdat) %>% summary

## RR
0.25 * c(quasi.hd(IV.RR, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End1, Freq = 4),
         quasi.hd(IV.RR, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End2, Freq = 4))

wrongshock_SW1 <- IV.SW$eps.ts %>% mutate(epsilon = epsilon*-1)
wrongshock_CVd1 <- non_mp_cv %>% select(time, demand_cv) %>% mutate(demand_cv = demand_cv*tempmod3$coefficients[1])
wrongshock_CVs1 <- non_mp_cv %>% select(time, supply_cv) %>% mutate(supply_cv = supply_cv*tempmod3$coefficients[2])

0.25 *c(decompose_Xi(IV.RR, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End1, Freq = 4, wrong_shocks = wrongshock_SW1),
        decompose_Xi(IV.RR, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End1, Freq = 4, wrong_shocks = wrongshock_CVd1),
        decompose_Xi(IV.RR, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End1, Freq = 4, wrong_shocks = wrongshock_CVs1)) #%>% sum

0.25 *c(decompose_Xi(IV.RR, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End2, Freq = 4, wrong_shocks = wrongshock_SW1),
        decompose_Xi(IV.RR, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End2, Freq = 4, wrong_shocks = wrongshock_CVd1),
        decompose_Xi(IV.RR, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End2, Freq = 4, wrong_shocks = wrongshock_CVs1))# %>% sum

## SZ
0.25 *c(quasi.hd(IV.SZ, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End1, Freq = 4),
        quasi.hd(IV.SZ, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End2, Freq = 4))

wrongshock_SW2 <- IV.SW$eps.ts  %>% mutate(epsilon = epsilon*-1)
wrongshock_CVd2 <- non_mp_cv %>% select(time, demand_cv) %>% mutate(demand_cv = demand_cv*tempmod4$coefficients[1])
wrongshock_CVs2 <- non_mp_cv %>% select(time, supply_cv) %>% mutate(supply_cv = supply_cv*tempmod4$coefficients[2])

0.25 *c(decompose_Xi(IV.SZ, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End1, Freq = 4, wrong_shocks = wrongshock_SW2),
        decompose_Xi(IV.SZ, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End1, Freq = 4, wrong_shocks = wrongshock_CVd2),
        decompose_Xi(IV.SZ, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End1, Freq = 4, wrong_shocks = wrongshock_CVs2))# %>% sum

0.25 *c(decompose_Xi(IV.SZ, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End2, Freq = 4, wrong_shocks = wrongshock_SW2),
        decompose_Xi(IV.SZ, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End2, Freq = 4, wrong_shocks = wrongshock_CVd2),
        decompose_Xi(IV.SZ, series = 2, Partial = 3, Epsname = c("d", "s", "mp", "f"), Start = Start, End = End2, Freq = 4, wrong_shocks = wrongshock_CVs2))# %>% sum



