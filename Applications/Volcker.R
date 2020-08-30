Start = 1979.75
End1 = 1982.5
End2 = 1983.5

(T0A = Start - (End1 - Start) - 0.25)
(T0B = Start - (End2 - Start) - 0.25)
(T0C = T0A - (End1 - Start) - 0.25)
(TTA = End1 + (End1 - Start) + 0.25)
(TTB = End2 + (End2 - Start) + 0.25)


Volcker.fct <- function(x, series = 2, Par = NULL, Epsname = c("d", "s", "mp"), Start, End){
        if (!is.null(Par)){
       erg <- quasi.hd(x, series, Partial = Par, Epsname = Epsname, Start = Start, End = End, Freq = 4) 
        } else {
       erg <- c(quasi.hd(x, series, Partial = 1, Epsname = Epsname, Start = Start, End = End, Freq = 4),
                quasi.hd(x, series, Partial = 2, Epsname = Epsname, Start = Start, End = End, Freq = 4),
                quasi.hd(x, series, Partial = 3, Epsname = Epsname, Start = Start, End = End, Freq = 4)) %>% t %>% t                
        }
        erg
}



# Pi ----------------------------------------------------------------------
# CV 
pi_RefA_CV <- Volcker.fct(CV, series = 2, Start = Start, End = End1) * 0.25
pi_RefB_CV <- Volcker.fct(CV, series = 2, Start = Start, End = End2) * 0.25
pi_anteA_CV <- Volcker.fct(CV, series = 2, Start = T0A, End = Start- 0.25) * 0.25
pi_anteB_CV <- Volcker.fct(CV, series = 2, Start = T0B, End = Start- 0.25) * 0.25
pi_anteC_CV <- Volcker.fct(CV, series = 2, Start = T0C, End = T0A- 0.25) * 0.25
pi_postA_CV <- Volcker.fct(CV, series = 2, Start = End1+ 0.25, End = TTA) * 0.25
pi_postB_CV <- Volcker.fct(CV, series = 2, Start = End2+ 0.25, End = TTB) * 0.25

# SW
pi_RefA_SW <- Volcker.fct(IV.SW, Par = 3, Start = Start, End = End1) * 0.25
pi_RefB_SW <- Volcker.fct(IV.SW, Par = 3, Start = Start, End = End2) * 0.25
pi_anteA_SW <- Volcker.fct(IV.SW, Par = 3, Start = T0A, End = Start- 0.25) * 0.25
pi_anteB_SW <- Volcker.fct(IV.SW, Par = 3, Start = T0B, End = Start- 0.25) * 0.25
pi_anteC_SW <- Volcker.fct(IV.SW, Par = 3, Start = T0C, End = T0A- 0.25) * 0.25
pi_postA_SW <- Volcker.fct(IV.SW, Par = 3, Start = End1+ 0.25, End = TTA) * 0.25
pi_postB_SW <- Volcker.fct(IV.SW, Par = 3, Start = End2+ 0.25, End = TTB) * 0.25


# output ------------------------------------------------------------------
# CV
x_RefA_CV <- Volcker.fct(CV, series = 1, Start = Start, End = End1)  
x_RefB_CV <- Volcker.fct(CV, series = 1, Start = Start, End = End2)  
x_anteA_CV <- Volcker.fct(CV, series = 1, Start = T0A, End = Start- 0.25)  
x_anteB_CV <- Volcker.fct(CV, series = 1, Start = T0B, End = Start- 0.25)  
x_anteC_CV <- Volcker.fct(CV, series = 1, Start = T0C, End = T0A- 0.25)  
x_postA_CV <- Volcker.fct(CV, series = 1, Start = End1+ 0.25, End = TTA)  
x_postB_CV <- Volcker.fct(CV, series = 1, Start = End2+ 0.25, End = TTB)  

# SW
x_RefA_SW <- Volcker.fct(IV.SW, series = 1, Par = 3, Start = Start, End = End1)  
x_RefB_SW <- Volcker.fct(IV.SW, series = 1, Par = 3, Start = Start, End = End2)  
x_anteA_SW <- Volcker.fct(IV.SW, series = 1, Par = 3, Start = T0A, End = Start- 0.25)  
x_anteB_SW <- Volcker.fct(IV.SW, series = 1, Par = 3, Start = T0B, End = Start- 0.25)  
x_anteC_SW <- Volcker.fct(IV.SW, series = 1, Par = 3, Start = T0C, End = T0A- 0.25)  
x_postA_SW <- Volcker.fct(IV.SW, series = 1, Par = 3, Start = End1+ 0.25, End = TTA)  
x_postB_SW <- Volcker.fct(IV.SW, series = 1, Par = 3, Start = End2+ 0.25, End = TTB)  




# REST --------------------------------------------------------------------

temphd <- hd(CV, series = 2)
temp <- temphd$hidec[,"Cumulative effect of flow  GBR1 shock on  pi"] %>% ts2df(Varname = "VP")
temp <- temp %>% dplyr::filter(time >= Start & time <= End1) %>% select(-1)
temp %>% sum()*0.25

tempcf <- cf(CV, series = 2)

# pi
USA_Tri_raw %>% filter(Year >= Start & Year <= End1) %>% select(pi) 
temp <- irf(CV, n.ahead = (End2-Start)*4+1)
temp$irf[,(3*(2-1)+1+1)]


Tn = 1966
T1 = 1979.5
T2 = 1968
T3 = 1983

Xi_before <- 0.25 *c(quasi.hd(CV, series = 2, Partial = 1, Epsname = c("d", "s", "mp"), Start = Tn, End = T1-0.25, Freq = 4),
                     quasi.hd(CV, series = 2, Partial = 2, Epsname = c("d", "s", "mp"), Start = Tn, End = T1-0.25, Freq = 4),
                     quasi.hd(CV, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Tn, End = T1-0.25, Freq = 4)) %>% t %>% t 


Xi_all <- 0.25 *c(quasi.hd(CV, series = 2, Partial = 1, Epsname = c("d", "s", "mp"), Start = Tn, End = T2, Freq = 4),
                  quasi.hd(CV, series = 2, Partial = 2, Epsname = c("d", "s", "mp"), Start = Tn, End = T2, Freq = 4),
                  quasi.hd(CV, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Tn, End = T2, Freq = 4)) %>% t %>% t

Xi_VP <- 0.25 *c(quasi.hd(CV, series = 2, Partial = 1, Epsname = c("d", "s", "mp"), Start = T1, End = T2, Freq = 4),
        quasi.hd(CV, series = 2, Partial = 2, Epsname = c("d", "s", "mp"), Start = T1, End = T2, Freq = 4),
        quasi.hd(CV, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = T1, End = T2, Freq = 4)) %>% t %>% t


0.25 *c(quasi.hd(CV, series = 2, Partial = 1, Epsname = c("d", "s", "mp"), Start = Start, End = End1, Freq = 4),
        quasi.hd(CV, series = 2, Partial = 2, Epsname = c("d", "s", "mp"), Start = Start, End = End1, Freq = 4),
        quasi.hd(CV, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End1, Freq = 4)) %>% t %>% t

0.25 *c(quasi.hd(CV, series = 2, Partial = 1, Epsname = c("d", "s", "mp"), Start = Start, End = End2, Freq = 4),
        quasi.hd(CV, series = 2, Partial = 2, Epsname = c("d", "s", "mp"), Start = Start, End = End2, Freq = 4),
        quasi.hd(CV, series = 2, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End2, Freq = 4)) %>% t %>% t 

# gbr1
c(quasi.hd(CV, series = 3, Partial = 1, Epsname = c("d", "s", "mp"), Start = Start, End = End1, Freq = 4),
        quasi.hd(CV, series = 3, Partial = 2, Epsname = c("d", "s", "mp"), Start = Start, End = End1, Freq = 4),
        quasi.hd(CV, series = 3, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End1, Freq = 4)) %>% t %>% t

c(quasi.hd(CV, series = 3, Partial = 1, Epsname = c("d", "s", "mp"), Start = Start, End = End2, Freq = 4),
        quasi.hd(CV, series = 3, Partial = 2, Epsname = c("d", "s", "mp"), Start = Start, End = End2, Freq = 4),
        quasi.hd(CV, series = 3, Partial = 3, Epsname = c("d", "s", "mp"), Start = Start, End = End2, Freq = 4)) %>% t %>% t



# New SVAR without VP -----------------------------------------------------
## some functions

# proxy svar
get.id.iv.base <- function(resid, resid_time, proxy){
        iv2ts <- function(x){
                erg <- data.frame(x %>% time %>% as.numeric(),  x)
                colnames(erg) <- c("time", "iv")
                return(erg)}
        proxy <- iv2ts(proxy)
        
        syn <- data.frame("time" = resid_time, "u" = resid) %>% left_join(proxy, by = "time") %>% tidyr::drop_na()
        K <- ncol(resid)
        u <- as.matrix(syn[,c(2:(K+1))])
        iv <- as.matrix(syn$iv)
        Tob <- nrow(u)
        
        Pi        <- solve(crossprod(u)/Tob) %*% (crossprod(u, iv)/Tob)
        phi2      <- (crossprod(iv, u)/Tob) %*% solve(crossprod(u)/Tob) %*% (crossprod(u, iv)/Tob)
        B_k       <- (crossprod(u, iv)/Tob) %*% (1/sqrt(phi2))
        s_errors  <- u %*% Pi %*% sqrt(phi2)^(-1)
        F_test    <- test.weak.IV(Tob = Tob, u = u, instruments = iv)
        
        erg <- list("B" = B_k,
                    "F_test" = F_test,
                    "epsilon" = s_errors,
                    "eps.ts" = data.frame("time" = syn$time, "epsilon" = s_errors))
        
        return(erg)
}



## reduced form parameters
USA_Tri_exVP <- USA_Tri_raw %>% filter(Year < Start | Year > End2)
var4_exVP <- USA_Tri_exVP[,-1] %>% VAR( p = 4, type = "const")
Beta_exVP <- t(Bcoef(var4_exVP))

# new residuals
u_exVP <- var4$datamat[,c(1:3)] - as.matrix(var4$datamat[,-c(1:3)]) %*% Beta_exVP

IV.SW_exVP <- get.id.iv.base(u_exVP, time(var4$y)[-c(1:4)], SW)

B_SW_exVP_stand <- IV.SW_exVP$B / IV.SW_exVP$B[3,1] * 0.25
Theta.temp <- get.my.irf(Simple = T, Bcoef = Beta_exVP, Bmat = B_SW_exVP_stand, Step = 15, Lag = 4, Integrate = F)

par(mfrow = c(3,1))
for (i in 1:3) {
plot(Theta.temp[i,3,], xlab = "Horizont", ylab = "Response", type = "l")
abline(h = 0, lty = 2, col = 2)
}
par(mfrow = c(1,1))

## new and old SW shocks
mp_SW_excum_VP <- cbind(IV.SW$eps.ts, IV.SW_exVP$epsilon)
colnames(mp_SW_excum_VP) <- c("time", "cum_VP", "ex_VP")

mp_SW_excum_VP %>% reshape2::melt(id = "time") %>% ggplot(aes(x = time, y = value, group = variable)) +
        geom_line(aes(linetype = variable), size = 0.7) + 
        scale_x_continuous(breaks = seq(mp_SW_excum_VP$time[1] %>% ceiling(),mp_SW_excum_VP$time[nrow(mp_SW_excum_VP)], 2)) + 
        xlab("") + ylab(" ")  + theme_bw() + theme(legend.title = element_blank()) + 
        geom_hline(yintercept = c(-1, 1), linetype = "dashed", alpha = 0.3) + 
        annotate("rect", xmin = 1979.5, xmax = 1983, ymin = -Inf, ymax = Inf, alpha = 0.3)
ggsave("../Plots/mp_SW_excum_VP.pdf", plot = last_plot(), width = 12, height = 6)

