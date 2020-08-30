
# time series --> data frame ----------------------------------------------
iv2ts <- function(x, varname){
  erg <- data.frame(x %>% time %>% as.numeric(),  x)
  colnames(erg) <- c("time", varname)
  return(erg)}
SW.ts <- iv2ts(SW, "SW")
RR.ts <- iv2ts(RR, "RR")
SZ.ts <- iv2ts(SZ, "SZ")

# Restriction for one KxK matrix
(C1 <- matrix(c(rep(c(0,0,0,1),3), rep(0,4)), nrow = 4, ncol = 4, byrow = T))
(C2 <- matrix(c(rep(c(0,0,0,0),3), 1,1,1,0), nrow = 4, ncol = 4, byrow = T))
# Granger Causality -------------------------------------------------------
# SW
USA_Tri_SW <- data.frame("time" = USA_Tri %>% time %>% as.numeric, USA_Tri) %>% left_join(SW.ts, by = "time") %>% tidyr::drop_na()
USA_Tri_SW <- ts(USA_Tri_SW %>% select(-1), frequency = 4, start = c(floor(USA_Tri_SW$time[1]), (USA_Tri_SW$time[1] %% 1) / 0.25 + 1 ))
VARselect(USA_Tri_SW, lag.max = 8)
var3_SW <- VAR(USA_Tri_SW, type = "const", p = 3)

# own test
test1 <- test.GC(var3_SW, C1)$p_val$p_F
test2 <- test.GC(var3_SW, C2)$p_val$p_F

# vars pkg
causality(var3_SW, cause = "SW")$Granger
causality(var3_SW, cause = c("x", "pi", "GBR1"))$Granger


# RR
USA_Tri_RR <- data.frame("time" = USA_Tri %>% time %>% as.numeric, USA_Tri) %>% left_join(RR.ts, by = "time") %>% tidyr::drop_na()
USA_Tri_RR <- ts(USA_Tri_RR %>% select(-1), frequency = 4, start = c(floor(USA_Tri_RR$time[1]), (USA_Tri_RR$time[1] %% 1) / 0.25 + 1 ))
VARselect(USA_Tri_RR, lag.max = 8)
var3_RR <- VAR(USA_Tri_RR, type = "const", p = 3)

test1 <- c(test1, test.GC(var3_RR, C1)$p_val$p_F)
test2 <- c(test2, test.GC(var3_RR, C2)$p_val$p_F)

causality(var3_RR, cause = "RR")$Granger
causality(var3_RR, cause = c("x", "pi", "GBR1"))$Granger

# SZ
USA_Tri_SZ <- data.frame("time" = USA_Tri %>% time %>% as.numeric, USA_Tri) %>% left_join(SZ.ts, by = "time") %>% tidyr::drop_na()
USA_Tri_SZ <- ts(USA_Tri_SZ %>% select(-1), frequency = 4, start = c(floor(USA_Tri_SZ$time[1]), (USA_Tri_SZ$time[1] %% 1) / 0.25 + 1 ))
VARselect(USA_Tri_SZ, lag.max = 8)
var3_SZ <- VAR(USA_Tri_SZ, p = 3, type = "const")

test1 <- c(test1, test.GC(var3_SZ, C1)$p_val$p_F)
test2 <- c(test2, test.GC(var3_SZ, C2)$p_val$p_F)

causality(var3_SZ, cause = "SZ")$Granger
causality(var3_SZ, cause = c("x", "pi", "GBR1"))$Granger



# FAVAR -------------------------------------------------------------------
r <- 2
# Restriction for one KxK matrix
(C3 <- matrix(c(rep(c(0,0,0,rep(0,r),1),3+r), rep(0,4+r)), nrow = 4+r, ncol = 4+r, byrow = T))
(C4 <- matrix(c(rep(c(0,0,0,0,rep(0,r)),3+r), rep(1,3+r),0), nrow = 4+r, ncol = 4+r, byrow = T))


SDFM <- FALSE
source("Fa.R")
colnames(USA_FAVAR) <- c(colnames(USA_Tri), 1:r)
# SW
USA_Tri_SW_Fa <- data.frame("time" = USA_FAVAR %>% time %>% as.numeric, USA_FAVAR) %>% left_join(SW.ts, by = "time") %>% tidyr::drop_na()
USA_Tri_SW_Fa <- ts(USA_Tri_SW_Fa %>% select(-1), frequency = 4, start = c(floor(USA_Tri_SW_Fa$time[1]), (USA_Tri_SW_Fa$time[1] %% 1) / 0.25 + 1 ))
VARselect(USA_Tri_SW_Fa, lag.max = 8)
var3_SW_Fa <- VAR(USA_Tri_SW_Fa, type = "const", p = 3)

test1 <- c(test1, test.GC(var3_SW_Fa, C3)$p_val$p_F)
test2 <- c(test2, test.GC(var3_SW_Fa, C4)$p_val$p_F)


# RR
USA_Tri_RR_Fa <- data.frame("time" = USA_FAVAR %>% time %>% as.numeric, USA_FAVAR) %>% left_join(RR.ts, by = "time") %>% tidyr::drop_na()
USA_Tri_RR_Fa <- ts(USA_Tri_RR_Fa %>% select(-1), frequency = 4, start = c(floor(USA_Tri_RR_Fa$time[1]), (USA_Tri_RR_Fa$time[1] %% 1) / 0.25 + 1 ))
VARselect(USA_Tri_RR_Fa, lag.max = 8)
var3_RR_Fa <- VAR(USA_Tri_RR_Fa, type = "const", p = 3)

test1 <- c(test1, test.GC(var3_RR_Fa, C3)$p_val$p_F)
test2 <- c(test2, test.GC(var3_RR_Fa, C4)$p_val$p_F)


# SZ
USA_Tri_SZ_Fa <- data.frame("time" = USA_FAVAR %>% time %>% as.numeric, USA_FAVAR) %>% left_join(SZ.ts, by = "time") %>% tidyr::drop_na()
USA_Tri_SZ_Fa <- ts(USA_Tri_SZ_Fa %>% select(-1), frequency = 4, start = c(floor(USA_Tri_SZ_Fa$time[1]), (USA_Tri_SZ_Fa$time[1] %% 1) / 0.25 + 1 ))
VARselect(USA_Tri_SZ_Fa, lag.max = 8)
var3_SZ_Fa <- VAR(USA_Tri_SZ_Fa, type = "const", p = 3)

test1 <- c(test1, test.GC(var3_SZ_Fa, C3)$p_val$p_F)
test2 <- c(test2, test.GC(var3_SZ_Fa, C4)$p_val$p_F)



GC_tab <- rbind(test1, test2)
colnames(GC_tab) <- c("SW", "RR", "SZ", "SW_fa", "RR_fa", "SZ_fa")
GC_tab %>% kableExtra::kable(format = "latex", booktabs = T, digits = 3, linesep = "")

# replace p_val$p_F by TestStat$lambda_F



# Pros --------------------------------------------------------------------
# RR
USA_Tri_RR_pro <- data.frame("time" = USA_Tri %>% time %>% as.numeric, USA_Tri) %>% left_join(iv2ts(RR_pro , "RR_pro"), by = "time") %>% tidyr::drop_na()
USA_Tri_RR_pro <- ts(USA_Tri_RR_pro %>% select(-1), frequency = 4, start = c(floor(USA_Tri_RR_pro$time[1]), (USA_Tri_RR_pro$time[1] %% 1) / 0.25 + 1 ))
VARselect(USA_Tri_RR_pro, lag.max = 8)
var3_RR_pro <- VAR(USA_Tri_RR_pro, type = "const", p = 3)

test.GC(var3_RR_pro, C1)
test.GC(var3_RR_pro, C2)

#test1 <- c(test1, test.GC(var3_RR_pro, C1)$p_val$p_F)
#test2 <- c(test2, test.GC(var3_RR_pro, C2)$p_val$p_F)



# SZ
USA_Tri_SZ_pro <- data.frame("time" = USA_Tri %>% time %>% as.numeric, USA_Tri) %>% left_join(iv2ts(SZ_pro , "SZ_pro"), by = "time") %>% tidyr::drop_na()
USA_Tri_SZ_pro <- ts(USA_Tri_SZ_pro %>% select(-1), frequency = 4, start = c(floor(USA_Tri_SZ_pro$time[1]), (USA_Tri_SZ_pro$time[1] %% 1) / 0.25 + 1 ))
VARselect(USA_Tri_SZ_pro, lag.max = 8)
var3_SZ_pro <- VAR(USA_Tri_SZ_pro, p = 3, type = "const")

test.GC(var3_SZ_pro, C1)
test.GC(var3_SZ_pro, C2)


# FAVAR 
r <- 3
# Restriction for one KxK matrix
(C3 <- matrix(c(rep(c(0,0,0,rep(0,r),1),3+r), rep(0,4+r)), nrow = 4+r, ncol = 4+r, byrow = T))
(C4 <- matrix(c(rep(c(0,0,0,0,rep(0,r)),3+r), rep(1,3+r),0), nrow = 4+r, ncol = 4+r, byrow = T))

SDFM <- FALSE
source("Fa.R")
colnames(USA_FAVAR) <- c(colnames(USA_Tri), 1:r)

# RR
USA_Tri_RR_Fa_pro <- data.frame("time" = USA_FAVAR %>% time %>% as.numeric, USA_FAVAR) %>% left_join(iv2ts(RR_pro , "RR_pro"), by = "time") %>% tidyr::drop_na()
USA_Tri_RR_Fa_pro <- ts(USA_Tri_RR_Fa_pro %>% select(-1), frequency = 4, start = c(floor(USA_Tri_RR_Fa_pro$time[1]), (USA_Tri_RR_Fa_pro$time[1] %% 1) / 0.25 + 1 ))
VARselect(USA_Tri_RR_Fa_pro, lag.max = 8)
var3_RR_Fa_pro <- VAR(USA_Tri_RR_Fa_pro, type = "const", p = 3)

test.GC(var3_RR_Fa_pro, C3)
test.GC(var3_RR_Fa_pro, C4)


# SZ
USA_Tri_SZ_Fa_pro <- data.frame("time" = USA_FAVAR %>% time %>% as.numeric, USA_FAVAR) %>% left_join(iv2ts(SZ_pro , "SZ_pro"), by = "time") %>% tidyr::drop_na()
USA_Tri_SZ_Fa_pro <- ts(USA_Tri_SZ_Fa_pro %>% select(-1), frequency = 4, start = c(floor(USA_Tri_SZ_Fa_pro$time[1]), (USA_Tri_SZ_Fa_pro$time[1] %% 1) / 0.25 + 1 ))
VARselect(USA_Tri_SZ_Fa_pro, lag.max = 8)
var3_SZ_Fa_pro <- VAR(USA_Tri_SZ_Fa_pro, type = "const", p = 3)

test.GC(var3_SZ_Fa_pro, C3)
test.GC(var3_SZ_Fa_pro, C4)


# gute schleife -----------------------------------------------------------
# SW
USA_Tri_SW <- data.frame("time" = USA_Tri %>% time %>% as.numeric, USA_Tri) %>% left_join(SW.ts, by = "time") %>% tidyr::drop_na()
USA_Tri_SW <- ts(USA_Tri_SW %>% select(-1), frequency = 4, start = c(floor(USA_Tri_SW$time[1]), (USA_Tri_SW$time[1] %% 1) / 0.25 + 1 ))
var3_SW <- VAR(USA_Tri_SW, type = "const", p = 4)

tab_SW <- c(
  test.GC(var3_SW, C1)$TestStat$lambda_F,
  test.GC(var3_SW, C1)$p_val$p_F)

SDFM <- FALSE
for (r in 1:3) {
  C3 <- matrix(c(rep(c(0,0,0,rep(0,r),1),3+r), rep(0,4+r)), nrow = 4+r, ncol = 4+r, byrow = T)
  source("Fa.R")
  colnames(USA_FAVAR) <- c(colnames(USA_Tri), 1:r)
  USA_Tri_SW_Fa <- data.frame("time" = USA_FAVAR %>% time %>% as.numeric, USA_FAVAR) %>% left_join(SW.ts, by = "time") %>% tidyr::drop_na()
  USA_Tri_SW_Fa <- ts(USA_Tri_SW_Fa %>% select(-1), frequency = 4, start = c(floor(USA_Tri_SW_Fa$time[1]), (USA_Tri_SW_Fa$time[1] %% 1) / 0.25 + 1 ))
  var3_SW_Fa <- VAR(USA_Tri_SW_Fa, type = "const", p = 4)
  tab_SW <- rbind(tab_SW, c(
    test.GC(var3_SW_Fa, C3)$TestStat$lambda_F,
    test.GC(var3_SW_Fa, C3)$p_val$p_F))
}


# RR
USA_Tri_RR <- data.frame("time" = USA_Tri %>% time %>% as.numeric, USA_Tri) %>% left_join(RR.ts, by = "time") %>% tidyr::drop_na()
USA_Tri_RR <- ts(USA_Tri_RR %>% select(-1), frequency = 4, start = c(floor(USA_Tri_RR$time[1]), (USA_Tri_RR$time[1] %% 1) / 0.25 + 1 ))
var3_RR <- VAR(USA_Tri_RR, type = "const", p = 4)

tab_RR <- c(
  test.GC(var3_RR, C1)$TestStat$lambda_F,
  test.GC(var3_RR, C1)$p_val$p_F)

for (r in 1:3) {
  C3 <- matrix(c(rep(c(0,0,0,rep(0,r),1),3+r), rep(0,4+r)), nrow = 4+r, ncol = 4+r, byrow = T)
  source("Fa.R")
  colnames(USA_FAVAR) <- c(colnames(USA_Tri), 1:r)
  USA_Tri_RR_Fa <- data.frame("time" = USA_FAVAR %>% time %>% as.numeric, USA_FAVAR) %>% left_join(RR.ts, by = "time") %>% tidyr::drop_na()
  USA_Tri_RR_Fa <- ts(USA_Tri_RR_Fa %>% select(-1), frequency = 4, start = c(floor(USA_Tri_RR_Fa$time[1]), (USA_Tri_RR_Fa$time[1] %% 1) / 0.25 + 1 ))
  var3_RR_Fa <- VAR(USA_Tri_RR_Fa, type = "const", p = 4)
  tab_RR <- rbind(tab_RR, c(
    test.GC(var3_RR_Fa, C3)$TestStat$lambda_F,
    test.GC(var3_RR_Fa, C3)$p_val$p_F))
}

# SZ
USA_Tri_SZ <- data.frame("time" = USA_Tri %>% time %>% as.numeric, USA_Tri) %>% left_join(SZ.ts, by = "time") %>% tidyr::drop_na()
USA_Tri_SZ <- ts(USA_Tri_SZ %>% select(-1), frequency = 4, start = c(floor(USA_Tri_SZ$time[1]), (USA_Tri_SZ$time[1] %% 1) / 0.25 + 1 ))
var3_SZ <- VAR(USA_Tri_SZ, p = 4, type = "const")

tab_SZ <- c(
  test.GC(var3_SZ, C1)$TestStat$lambda_F,
  test.GC(var3_SZ, C1)$p_val$p_F)

for (r in 1:3) {
  C3 <- matrix(c(rep(c(0,0,0,rep(0,r),1),3+r), rep(0,4+r)), nrow = 4+r, ncol = 4+r, byrow = T)
  source("Fa.R")
  colnames(USA_FAVAR) <- c(colnames(USA_Tri), 1:r)
  USA_Tri_SZ_Fa <- data.frame("time" = USA_FAVAR %>% time %>% as.numeric, USA_FAVAR) %>% left_join(SZ.ts, by = "time") %>% tidyr::drop_na()
  USA_Tri_SZ_Fa <- ts(USA_Tri_SZ_Fa %>% select(-1), frequency = 4, start = c(floor(USA_Tri_SZ_Fa$time[1]), (USA_Tri_SZ_Fa$time[1] %% 1) / 0.25 + 1 ))
  var3_SZ_Fa <- VAR(USA_Tri_SZ_Fa, type = "const", p = 4)
  tab_SZ <- rbind(tab_SZ, c(
    test.GC(var3_SZ_Fa, C3)$TestStat$lambda_F,
    test.GC(var3_SZ_Fa, C3)$p_val$p_F))
}

# RR_pro
# RR
USA_Tri_RR_pro <- data.frame("time" = USA_Tri %>% time %>% as.numeric, USA_Tri) %>% left_join(iv2ts(RR_pro , "RR_pro"), by = "time") %>% tidyr::drop_na()
USA_Tri_RR_pro <- ts(USA_Tri_RR_pro %>% select(-1), frequency = 4, start = c(floor(USA_Tri_RR_pro$time[1]), (USA_Tri_RR_pro$time[1] %% 1) / 0.25 + 1 ))
var3_RR_pro <- VAR(USA_Tri_RR_pro, type = "const", p = 4)

tab_RR_pro <- c(
  test.GC(var3_RR_pro, C1)$TestStat$lambda_F,
  test.GC(var3_RR_pro, C1)$p_val$p_F)


for (r in 1:3) {
  C3 <- matrix(c(rep(c(0,0,0,rep(0,r),1),3+r), rep(0,4+r)), nrow = 4+r, ncol = 4+r, byrow = T)
  source("Fa.R")
  colnames(USA_FAVAR) <- c(colnames(USA_Tri), 1:r)
  USA_Tri_RR_Fa_pro <- data.frame("time" = USA_FAVAR %>% time %>% as.numeric, USA_FAVAR) %>% left_join(iv2ts(RR_pro , "RR_pro"), by = "time") %>% tidyr::drop_na()
  USA_Tri_RR_Fa_pro <- ts(USA_Tri_RR_Fa_pro %>% select(-1), frequency = 4, start = c(floor(USA_Tri_RR_Fa_pro$time[1]), (USA_Tri_RR_Fa_pro$time[1] %% 1) / 0.25 + 1 ))
  var3_RR_Fa_pro <- VAR(USA_Tri_RR_Fa_pro, type = "const", p = 4)
  tab_RR_pro <- rbind(tab_RR_pro, c(
    test.GC(var3_RR_Fa_pro, C3)$TestStat$lambda_F,
    test.GC(var3_RR_Fa_pro, C3)$p_val$p_F))
}

# SZ_pro
USA_Tri_SZ_pro <- data.frame("time" = USA_Tri %>% time %>% as.numeric, USA_Tri) %>% left_join(iv2ts(SZ_pro , "SZ_pro"), by = "time") %>% tidyr::drop_na()
USA_Tri_SZ_pro <- ts(USA_Tri_SZ_pro %>% select(-1), frequency = 4, start = c(floor(USA_Tri_SZ_pro$time[1]), (USA_Tri_SZ_pro$time[1] %% 1) / 0.25 + 1 ))
var3_SZ_pro <- VAR(USA_Tri_SZ_pro, p = 4, type = "const")

tab_SZ_pro <- c(
  test.GC(var3_SZ_pro, C1)$TestStat$lambda_F,
  test.GC(var3_SZ_pro, C1)$p_val$p_F)

for (r in 1:3) {
  C3 <- matrix(c(rep(c(0,0,0,rep(0,r),1),3+r), rep(0,4+r)), nrow = 4+r, ncol = 4+r, byrow = T)
  source("Fa.R")
  colnames(USA_FAVAR) <- c(colnames(USA_Tri), 1:r)
  USA_Tri_SZ_Fa_pro <- data.frame("time" = USA_FAVAR %>% time %>% as.numeric, USA_FAVAR) %>% left_join(iv2ts(SZ_pro , "SZ_pro"), by = "time") %>% tidyr::drop_na()
  USA_Tri_SZ_Fa_pro <- ts(USA_Tri_SZ_Fa_pro %>% select(-1), frequency = 4, start = c(floor(USA_Tri_SZ_Fa_pro$time[1]), (USA_Tri_SZ_Fa_pro$time[1] %% 1) / 0.25 + 1 ))
  var3_SZ_Fa_pro <- VAR(USA_Tri_SZ_Fa_pro, type = "const", p = 4)
  tab_SZ_pro <- rbind(tab_SZ_pro, c(
    test.GC(var3_SZ_Fa_pro, C3)$TestStat$lambda_F,
    test.GC(var3_SZ_Fa_pro, C3)$p_val$p_F))
}

tab <- cbind(tab_SW, tab_RR, tab_SZ, tab_RR_pro, tab_SZ_pro)
rownames(tab) <- 0:3
tab %>% kableExtra::kable(format = "latex", booktabs = T, digits = 3, linesep = "")



