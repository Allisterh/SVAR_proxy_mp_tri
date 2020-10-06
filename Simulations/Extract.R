# Load functions ----------------------------------------------------------
rm(list = ls())
Functions <- devtools::as.package("../Functions")
devtools::load_all("../Functions")

IDM   <- c("SR", "CV", "dCov", "CvM", "IV", "IV_W", "IV_P")
Lab   <- c("SR", "CV", "dCov", "CvM", "IV", "IV-W", "IV-E")
Dist  <- c("L1", "L2a", "L2b", "L3a", "L3b")

# Graphics ----------------------------------------------------------------
Tobs <- 240
for (i in 1:5) {
  aoa_temp <- list()
  for (j in 1:length(IDM)) {
    Path_In <- paste0("../Server/", "T", Tobs, "/", Dist[i], "/return/", IDM[j], ".rds")
    temp <- readRDS(Path_In)
    if (IDM[j] == "SR"){
      aoa_temp[[j]] <- temp$M
    } else {
      aoa_temp[[j]] <- temp$aoa
    }
  }
  plot_temp <- plot.simu.multi(aoa_temp, Label = Lab, shape_manual = c(5,1,2, 15:17,3), level_manual = c("IV", "IV-W", "IV-E", "CV", "CvM", "dCov", "SR"))
  Path_Out <- paste0("../Plots/")
  ggsave(paste0(Path_Out, "T", Tobs, "_", Dist[i], ".pdf"), plot = plot_temp, width = 12, height = 8)
  rm(temp, aoa_temp, plot_temp)
}

# Tables ------------------------------------------------------------------
TT <- c(120, 240, 480)

## AoA tables
aoa_s <- matrix(NA, nrow = length(IDM)*5, ncol = 12)
aoa_l <- matrix(NA, nrow = length(IDM)*5, ncol = 12)
for (i in 1:5) {
  for (t in 1:3) {
    aoa_temp <- list()
    for (j in 1:length(IDM)) {
      Path_In <- paste0("../Server/", "T", TT[t], "/", Dist[i], "/return/", IDM[j], ".rds")
      temp <- readRDS(Path_In)
      if (IDM[j] == "SR"){
        aoa_temp[[j]] <- temp$M
      } else {
        aoa_temp[[j]] <- temp$aoa
      }
    }
    tab_temp <- tabulate.simu.aoa(aoa_temp, Label = Lab)
    aoa_s[c((length(IDM)*(i-1)+1):(i*length(IDM))),c((4*(t-1)+1):(t*4))] <- tab_temp$short
    aoa_l[c((length(IDM)*(i-1)+1):(i*length(IDM))),c((4*(t-1)+1):(t*4))] <- tab_temp$long
  }
  rm(temp, aoa_temp, tab_temp)
}
rownames(aoa_s) <- rownames(aoa_l) <- rep(Lab, 5)  

kableExtra::kable(aoa_s, format = "latex", col.names = 1:12, booktabs = T, digits = 3, linesep = "") %>% kableExtra::row_spec(c(1:5)*length(IDM), hline_after = T)
kableExtra::kable(aoa_l, format = "latex", col.names = 1:12, booktabs = T, digits = 3, linesep = "") %>% kableExtra::row_spec(c(1:5)*length(IDM), hline_after = T)

## UMP table
ump_tab <- matrix(NA, nrow = length(IDM), ncol = 15)
for (i in 1:5) {
  for (t in 1:3) {
    ump_temp <- list()
    for (j in 1:length(IDM)) {
      Path_In <- paste0("../Server/", "T", TT[t], "/", Dist[i], "/return/", IDM[j], ".rds")
      ump_temp[[j]] <- readRDS(Path_In)
    }
    ump_tab[, ((i-1)*3+t)] <- tabulate.simu.ump(ump_temp, Lab)$ump
  }
  rm(ump_temp)
}
rownames(ump_tab) <- Lab
kableExtra::kable(ump_tab, format = "latex", col.names = 1:15, booktabs = T, digits = 3, linesep = "")

## RMSE table
rmse_tab <- matrix(NA, nrow = length(IDM)*5, ncol = 12)

for (i in 1:5) {
  for (t in 1:3) {
    rmse_temp <- list()
    for (j in 1:length(IDM)) {
      Path_In <- paste0("../Server/", "T", TT[t], "/", Dist[i], "/return/", IDM[j], ".rds")
      rmse_temp[[j]] <- readRDS(Path_In)
    }
    rmse_tab[c((length(IDM)*(i-1)+1):(i*length(IDM))),c((4*(t-1)+1):(t*4))] <- tabulate.simu.ump(rmse_temp, Lab)$rmse
  }
  rm(rmse_temp)
}
rownames(rmse_tab) <- rep(Lab, 5)
kableExtra::kable(rmse_tab, format = "latex", col.names = 1:12, booktabs = T, digits = 3, linesep = "") %>% kableExtra::row_spec(c(1:5)*length(IDM), hline_after = T)


## RMSE plot
rmse_plot_temp <- rmse_tab[-7*c(1:5),c(-1, -5, -9)] %>% as.data.frame()

Dist_lab <- c("L[1]", "L[2]^a", "L[2]^b", "L[3]^a", "L[3]^b")

rmse_plot <- data.frame("T" = rep(TT,each = 3), "variable" = c("x", "pi", "r"), "value" = unlist(rmse_plot_temp[1,]), "Methods" = Lab[1])
for (m in 2:6) {
  rmse_temp <- data.frame("T" = rep(TT,each = 3), "variable" = c("x", "pi", "r"), "value" = unlist(rmse_plot_temp[m,]), "Methods" = Lab[m])
  rmse_plot <- rbind(rmse_plot, rmse_temp)
}
rmse_plot <- rmse_plot %>% add_column(case = "L[1]")

for (d in 2:5) {
  rmse_plot_new <- data.frame("T" = rep(TT,each = 3), "variable" = c("x", "pi", "r"), "value" = unlist(rmse_plot_temp[(6*(d-1)+1),]), "Methods" = Lab[1])
  for (m in 2:6) {
    rmse_temp <- data.frame("T" = rep(TT,each = 3), "variable" = c("x", "pi", "r"), "value" = unlist(rmse_plot_temp[(6*(d-1)+m),]), "Methods" = Lab[m])
    rmse_plot_new <- rbind(rmse_plot_new, rmse_temp)
  }
  rmse_plot_new <- rmse_plot_new %>% add_column(case = Dist_lab[d])
  rmse_plot <- rbind(rmse_plot, rmse_plot_new)
}

library(latex2exp)
levels(rmse_plot$variable) <- c(x = TeX("$b_{1,3}$"), pi = TeX("$b_{2,3}$"), r = TeX("$b_{3,3}$"))

rmse_plot$T <- as.character(rmse_plot$T)

ggplot(data = rmse_plot, aes(x = T, y = value, group = Methods)) + geom_line(aes(linetype = Methods, color = Methods)) +
  geom_point(aes(shape = Methods, color = Methods), size=3) +
  facet_grid(variable~case, scales = "free_y", labeller = label_parsed) +
  xlab(" ") + ylab(" ")  + theme_bw() + scale_shape_manual(values = c(15:17, 5, 2:3))

ggsave("../Plots/RMSE.pdf", plot = last_plot(), width = 12, height = 8)


# RMSE plot (new) ---------------------------------------------------------
rmse_plot_temp <- rmse_tab[-7*c(1:5),] %>% as.data.frame()

Dist_lab <- c("L[1]", "L[2]^a", "L[2]^b", "L[3]^a", "L[3]^b")

rmse_plot <- data.frame("T" = rep(TT,each = 4), "variable" = c("mean", "x", "pi", "r"), "value" = unlist(rmse_plot_temp[1,]), "Methods" = Lab[1])
for (m in 2:6) {
  rmse_temp <- data.frame("T" = rep(TT,each =  4), "variable" = c("mean", "x", "pi", "r"), "value" = unlist(rmse_plot_temp[m,]), "Methods" = Lab[m])
  rmse_plot <- rbind(rmse_plot, rmse_temp)
}
rmse_plot <- rmse_plot %>% add_column(case = "L[1]")

for (d in 2:5) {
  rmse_plot_new <- data.frame("T" = rep(TT,each = 4), "variable" = c("mean","x", "pi", "r"), "value" = unlist(rmse_plot_temp[(6*(d-1)+1),]), "Methods" = Lab[1])
  for (m in 2:6) {
    rmse_temp <- data.frame("T" = rep(TT,each = 4), "variable" = c("mean","x", "pi", "r"), "value" = unlist(rmse_plot_temp[(6*(d-1)+m),]), "Methods" = Lab[m])
    rmse_plot_new <- rbind(rmse_plot_new, rmse_temp)
  }
  rmse_plot_new <- rmse_plot_new %>% add_column(case = Dist_lab[d])
  rmse_plot <- rbind(rmse_plot, rmse_plot_new)
}

library(latex2exp)
levels(rmse_plot$variable) <- c(mean = TeX("$\\bar{MSE}$"), x = TeX("$b_{1,3}$"), pi = TeX("$b_{2,3}$"), r = TeX("$b_{3,3}$"))
rmse_plot$case %>% factor()

rmse_plot$T <- as.character(rmse_plot$T)

ggplot(data = rmse_plot, aes(x = T, y = value, group = Methods)) + geom_line(aes(linetype = Methods, color = Methods)) +
  geom_point(aes(shape = Methods, color = Methods), size=3) +
  facet_grid(variable~case, scales = "free_y", labeller = label_parsed) +
  xlab(" ") + ylab(" ")  + theme_bw() + scale_shape_manual(values = c(15:17, 5, 2:3))

ggsave("../Plots/RMSE_new.pdf", plot = last_plot(), width = 12, height = 8)

# UMPS Plot ---------------------------------------------------------------
ump_plot_int <- as.data.frame(ump_tab)
ump_plot_int["SR",] <- NA
ump_plot_int <- t(ump_plot_int) %>% as.data.frame()
ump_plot_int <- cbind("Dist" = rep(Dist, each = 3), ump_plot_int)  %>% reshape2::melt(id = "Dist")
ump_plot_int <- ump_plot_int %>% add_column(Tob = rep(TT, 5*7))
ump_plot_int <- ump_plot_int[,c(1,2,4,3)]
ump_plot_int <- ump_plot_int %>% add_column(measure = "Frequency")

mean_rmse <- rmse_tab[,c(1,5,9)] %>% as.data.frame()
rmse_plot_int <- mean_rmse[1:7,]
for (i in 2:5) {
  rmse_plot_int <- cbind(rmse_plot_int, mean_rmse[(7*(i-1)+1):(7*i),])
}
rownames(rmse_plot_int) <- Lab
#rmse_plot_int["IV-P",] <- NA
rmse_plot_int <- t(rmse_plot_int) %>% as.data.frame()
rmse_plot_int <- cbind("Dist" = rep(Dist, each = 3), rmse_plot_int)  %>% reshape2::melt(id = "Dist")
rmse_plot_int <- rmse_plot_int %>% add_column(Tob = rep(TT, 5*7))
#rmse_plot_int <- rmse_plot_int[,c(1,2,4,3)]
rmse_plot_int <- rmse_plot_int %>% add_column(measure = "RMSE")

UMPS_plot <- rbind(ump_plot_int, rmse_plot_int)


library(latex2exp)
levels(UMPS_plot$Dist) <- c(TeX("$\\textbf{L}_I$"), TeX("$\\textbf{L}^{(a)}_{II}$"), 
                            TeX("$\\textbf{L}^{(b)}_{II}$"), TeX("$\\textbf{L}^{(a)}_{III}$"),
                            TeX("$\\textbf{L}^{(b)}_{III}$"))

colnames(UMPS_plot)[2] <- "Methods"

UMPS_plot$measure<- UMPS_plot$measure %>% as.factor()
levels(UMPS_plot$measure)[2] <- TeX("$\\bar{RMSE}$")

UMPS_plot$Methods <- UMPS_plot$Methods %>% as.character()
UMPS_plot$Tob <- UMPS_plot$Tob %>% as.character()

UMPS_plot <- UMPS_plot %>% mutate(Methods = replace(Methods, Methods == "IV-P", "IV-E"))

UMPS_plot$Methods <- factor(  UMPS_plot$Methods, levels = c("IV", "IV-W", "IV-E", "CV", "CvM", "dCov", "SR"))


ggplot(data = UMPS_plot, aes(x = Tob, y = value, group = Methods)) + geom_line(aes(linetype = Methods, color = Methods)) +
  geom_point(aes(shape = Methods, color = Methods), size=3) +
  facet_grid(measure~Dist, scales = "free_y", labeller = label_parsed) +
  xlab(" T ") + ylab(" ")  + theme_bw() + scale_shape_manual(values = c(15:17, 5, 1:3)) +
  scale_shape_manual(values = c(5,1,2, 15:17,3))

ggsave("../Plots/UMPS_all.pdf", plot = last_plot(), width = 12, height = 7)



