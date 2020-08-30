# Load functions ----------------------------------------------------------
rm(list = ls())
graphics.off()
Functions <- devtools::as.package("../Functions")
devtools::load_all("../Functions")


# AOA ---------------------------------------------------------------------
AOA <- readRDS("Replica/DGP1.AOA.rds")


# Baseline distribution scenario ------------------------------------------
dist  <- "Gaussian"


# Read simulation results -------------------------------------------------
path  <- paste0("./Replica/", dist, "/results/")

SR    <- readRDS(paste0(path, "simu.sign.rds"))
CV    <- readRDS(paste0(path, "simu.cv.1.rds"))
dCov  <- readRDS(paste0(path, "simu.dcov.rds"))
CvM   <- readRDS(paste0(path, "simu.cvm.rds"))
IV    <- readRDS(paste0(path, "simu.iv.rds"))
IV_W  <- readRDS(paste0(path, "simu.wiv.rds"))
IV_P  <- readRDS(paste0(path, "simu.iiv.rds"))


# Plot --------------------------------------------------------------------
plot.simu.multi(list(SR$M, SR$MT, CV, dCov, CvM, IV, IV_W, IV_P), 
                Label = c("SR", "SR-MT", "CV", "dCov", "CvM", "IV", "IV-W", "IV-P"), 
                shape_manual = c(15:17, 5, 1:4))

# Other distribution scenarios? -------------------------------------------
dist  <- "Chisq"        # Run code line 16-31
dist  <- "t"            # Run code line 16-31
dist  <- "mixed_a"      # Run code line 16-31
dist  <- "mixed_b"      # Run code line 16-31


# Large sample size? ------------------------------------------------------
dist  <- "large_sample" # Run code line 16-31
