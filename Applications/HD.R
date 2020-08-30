
Epsname = c("d", "s", "mp")

# pi ----------------------------------------------------------------------
hd_pi_SR <- hd.mod(SR, series = 2, Partial = NULL, Epsname)
hd_pi_SR_mp <- hd.mod(SR, series = 2, Partial = 3, Epsname)
hd_pi_CV <- hd.mod(CV, series = 2, Partial = NULL, Epsname)
hd_pi_CV_mp <- hd.mod(CV, series = 2, Partial = 3, Epsname)
hd_pi_dCov <- hd.mod(dCov, series = 2, Partial = NULL, Epsname)
hd_pi_dCov_mp <- hd.mod(dCov, series = 2, Partial = 3, Epsname)
hd_pi_SW_mp <- hd.mod(IV.SW, series = 2, Partial = 3, Epsname)
hd_pi_RR_mp <- hd.mod(IV.RR, series = 2, Partial = 3, Epsname)
hd_pi_SZ_mp <- hd.mod(IV.SZ, series = 2, Partial = 3, Epsname)


hd_pi_SR %>% plot.hd.mod()
hd_pi_CV %>% plot.hd.mod()
hd_pi_dCov %>% plot.hd.mod()
hd_pi_SW_mp %>% plot.hd.mod()
hd_pi_RR_mp %>% plot.hd.mod()
hd_pi_SZ_mp %>% plot.hd.mod()

## combined plot
hd_pi_mp <- hd_pi_SR_mp %>% left_join(hd_pi_dCov_mp, by = "t") %>% 
  left_join(hd_pi_CV_mp, by = "t") %>% left_join(hd_pi_SW_mp, by = "t") %>% 
  left_join(hd_pi_SZ_mp, by = "t") %>% left_join(hd_pi_RR_mp, by = "t") %>% tidyr::drop_na()

hd_pi_mp %<>% dplyr::filter(t >= 1970 & t <= 1985) 

colnames(hd_pi_mp)[-1] <- c("SR", "dCov", "CV", "SW", "SZ", "RR")

recessions_start <- c(1969.75, 1973.75, 1980)
recessions_end <- c(1970.75, 1975, 1982.75)

hd_pi_mp_plot <- hd_pi_mp %>% reshape2::melt(id = "t") %>% ggplot(aes(x = t, y = value)) +
  geom_line() + facet_wrap(~variable, ncol = 3) + xlab("") + ylab("") + 
  theme_bw() + geom_hline(yintercept = 0, color = "red", alpha = 0.5) +
  annotate("rect", xmin = recessions_start, xmax = recessions_end, ymin = -Inf, ymax = Inf, alpha = 0.3)


# x ----------------------------------------------------------------------
hd_x_SR <- hd.mod(SR, series = 1, Partial = NULL, Epsname)
hd_x_SR_mp <- hd.mod(SR, series = 1, Partial = 3, Epsname)
hd_x_CV <- hd.mod(CV, series = 1, Partial = NULL, Epsname)
hd_x_CV_mp <- hd.mod(CV, series = 1, Partial = 3, Epsname)
hd_x_dCov <- hd.mod(dCov, series = 1, Partial = NULL, Epsname)
hd_x_dCov_mp <- hd.mod(dCov, series = 1, Partial = 3, Epsname)
hd_x_SW_mp <- hd.mod(IV.SW, series = 1, Partial = 3, Epsname)
hd_x_RR_mp <- hd.mod(IV.RR, series = 1, Partial = 3, Epsname)
hd_x_SZ_mp <- hd.mod(IV.SZ, series = 1, Partial = 3, Epsname)

hd_x_SR %>% plot.hd.mod()
hd_x_CV %>% plot.hd.mod()
hd_x_dCov %>% plot.hd.mod()
hd_x_SW_mp %>% plot.hd.mod()
hd_x_RR_mp %>% plot.hd.mod()
hd_x_SZ_mp %>% plot.hd.mod()

## combined plot
hd_x_mp <- hd_x_SR_mp %>% left_join(hd_x_dCov_mp, by = "t") %>% 
  left_join(hd_x_CV_mp, by = "t") %>% left_join(hd_x_SW_mp, by = "t") %>% 
  left_join(hd_x_SZ_mp, by = "t") %>% left_join(hd_x_RR_mp, by = "t") %>% tidyr::drop_na()

hd_x_mp %<>% dplyr::filter(t >= 1970 & t <= 1985) 

colnames(hd_x_mp)[-1] <- c("SR", "dCov", "CV", "SW", "SZ", "RR")

hd_x_mp %>% reshape2::melt(id = "t") %>% ggplot(aes(x = t, y = value)) +
  geom_line() + facet_wrap(~variable, ncol = 3) + xlab("") + ylab("") + 
  theme_bw() + geom_hline(yintercept = 0, color = "red", alpha = 0.5)



