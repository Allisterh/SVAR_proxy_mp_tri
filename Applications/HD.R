
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
hd_pi_RR_mp_pro <- hd.mod(IV.RR_pro, series = 2, Partial = 3, Epsname)
hd_pi_SZ_mp_pro <- hd.mod(IV.SZ_pro, series = 2, Partial = 3, Epsname)

hd_pi_SR %>% plot.hd.mod()
hd_pi_CV %>% plot.hd.mod()
hd_pi_dCov %>% plot.hd.mod()
hd_pi_SW_mp %>% plot.hd.mod()
hd_pi_RR_mp %>% plot.hd.mod()
hd_pi_SZ_mp %>% plot.hd.mod()
hd_pi_RR_mp_pro %>% plot.hd.mod()
hd_pi_SZ_mp_pro %>% plot.hd.mod()

## combined plot
hd_pi_mp <- hd_pi_SR_mp %>% left_join(hd_pi_dCov_mp, by = "t") %>% 
  left_join(hd_pi_CV_mp, by = "t") %>% left_join(hd_pi_SW_mp, by = "t") %>% 
  left_join(hd_pi_SZ_mp, by = "t") %>% left_join(hd_pi_RR_mp, by = "t") %>%
  left_join(hd_pi_SZ_mp_pro, by = "t") %>% left_join(hd_pi_RR_mp_pro, by = "t") %>% tidyr::drop_na()

hd_pi_mp %<>% dplyr::filter(t >= 1973 & t <= 1983) 

colnames(hd_pi_mp)[-1] <- c("SR", "dCov", "CV", "SW", "SZ", "RR", "SZ_pro", "RR_pro")

recessions_start <- c( 1973.75, 1980,1981.5)
recessions_end <- c( 1975, 1980.5,1982.75)

hd_pi_mp %<>% reshape2::melt(id = "t")

## subsets
hd_pi_mp %<>% add_column(lab = "old") %>% mutate(lab = ifelse(variable == "SZ_pro" | variable == "RR_pro", "new", "old"))
hd_pi_mp[247:287,2] <- "SZ"
hd_pi_mp[288:328,2] <- "RR"

hd_pi_mp$lab %>% factor()

hd_pi_mp_plot <- hd_pi_mp %>% ggplot(aes(x = t, y = value)) + facet_wrap(~variable, ncol = 3) +
  theme_bw() + geom_hline(yintercept = 0, color = "red", alpha = 0.5) +
  scale_x_continuous(breaks = seq(hd_pi_mp$t[1] %>% ceiling(),hd_pi_mp$t[nrow(hd_pi_mp)], by = 2)) +  
  annotate("rect", xmin = recessions_start, xmax = recessions_end, ymin = -Inf, ymax = Inf, alpha = 0.3) +
  geom_line(data = subset(hd_pi_mp, lab == "old"), linetype = "solid") +
  geom_line(data = subset(hd_pi_mp, lab == "new"), linetype = "dashed") +
  xlab("") + ylab("")

ggsave("../Plots/hd_pi_mp_pro.pdf", plot = hd_pi_mp_plot, width = 12, height = 6)


## CV HD
hd_pi_CV_df <- hd_pi_CV %>% dplyr::filter(t >= 1973 & t <= 1983) %>% dplyr::select(c(1, 4, 5))
colnames(hd_pi_CV_df)[c(2,3)] <- c("d", "s")
hd_pi_CV_df %<>% reshape2::melt(id = "t")

hd_pi_CV_plot <- hd_pi_CV_df %>% ggplot(aes(x = t, y = value)) + facet_wrap(~variable, ncol = 1) +
  theme_bw() + geom_hline(yintercept = 0, color = "red", alpha = 0.5) +
  scale_x_continuous(breaks = seq(hd_pi_CV_df$t[1] %>% ceiling(),hd_pi_CV_df$t[nrow(hd_pi_CV_df)], by = 2)) +  
  annotate("rect", xmin = recessions_start, xmax = recessions_end, ymin = -Inf, ymax = Inf, alpha = 0.3) +
  geom_line() +  xlab("") + ylab("")

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

hd_x_mp %<>% dplyr::filter(t >= 1973 & t <= 1983) 

colnames(hd_x_mp)[-1] <- c("SR", "dCov", "CV", "SW", "SZ", "RR")

hd_x_mp %>% reshape2::melt(id = "t") %>% ggplot(aes(x = t, y = value)) +
  geom_line() + facet_wrap(~variable, ncol = 3) + xlab("") + ylab("") + 
  theme_bw() + geom_hline(yintercept = 0, color = "red", alpha = 0.5) +
  scale_x_continuous(breaks = seq(hd_x_mp$t[1] %>% ceiling(),hd_x_mp$t[nrow(hd_x_mp)], by = 2)) +  
  annotate("rect", xmin = recessions_start, xmax = recessions_end, ymin = -Inf, ymax = Inf, alpha = 0.3)



