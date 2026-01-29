
# posterior histograms
# of state transition probabilities
# using BUGS output


data("jags_output")
R2jags::attach.jags(mm1)


################
# base version #
################

# with ICD
x11()
par(mfrow = c(2,3))
# obs
hist(lambda.0[, 1, 1], breaks = 20, xlab = "Prob", main = "ICD -> ICD")
abline(v = mean(lambda.0[, 1, 1]), col = "red")
hist(lambda.0[, 1, 2], breaks = 20, xlab = "Prob", main = "ICD -> Shock")
abline(v = mean(lambda.0[, 1, 2]), col = "red")
hist(lambda.0[, 1, 3], breaks = 20, xlab = "Prob", main = "ICD -> All-cause death")
abline(v = mean(lambda.0[, 1, 3]), col = "red")
# > 6%
hist(lambda.1[, 1, 1], breaks = 20, xlab = "Prob", main = "ICD -> ICD")
abline(v = mean(lambda.1[, 1, 1]), col = "red")
hist(lambda.1[, 1, 2], breaks = 20, xlab = "Prob", main = "ICD -> Shock")
abline(v = mean(lambda.1[, 1, 2]), col = "red")
hist(lambda.1[, 1, 3], breaks = 20, xlab = "Prob", main = "ICD -> All-cause death")
abline(v = mean(lambda.1[, 1, 3]), col = "red")

# low risk
x11()
# obs
par(mfrow = c(2,3))
hist(lambda.0[, 4, 4], breaks = 20, xlab = "Prob", main = "Low risk -> Low risk")
abline(v = mean(lambda.0[, 4, 4]), col = "red")
hist(lambda.0[, 4, 5], breaks = 20, xlab = "Prob", main = "Low risk -> SCD")
abline(v = mean(lambda.0[, 4, 5]), col = "red")
hist(lambda.0[, 4, 6], breaks = 20, xlab = "Prob", main = "Low risk -> All-cause death")
abline(v = mean(lambda.0[, 4, 6]), col = "red")
# > 6%
hist(lambda.1[, 4, 4], breaks = 20, xlab = "Prob", main = "Low risk -> Low risk")
abline(v = mean(lambda.1[, 4, 4]), col = "red")
hist(lambda.1[, 4, 5], breaks = 20, xlab = "Prob", main = "Low risk -> SCD")
abline(v = mean(lambda.1[, 4, 5]), col = "red")
hist(lambda.1[, 4, 6], breaks = 20, xlab = "Prob", main = "Low risk -> All-cause death")
abline(v = mean(lambda.1[, 4, 6]), col = "red")


###################
# ggplot2 version #
###################

library(reshape2)
library(ggplot2)
library(dplyr)

## no ICD
# lambda: [1:1000, 1:6, 1:6]
pobs <- as_tibble(melt(lambda.0[, 4, ]))
prisk6 <- as_tibble(melt(lambda.1[, 4, ]))
prisk4 <- as_tibble(melt(lambda.2[, 4, ]))

pobs <- cbind(pobs, strat = "obs")
prisk6 <- cbind(prisk6, strat = "risk6")
prisk4 <- cbind(prisk4, strat = "risk4")

plot_dat_no <-
  rbind.data.frame(pobs, prisk6, prisk4) %>%
  dplyr::rename(sim = Var1,
                state = Var2) %>%
  filter(state %in% c(4,5,6)) %>%
  mutate(state =
           ifelse(state == 4, "HCM",
                  ifelse(state == 5, "SCD/shock",
                         "All-cause mortality")))

ggplot(plot_dat_no, aes(x = value, y = ..density..)) +
  facet_grid(strat ~ state, scales = "free") +
  # geom_histogram() +
  geom_density() +
  ggtitle("No ICD")

ggsave(filename = "images/post_hist_noICD.png")

## with ICD
pobs <- as_tibble(melt(lambda.0[, 1, ]))
prisk6 <- as_tibble(melt(lambda.1[, 1, ]))
prisk4 <- as_tibble(melt(lambda.2[, 1, ]))

pobs <- cbind(pobs, strat = "obs")
prisk6 <- cbind(prisk6, strat = "risk6")
prisk4 <- cbind(prisk4, strat = "risk4")

plot_dat_icd <-
  rbind.data.frame(pobs, prisk6, prisk4) %>%
  filter(Var2 %in% c(1,2,3)) %>%
  dplyr::rename(sim = Var1,
                state = Var2) %>%
  mutate(state =
           ifelse(state == 1, "HCM",
                  ifelse(state == 2, "SCD/shock",
                         "All-cause mortality")))

ggplot(plot_dat_icd, aes(x = value, y = ..density..)) +
  facet_grid(strat ~ state, scales = "free") +
  # geom_histogram() +
  geom_density() +
  ggtitle("With ICD")

ggsave(filename = "images/post_hist_withICD.png")

## combined plot

plot_dat_no <- cbind(plot_dat_no, icd = FALSE)
plot_dat_icd <- cbind(plot_dat_icd, icd = TRUE)

plot_dat <-
  rbind.data.frame(plot_dat_icd,
                   plot_dat_no) %>%
  mutate(strat = ifelse(strat == "obs",
                         "Observed",
                        ifelse(strat == "risk6",
                               "Cox 5-year SCD risk >6%",
                               "Cox 5-year SCD risk >4%"))) %>%
  rename(Probability = value,
         ICD = icd)

ggplot(plot_dat, aes(x = Probability, y = ..density.., group = ICD, col = ICD)) +
  facet_grid(strat ~ state, scales = "free") +
  geom_density(size = 1.5) +
  ggtitle("") +
  theme_bw()

ggsave(filename = "images/post_hist.png", dpi = 640, width = 10, height = 10)


#
library(ggdist)

plot_dat %>%
  ggplot(aes(x = Probability, fill = ICD, group = ICD)) +
  facet_grid(strat ~ state, scales = "free") +
  stat_slab(alpha = 0.3) +
  stat_pointinterval(position = position_dodge(width = .4, preserve = "single")) +
  theme_bw() +
  ylab("Density") +
  theme(text = element_text(size = 20),
        legend.position="top",
        legend.title =element_blank()) +
  scale_fill_discrete(labels = c('No ICD', 'ICD'))

ggsave(filename = "images/post_hist_ggdist.png", dpi = 640, width = 12, height = 12)

