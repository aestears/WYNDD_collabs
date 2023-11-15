# Different figures and data exploration here
# Note I have not updated counts beyond 2021

# Setup ----

library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(lme4)
library(viridis)
library(mgcv)
library(ggpubr)
library(fitdistrplus)

theme_set(theme_bw(base_size = 14))
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Input files ----
census <- read.csv("data/creek.counts.csv", header = TRUE) %>% 
  pivot_longer(cols = crow:unnamed, names_to = "creek", values_to = "count")
tots <- census %>% 
  group_by(year) %>% 
  summarise(creek = "total",
            count = sum(count, na.rm = FALSE))
census <- bind_rows(census, tots)
census$creek <- factor(census$creek, levels = c("total", "crow", "diamond", "unnamed"),
                         labels = c("FEWAFB Total", "Crow creek", "Diamond creek", "Unnamed creek"))

# Milestones of ESA timeline for plot annotations
milestones <- data.frame(year = c(2000, 2005, 2011, 2019),
                         event = c("2000\nListed as\nthreatened",
                                   "2005\nCritical habitat\ndesignated", "2011\nStatus\nreview",
                                   "2019\nRemoved\nfrom listing"),
                         count = 500, creek = "total")
milestones$creek <- factor(milestones$creek, levels = c("total", "crow", "diamond", "unnamed"),
                           labels = c("WAFB Total", "Crow creek", "Diamond creek", "Unnamed creek"))

# Annual growth rates: log(lambda) and counts by creek ----
census <- census %>% 
  filter(year > 1987) %>% 
  group_by(creek) %>% 
  arrange(year) %>% 
  mutate(lead1count = lead(count, n = 1L),
         gr = lead1count/count,
         loglam = log(lead1count/count))

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

meanlam <- census %>% 
  filter(year < 2021) %>% 
  summarise(arithlam = mean(loglam, na.rm = TRUE),
            geomlam = log(gm_mean(gr, na.rm = TRUE)))

plt1 <- ggplot(census, aes(x = year, y = count, group = creek)) +
  scale_color_manual(name = NULL, values = "black") +
  scale_y_continuous(limits=c(0,NA), expand = expansion(mult = c(0, .1))) +
  geom_point(alpha = .5) +
  geom_path() +
  facet_wrap(~creek, nrow = 1, scales = "free_y") +
  # geom_text(data = milestones, aes(label = event), inherit.aes = TRUE, size = 3) +
  ylab("Flowering plant count") +
  xlab("Year") +
  theme(axis.text.y = element_text(angle = 90, vjust = .5, hjust=.5))
plt1

plt2 <- ggplot(census, aes(x = year, y = loglam, group = creek)) +
  scale_color_manual(name = NULL, values = "black") +
  scale_y_continuous(limits=c(-2.3,2)) +
  geom_point(alpha = .5) +
  geom_path() +
  # geom_hline(yintercept = 0, linetype = "dotted") +
  facet_wrap(~creek, nrow = 1, scales = "free_y") +
  geom_hline(data = meanlam, aes(yintercept = arithlam, group = creek), linetype = "dashed") +
  # geom_text(data = milestones, aes(label = event), inherit.aes = TRUE, size = 3) +
  ylab("log Lambda") +
  xlab("Year") +
  theme(axis.text.y = element_text(angle = 90, vjust = .5, hjust=.5))

plt2

plt3 <- ggarrange(plt1, plt2, nrow = 2, ncol = 1, heights = c(1, 1))
plt3
ggsave(plt3, filename = "plots/creekcounts_revised.png", device = "png", width = 7, height = 10, units = "in")


# Nonlinear trends with GAM ----

plt1 <- ggplot(census %>% filter(creek == "FEWAFB Total"), aes(x = year, y = count, group = creek, color = creek)) +
  geom_smooth(method = "gam", se = FALSE, formula = y ~ s(x, bs = "cs", k = 10), 
              method.args = list(family = poisson(link = "log")),  alpha = .5) +
  scale_color_manual(name = NULL, values = "black") +
  # geom_ribbon(aes(x = year, ymin = lwr, ymax = upr), alpha = .2) +
  geom_point(alpha = .5) +
  coord_cartesian(ylim = c(0, NA)) +
  # facet_wrap(~creek, ncol = 1, scales = "free_y") +
  # geom_text(data = milestones, aes(label = event), inherit.aes = TRUE, size = 3) +
  ylab("Flower census count") +
  xlab("Year") +
  theme(legend.position = c(.17, .83))

plt1

plt2 <- ggplot(census %>% filter(creek != "FEWAFB Total") %>% droplevels(), aes(x = year, y = count, group = creek, color = creek)) +
  geom_smooth(method = "gam", se = FALSE, formula = y ~ s(x, bs = "cs", k = 10), 
              method.args = list(family = poisson(link = "log")),  alpha = .5) +
  scale_color_brewer(name = NULL, palette = "Dark2") +
  # geom_ribbon(aes(x = year, ymin = lwr, ymax = upr), alpha = .2) +
  geom_point(alpha = .75) +
  coord_cartesian(ylim = c(0, NA)) +
  # facet_wrap(~creek, ncol = 1) +
  ylab("Flower census count") +
  xlab("Year") +
  theme(legend.position = c(.17, .83))
plt2

plt3 <- ggarrange(plt1, plt2, nrow = 2, ncol = 1, heights = c(1, 1), labels = c('A', 'B'))
plt3

ggsave(plt3, filename = "plots/creekcounts_revised.png", device = "png", width = 7, height = 10, units = "in")



# Rosette Ratios ----
# I tried and failed to make Bayesian population models that would estimate
# the total population of flowers + rosettes, allowing for annual variation 
# in the proportion of plants flowering. Something like a state-space model but
# with a binomial distribution for the flowering proportion.
# Below is my attempt to find a distribution for the flower:rosette ratio 
# for an informative prior.
# Data comes from Bonnie and from the Floyd & Ranker 1995 paper

dat <- read.csv("data/allrosetteratios.csv", header = TRUE) %>% 
  filter(segment != "average") %>% 
  # filter(year <= 1995) %>%
  mutate(id_ycs = paste(year, creek, segment, sep = "_"),
         id_yc = paste(year, creek, sep = "_"),
         id_cs = paste(creek, segment, sep = "_"),
         total = rosettes + flowers) %>% 
  filter(complete.cases(.)) %>% 
  group_by(collector, year, creek, segment, id_ycs, id_yc, id_cs) %>% 
  summarise(rosettes = sum(rosettes),
            flowers = sum(flowers),
            total = rosettes + flowers)

mod <- glmer(cbind(flowers, total) ~ 1 + (1|year) + (1|id_yc) + (1|id_cs) + (1|id_ycs),
             family = binomial, data = dat)
mod1 <- glmer(cbind(flowers, total) ~ 1 + (1|year) + (1|id_cs),
              family = binomial, data = dat)
mod2 <- glmer(cbind(flowers, total) ~ collector + (1|id_ycs),
              family = binomial, data = dat)
# creeks not significantly different, but Floyd and Ranker not representative
mod3 <- glmer(cbind(flowers, total) ~ collector + creek + (1|year) + (1|id_ycs),
              family = binomial, data = dat)
# mod1 <- glmer(cbind(flowers, total) ~ log(total + 1) + (1|year) + (1|creek) + (1|id_ycs),
#              family = binomial, data = dat)
AIC(mod, mod1, mod2, mod3)

dat$pred <- predict(mod, dat)

# Counts vs rosette ratios ----
flowprop <- coef(mod)$id_yc
propdf <- data.frame(id = row.names(flowprop), perc = boot::inv.logit(flowprop[,1])) %>% 
  mutate(ratio = (1-perc)/perc,
         year = as.numeric(str_split_fixed(id, "_", 2)[,1]),
         creek = tolower(str_split_fixed(id, "_", 2)[,2])) %>% 
  filter(creek != "soapstone")

propdf$creek <- factor(propdf$creek, levels = c("crow", "diamond", "unnamed"),
                       labels = c("Crow creek", "Diamond creek", "Unnamed creek"))

census <- left_join(census, propdf)


plt <- ggplot(census, aes(x = year, y = count, group = creek)) +
  geom_smooth(method = "gam", se = FALSE, formula = y ~ s(x, bs = "cs", k = 10), 
              method.args = list(family = poisson(link = "log")),  alpha = .5) +
  geom_point(aes(color = perc), size = 3) +
  scale_color_viridis() +
  coord_cartesian(ylim = c(0, NA)) +
  facet_wrap(~creek, ncol = 1, scales = "free") +
  ylab("Flower census count") +
  xlab("Year")
plt

# lognormal good enough, easy to use with JAGS code with logit link functions and normal priors

vals <- boot::inv.logit(dat$pred)
descdist(vals, boot = 1000)
fb <- fitdist(vals, "beta")
fg <- fitdist(vals, "gamma")
fln <- fitdist(vals, "lnorm")
denscomp(list(fb, fg, fln))
summary(fln)



