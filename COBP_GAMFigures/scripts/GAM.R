#///////////////////////////
# cript to generate a figure showing Colorado butterfly plant subpopulation 
# change over time at FE Warren Airforce Base, with trendlines from GAM analysis
# Alice Stears
# 18 May 2023
#///////////////////////////

# Load packages -----------------------------------------------------------
library(tidyverse)
library(readxl)


# Verbiage from previous report: We fit separate GAMs for the annual counts at
# each creek, creek segment, and the full population as smoothed functions of
# year.  All models used a log link function and a negative-binomial error
# structure to account for over-dispersion in the count data.  We selected the
# smoothing parameter (k) by refitting models with increasing values of k until
# the effective degrees of freedom (edf) remained stable (Wood et al. 2011), and
# assessed goodness of fit and serial auto-correlation by examining residual
# plots.  To identify periods of change within each curvilinear trend, we
# estimated the first derivative of the trend line with 95% point-wise
# confidence intervals (CI) and interpreted periods when CI did not overlap zero
# as evidence that the slope was significantly different from zero (Curtis and
# Simpson 2014).  We fit GAMs using restricted maximum likelihood functions with
# the mgcv package (version 1.8-28; Wood 2011) in the statistical software
# program R (version 3.6.1; R Development Core Team 2014).  While our analyses
# were conducted piece-wise, we suggest future efforts consider using
# hierarchical GAMs to explore non-linear variation at the nested spatial scales
# of creeks and segments in a single model of the full dataset.


# Load data ---------------------------------------------------------------
## load creek-level data
dat_creeks <- read_excel("./COBP_GAMFigures/data/COBP_wafb_2022.xlsx", sheet =1) %>% 
  rename(`Crow Creek` = `Crow Cr`, `Diamond Creek` = `Diamond Cr`, `Unnamed Creek` = `Unnamed Cr`, Total =  `WAFB (Total)`, Notes =  ...6 ) %>% 
  filter(Year != "av") %>% 
  mutate(CrowCreek = as.integer(`Crow Creek`),
         DiamondCreek = as.integer(`Diamond Creek`), 
         UnnamedCreek = as.integer(`Unnamed Creek`), 
         Total = as.integer(Total)) 
# replace asterisks by 2007 and 2008 (not sure what they mean?)
dat_creeks[dat_creeks$Year %in% c("2007*", "2008*"),"Year"] <- c("2007", "2008")

dat_creeks <- dat_creeks %>% 
  mutate(
    Year = as.integer(Year)
  ) %>% 
  pivot_longer(
    cols = c(`Crow Creek`, `Diamond Creek`, `Unnamed Creek`, Total), 
    names_to = "Creek",
    values_to = "popSize"
  ) %>% 
  mutate(Creek = as.factor(Creek))

dat_creeks <- dat_creeks[dat_creeks$Creek != "Total",]

# remove NAs 
dat_creeks <- dat_creeks[!is.na(dat_creeks$popSize),]

## load creek segment-level data
dat_seg <- read_excel("./COBP_GAMFigures/data/COBP_wafb_2022.xlsx", sheet =2)
names(dat_seg) <- c("Segment", paste(1989:2022)) 
dat_seg <- dat_seg %>% 
  pivot_longer(cols = paste(1989:2022), 
               names_to = "Year", 
               values_to = "popSize") %>% 
  mutate(Year = as.integer(Year), 
         popSize = as.integer(popSize))
dat_seg <- dat_seg[is.na(dat_seg$Segment)==FALSE,]

dat_seg <- dat_seg[!(dat_seg$Segment %in% c("Crow", "Diamond", "Unnamed", "Total")),]

# add a column for creek type 
dat_seg$Creek <- NA
dat_seg[dat_seg$Segment %in% c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"), "Creek"] <- "Crow Creek"
dat_seg[dat_seg$Segment %in% c("D1", "D2", "D3", "D4", "D5"), "Creek"] <- "Diamond Creek"
dat_seg[dat_seg$Segment %in% c("U1", "U2"), "Creek"] <- "Unnamed Creek"

# remove the year(s) w/ NA
dat_seg <- dat_seg[!is.na(dat_seg$popSize),]

dat_seg$Creek <- as.factor(dat_seg$Creek)
dat_seg$Segment <- as.factor(dat_seg$Segment)

# Fit hierarchical GAM --------------------------------------------------------------
# based on advice from Pedersen et. al, 2019
## start by fitting all types of models w/ global vs. non-global smoothers as well as same and different amounts of wigliness

## model w/ global smoother and same wigliness (model "G" from the paper)
# with a random effect for creek 
mod_G <- gam(popSize ~ s(Year, k = 10, bs = "tp") + 
               s(Creek, bs = "re") ,
             data = dat_creeks, method = "REML", family = "nb")

# check model
gam.check(mod_G)
gratia::draw(mod_G)

# setup prediction data
mod_G_pred <- with(dat_seg,
                   expand.grid(Year = unique(dat_seg$Year),
                               Creek = levels(dat_seg$Creek)
                   ))

#  make the prediction, add this and a column of standard errors to the prediction
# data.frame. Predictions are on the log scale.
mod_G_pred <- cbind(mod_G_pred,
                    predict(mod_G, 
                            mod_G_pred, 
                            se.fit=TRUE, 
                            type="response"))

# make the plot. Note here the use of the exp() function to back-transform the
# predictions (which are for log-uptake) to the original scale
plot_G <- ggplot() +
  facet_wrap(~Creek) +
  geom_ribbon(aes(ymin=(fit - 2*se.fit), ymax=(fit + 2*se.fit), x=Year),
              data=mod_G_pred,
              alpha=0.3,
              inherit.aes=FALSE) +
  geom_line(aes(x = Year, y=fit), data=mod_G_pred) +
  geom_point(data=dat_creeks, aes(x=Year, y=popSize)) +
  labs(x="Year",
       y="Population Size")


## A single common smoother plus group-level smoothers that have the same wiggliness (model "GS" from the paper)
mod_GS <- gam(popSize ~ s(Year, k = 10, m = 2, bs = "tp") + 
                s(Year, Creek, k = 5, bs = "fs", m = 2),
             data = dat_creeks, method = "REML", family = "nb")

gam.check(mod_GS)
gratia::draw(mod_GS)


mod_GS_pred <- predict(mod_GS, se.fit=TRUE, type = "response")

mod_GS_pred <- transform(dat_creeks, 
                 modGS = mod_GS_pred$fit, 
                 modGS_se = mod_GS_pred$se.fit)

plot_GS <- ggplot() +
  facet_wrap(~Creek) +
  geom_ribbon(aes(ymin=(modGS - 2*modGS_se), ymax=(modGS + 2*modGS_se), x=Year),
              data=mod_GS_pred,
              alpha=0.3,
              inherit.aes=FALSE) +
  geom_line(aes(x = Year, y=modGS), data=mod_GS_pred) +
  geom_point(data=dat_creeks, aes(x=Year, y=popSize)) +
  labs(x="Year",
       y="Population Size")

## A single common smoother plus group-level smoothers with differing wiggliness (Model GI in the paper)
mod_GI <- gam(popSize ~ s(Year, k = 10, m = 2) + 
                s(Year, by = Creek, k = 5, m = 1, bs = "tp") +
                s(Creek, bs = "re", k = 5),
              data = dat_creeks, method = "REML", family = "nb")

gam.check(mod_GI)
gratia::draw(mod_GI)


mod_GI_pred <- predict(mod_GI, se.fit=TRUE, type = "response")

mod_GI_pred <- transform(dat_creeks, 
                         modGI = mod_GI_pred$fit, 
                         modGI_se = mod_GI_pred$se.fit)

plot_GI <- ggplot() +
  facet_wrap(~Creek) +
  geom_ribbon(aes(ymin=(modGI - 2*modGI_se), ymax=(modGI + 2*modGI_se), x=Year),
              data=mod_GI_pred,
              alpha=0.3,
              inherit.aes=FALSE) +
  geom_line(aes(x = Year, y=modGI), data=mod_GI_pred) +
  geom_point(data=dat_creeks, aes(x=Year, y=popSize)) +
  labs(x="Year",
       y="Population Size")

## Model S (shared smoothers) is model GS without the global smoother term; 
mod_S <- gam(popSize ~  
                s(Year, Creek, k = 10, bs = "fs", m = 2),
              data = dat_creeks, method = "REML", family = "nb")

gam.check(mod_S)
gratia::draw(mod_S)


mod_S_pred <- predict(mod_S, se.fit=TRUE, type = "response")

mod_S_pred <- transform(dat_creeks, 
                         modS = mod_S_pred$fit, 
                         modS_se = mod_S_pred$se.fit)

plot_S <- ggplot() +
  facet_wrap(~Creek) +
  geom_ribbon(aes(ymin=(modS - 2*modS_se), ymax=(modS + 2*modS_se), x=Year),
              data=mod_S_pred,
              alpha=0.3,
              inherit.aes=FALSE) +
  geom_line(aes(x = Year, y=modS), data=mod_S_pred) +
  geom_point(data=dat_creeks, aes(x=Year, y=popSize)) +
  labs(x="Year",
       y="Population Size")

## Model I is model GI without the first term
mod_I <- gam(popSize ~  
                s(Year, by = Creek, k = 10, m = 1, bs = "tp") +
                s(Creek, bs = "re", k = 5),
              data = dat_creeks, method = "REML", family = "nb")

gam.check(mod_I)
gratia::draw(mod_I)


mod_I_pred <- predict(mod_I, se.fit=TRUE, type = "response")

mod_I_pred <- transform(dat_creeks, 
                        modI = mod_I_pred$fit, 
                        modI_se = mod_I_pred$se.fit)

plot_I <- ggplot() +
  facet_wrap(~Creek) +
  geom_ribbon(aes(ymin=(modI - 2*modI_se), ymax=(modI + 2*modI_se), x=Year),
              data=mod_I_pred,
              alpha=0.3,
              inherit.aes=FALSE) +
  geom_line(aes(x = Year, y=modI), data=mod_I_pred) +
  geom_point(data=dat_creeks, aes(x=Year, y=popSize)) +
  labs(x="Year",
       y="Population Size")

## compare model types
cowplot::plot_grid(plot_G, plot_GI, plot_GS, plot_S, plot_I)

AIC(mod_G, mod_S, mod_I, mod_GS, mod_GI)
## according to AIC, model GS is the best


# Make models by creek ----------------------------------------------------
## Crow Creek
dat_Crow <- dat_creeks[dat_creeks$Creek == "Crow Creek",]

# make model with gamm() and calculate confidence intervals
mod_Crow <- gam(popSize ~  
                   s(Year, k = 15),
                 data = dat_Crow, method = "REML", family = "poisson")
# check for overdispersion
sum(residuals(mod_Crow, type = "pearson")^2) / df.residual(mod_Crow)
# try again with negative binomial model
mod_Crow <- gam(popSize ~  
                   s(Year, k = 15),
                 data = dat_Crow, method = "REML", family = "nb")
# check for overdispersion
sum(residuals(mod_Crow, type = "pearson")^2) / df.residual(mod_Crow)
# now underdispersed, but better? 


# waaay overdispersed, try a negative binomial 
want <- seq(1, nrow(dat_Crow), length.out = 36)

pdat_Crow <- with(dat_Crow,
             data.frame(Year = Year[want]))

p2 <- predict(mod_Crow, newdata = pdat_Crow, type = "response", se.fit = TRUE)
pdat_Crow <- transform(pdat_Crow, p2 = p2$fit, se2 = p2$se.fit)

df.res <- df.residual(mod_Crow)
crit.t <- qt(0.025, df.res, lower.tail = FALSE)
pdat_Crow <- transform(pdat_Crow,
                  upper = p2 + (crit.t * se2),
                  lower = p2 - (crit.t * se2))


# get functions from github
# tmpf <- tempfile()
# download.file("https://gist.github.com/gavinsimpson/e73f011fdaaab4bb5a30/raw/82118ee30c9ef1254795d2ec6d356a664cc138ab/Deriv.R",
#               tmpf, method = "wget")
# source(tmpf)
# ls()

# estimate second derivatives
Term <- "Year"
m2.d <- Deriv(mod_Crow, n = 36)
m2.dci <- confint(m2.d, term = Term)
m2.dsig <- signifD(pdat_Crow$p2, d = m2.d[[Term]]$deriv,
                     +                    m2.dci[[Term]]$upper, m2.dci[[Term]]$lower)
# add significant trend data to the pdat_Crow data.frame
pdat_Crow$incr_sig <- unlist(m2.dsig$incr)
pdat_Crow$decr_sig <- unlist(m2.dsig$decr)

## Unnamed Creek
  dat_Unnamed <- dat_creeks[dat_creeks$Creek == "Unnamed Creek",]
  
  # make model with gamm() and calculate confidence intervals
  mod_Unnamed <- gam(popSize ~  
                     s(Year, k = 15),
                   data = dat_Unnamed, method = "REML", family = "poisson")
  # check for overdispersion
  sum(residuals(mod_Unnamed, type = "pearson")^2) / df.residual(mod_Unnamed)
  # try again with negative binomial model
  mod_Unnamed <- gam(popSize ~  
                     s(Year, k = 15),
                   data = dat_Unnamed, method = "REML", family = "nb")
  # check for overdispersion
  sum(residuals(mod_Unnamed, type = "pearson")^2) / df.residual(mod_Unnamed)
  # now underdispersed, but better? 
  
  
  want <- seq(1, nrow(dat_Unnamed), length.out = 36)
  
  pdat_Unnamed <- with(dat_Unnamed,
               data.frame(Year = Year[want]))
  
  p2 <- predict(mod_Unnamed, newdata = pdat_Unnamed, type = "response", se.fit = TRUE)
  pdat_Unnamed <- transform(pdat_Unnamed, p2 = p2$fit, se2 = p2$se.fit)
  
  df.res <- df.residual(mod_Unnamed)
  crit.t <- qt(0.025, df.res, lower.tail = FALSE)
  pdat_Unnamed <- transform(pdat_Unnamed,
                    upper = p2 + (crit.t * se2),
                    lower = p2 - (crit.t * se2))

  # estimate second derivatives
  Term <- "Year"
  m2.d <- Deriv(mod_Unnamed, n = 36)
  m2.dci <- confint(m2.d, term = Term)
  m2.dsig <- signifD(pdat_Unnamed$p2, d = m2.d[[Term]]$deriv,
                     +                    m2.dci[[Term]]$upper, m2.dci[[Term]]$lower)
  # add significant trend data to the pdat_Unnamed data.frame
  pdat_Unnamed$incr_sig <- unlist(m2.dsig$incr)
  pdat_Unnamed$decr_sig <- unlist(m2.dsig$decr)
  
## Diamond Creek
  dat_Diamond <- dat_creeks[dat_creeks$Creek == "Diamond Creek",]
  
  # make model with gamm() and calculate confidence intervals
  mod_Diamond <- gam(popSize ~  
                     s(Year, k = 15),
                   data = dat_Diamond, method = "REML", family = "poisson")
  # check for overdispersion
  sum(residuals(mod_Diamond, type = "pearson")^2) / df.residual(mod_Diamond)
  # try again with negative binomial model
  mod_Diamond <- gam(popSize ~  
                     s(Year, k = 15),
                   data = dat_Diamond, method = "REML", family = "nb",
                   niterPQL = 50)
  # check for overdispersion
  sum(residuals(mod_Diamond, type = "pearson")^2) / df.residual(mod_Diamond)
  # now underdispersed, but better? 
  
  # waaay overdispersed, try a negative binomial 
  want <- seq(1, nrow(dat_Diamond), length.out = 36)
  
  pdat_Diamond <- with(dat_Diamond,
               data.frame(Year = Year[want]))
  
  p2 <- predict(mod_Diamond, newdata = pdat_Diamond, type = "response", se.fit = TRUE)
  pdat_Diamond <- transform(pdat_Diamond, p2 = p2$fit, se2 = p2$se.fit)
  
  df.res <- df.residual(mod_Diamond)
  crit.t <- qt(0.025, df.res, lower.tail = FALSE)
  pdat_Diamond <- transform(pdat_Diamond,
                    upper = p2 + (crit.t * se2),
                    lower = p2 - (crit.t * se2))

  # estimate second derivatives
  Term <- "Year"
  m2.d <- Deriv(mod_Diamond, n = 36)
  m2.dci <- confint(m2.d, term = Term)
  m2.dsig <- signifD(pdat_Diamond$p2, d = m2.d[[Term]]$deriv,
                     +                    m2.dci[[Term]]$upper, m2.dci[[Term]]$lower)
  # add significant trend data to the pdat_Diamond data.frame
  pdat_Diamond$incr_sig <- unlist(m2.dsig$incr)
  pdat_Diamond$decr_sig <- unlist(m2.dsig$decr)

## put all of the data into on DF
  # add creek name
  pdat_Crow$Creek <-  "Crow Creek"
  pdat_Diamond$Creek <-  "Diamond Creek"
  pdat_Unnamed$Creek <- "Unnamed Creek"
  # add % deviance explained
  pdat_Crow$pDev <- summary(mod_Crow)$dev.expl * 100
  pdat_Diamond$pDev <- summary(mod_Diamond)$dev.expl * 100
  pdat_Unnamed$pDev <- summary(mod_Unnamed)$dev.expl * 100

predDat <- rbind(pdat_Crow, pdat_Diamond, pdat_Unnamed)
predDat$Creek <- as.factor(predDat$Creek)

# make a plot with the significant increases or decreases
byCreek_figure <- ggplot() +
    geom_ribbon(aes(ymin=lower, ymax=upper, x=Year),
                data=predDat,
                alpha=0.3,
                inherit.aes=FALSE) +
    geom_point(data=dat_creeks, aes(x=Year, y=popSize)) +
    labs(x="Year",
         y="Number of Individuals") + 
    geom_line(aes(x = Year, y = incr_sig), data = predDat, col = "royalblue", lwd = 1.5, alpha = .9) + 
    geom_line(aes(x = Year, y = decr_sig), data = predDat, col = "tomato", lwd = 1.5) +
    geom_line(aes(x = Year, y=p2), data=predDat) + 
    facet_wrap(.~ Creek, ncol = 1) +#, scales = "free_y") +
    geom_text(aes(x = 2021, y = 7500, label = paste0("%d = ",round(pDev,1))), data = predDat, size = 3) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgrey"),
          strip.text = element_text(size = 10, face = "bold"))
  #geom_smooth(aes(x = Year, y = popSize), data = dat_Diamond, method = "gam", method.args = list(family = "nb"))
  
# save to file
  ggsave(filename = "byCreek_GAM_figure.pdf", plot = byCreek_figure, device = "pdf", path = "./COBP_GAMFigures/")


# Whole-Population Model --------------------------------------------------
dat_all <- dat_creeks %>%  
    dplyr::select(c(Year, Creek, popSize)) %>% 
    group_by(Year) %>% 
    summarize(popSize = sum(popSize))
    

# make model with gamm() and calculate confidence intervals
mod_all <- gam(popSize ~  
                   s(Year, k = 15),
                 data = dat_all, method = "REML", family = "poisson")
# check for overdispersion
sum(residuals(mod_all, type = "pearson")^2) / df.residual(mod_all)
# is very overdispersed, so use nb 
mod_all <- gam(popSize ~  
                 s(Year, k = 15),
               data = dat_all, method = "REML", family = "nb")
# check for overdispersion
sum(residuals(mod_all, type = "pearson")^2) / df.residual(mod_all)

# waaay overdispersed, try a negative binomial 
want <- seq(1, nrow(dat_all), length.out = 36)

pdat_all <- with(dat_all,
             data.frame(Year = Year[want]))

p2 <- predict(mod_all, newdata = pdat_all, type = "response", se.fit = TRUE)
pdat_all <- transform(pdat_all, p2 = p2$fit, se2 = p2$se.fit)

df.res <- df.residual(mod_all)
crit.t <- qt(0.025, df.res, lower.tail = FALSE)
pdat_all <- transform(pdat_all,
                  upper = p2 + (crit.t * se2),
                  lower = p2 - (crit.t * se2))

# estimate second derivatives
Term <- "Year"
m2.d <- Deriv(mod_all, n = 36)
m2.dci <- confint(m2.d, term = Term)
m2.dsig <- signifD(pdat_all$p2, d = m2.d[[Term]]$deriv,
                     +                    m2.dci[[Term]]$upper, m2.dci[[Term]]$lower)
# add significant trend data to the pdat_all data.frame
pdat_all$incr_sig <- unlist(m2.dsig$incr)
pdat_all$decr_sig <- unlist(m2.dsig$decr)
pdat_all$pDev <- summary(mod_all)$dev.expl * 100

allCreeks_figure <- ggplot() +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=Year),
              data=pdat_all,
              alpha=0.3,
              inherit.aes=FALSE) +
  geom_point(data=dat_all, aes(x=Year, y=popSize)) +
  labs(x="Year",
       y="Number of Individuals") + 
  geom_line(aes(x = Year, y = incr_sig), data = pdat_all, col = "royalblue", lwd = 1.5, alpha = .9) + 
  geom_line(aes(x = Year, y = decr_sig), data = pdat_all, col = "tomato", lwd = 1.5) +
  geom_line(aes(x = Year, y=p2), data=pdat_all) + 
  geom_text(aes(x = 2021, y = 15000, label = paste0("%d = ",round(pDev,1))), data = pdat_all, size = 3) +
  theme_minimal() +
  theme(strip.background = element_rect(fill = "lightgrey"),
        strip.text = element_text(size = 10, face = "bold"))
  
# save to file
ggsave(filename = "allCreeks_GAM_figure.pdf", plot = allCreeks_figure, 
       device = "pdf", path = "./COBP_GAMFigures/",
       height = 3, width = 6)
