#///////////////////////////
# Script to generate a figure showing Colorado butterfly plant subpopulation 
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
  rename(CrowCreek = `Crow Cr`, DiamondCreek = `Diamond Cr`, UnnamedCreek = `Unnamed Cr`, Total =  `WAFB (Total)`, Notes =  ...6 ) %>% 
  filter(Year != "av") %>% 
  mutate(CrowCreek = as.integer(CrowCreek),
         DiamondCreek = as.integer(DiamondCreek), 
         UnnamedCreek = as.integer(UnnamedCreek), 
         Total = as.integer(Total)) 
# replace asterisks by 2007 and 2008 (not sure what they mean?)
dat_creeks[dat_creeks$Year %in% c("2007*", "2008*"),"Year"] <- c("2007", "2008")

dat_creeks <- dat_creeks %>% 
  mutate(
    Year = as.integer(Year)
  ) %>% 
  pivot_longer(
    cols = c(CrowCreek, DiamondCreek, UnnamedCreek, Total), 
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
dat_seg[dat_seg$Segment %in% c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"), "Creek"] <- "CrowCreek"
dat_seg[dat_seg$Segment %in% c("D1", "D2", "D3", "D4", "D5"), "Creek"] <- "DiamondCreek"
dat_seg[dat_seg$Segment %in% c("U1", "U2"), "Creek"] <- "UnnamedCreek"

# remove the year(s) w/ NA
dat_seg <- dat_seg[!is.na(dat_seg$popSize),]

dat_seg$Creek <- as.factor(dat_seg$Creek)
dat_seg$Segment <- as.factor(dat_seg$Segment)

# Fit hierarchical GAM --------------------------------------------------------------
# based on advice from Pedersen et. al, 2019
## start by fitting all types of models w/ global vs. non-global smoothers as well as same and different amounts of wigliness

## model w/ global smoother and same wigliness (model "G" from the paper)
# with a random effect for creek 
mod_G <- gam(popSize ~ s(Year, k = 15, bs = "tp") + 
               s(Creek, bs = "re") ,
             data = dat_creeks, method = "REML", family = "poisson")

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
ggplot() +
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
mod_GS <- gam(popSize ~ s(Year, k = 15, m = 2) + 
                s(Year, Creek, k = 5, bs = "fs", m = 2),
             data = dat_creeks, method = "REML", family = "poisson")


mod_GS <- gam(log(uptake) ~ s(log(conc), k=5, m=2) + 
                   s(log(conc), Plant_uo, k=5,  bs="fs", m=2),
                 data=CO2, method="REML")
gam.check(mod_GS)
gratia::draw(mod_GS)


mod_GS_pred <- predict(mod_GS, se.fit=TRUE, type = "response")

mod_GS_pred <- transform(dat_creeks, 
                 modGS = mod_GS_pred$fit, 
                 modGS_se = mod_GS_pred$se.fit)

ggplot() +
  facet_wrap(~Creek) +
  geom_ribbon(aes(ymin=(modGS - 2*modGS_se), ymax=(modGS + 2*modGS_se), x=Year),
              data=mod_GS_pred,
              alpha=0.3,
              inherit.aes=FALSE) +
  geom_line(aes(x = Year, y=modGS), data=mod_GS_pred) +
  geom_point(data=dat_creeks, aes(x=Year, y=popSize)) +
  labs(x="Year",
       y="Population Size")

ggplot(data=CO2, aes(x=conc, y=uptake, group=Plant_uo)) +
  facet_wrap(~Plant_uo) +
  geom_ribbon(aes(ymin=exp(modGS-2*modGS_se),
                  ymax=exp(modGS+2*modGS_se)), alpha=0.25) +
  geom_line(aes(y=exp(modGS))) +
  geom_point() +
  labs(x=expression(CO[2] ~ concentration ~ (mL ~ L^{-1})),
       y=expression(CO[2] ~ uptake ~ (mu*mol ~ m^{-2})))
