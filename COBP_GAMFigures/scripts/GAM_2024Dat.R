#///////////////////////////
# cript to generate a figure showing Colorado butterfly plant subpopulation 
# change over time at FE Warren Airforce Base, with trendlines from GAM analysis
# Alice Stears
# 10 October 2023
#///////////////////////////

# Load packages -----------------------------------------------------------
library(tidyverse)
library(readxl)
library(mgcv)

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
dat_creeks <- read_excel("~/Documents/Dropbox_static/Work/collab_papers/WYNDD_collabs/COBP_GAMFigures/data/COBP_wafb_2024data.xlsx", sheet =1) %>% 
  rename(`Crow Creek` = `Crow Cr`, `Diamond Creek` = `Diamond Cr`, `Unnamed Creek` = `Unnamed Cr`, Total =  `WAFB (Total)`, Notes =  ...6 ) %>% 
  filter(Year != "Average") %>% 
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
dat_seg <- read_excel("~/Documents/Dropbox_static/Work/collab_papers/WYNDD_collabs/COBP_GAMFigures/data/COBP_wafb_2024data.xlsx", sheet =2)
names(dat_seg) <- c("Segment", paste(1989:2024)) 
dat_seg <- dat_seg %>% 
  pivot_longer(cols = paste(1989:2024), 
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
                  s(Year, k = 5),
                data = dat_Crow, method = "REML", family = "poisson")
# check for overdispersion
sum(residuals(mod_Crow, type = "pearson")^2) / df.residual(mod_Crow)
# not overdispersed (a bit underdispersed!)

## test different k values (max is (n/2) = 18)
mod_Crow # k = 5, edf = 3.68 , REML = 272.5435
# increase k 
mod_Crow <- gam(popSize ~  s(Year, k = 6), data = dat_Crow, method = "REML", family = "nb")
mod_Crow # k = 6, edf = 3.03 , REML = 272.6241
mod_Crow <- gam(popSize ~  s(Year, k = 7), data = dat_Crow, method = "REML", family = "nb")
mod_Crow # k = 7, edf = 3.79, REML = 272.4834
mod_Crow <- gam(popSize ~  s(Year, k = 8), data = dat_Crow, method = "REML", family = "nb")
mod_Crow # k = 8, edf = 3.74, REML = 272.5202
mod_Crow <- gam(popSize ~  s(Year, k = 9), data = dat_Crow, method = "REML", family = "nb")
mod_Crow # k = 8, edf = 3.68, REML = 272.5435

# use k=7
mod_Crow <- gam(popSize ~  s(Year, k = 7), data = dat_Crow, method = "REML", family = "nb")

## calculate significant trends
want <- seq(1, nrow(dat_Crow), length.out = 37)

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
tmpf <- tempfile()
download.file("https://gist.github.com/gavinsimpson/e73f011fdaaab4bb5a30/raw/82118ee30c9ef1254795d2ec6d356a664cc138ab/Deriv.R",
              tmpf, method = "wget")
source(tmpf)
ls()

# estimate second derivatives
Term <- "Year"
m2.d <- Deriv(mod_Crow, n = 37)
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
                     s(Year, k = 5),
                   data = dat_Unnamed, method = "REML", family = "poisson")
# check for overdispersion
sum(residuals(mod_Unnamed, type = "pearson")^2) / df.residual(mod_Unnamed)
# try again with negative binomial model
mod_Unnamed <- gam(popSize ~  
                     s(Year, k = 5),
                   data = dat_Unnamed, method = "REML", family = "nb")
# check for overdispersion
sum(residuals(mod_Unnamed, type = "pearson")^2) / df.residual(mod_Unnamed)
# now underdispersed, but better? 


## test different k values (max is (n/2) = 18)
mod_Unnamed # k = 5, edf = 3.21 , REML = 309.56
# increase k 
mod_Unnamed <- gam(popSize ~  s(Year, k = 6), data = dat_Unnamed, method = "REML", family = "nb")
mod_Unnamed # k = 6, edf = 3.24 , REML = 309.56
mod_Unnamed <- gam(popSize ~  s(Year, k = 7), data = dat_Unnamed, method = "REML", family = "nb")
mod_Unnamed # k = 7, edf = 3.24, REML = 309.57
mod_Unnamed <- gam(popSize ~  s(Year, k = 8), data = dat_Unnamed, method = "REML", family = "nb")
mod_Unnamed # k = 8, edf = 3.25, REML = 309.58

# use k=6
mod_Unnamed <- gam(popSize ~  s(Year, k = 6), data = dat_Unnamed, method = "REML", family = "nb")

## calculate significant trends

want <- seq(1, nrow(dat_Unnamed), length.out = 37)

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
m2.d <- Deriv(mod_Unnamed, n = 37)
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
                     s(Year, k = 5),
                   data = dat_Diamond, method = "REML", family = "poisson")
# check for overdispersion
sum(residuals(mod_Diamond, type = "pearson")^2) / df.residual(mod_Diamond)
# try again with negative binomial model
mod_Diamond <- gam(popSize ~  
                     s(Year, k = 5),
                   data = dat_Diamond, method = "REML", family = "nb",
                   niterPQL = 50)
# check for overdispersion
sum(residuals(mod_Diamond, type = "pearson")^2) / df.residual(mod_Diamond)
# now underdispersed, but better? 


## test different k values (max is (n/2) = 18)
mod_Diamond # k = 5, edf = 2.14 , REML = 333.51
# increase k 
mod_Diamond <- gam(popSize ~  s(Year, k = 6), data = dat_Diamond, method = "REML", family = "nb")
mod_Diamond # k = 6, edf = 5.43 , REML = 332.32
mod_Diamond <- gam(popSize ~  s(Year, k = 7), data = dat_Diamond, method = "REML", family = "nb")
mod_Diamond # k = 7, edf = 5.67, REML = 332.61
mod_Diamond <- gam(popSize ~  s(Year, k = 8), data = dat_Diamond, method = "REML", family = "nb")
mod_Diamond # k = 8, edf = 5.96, REML = 332.54
mod_Diamond <- gam(popSize ~  s(Year, k = 9), data = dat_Diamond, method = "REML", family = "nb")
mod_Diamond # k = 9, edf = 6.03, REML = 3332.62
mod_Diamond <- gam(popSize ~  s(Year, k = 10), data = dat_Diamond, method = "REML", family = "nb")
mod_Diamond # k = 10, edf = 6.25, REML = 332.55
mod_Diamond <- gam(popSize ~  s(Year, k = 11), data = dat_Diamond, method = "REML", family = "nb")
mod_Diamond # k = 11, edf = 6.3, REML = 332.54
mod_Diamond <- gam(popSize ~  s(Year, k = 12), data = dat_Diamond, method = "REML", family = "nb")
mod_Diamond # k = 12, edf = 6.35, REML = 332.54
mod_Diamond <- gam(popSize ~  s(Year, k = 13), data = dat_Diamond, method = "REML", family = "nb")
mod_Diamond # k = 13, edf = 6.61, REML = 332.54
mod_Diamond <- gam(popSize ~  s(Year, k = 14), data = dat_Diamond, method = "REML", family = "nb")
mod_Diamond # k = 14, edf = 6.69, REML = 332.54
mod_Diamond <- gam(popSize ~  s(Year, k = 15), data = dat_Diamond, method = "REML", family = "nb")
mod_Diamond # k = 15, edf = 6.69, REML = 332.44

# use k=15
mod_Diamond <- gam(popSize ~  s(Year, k = 10), data = dat_Diamond, method = "REML", family = "nb")

## calculate significant trends

want <- seq(1, nrow(dat_Diamond), length.out = 37)

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
m2.d <- Deriv(mod_Diamond, n = 37)
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

predDat <- predDat %>% 
  left_join(dat_creeks[,c("Year", "Creek", "popSize")], by = c("Year", "Creek"))

temp <- predDat %>% 
  group_by(Creek) %>% 
  summarize(maxPop = max(popSize),
            maxInt = max(upper))

temp$maxY <- apply(temp, MARGIN = 1, FUN = function(x)
  as.integer(max(x[2:3]))
)

predDat <- predDat %>% 
  left_join(temp[,c("Creek", "maxY")], by = "Creek")

# make a plot with the significant increases or decreases
byCreek_figure <- ggplot() +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=Year),
              data=predDat,
              alpha=0.3,
              inherit.aes=FALSE) +
  geom_point(data=dat_creeks, aes(x=Year, y=popSize)) +
  labs(x="Year",
       y="Number of Individuals") + 
  geom_line(aes(x = Year, y = incr_sig), data = predDat, col = "tomato", lwd = 1.5, alpha = .9) + 
  geom_line(aes(x = Year, y = decr_sig), data = predDat, col = "royalblue", lwd = 1.5) +
  geom_line(aes(x = Year, y=p2), data=predDat) + 
  facet_wrap(.~ Creek, ncol = 1, scales = "free_y") +
  geom_text(aes(x = 2021, y = maxY*.95, label = paste0("%d = ",round(pDev,1))), data = predDat, size = 3) +
  theme_minimal() +
  ylim(c(0,NA)) +
  theme(strip.background = element_rect(fill = c("lightgrey"), colour = NA),
        strip.text = element_text(size = 10, face = "bold", hjust = -.005))
#geom_smooth(aes(x = Year, y = popSize), data = dat_Diamond, method = "gam", method.args = list(family = "nb"))

# save to file
ggsave(filename = "byCreek_GAM_figure_2024.pdf", plot = byCreek_figure, device = "pdf", path = "./COBP_GAMFigures/COBP_GAMFigures_2024")


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
                 s(Year, k = 10),
               data = dat_all, method = "REML", family = "nb")
# check for overdispersion
sum(residuals(mod_all, type = "pearson")^2) / df.residual(mod_all)

## test different k values (max is (n/2) = 18)
mod_all # k = 10, edf = 2 , REML = 348.10
# increase k 
mod_all <- gam(popSize ~  s(Year, k = 11), data = dat_all, method = "REML", family = "nb")
mod_all # k = 11, edf = 5.47 , REML = 348.03
mod_all <- gam(popSize ~  s(Year, k = 13), data = dat_all, method = "REML", family = "nb")
mod_all # k = 13, edf = 5.67, REML = 347.98
mod_all <- gam(popSize ~  s(Year, k = 15), data = dat_all, method = "REML", family = "nb")
mod_all # k = 15, edf =  5.71, REML = 3347.98
mod_all <- gam(popSize ~  s(Year, k = 16), data = dat_all, method = "REML", family = "nb")
mod_all # k = 8, edf = 5.71, REML = 347.99

# use k=15
mod_all <- gam(popSize ~  s(Year, k = 15), data = dat_all, method = "REML", family = "nb")


want <- seq(1, nrow(dat_all), length.out = 37)

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
m2.d <- Deriv(mod_all, n = 37)
m2.dci <- confint(m2.d, term = Term)
m2.dsig <- signifD(pdat_all$p2, d = m2.d[[Term]]$deriv,
                   +                    m2.dci[[Term]]$upper, m2.dci[[Term]]$lower)
# add significant trend data to the pdat_all data.frame
pdat_all$incr_sig <- unlist(m2.dsig$incr)
pdat_all$decr_sig <- unlist(m2.dsig$decr)
pdat_all$pDev <- summary(mod_all)$dev.expl * 100


# add information for "Creek"
pdat_all$Creek <- "All Creeks"
pdat_all <- left_join(pdat_all, dat_all, by = "Year")


pdat_all$maxY <- max(pdat_all$popSize)
pdat_all <- pdat_all[,names(predDat)]

allCreeks_figure <- ggplot() +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=Year),
              data=pdat_all,
              alpha=0.3,
              inherit.aes=FALSE) +
  geom_point(data=dat_all, aes(x=Year, y=popSize)) +
  labs(x="Year",
       y="Number of Individuals") + 
  geom_line(aes(x = Year, y = incr_sig), data = pdat_all, col = "tomato", lwd = 1.5, alpha = .9) + 
  geom_line(aes(x = Year, y = decr_sig), data = pdat_all, col = "royalblue", lwd = 1.5) +
  geom_line(aes(x = Year, y=p2), data=pdat_all) + 
  geom_text(aes(x = 2021, y = 15000, label = paste0("%d = ",round(pDev,1))), data = pdat_all, size = 3) +
  theme_minimal() +
  theme(strip.background = element_rect(fill = "lightgrey"),
        strip.text = element_text(size = 10, face = "bold"))

# save to file
ggsave(filename = "allCreeks_GAM_figure.pdf", plot = allCreeks_figure, 
       device = "pdf", path = "~/Documents/Dropbox_static/Work/collab_papers/WYNDD_collabs/COBP_GAMFigures/COBP_GAMFigures_2024/",
       height = 3, width = 6)


# figure for both all data and by creek -----------------------------------

#add data together
preds <- rbind(predDat, pdat_all)
# order the factor of "creek" so that "all creeks" is on the top
preds$Creek <- factor(preds$Creek, levels = c("All Creeks", "Crow Creek", "Diamond Creek", "Unnamed Creek"), ordered = TRUE)

# make figure for creeks
# make labelling function
creek_label = c("All Creeks" = "A) All Creeks",
                "Crow Creek" = "B) Crow Creek", 
                "Diamond Creek" = "C) Diamond Creek", 
                "Unnamed Creek" = "D) Unnamed Creek")

byCreeks_partFig <- ggplot(data = preds[preds$Creek != "All Creeks",]) +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=Year),
              alpha=0.3,
              inherit.aes=FALSE) +
  geom_point( aes(x=Year, y=popSize)) +
  labs(x="Year",
       y="Number of Individuals") + 
  geom_line(aes(x = Year, y = incr_sig), col = "tomato", lwd = 1.5, alpha = .9) + 
  geom_line(aes(x = Year, y = decr_sig), col = "royalblue", lwd = 1.5) +
  geom_line(aes(x = Year, y=p2)) + 
  facet_wrap(.~ Creek, ncol = 1, scales = "free_y", labeller = labeller(Creek = creek_label)) +
  geom_text(aes(x = 2021, y = maxY*.95, label = paste0("%d = ",round(pDev,1))), size = 3) +
  theme_minimal() +
  ylim(c(0,NA)) +
  xlim(c(1986, 2022)) + 
  theme(strip.background = element_rect(fill = c("lightgrey"), colour = NA),
        strip.text = element_text(size = 10, face = "bold", hjust = -.005))
#geom_smooth(aes(x = Year, y = popSize), data = dat_Diamond, method = "gam", method.args = list(family = "nb"))

allCreeks_partFig <- ggplot(data = preds[preds$Creek == "All Creeks",]) +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=Year),
              alpha=0.3,
              inherit.aes=FALSE) +
  geom_point( aes(x=Year, y=popSize)) +
  labs(x="Year",
       y="Number of Individuals") + 
  geom_line(aes(x = Year, y = incr_sig), col = "tomato", lwd = 1.5, alpha = .9) + 
  geom_line(aes(x = Year, y = decr_sig), col = "royalblue", lwd = 1.5) +
  geom_line(aes(x = Year, y=p2)) + 
  facet_wrap(.~ Creek, ncol = 1, scales = "free_y", labeller = labeller(Creek = creek_label)) +
  geom_text(aes(x = 2021, y = maxY*.95, label = paste0("%d = ",round(pDev,1))), size = 3) +
  theme_minimal() +
  ylim(c(0,NA)) +
  xlim(c(1986, 2022)) + 
  theme(strip.background = element_rect(fill = c("lightgrey"), colour = NA),
        strip.text = element_text(size = 10, face = "bold", hjust = -.005))


allCreekByCreek_Fig <- cowplot::plot_grid(allCreeks_partFig, byCreeks_partFig, ncol = 1, rel_heights = c(.3,.6))

# save figure
ggsave(filename = "allDat_byCreeks_GAM_figure.pdf", plot =  allCreekByCreek_Fig, 
       device = "pdf", path = "~/Documents/Dropbox_static/Work/collab_papers/WYNDD_collabs/COBP_GAMFigures/COBP_GAMFigures_2024/",
       height = 7, width = 6)

# Segment-level figures ---------------------------------------------------

for (i in 1:length(unique(dat_seg$Segment))) {
  # get the name of the ith segment
  seg_i <- unique(dat_seg$Segment)[i]
  # get data for that model
  dat_i <- dat_seg[dat_seg$Segment ==  seg_i,]
  
  ## chose the appropriate k 
  # make a vector of possible k values
  k_poss <- as.array((3:(nrow(dat_i)/2)))
  # fit the gam w/ each k
  modResults <- apply(X = k_poss, MARGIN = 1, FUN = function(x) {
    # run the model w/ the current value of "k"
    modTemp <- gam(popSize ~ s(Year, k = x),
                   data = dat_i, method = "REML", family = "nb")
    # save the relevant values from that run
    c("K" = x, "EDF" = sum(modTemp$edf)-1, "k_index" = k.check(modTemp)[3])
  }
  )
  
  # get the k for which the EDF stabilizes (the most common number in the vector of EDFs?)
  commonValues <- as.data.frame(table(round(modResults[2,],1)))
  if (nrow(commonValues)==1) {
    commonValues <- commonValues[order(commonValues$Freq, decreasing = TRUE),][1,]
  } else {
    commonValues <- commonValues[order(commonValues$Freq, decreasing = TRUE),][1:2,] 
  }
  # get the k's with the "most common" EDFs 
  bestEDF <- as.data.frame(modResults[,which(round(modResults[2,],1) %in% as.numeric(as.character(commonValues$Var1)))])
  # get the k with the smallest k index
  bestK <- bestEDF[1,which.max(bestEDF[3,])]
  
  # make the model
  mod_i <- gam(popSize ~  
                 s(Year, k = bestK),
               data = dat_i, method = "REML", family = "nb")
  
  # get the predictions and ses
  want <- seq(1, nrow(dat_i), length.out = nrow(dat_i))
  
  pdat_i <- with(dat_i,
                 data.frame(Year = Year[want]))
  
  p2 <- predict(mod_i, newdata = pdat_i, type = "response", se.fit = TRUE)
  pdat_i <- transform(pdat_i, p2 = p2$fit, se2 = p2$se.fit)
  
  df.res <- df.residual(mod_i)
  crit.t <- qt(0.025, df.res, lower.tail = FALSE)
  pdat_i <- transform(pdat_i,
                      upper = p2 + (crit.t * se2),
                      lower = p2 - (crit.t * se2))
  
  # estimate second derivatives
  Term <- "Year"
  m2.d <- Deriv(mod_i, n = nrow(pdat_i))
  m2.dci <- confint(m2.d, term = Term)
  m2.dsig <- signifD(pdat_i$p2, d = m2.d[[Term]]$deriv,
                     +                    m2.dci[[Term]]$upper, m2.dci[[Term]]$lower)
  
  # add significant trend data to the pdat_i data.frame
  pdat_i$incr_sig <- unlist(m2.dsig$incr)
  pdat_i$decr_sig <- unlist(m2.dsig$decr)
  pdat_i$pDev <- summary(mod_i)$dev.expl * 100
  pdat_i$Segment <- seg_i
  pdat_i$k_i <- bestK
  # save model terms
  pdat_i$edf <- summary(mod_i)$edf
  pdat_i$pVal <- summary(mod_i)$s.table[4]
  
  ## save the data
  if (i == 1) {
    pdat_allSegs <- pdat_i
  } else {
    pdat_allSegs <- rbind(pdat_allSegs, pdat_i)
  }
}

## add actual data to the same data.frame
pdat_allSegsNew <- pdat_allSegs %>% 
  left_join(dat_seg, by = c("Year", "Segment"))

temp <- pdat_allSegsNew %>% 
  group_by(Segment) %>% 
  summarize(maxPop = max(popSize),
            maxInt = max(upper))

temp$maxY <- apply(temp, MARGIN = 1, FUN = function(x)
  as.integer(max(x[2:3]))
)

pdat_allSegs <- pdat_allSegs %>% 
  left_join(temp[,c("Segment", "maxY")], by = "Segment")

# make lower bound of confidence intervals 0
pdat_allSegs[pdat_allSegs$lower <0,"lower"] <- 0
## put as a figure
# make a plot with the significant increases or decreases
bySegment_figure <- 
  ggplot() +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=Year),
              data=pdat_allSegs,
              alpha=0.3,
              inherit.aes=FALSE) +
  geom_point(data=dat_seg, aes(x=Year, y=popSize)) +
  labs(x="Year",
       y="Number of Individuals") + 
  geom_line(aes(x = Year, y = incr_sig), data = pdat_allSegs, col = "tomato", lwd = 1.5, alpha = .9) + 
  geom_line(aes(x = Year, y = decr_sig), data = pdat_allSegs, col = "royalblue", lwd = 1.5) +
  geom_line(aes(x = Year, y=p2), data=pdat_allSegs) + 
  facet_wrap(.~ Segment, scales = "free_y") +
  geom_text(aes(x = 2018, y = maxY*.95, label = paste0("%d = ",round(pDev,0))), 
            data = pdat_allSegs, size = 3) +
  theme_minimal() +
  theme(strip.background = element_rect(fill = c("lightgrey"), colour = NA),
        strip.text = element_text(size = 10, face = "bold", hjust = -.005)) +
  ylim(c(0,NA))

# save to file
ggsave(filename = "bySegment_GAM_figure.pdf", plot = bySegment_figure, 
       device = "pdf", path = "~/Documents/Dropbox_static/Work/collab_papers/WYNDD_collabs/COBP_GAMFigures/COBP_GAMFigures_2024/",
       height = 8, width = 8)


# Log-lambda analysis -----------------------------------------------------
# generate plots showing change in population growth rate at each creek
# data: population = dat_all; creek = dat_creek
# for the entire population

lambda_all <- dat_all %>% 
  mutate(popSize_tMinus1 = as.numeric(dplyr::lag(dat_all$popSize,n = 1))) %>% 
  rowwise() %>% 
  mutate(
         lambda = popSize/popSize_tMinus1,
         logLam = log(lambda)) %>% 
  mutate(Creek = "All Creeks")

# by creek
temp <- dat_creeks %>% 
  select(Year, CrowCreek, DiamondCreek, UnnamedCreek) %>% 
  group_by(Year) %>% 
  summarize(CrowCreek_tplus1 = mean(CrowCreek),
            DiamondCreek_tplus1 = mean(DiamondCreek),
            UnnamedCreek_tplus1 = mean(UnnamedCreek)) 

lambda_creek <- temp %>% 
  lag() %>% 
  bind_cols(temp) %>%
  dplyr::select(2:8) %>% 
  rename("CrowCreek_tplus1" = "CrowCreek_tplus1...2", 
         "DiamondCreek_tplus1" = "DiamondCreek_tplus1...3",
         "UnnamedCreek_tplus1" = "UnnamedCreek_tplus1...4", 
         "Year" = "Year...5",
         "CrowCreek_t" = "CrowCreek_tplus1...6",
         "DiamondCreek_t" = "DiamondCreek_tplus1...7",
         "UnnamedCreek_t" = "UnnamedCreek_tplus1...8") %>% 
  mutate(CrowCreek_lambda = CrowCreek_tplus1/CrowCreek_t,
         DiamondCreek_lambda = DiamondCreek_tplus1/DiamondCreek_t,
         UnnamedCreek_lambda = UnnamedCreek_tplus1/UnnamedCreek_t,
         CrowCreek_logLam = log(CrowCreek_lambda),
         DiamondCreek_logLam = log(DiamondCreek_lambda),
         UnnamedCreek_logLam = log(UnnamedCreek_lambda)) %>% 
  pivot_longer(cols = c(1:3,5:13),
               names_to = c("Creek", "parameter"),
               values_to = "values",
               names_sep = "_")

lambda_creek <- lambda_creek %>% 
  mutate(Creek = str_glue("{name} Creek",
                          name = str_split(lambda_creek$Creek, "Creek", simplify = TRUE)[,1]))

creek_label = c("All Creeks" = "A) All Creeks",
                "Crow Creek" = "B) Crow Creek", 
                "Diamond Creek" = "C) Diamond Creek", 
                "Unnamed Creek" = "D) Unnamed Creek")
# plot the values over time for each creek 
creekLam_plot <- ggplot() +
  geom_hline(aes(), yintercept = 0, col = "grey") +
  geom_line(aes(x = Year, y = values), data = lambda_creek[lambda_creek$parameter %in% c("logLam"),]) +
  facet_wrap(.~Creek, ncol = 1, labeller = labeller(Creek = creek_label)) +
  theme_minimal() +
  theme(strip.background = element_rect(fill = c("lightgrey"), colour = NA),
        strip.text = element_text(size = 10, face = "bold", hjust = -.005)) +
  ylab("log(Lambda)") +
  xlim(c(1988,2022))

allLam_plot <- 
  ggplot() +
  geom_hline(aes(), yintercept = 0, col = "grey") +
  geom_line(aes(x = Year, y = logLam), data = lambda_all) +
  facet_wrap(.~Creek, ncol = 1, labeller = labeller(Creek = creek_label)) +
  theme_minimal() +
  theme(strip.background = element_rect(fill = c("lightgrey"), colour = NA),
        strip.text = element_text(size = 10, face = "bold", hjust = -.005)) +
  ylab("log(Lambda)") +
  xlim(c(1988,2022)) +
  ylim(c(-1.5,2))

lambdaPlots <- cowplot::plot_grid(allLam_plot, creekLam_plot, ncol = 1, rel_heights = c(.25,.75))
# save to file
ggsave(filename = "logLambda_figure.pdf", plot = lambdaPlots, 
       device = "pdf", path = "~/Documents/Dropbox_static/Work/collab_papers/WYNDD_collabs/COBP_GAMFigures/COBP_GAMFigures_2024/",
       height =6, width = 5)

