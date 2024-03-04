# Analysis for the April 2021 PVA Update
# Code from Tyson Wepperich, updated by Alice Stears Nov. 2023


# 1. Annual growth rates & creek carrying capacity 
# 2. Population models with site/year random intercepts and density-dependence
# 3. Climate windows with 'climwin' package
# 4. Dynamic shift detector
# 5. Simulations: creek vs segment, especially for weather effects and future resilience

# Setup ----
library(dplyr)
library(tidyr)
library(lme4)
library(climwin)
library(lubridate)
library(cowplot)
library(sjPlot)
library(minpack.lm)
library(Matrix)

'%!in%' <- function(x,y)!('%in%'(x,y))
theme_set(theme_bw(base_size = 14))
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
scale_colour_discrete <- function(...) {
  scale_colour_brewer(..., type = "qual", palette = "Dark2")
}
source('./COBP_PVA/dynamic_shift_detector.R') # from https://github.com/cbahlai/dynamic_shift_detector

creeks <- read.csv("./COBP_PVA/data/creek.counts.csv")
segments <- read.csv("./COBP_PVA/data/segment.counts.csv")
segdat <- segments %>% 
  gather(segment, flower.count, C1:U2) %>% 
  filter(segment %!in% c("C7", "C8"))

census <- read.csv("./COBP_PVA/data/creek.counts.csv") %>% 
  filter(year > 1987) %>% 
  gather(segment, flower.count, CrowCreek:UnnamedCreek) %>% 
  bind_rows(segdat) %>% 
  select(-Total)

# 1. Annual growth rates ----
popdat <- census %>% 
  group_by(segment) %>% 
  arrange(year) %>% 
  mutate(N = flower.count,
         Nt1 = lag(N, 1),
         loglam = log((N + 1)/(Nt1 + 1)),
         zyear = year - median(year),
         yearfact = as.factor(year),
         logNt1 = log(Nt1 + 1),
         date = paste0("01/08/", year)) %>% 
  filter(segment %!in% c("CrowCreek", "DiamondCreek", "UnnamedCreek")) %>% 
  select(-Notes)

# 2. Population models ----
# Gompertz (mod2) better than Ricker (mod1), mod3 (varying int/slope singular fit)
# Ricker: K = -a/b
# Gompertz: K = exp(-a/b)
# basic model w/ no predictor 
mod0 <- lme4::lmer(loglam ~ 1 + (1|segment) + (1|yearfact), 
            data = popdat, REML = FALSE)
# Ricker model
mod1 <- lmer(loglam ~ Nt1  + (1|segment) + (1|yearfact), 
            data = popdat)
#Gompertz model
mod2 <- lmer(loglam ~ logNt1  + (1|segment) + (1|yearfact), 
            data = popdat)
#Gompertz model w/ random slope for log(pop size in t-1) in each segment
mod3 <- lmer(loglam ~ logNt1  + (0 + logNt1|segment) + (1|yearfact),
             data = popdat) # not better than simpler, but segment-specific K is nice
# not sure what model type this is? ... didn't converge
mod4 <- lmer(loglam ~ zyear + logNt1  + (1 + zyear |segment) + (1|yearfact), 
             data = popdat)
AIC(mod0, mod1, mod2, mod3, mod4)
# model 2 is the best 

# Carrying capacity estimates for all years together
df <- data.frame(coef(mod2)$segment) %>% 
  mutate(K = exp(-X.Intercept. / logNt1),
         segment = row.names(coef(mod2)$segment),
         creek = c(rep("crow", 6), rep("diamond", 5), rep("unnamed", 2)))

moddat <- popdat %>% 
  ungroup() %>% 
  na.omit %>% 
  dplyr::mutate(loglam_pred = predict(mod2),
         loglam_res = residuals(mod2),
         date = paste0("01-08-", year))


# 3. Climate windows ----

# Load monthly climate data
# from https://wrcc.dri.edu/cgi-bin/cliMAIN.pl?wy1675
temp <- read.csv("./COBP_PVA/data/temp.cheyenne.csv", skip = 1, header = TRUE) %>% 
  tidyr::gather(key = month, value = temp, jan:dec) %>% 
  mutate(date = lubridate::ymd(paste(year, month, 1))) %>% 
  dplyr::select(date, temp) %>% 
  filter(year(date)>=1985) %>% 
  mutate(date = strftime(date, format = "%d/%m/%Y"))

prec <- read.csv("./COBP_PVA/data/precip.cheyenne.csv", skip = 1, header = TRUE) %>% 
  tidyr::gather(key = month, value = prec, jan:dec) %>% 
  mutate(date = lubridate::ymd(paste(year, month, 1))) %>% 
  dplyr::select(date, prec) %>% 
  filter(year(date)>=1985) %>% 
  mutate(date = strftime(date, format = "%d/%m/%Y"))

snow <- read.csv("./COBP_PVA/data/snow.cheyenne.csv", skip = 1, header = TRUE) %>% 
  tidyr::gather(key = month, value = snow, jul:jun) %>% 
  mutate(year_split = ifelse(month %in% c("jul", "aug", "sep", "oct", "nov", "dec"),
                             substr(year, 1, 4),
                             ifelse(substr(year,1,4) == "1999", 2000,
                                    paste0(substr(year, 1, 2), substr(year,6,7)))),
         date = lubridate::ymd(paste(year_split, month, 1))) %>% 
  dplyr::select(date, snow) %>% 
  filter(year(date)>=1985) %>% 
  mutate(date = strftime(date, format = "%d/%m/%Y"),
         zsnow = log(snow + 1))

# water only there for 1993-2018
# from https://waterdata.usgs.gov/wy/nwis/inventory/?site_no=06755960
water <- read.csv("./COBP_PVA/data/crow.flow.csv", header = TRUE, fileEncoding="latin1") %>% 
  tidyr::gather(key = month, value = flow, jan:dec) %>% 
  mutate(date = lubridate::ymd(paste(year, month, 1))) %>% 
  dplyr::select(date, flow) %>% 
  filter(year(date)>=1985) %>% 
  mutate(date = strftime(date, format = "%d/%m/%Y"),
         zflow = log(flow + 1))

# from https://www1.ncdc.noaa.gov/pub/data/cirs/climdiv/climdiv-pdsidv-v1.0.0-20231106
testname <- read.csv("./COBP_PVA/data/precip.cheyenne.csv", skip = 1, header = TRUE)
pdsi <- read.delim("./COBP_PVA/data/pdsi.txt", header = FALSE, sep = "")
names(pdsi) <- names(testname)[1:13]
pdsi$year <- gsub("480805", replacement = "", x = pdsi$year)
pdsi <- pdsi %>% 
  tidyr::gather(key = month, value = pdsi, jan:dec) %>% 
  mutate(date = lubridate::ymd(paste(year, month, 1))) %>% 
  dplyr::select(date, pdsi) %>% 
  filter(year(date)>=1985) %>%
  filter(!(year(date) == 2023 & month(date) %in% c(11, 12))) %>% 
  mutate(date = strftime(date, format = "%d/%m/%Y"))

# Pair environmental variables with appropriate census dataset
# Run climate window analysis (climwin package)
bestvars <- list()

# Discharge ----
popdat2 <- popdat %>% 
  ungroup() %>% 
  na.omit %>% 
  filter(year >= 1996 & year <= 2018) %>% # for water flow
  dplyr::mutate(climate = 1) # for interactions

# model with DD and annual REs? or model residuals?
outwin <- slidingwin(
                      exclude = c(3, 0),
                      xvar = list(Prec = water$zflow),
                      cdate = water$date,
                      bdate = popdat2$date,
                      baseline = lmer(loglam ~ logNt1  + (1 |segment), 
                                      data = popdat2),
                      cinterval = "month",
                      cmissing = "method2",
                      range = c(36, 1),
                      type = "absolute", refday = c(01, 08),
                      stat = "mean",
                      func = c("quad", "lin"))

# can compare to randomized repeats for pvalues
outwin$combos
outwin[[1]]$BestModel
head(outwin[[1]]$Dataset)
plt1 <- plotdelta(dataset = outwin[[1]]$Dataset)
plotweights(dataset = outwin[[1]]$Dataset)
plotbetas(dataset = outwin[[1]]$Dataset)
plotwin(dataset = outwin[[1]]$Dataset)

outdat <- outwin[[1]]$BestModelData %>% 
  mutate(pred = predict(outwin[[1]]$BestModel),
         year = popdat2$year,
         date = popdat2$date,
         climvar = "zflow")
bestvars[[1]] <- outdat

# # Commented out redundant code for model predictions, geom_smooth was sufficient

# newdata <- data.frame(expand.grid(segment = unique(outdat$segment),
#   climate = seq(min(outdat$climate), max(outdat$climate), length.out = 100),
# logNt1 = mean(outdat$logNt1))) 
# # with re.form = NA, no need to have all segments in newdata
# newdata$pred <-  predict(outwin[[1]]$BestModel, newdata = newdata, re.form = NA)
#                       
# 
# # need newdata for prediction at mean logNt1 and range of climate vars to get best line
# ggplot(outdat, aes(x = climate, y = yvar, color = year, group = segment)) +
#   geom_point() +
#   scale_color_continuous() +
#   # geom_smooth(method = "lm", formula = y ~ poly(x, 2), alpha = .05, color = "light gray")
#   geom_line(data = newdata, aes(x = climate, y = pred, group = segment), inherit.aes = FALSE)

plt2 <- ggplot(outdat, aes(x = climate, y = yvar, color = year)) +
  geom_point() +
  scale_color_continuous() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), alpha = .5, color = "black") +
  xlab("Stream discharge (23 - 7 months prior)") +
  ylab("Annual log population growth rate (flower)") +
  climwin:::theme_climwin() +
  theme(legend.position = c(.85, .2)) +
  ggtitle("Best fit model")

# Figure 2
cowplot::plot_grid(plt1, plt2, labels = c('A', 'B'), label_size = 12)


# Temperature ----
popdat2 <- popdat %>% 
  ungroup() %>% 
  na.omit %>% 
  # filter(year >= 1996 & year <= 2018) %>% # for water flow
  dplyr::mutate(climate = 1) # for interactions

# model with DD and annual REs? or model residuals?
outwin <- slidingwin(
  exclude = c(3, 0),
  xvar = list(Prec = temp$temp),
  cdate = temp$date,
  bdate = popdat2$date,
  baseline = lmer(loglam ~ logNt1  + (1 |segment), 
                  data = popdat2),
  cinterval = "month",
  cmissing = "method2",
  range = c(36, 1),
  type = "absolute", refday = c(01, 08),
  stat = "mean",
  func = c("lin","quad"))

# can compare to randomized repeats for pvalues
outwin$combos
outwin[[1]]$BestModel
head(outwin[[1]]$Dataset)
plt1 <- plotdelta(dataset = outwin[[1]]$Dataset)
plotweights(dataset = outwin[[1]]$Dataset)
plotbetas(dataset = outwin[[1]]$Dataset)
plotwin(dataset = outwin[[1]]$Dataset)

outdat <- outwin[[1]]$BestModelData %>% 
  mutate(pred = predict(outwin[[1]]$BestModel),
         year = popdat2$year,
         date = popdat2$date,
         climvar = "temp")
bestvars[[2]] <- outdat

# # Commented out redundant code for model predictions, geom_smooth was sufficient

# newdata <- data.frame(expand.grid(segment = unique(outdat$segment),
#                                   climate = seq(min(outdat$climate), max(outdat$climate), length.out = 100),
#                                   logNt1 = mean(outdat$logNt1))) 
# # with re.form = NA, no need to have all segments in newdata
# newdata$pred <-  predict(outwin[[1]]$BestModel, newdata = newdata, re.form = NA)
# 
# 
# # need newdata for prediction at mean logNt1 and range of climate vars to get best line
# ggplot(outdat, aes(x = climate, y = yvar, color = year, group = segment)) +
#   geom_point() +
#   scale_color_continuous() +
#   # geom_smooth(method = "lm", formula = y ~ poly(x, 2), alpha = .05, color = "light gray")
#   geom_line(data = newdata, aes(x = climate, y = pred, group = segment), inherit.aes = FALSE)

plt2 <- ggplot(outdat, aes(x = climate, y = yvar, color = year)) +
  geom_point() +
  scale_color_continuous() +
  geom_smooth(method = "lm", alpha = .5, color = "black") +
  xlab("Temperature (23 - 19 months prior)") +
  ylab("Annual log population growth rate (flower)") +
  climwin:::theme_climwin() +
  theme(legend.position = c(.1, .2)) +
  ggtitle("Best fit model")

# Figure 3
cowplot::plot_grid(plt1, plt2, labels = c('A', 'B'), label_size = 12)
# 875x450 copy plot


# Precipitation ----
popdat2 <- popdat %>% 
  ungroup() %>% 
  na.omit %>% 
  # filter(year >= 1996 & year <= 2018) %>% # for water flow
  dplyr::mutate(climate = 1) # for interactions

# model with DD and annual REs? or model residuals?
outwin <- slidingwin(
  exclude = c(3, 0),
  xvar = list(Prec = prec$prec),
  cdate = prec$date,
  bdate = popdat2$date,
  baseline = lmer(loglam ~ logNt1  + (1 |segment), 
                  data = popdat2),
  cinterval = "month",
  cmissing = "method2",
  range = c(36, 1),
  type = "absolute", refday = c(01, 08),
  stat = "mean",
  func = c("quad", "lin"))

# can compare to randomized repeats for pvalues
outwin$combos
outwin[[1]]$BestModel
head(outwin[[1]]$Dataset)
plt1 <- plotdelta(dataset = outwin[[1]]$Dataset)
plotweights(dataset = outwin[[1]]$Dataset)
plotbetas(dataset = outwin[[1]]$Dataset)
plotwin(dataset = outwin[[1]]$Dataset)
plotbest(dataset = outwin[[1]]$Dataset,
         bestmodel = outwin[[1]]$BestModel, 
         bestmodeldata = outwin[[1]]$BestModelData)

outdat <- outwin[[1]]$BestModelData %>% 
  mutate(pred = predict(outwin[[1]]$BestModel),
         year = popdat2$year,
         date = popdat2$date,         
         climvar = "prec")
bestvars[[3]] <- outdat

# # Commented out redundant code for model predictions, geom_smooth was sufficient

# newdata <- data.frame(expand.grid(segment = unique(outdat$segment),
#                                   climate = seq(min(outdat$climate), max(outdat$climate), length.out = 100),
#                                   logNt1 = mean(outdat$logNt1)))
# # with re.form = NA, no need to have all segments in newdata
# newdata$pred <-  predict(outwin[[1]]$BestModel, newdata = newdata, re.form = NA)
# 
# 
# # need newdata for prediction at mean logNt1 and range of climate vars to get best line
# ggplot(outdat, aes(x = climate, y = yvar, color = year, group = segment)) +
#   geom_point() +
#   scale_color_continuous() +
#   # geom_smooth(method = "lm", formula = y ~ poly(x, 2), alpha = .05, color = "light gray")
#   geom_line(data = newdata, aes(x = climate, y = pred, group = segment), inherit.aes = FALSE)

plt2 <- ggplot(outdat, aes(x = climate, y = yvar, color = year)) +
  geom_point() +
  scale_color_continuous() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), alpha = .5, color = "black") +
  xlab("Precipitation (24 - 2 months prior)") +
  ylab("Annual log population growth rate (flower)") +
  climwin:::theme_climwin() +
  theme(legend.position = c(.82, .2)) +
  ggtitle("Best fit model")

# Figure 4
cowplot::plot_grid(plt1, plt2, labels = c('A', 'B'), label_size = 12)

# Snow ----
popdat2 <- popdat %>% 
  ungroup() %>% 
  na.omit %>% 
  # filter(year >= 1996 & year <= 2018) %>% # for water flow
  dplyr::mutate(climate = 1) # for interactions

# model with DD and annual REs? or model residuals?
outwin <- slidingwin(
  exclude = c(3, 0),
  xvar = list(Prec = snow$zsnow),
  cdate = snow$date,
  bdate = popdat2$date,
  baseline = lmer(loglam ~ logNt1  + (1 |segment), 
                  data = popdat2),
  cinterval = "month",
  cmissing = "method2",
  range = c(36, 1),
  type = "absolute", refday = c(01, 08),
  stat = "mean",
  func = c("quad", "lin"))

# can compare to randomized repeats for pvalues
outwin$combos
outwin[[1]]$BestModel
head(outwin[[1]]$Dataset)
plt1 <- plotdelta(dataset = outwin[[1]]$Dataset)
plotweights(dataset = outwin[[1]]$Dataset)
plotbetas(dataset = outwin[[1]]$Dataset)
plotwin(dataset = outwin[[1]]$Dataset)
plotbest(dataset = outwin[[1]]$Dataset,
         bestmodel = outwin[[1]]$BestModel, 
         bestmodeldata = outwin[[1]]$BestModelData)

outdat <- outwin[[1]]$BestModelData %>% 
  mutate(pred = predict(outwin[[1]]$BestModel),
         year = popdat2$year,
         date = popdat2$date,
         climvar = "zsnow")

bestvars[[4]] <- outdat

# # Commented out redundant code for model predictions, geom_smooth was sufficient

# newdata <- data.frame(expand.grid(segment = unique(outdat$segment),
#                                   climate = seq(min(outdat$climate), max(outdat$climate), length.out = 100),
#                                   logNt1 = mean(outdat$logNt1))) 
# # with re.form = NA, no need to have all segments in newdata
# newdata$pred <-  predict(outwin[[1]]$BestModel, newdata = newdata, re.form = NA)
# 
# 
# # need newdata for prediction at mean logNt1 and range of climate vars to get best line
# ggplot(outdat, aes(x = climate, y = yvar, color = year, group = segment)) +
#   geom_point() +
#   scale_color_continuous() +
#   # geom_smooth(method = "lm", formula = y ~ poly(x, 2), alpha = .05, color = "light gray")
#   geom_line(data = newdata, aes(x = climate, y = pred, group = segment), inherit.aes = FALSE)

plt2 <- ggplot(outdat, aes(x = climate, y = yvar, color = year)) +
  geom_point() +
  scale_color_continuous() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), alpha = .5, color = "black") +
  xlab("Snowfall (4 - 1 months prior)") +
  ylab("Annual log population growth rate (flower)") +
  theme(legend.position = c(.9, .2)) +
  climwin:::theme_climwin() +
  ggtitle("Best fit model")

# Figure 5
cowplot::plot_grid(plt1, plt2, labels = c('A', 'B'), label_size = 12)


# PDSI ----
popdat2 <- popdat %>% 
  ungroup() %>% 
  na.omit %>% 
  # filter(year >= 1996 & year <= 2018) %>% # for water flow
  dplyr::mutate(climate = 1) # for interactions

# model with DD and annual REs? or model residuals?
outwin <- slidingwin(
  exclude = c(3, 0),
  xvar = list(Prec = pdsi$pdsi),
  cdate = pdsi$date,
  bdate = popdat2$date,
  baseline = lmer(loglam ~ logNt1  + (1 |segment), 
                  data = popdat2),
  cinterval = "month",
  cmissing = "method2",
  range = c(36, 1),
  type = "absolute", refday = c(01, 08),
  stat = "mean",
  func = c("lin", "quad"))

# can compare to randomized repeats for pvalues
outwin$combos
outwin[[1]]$BestModel
head(outwin[[1]]$Dataset)
plt1 <- plotdelta(dataset = outwin[[1]]$Dataset)
plotweights(dataset = outwin[[1]]$Dataset)
plotbetas(dataset = outwin[[1]]$Dataset)
plotwin(dataset = outwin[[1]]$Dataset)
plotbest(dataset = outwin[[1]]$Dataset,
         bestmodel = outwin[[1]]$BestModel, 
         bestmodeldata = outwin[[1]]$BestModelData)

outdat <- outwin[[1]]$BestModelData %>% 
  mutate(pred = predict(outwin[[1]]$BestModel),
         year = popdat2$year,
         date = popdat2$date,
         climvar = "pdsi")
bestvars[[5]] <- outdat

# # Commented out redundant code for model predictions, geom_smooth was sufficient

# newdata <- data.frame(expand.grid(segment = unique(outdat$segment),
#                                   climate = seq(min(outdat$climate), max(outdat$climate), length.out = 100),
#                                   logNt1 = mean(outdat$logNt1))) 
# # with re.form = NA, no need to have all segments in newdata
# newdata$pred <-  predict(outwin[[1]]$BestModel, newdata = newdata, re.form = NA)
# 
# 
# # need newdata for prediction at mean logNt1 and range of climate vars to get best line
# ggplot(outdat, aes(x = climate, y = yvar, color = year, group = segment)) +
#   geom_point() +
#   scale_color_continuous() +
#   # geom_smooth(method = "lm", formula = y ~ poly(x, 2), alpha = .05, color = "light gray")
#   geom_line(data = newdata, aes(x = climate, y = pred, group = segment), inherit.aes = FALSE)

plt2 <- ggplot(outdat, aes(x = climate, y = yvar, color = year)) +
  geom_point() +
  scale_color_continuous() +
  geom_smooth(method = "lm", alpha = .5, color = "black") +
  xlab("PDSI (36-1 months prior)") +
  ylab("Annual log population growth rate (flower)") +
  climwin:::theme_climwin() +
  theme(legend.position = c(.15, .2))+
  ggtitle("Best fit model")

#Figure 6
cowplot::plot_grid(plt1, plt2, labels = c('A', 'B'), label_size = 12)

# Correlations ----
envdat <- bind_rows(bestvars)
# write.csv(envdat, "data/bestwindowvars.csv", row.names = FALSE)
# envdat <- read.csv("data/bestwindowvars.csv", header = TRUE)
# 
# envdat <- envdat %>% 
#   dplyr::select(yvar, logNt1, segment, climate, year, climvar) %>% 
#   tidyr::pivot_wider(names_from = climvar, values_from = climate) %>% 
#   filter(complete.cases(.))
# 
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
# # with zflow, limited years bc missing data
# pairs(envdat[,5:9], lower.panel = panel.smooth, upper.panel = panel.cor,
#       gap=0, row1attop=FALSE)

# removing zflow and using all years
envdat <- envdat %>% 
  filter(climvar != "zflow") %>% 
  dplyr::select(yvar, logNt1, segment, climate, year, climvar) %>% 
  tidyr::pivot_wider(names_from = climvar, values_from = climate) %>% 
  mutate(yearfact = as.character(year),
         zyear = year - 2005,
         ztemp = scale(temp)[,1],
         zprec = scale(prec)[,1],
         zprec2 = poly(zprec,2)[,2],
         zsnow2 = poly(zsnow,2)[,2])
pairs(envdat[,c(7,8,11,12,13)], lower.panel = panel.smooth, upper.panel = panel.cor,
      gap=0, row1attop=FALSE) 

# Table 1: Mixed effects models  ----
# Optimal climate windows for variables
# Full model
# Model 2 in Table 1

mod2 <- lmer(yvar ~ logNt1  + 
               ztemp +
               zprec +
               zprec2 +
               zsnow +
               # zsnow2 +
               #pdsi +
               (1 |segment) + (1|yearfact), 
             data = envdat)
summary(mod2)
 
# Model 1 in Table 1                   
mod4 <- lmer(yvar ~ logNt1  + 
               (1 |segment) + (1|yearfact), 
             data = envdat)
summary(mod4)

stargazer(mod4, mod2, type = "text", ci = TRUE, p.auto = TRUE, report = c("vcsp"))
sjPlot::tab_model(mod4, mod2)

# # other varying slopes lead to fitting problems or singular fit
# # precip probably most important/robust anyways
# mod3 <- lmer(yvar ~ logNt1  + 
#                ztemp +
#                zprec +
#                zprec2 +
#                zsnow +
#                (1 + zprec + zprec2 + ztemp|segment) + (1|yearfact), 
#              data = envdat)
# summary(mod3)
# AIC(mod2, mod3)

df <- data.frame(coef(mod2)$segment) 
df <- df %>% 
  mutate(K = exp(-X.Intercept. / logNt1),
         segment = row.names(coef(mod2)$segment))

# Main effect of precip with lines for each segments varying slopes (messy!)

newdata <- data.frame(expand.grid(segment = unique(envdat$segment),
                                  logNt1 = c(2, 4, 6, 8),
                                  ztemp = mean(envdat$ztemp),
                                  zprec = seq(min(envdat$zprec), max(envdat$zprec), length.out = 100),
                                  zprec2 = seq(min(envdat$zprec), max(envdat$zprec), length.out = 100)^2,
                                  zsnow = mean(envdat$zsnow),
                                  pdsi = mean(envdat$pdsi))) %>% 
  mutate(zprec2 = zprec^2)
newdata$pred <- predict(mod2, newdata = newdata, re.form = ~ (1 |segment))

#? not sure what this is for? 
plt5 <- plot_model(mod2, type = "pred", terms = "zprec",
                   axis.title = c("Mean precipitation in inches (24-2 months prior)", "Log(annual population growth rate)"),  
                   title = "") +
  geom_line(data = newdata, aes(x = zprec, y = pred, group = segment), color = "light blue", inherit.aes = FALSE)

ggplot(newdata, aes(x = zprec, y = pred, group = segment)) +
  geom_line() +
  facet_wrap(~logNt1)

# marginal effects for 4 fixed effects
# Didn't use this in the report. Would need to edit to include quadratic precip effect.
plt1 <- plot_model(mod2, type = "pred", terms = "logNt1", 
           axis.title = c("Log(count prior year)", "Log(annual population growth rate)"), 
           title = "", axis.lim = c(-3, 3))
plt2 <- plot_model(mod2, type = "pred", terms = "ztemp",
                   axis.title = c("Mean temperature in F (23-19 months prior)", "Log(annual population growth rate)"), 
                   title = "", axis.lim = c(-3, 3))
plt3 <- plot_model(mod2, type = "pred", terms = "zprec",
                   axis.title = c("Mean precipitation in inches (24-2 months prior)", "Log(annual population growth rate)"),  
                   title = "", axis.lim = c(-3, 3))
plt4 <- plot_model(mod2, type = "pred", terms = "zsnow",
                   axis.title = c("Log(mean snowfall in inches) (4-1 months prior)", "Log(annual population growth rate)"), 
                   title = "", axis.lim = c(-3, 3))
cowplot::plot_grid(plt1, plt2, plt3, plt4, labels = c('A', 'B', "C", "D"), label_size = 12)


# 4. Dynamic shift detector ----
# Gompertz (modified the Ricker model in dynamic_shift_detector.R)  

out <- census %>% 
  select(-Notes) %>% 
  group_by(segment) %>%
  nest()

# Warning! takes >30 minutes to run algorithm on all segments
# Load saved results below
breaks <- lapply(out$data, function(x) DSdetector(data = data.frame(x), criterion = "AICc"))

# names(breaks) <- unique(census$segment)
#saveRDS(breaks, "DSD_segment_naive_gompertz.rds")

#breaks <- readRDS("./COBP_PVA/DSD_segment_naive_gompertz.rds")

# calcute K for each model 
breaks_results <- lapply(1:length(breaks), FUN = function(x) dplyr::mutate(breaks[[x]][[3]], segment = x)) %>% 
  bind_rows() %>% 
  rename(a = r, b = k, a_se = rse, b_se = kse) %>% # update names for gompertz model 
  mutate(K = exp(-a/b)) # calculate K 



# 5. Simulations ----

# simulate future by decade (simplistic, maybe MACAv2 could be source of more realistic scenarios)
# divide into two 8 year segments (and two 9 year segment) (33 years of data)
future_env <- envdat %>% 
  select(year, temp, prec, zsnow, pdsi, ztemp, zprec, zprec2) %>% 
  #filter(year != 2020) %>% 
  group_by(year) %>% 
  slice(1) %>% 
  mutate(decade = case_when(year <= 1998 ~ "90-98",
                            year >= 1999 & year < 2007 ~ "99-06",
                            year >= 2007 & year < 2016 ~ "07-15",
                            year >= 2016 ~ "16-23"))

future_summ <- future_env %>% 
  group_by(decade) %>% 
  summarise_all(mean)

# 16 different 20-year scenarios
scens <- expand.grid(unique(future_env$decade), unique(future_env$decade)) %>% 
  mutate(scenario = row_number()) #%>% 
  # filter(scenario == 1)
  # filter(scenario %in% c(1, 2, 3, 5, 6, 9))
nsim <- 100

segcounts <- census %>% filter(segment %in% envdat$segment, year == 2020)

# loop to get model coefficients for end of monitoring period
# ignore errors resulting from assigning values past the end of the time period.
# Warning: takes a while to run nsim x scenario!
# Load saved simulation results after loop
outseg <- list()
outcreek <- list()
for (sn in unique(scens$scenario)){
  for(sim in 1:nsim){
    
    weather1 <- future_env %>% filter(decade == scens$Var1[which(scens$scenario == sn)])
    weather2 <- future_env %>% filter(decade == scens$Var2[which(scens$scenario == sn)])
    weather <- bind_rows(weather1, weather2) %>% 
      ungroup() %>% 
      mutate(#year = 2021:2040,
              year = sample(c(2023:2040), size = (nrow(weather1) + nrow(weather2)), replace = FALSE), # doesn't do autocorrelated years
             sim = sim,
             scenario = sn,
             annvar = rnorm((nrow(weather1) + nrow(weather2)), mean = 0, sd = .36)) %>% # from mod2 random intercept
      arrange(year)
    
    moddat <- expand_grid(weather, segment = unique(envdat$segment)) %>% 
      mutate(flower.count = NA,
             logNt1 = NA,
             loglam = NA)
    
    # setup counts from 2020
    moddat$logNt1[moddat$year == min(moddat$year)] <- log(segcounts$flower.count + 1)
    
    for (nr in 1:nrow(moddat)){
      moddat$loglam[nr] <- predict(mod2, newdata = moddat[nr,], re.form = ~(1|segment)) + moddat$annvar[nr] + rnorm(1, 0, 0.721) # from mod2 random intercept
      mu <- exp(moddat$logNt1[nr]) * exp(moddat$loglam[nr]) 
      moddat$flower.count[nr] <- rpois(n = 1, lambda = mu)
      try(moddat$logNt1[nr+13] <- log(moddat$flower.count[nr] + 1))
    } # note: we do see error messages, but the code still runs b/c it's wrapped in a "try" 
    
    outseg[[length(outseg)+1]] <- moddat
  }
}
segs <- bind_rows(outseg)
saveRDS(segs, "segpreds.rds")


segs <- readRDS("segpreds.rds")

## remove the "NA" flower count predictions from 'segs'??
segs <- segs[!is.na(segs$flower.count),]
creeks <- segs %>% 
  mutate(scenario = as.factor(scenario)) %>% 
  # mutate(scenario = case_when(scenario == 1 ~ "1990-1999",
  #                             scenario %in% c(2, 4) ~ "1990-2009",
  #                             scenario %in% c(3, 7) ~ "1990-1999, 2010-2019",
  #                             scenario %in% c(6, 8) ~ "2000-2019",
  #                             scenario == 5 ~ "2000-2009",
  #                             scenario == 9 ~ "2010-2019")) %>%
  ungroup() %>% 
  mutate(creek = stringr::str_sub(segment, start = 0, end = 1)) %>% 
  group_by(creek, year, scenario, sim) %>% 
  summarise(flower.count = sum(flower.count)) %>% 
  group_by(creek, year, scenario) %>% 
  summarise(pred_median = quantile(flower.count, probs = .5),
            pred_lo = quantile(flower.count, probs = .75),
            pred_hi = quantile(flower.count, probs = .25),
            pred_lower = quantile(flower.count, probs = .05),
            pred_higher = quantile(flower.count, probs = .95))
all <- segs %>% 
  mutate(scenario = as.factor(scenario)) %>% 
  # mutate(scenario = case_when(scenario == 1 ~ "1990-1999",
  #                             scenario %in% c(2, 4) ~ "1990-2009",
  #                             scenario %in% c(3, 7) ~ "1990-1999, 2010-2019",
  #                             scenario %in% c(6, 8) ~ "2000-2019",
  #                             scenario == 5 ~ "2000-2009",
  #                             scenario == 9 ~ "2010-2019")) %>%
  ungroup() %>%
  group_by(year, scenario, sim) %>% 
  summarise(flower.count = sum(flower.count)) %>% 
  group_by(year, scenario) %>% 
  summarise(pred_median = quantile(flower.count, probs = .5, na.rm = TRUE),
            pred_lo = quantile(flower.count, probs = .25),
            pred_hi = quantile(flower.count, probs = .75),
            pred_lower = quantile(flower.count, probs = .05),
            pred_higher = quantile(flower.count, probs = .95),
            creek = "Total")
creeks <- bind_rows(creeks, all)

creeks$creek <- factor(creeks$creek, levels = c("Total", "C", "D", "U"),
                       labels = c("Total", "Crow", "Diamond", "Unnamed"))

# Figure 8 ----
ggplot(creeks %>% filter(scenario %in% c("1", "6", "11", "16")), 
       aes(x = year, y = pred_median, group = scenario, color = scenario)) +
  facet_wrap(~creek, ncol = 2, scales = "free") +
  geom_point() +
  geom_path() +
  scale_colour_discrete(labels = c("1990-98", "1999-2006", "2007-15", "2016-23")) +
  # geom_smooth(method = "gam", se = FALSE) +
  geom_linerange(aes(x = year, ymin = pred_lo, ymax = pred_hi)) +
  scale_y_continuous(limits = c(0, NA)) +
  xlab("Year") +
  ylab("Forecasted flower count")


total <- census %>% 
  filter(year > 1987, segment %in% c("CrowCreek", "DiamondCreek", "UnnamedCreek")) %>%
  group_by(year) %>%
  summarise(flower.count = sum(flower.count),
            segment = "total")
census2 <- bind_rows(census, total) %>% 
  filter(segment %in% c("CrowCreek", "DiamondCreek", "UnnamedCreek", "total")) %>% 
  mutate(creek = as.character(segment),
         pred_median = flower.count) %>% 
  expand_grid(scens) %>% 
  mutate(scenario = as.factor(scenario)) %>% 
  select(-Var1, -Var2, -segment)
creeks2 <- creeks %>%
  mutate(creek = tolower(as.character(creek)))
census3 <- bind_rows(census2, creeks2) %>% 
  filter(scenario %in% c("1", "6", "11", "16"))
# rename creeks 
census3[census3$creek == "crow", "creek"] <- "CrowCreek"
census3[census3$creek == "unnamed", "creek"] <- "UnnamedCreek"
census3[census3$creek == "diamond", "creek"] <- "DiamondCreek"

census3$creek <- factor(census3$creek)
census3$creek <- relevel(census3$creek, ref = "total")

census3 <- census3 %>% 
  mutate(scenario = as.character(scenario)) %>% 
  mutate(scenario = case_when(year < 2024 ~ "observed",
                              TRUE ~ scenario))

# Figure 9 ----
census3$scenario <- factor(census3$scenario, levels = c(1, 6, 11, 16, "observed"))

ggplot(census3, aes(x = year, y = pred_median, group = scenario, color = scenario)) +
  geom_point() +
  geom_path() +
  scale_y_continuous(limits = c(0, NA)) +
  scale_colour_discrete(labels = c("1990-98", "1999-2006",  "2007-15", "2016-23", "Observed")) +
  # geom_linerange(aes(x = year, ymin = pred_lo, ymax = pred_hi)) +
  # geom_smooth(method = "gam", se = FALSE, formula = y ~ s(x, bs = "cs", k = 17),  alpha = .5) +  
  facet_wrap(~creek, ncol = 2, scales = "free") +
  xlab("Year") +
  ylab("Median forecasted flower count")



# segments only
# segs <- bind_rows(outseg)
creeks <- segs %>% 
  mutate(scenario = as.factor(scenario)) %>% 
  # mutate(scenario = case_when(scenario == 1 ~ "1990-1999",
  #                             scenario %in% c(2, 4) ~ "1990-2009",
  #                             scenario %in% c(3, 7) ~ "1990-1999, 2010-2019",
  #                             scenario %in% c(6, 8) ~ "2000-2019",
  #                             scenario == 5 ~ "2000-2009",
  #                             scenario == 9 ~ "2010-2019")) %>% 
  ungroup() %>% 
  group_by(segment, year, scenario) %>% 
  summarise(pred_median = quantile(flower.count, probs = .5),
            pred_lo = quantile(flower.count, probs = .25),
            pred_hi = quantile(flower.count, probs = .75),
            pred_lower = quantile(flower.count, probs = .05),
            pred_higher = quantile(flower.count, probs = .95))

segsumm <- segs %>% 
  mutate(scenario = as.factor(scenario)) %>% 
  filter(scenario %in% c("1", "6", "11", "16")) %>% 
  group_by(segment, scenario) %>% 
  summarise(count_median = quantile(flower.count, probs = .5),
            count_lo = quantile(flower.count, probs = .25),
            count_hi = quantile(flower.count, probs = .75),
            count_lower = quantile(flower.count, probs = .05),
            count_higher = quantile(flower.count, probs = .95))
segsumm$scenario <- factor(segsumm$scenario, labels = c("1990-98","1999-2006", "2007-15", "2016-23"))

obssumm <- census %>% 
  filter(segment %!in% c("CrowCreek", "DiamondCreek", "UnnamedCreek")) %>% 
  group_by(segment) %>% 
  summarise(scenario = "observed",
            count_median = quantile(flower.count, probs = .5),
            count_lo = quantile(flower.count, probs = .25),
            count_hi = quantile(flower.count, probs = .75),
            count_lower = quantile(flower.count, probs = .05),
            count_higher = quantile(flower.count, probs = .95)) 

segsumm2 <- bind_rows(segsumm, obssumm)

# likelihood below threshold?
probs <- obssumm %>% 
  rename_at(vars(count_median:count_higher), toupper) %>% 
  select(-scenario) %>% 
  left_join(segsumm) %>% 
  mutate(percmed = paste0(round(100 * count_median / COUNT_MEDIAN), "%")) %>% 
  filter(!is.na(scenario))

# Figure 11 ----
#total range of flower.counts by segments
fig11 <- ggplot(segsumm2, aes(x = segment, y = count_median, group = scenario, color = scenario)) +
  geom_point(position = position_dodge(width = .7), size = 2) +
  geom_linerange(aes(ymin = count_lo, ymax = count_hi), inherit.aes = TRUE, size = 2, position = position_dodge(width = .7)) +
  geom_linerange(aes(ymin = count_lower, ymax = count_higher), inherit.aes = TRUE, size = 1, position = position_dodge(width = .7)) +
  geom_text(probs, mapping = aes(x = segment, y = -500, group = scenario, color = scenario, label = percmed), 
            inherit.aes = TRUE, position = position_dodge(width = .7), size = 3) +
  ylab("Distribution of flower counts (Median, 50%, 95%)") +
  coord_flip()  #+
  #scale_colour_discrete(labels = c("1990-98", "1999-2006",  "2007-15", "2016-23", "Observed")) 

pdf(file = "./COBP_PVA/report_fig11.pdf", width = 8, height = 9.5)  
fig11
dev.off()

#annual segment forecasts
ggplot(creeks %>% filter(scenario %in% c("1", "5", "9")), aes(x = year, y = pred_median, group = scenario, color = scenario)) +
  geom_point() +
  geom_path() +
  # geom_smooth(method = "gam", se = FALSE) +
  geom_linerange(aes(x = year, ymin = pred_lo, ymax = pred_hi)) +
  scale_y_continuous(limits = c(0, NA)) +
  facet_wrap(~segment, ncol = 3, scales = "free") +
  xlab("Year") +
  ylab("Median forecasted flower count")



# same but for creek level

total <- census %>% 
  filter(year > 1987, segment %in% c("CrowCreek", "DiamondCreek", "UnnamedCreek")) %>%
  group_by(year) %>%
  summarise(flower.count = sum(flower.count),
            segment = "Total")
censumm <- census %>% 
  filter(year > 1987, segment %in% c("CrowCreek", "DiamondCreek", "UnnamedCreek")) %>%
  bind_rows(total) %>% 
  group_by(segment) %>% 
  summarise(scenario = "observed",
            pred_median = quantile(flower.count, probs = .5),
            pred_lo = quantile(flower.count, probs = .25),
            pred_hi = quantile(flower.count, probs = .75),
            pred_lower = quantile(flower.count, probs = .05),
            pred_higher = quantile(flower.count, probs = .95)) %>% 
  ungroup() %>% 
  mutate(creek = segment)

creeks <- segs %>% 
  mutate(scenario = as.factor(scenario)) %>% 
  ungroup() %>% 
  mutate(creek =stringr::str_sub(segment, start = 0, end = 1)) %>% 
  group_by(year, creek, scenario, sim) %>% 
  summarise(flower.count = sum(flower.count)) %>% 
  group_by(creek, scenario) %>% 
  summarise(pred_median = quantile(flower.count, probs = .5),
            pred_lo = quantile(flower.count, probs = .25),
            pred_hi = quantile(flower.count, probs = .75),
            pred_lower = quantile(flower.count, probs = .05),
            pred_higher = quantile(flower.count, probs = .95))
all <- segs %>% 
  mutate(scenario = as.factor(scenario)) %>% 
  ungroup() %>% 
  group_by(year, scenario, sim) %>% 
  summarise(flower.count = sum(flower.count)) %>% 
  group_by(scenario) %>% 
  summarise(pred_median = quantile(flower.count, probs = .5),
            pred_lo = quantile(flower.count, probs = .25),
            pred_hi = quantile(flower.count, probs = .75),
            pred_lower = quantile(flower.count, probs = .05),
            pred_higher = quantile(flower.count, probs = .95),
            creek = "Total")
creeks <- bind_rows(creeks, all)

creeks$creek <- factor(creeks$creek, levels = c("Total", "C", "D", "U"),
                       labels = c("Total", "CrowCreek", "DiamondCreek", "UnnamedCreek"))
segsumm2 <- creeks %>% 
  ungroup() %>% 
  #mutate(creek = tolower(as.character(creek))) %>% 
  bind_rows(censumm) %>% 
  filter(scenario %in% c("observed", "1", "6", "11", "16"))

# likelihood below threshold?
probs <- censumm %>% 
  rename_at(vars(pred_median:pred_higher), toupper) %>% 
  select(-scenario) %>% 
  left_join(segsumm2, by = "creek") %>% 
  mutate(percmed = paste0(round(100 * pred_median / PRED_MEDIAN), "%")) %>% 
  filter(scenario %in% c("1", "6", "11", "16")) %>% 
  bind_rows(data.frame(creek = c("CrowCreek", "DiamondCreek", "Total", "UnnamedCreek"),
                       scenario = "observed",
                       percmed = ""))

segsumm2$creek <- factor(segsumm2$creek)
segsumm2$creek <- relevel(segsumm2$creek, ref = "Total")

segsumm2$scenario <- factor(segsumm2$scenario, labels = c("1990-98", "1999-2006", "2007-15", "2016-23", "observed"))
probs$scenario <- factor(probs$scenario, labels = c("1990-98", "1999-2006", "2007-15", "2016-23", "observed"))
# Figure 10 ----
# total range of flower.counts by segments
ggplot(segsumm2, aes(x = creek, y = pred_median, group = scenario, color = scenario)) +
  geom_point(position = position_dodge(width = .5), size = 2) +
  geom_linerange(aes(ymin = pred_lo, ymax = pred_hi), inherit.aes = TRUE, size = 2, position = position_dodge(width = .5)) +
  geom_linerange(aes(ymin = pred_lower, ymax = pred_higher), inherit.aes = TRUE, size = 1, position = position_dodge(width = .5)) +
  geom_text(probs, mapping = aes(x = creek, y = -1000, group = scenario, color = scenario, label = percmed), 
            inherit.aes = TRUE, position = position_dodge(width = .5), size = 3) +
  #scale_colour_discrete(labels = c("1990s", "2000s", "2010s", "observed")) +
  ylab("Distribution of flower counts (Median, 50%, 95%)") +
  coord_flip()  





