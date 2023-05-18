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
  )

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


# Fit heirarchical GAM --------------------------------------------------------------
# based on advice from Pedersen et. al, 2019
# fit a model that is consitent with model "I" from this paper 
# i.e. Group-specific smoothers with different wiggliness. 


