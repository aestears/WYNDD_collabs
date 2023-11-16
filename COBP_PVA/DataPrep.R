# Preparing raw data for further analysis
# Alice Stears Nov. 2023

# load packages
library(tidyverse)
library(readxl)

# read in raw data (in .xslx format)
# for creeks
creeks <- read_xlsx("./COBP_PVA/data/2023_wafb_gaura_master.xlsx", sheet = 1) %>% 
  filter(Year != "Average") %>% # remove last row, which is "average" 
  mutate(Year = replace(Year, Year == "2007*", "2007")) %>% 
  mutate(Year = replace(Year, Year == "2008*", "2008")) %>% 
  mutate(Year = as.integer(Year)) %>% 
  rename("year" = "Year", "CrowCreek" = "Crow Cr", "DiamondCreek" = "Diamond Cr", "UnnamedCreek" = "Unnamed Cr",
         "Total" = "WAFB (Total)", "Notes" = "...6")

# for segments
segments <- read_xlsx("./COBP_PVA/data/2023_wafb_gaura_master.xlsx", sheet = 2) %>% 
  rename("Segment" = "...1" )  %>% 
  filter(!is.na(Segment)) %>% #remove rows that are redundant
  filter(!(Segment %in% c("Crow", "Diamond", "Unnamed", "Total"))) %>% 
  pivot_longer(cols = c('1989':'2023'), names_to = "year", values_to = "count") %>% 
  pivot_wider( names_from = Segment, values_from = count)

# for polygons
# make lookup table for segment id and polygon id
seg_poly <- data.frame("Polygon_seg" = unique(polygons$Segment),
                       "Segment_seg" = c("C1", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "D1", "D2", "D3", "D4", "D5", "U1", "U2"))

polygons <- read_xlsx("./COBP_PVA/data/2023_wafb_gaura_master.xlsx", sheet = 3) %>% 
  rename("Polygon" = "Polygon ID") %>% 
  filter(!is.na(Polygon)) %>% 
  filter(!(Polygon %in% c("16", "69", "71", "CI", "CII", "CIII", "CIV", "Crow", 
                          "CV", "CVI", "CVII", "CVIII", "DI",  "Diamond", "DII", 
                          "DIII", "DIV", "DV", "Total", "UI", "UII", "Unnamed"))) %>% 
  mutate("Segment" = str_extract(polygons$Polygon, pattern = regex("[a-z][-][:alnum:]*", TRUE))) %>% 
  select(-Date) %>% 
  left_join(seg_poly, by = c("Segment" = "Polygon_seg")) %>% 
  select(-Segment) %>% 
  rename("Segment" ="Segment_seg")
  
# save output
write.csv(creeks, "./COBP_PVA/data/creek.counts.csv", row.names = FALSE)
write.csv(segments, "./COBP_PVA/data/segment.counts.csv", row.names = FALSE)
write.csv(polygons, "./COBP_PVA/data/polygon.counts.csv", row.names = FALSE)
