#Labelling Code
library(tidyverse)
library(GeoPressureR)
library(raster)
library(plotly)
library(RColorBrewer)

ID <- "CB594"
pathname <- "C:/Users/Garrett Rhyne/Documents/SWWA_Connectivity/Analysis/SWWA_Pressure-main/data/0_PAM/CB594"
#Read in Data
pam_data <- pam_read(
  pathname = pathname,
  pressure_file = ".deg",
  light_file = "driftadj.lux",
  acceleration_file = NA)
#Setup Pam Data for TRAINSET
pam_data <- pam_classify(pam_data, min_duration = 5)
trainset_write(pam_data, pathname = pathname, filename = "CB594_lt_pres-labeled-v1.csv")
