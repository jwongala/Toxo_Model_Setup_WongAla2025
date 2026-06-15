#!/usr/bin/env R

# setwd('/nfs7/CEOAS/Ciannelli_Lab/wongalaj/TOOTS') # set working directory
# load libraries
library(arrow)
library(tidyverse)

ptm<-proc.time() # start timer

path_f<-"/nfs7/CEOAS/Ciannelli_Lab/wongalaj/TOOTS/toots_test_06282023_oocysts_model_2019.feather" # file to open
path_r<-"toots_test_06282023_oocysts_model_2019.RData"

## load in feather file

# setwd('/nfs7/CEOAS/Ciannelli_Lab/wongalaj/TOOTS') # set working directory

df<-arrow::read_feather(path_f)

toots<-df # rename it to toots (name I use when creating figures)

save(toots, file = path_r)

proc.time() - ptm #stop time

