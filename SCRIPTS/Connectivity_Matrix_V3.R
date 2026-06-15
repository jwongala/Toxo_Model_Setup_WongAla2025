### Jennifer Wong-Ala
### Purpose: Make connectivity matrices for ocean parcels TOOTS output


#***NOTE***: 
# - about number of particles released for each particle forcing run
# - less particles released for oocyst model forcing from OA_kane buffer
# - more particles released from the three other Oahu buffer areas (ss, sw, nw) using the oocyst model forcing 


##########################################################
### Clear workspace 

rm(list=ls())

ptm<- proc.time() # start timer


##########################################################
### Load input variables

file_con<-'TOOTS_NESTED_OAHUSS_2018_constant_particleforcing_01'
paste_filecon<-paste0(file_con, '.','RData')

file_oocyst<-'TOOTS_NESTED_OAHUSS_2018_oocysts_particleforcing_01'
paste_fileoocy<-paste0(file_oocyst, '.','RData')
 
wd<-"/Users/wongalaj/TOOTS/TOOTS_parcels/DATA/RData"

age_par<-10 # age in days for CM you want to make 
  # (10,20,30,60,90) 
  ## I checked for longer PLDs and all of the particles are gone. Would be better to do traj_end plot 
  
  # Max age for each year
  ## 2018 = 240 days
  ## 2019 = 365 days
  ## 2020 = 365 days
  ## 2021 = 365 days
  
sec_to_day<-86400 # seconds to day conversion

island<-9 # do numbers from 0-8 (each island is 1-8, and I'll use 0 to keep all islands)
  # All islands = 9
  # Maui Nui = 'MN'
  # niihau = 1
  # kauai = 2
  # oahu = 3
  # molokai = 4
  # lanai = 5
  # kahoolawe = 6
  # maui = 7
  # Hawaii/Big Island = 8

island_name<- 'oahu'
  # All islands = MHI
  # Maui nui = MN
  # Hawaii Island = BI
  # all other islands put full names

##########################################################
### Load libraries

library(tidyverse)
library(marmap)
library(ncdf4)
library(maps)
library(fields)
library(sf)
library(raster)
library(grDevices)
library(sgeostat)
library(mapdata)
library(rgeos)


##########################################################
### Load data


#### Load constant simulation data
setwd(wd)
load(paste_filecon)

# dim(toots)
# head(toots)
# class(toots)

toots_con<-data.frame(toots) # convert it to a non-tbl dataframe 
rm(toots)
head(toots_con)
class(toots_con) # now it's a dataframe

colnames(toots_con)<-c('traj_id', 'time', 'lat_dat', 'lon_dat', 'depth_dat', 'dU', 'dV', 'd2s', 'age_dat', 'release_dat', 'island_dat', 'obs', 'traj')


##########################################################

#### Load constant simulation data
setwd(wd)
load(paste_fileoocy)

# dim(toots)
# head(toots)
# class(toots)

toots_oocyst<-data.frame(toots) # convert it to a non-tbl dataframe 
rm(toots)
head(toots_oocyst)
class(toots_oocyst) # now it's a dataframe

colnames(toots_oocyst)<-c('traj_id', 'time', 'lat_dat', 'lon_dat', 'depth_dat', 'dU', 'dV', 'd2s', 'age_dat', 'release_dat', 'island_dat', 'obs', 'traj')

# head(toots)
# range(toots$age_dat, na.rm=T)


##########################################################

### load in buffer polygon points
setwd( "/Users/wongalaj/TOOTS/TOOTS_parcels/")

load('BIdpts_shpfile_polygon.RData')
load('MNdpts_shpfile_polygon.RData')
load('OAdpts_kane_shpfile_polygon.RData')
load('OAdpts_oahuss_shpfile_polygon.RData')
load('OAdpts_swoahu_shpfile_polygon.RData')
load('OAdpts_nwoahu_shpfile_polygon.RData')
load('KAdpts_shpfile_polygon.RData')
load('NIdpts_shpfile_polygon.RData')


##########################################################
### Subset particles with ages of 0, 10, 20, 30 days

# subset out for one island
if(island==9) {
  
  ## ALL ISLANDS
  ## constant
  age_df_con<-subset(toots_con, age_dat==age_par*sec_to_day) # particles that made it to 90 days
  age0_con<-subset(toots_con, age_dat<=0) # initial particles at each location
  
  ## oocysts model
  age_df_oocyst<-subset(toots_oocyst, age_dat==age_par*sec_to_day) # particles that made it to 90 days
  age0_oocyst<-subset(toots_oocyst, age_dat<=0) # initial particles at each location
  
} else if (island==2) {
  
  ## Kauai
  ## constant
  toots_con2<-subset(toots_con, island_dat == 2) 
  age_df_con<-subset(toots_con2, age_dat==age_par*sec_to_day) # particles that made it to 90 days
  age0_con<-subset(toots_con2, age_dat<=0) # initial particles at each location
  
  ## oocysts model
  toots_oocyst2<-subset(toots_oocyst, island_dat == 2) 
  age_df_oocyst<-subset(toots_oocyst2, age_dat==age_par*sec_to_day) # particles that made it to 90 days
  age0_oocyst<-subset(toots_oocyst2, age_dat<=0) # initial particles at each location
  
} else if (island==3) {
  
  ## Oahu
  ## constant
  toots_con2<-subset(toots_con, island_dat == 3) 
  age_df_con<-subset(toots_con2, age_dat==age_par*sec_to_day) # particles that made it to 90 days
  age0_con<-subset(toots_con2, age_dat<=0) # initial particles at each location
  
  ## oocysts model
  toots_oocyst2<-subset(toots_oocyst, island_dat == 3) 
  age_df_oocyst<-subset(toots_oocyst2, age_dat==age_par*sec_to_day) # particles that made it to 90 days
  age0_oocyst<-subset(toots_oocyst2, age_dat<=0) # initial particles at each location
  
} else if (island==4){
  
  ## Molokai
  ## constant
  toots_con2<-subset(toots_con, island_dat == 4) 
  age_df_con<-subset(toots_con2, age_dat==age_par*sec_to_day) # particles that made it to 90 days
  age0_con<-subset(toots_con2, age_dat<=0) # initial particles at each location
  
  ## oocysts model
  toots_oocyst2<-subset(toots_oocyst, island_dat == 4) 
  age_df_oocyst<-subset(toots_oocyst2, age_dat==age_par*sec_to_day) # particles that made it to 90 days
  age0_oocyst<-subset(toots_oocyst2, age_dat<=0) # initial particles at each location
  
} else if (island==5){
  
  ## Lanai
  ## constant
  toots_con2<-subset(toots_con, island_dat == 5) 
  age_df_con<-subset(toots_con2, age_dat==age_par*sec_to_day) # particles that made it to 90 days
  age0_con<-subset(toots_con2, age_dat<=0) # initial particles at each location
  
  ## oocysts model
  toots_oocyst2<-subset(toots_oocyst, island_dat == 5) 
  age_df_oocyst<-subset(toots_oocyst2, age_dat==age_par*sec_to_day) # particles that made it to 90 days
  age0_oocyst<-subset(toots_oocyst2, age_dat<=0) # initial particles at each location
  
} else if (island==6){
  
  ## Kahoolawe
  ## constant
  toots_con2<-subset(toots_con, island_dat == 6) 
  age_df_con<-subset(toots_con2, age_dat==age_par*sec_to_day) # particles that made it to 90 days
  age0_con<-subset(toots_con2, age_dat<=0) # initial particles at each location
  
  ## oocysts model
  toots_oocyst2<-subset(toots_oocyst, island_dat == 6) 
  age_df_oocyst<-subset(toots_oocyst2, age_dat==age_par*sec_to_day) # particles that made it to 90 days
  age0_oocyst<-subset(toots_oocyst2, age_dat<=0) # initial particles at each location
  
} else if (island ==7) {
  
  ## Maui
  ## constant
  toots_con2<-subset(toots_con, island_dat == 7) 
  age_df_con<-subset(toots_con2, age_dat==age_par*sec_to_day) # particles that made it to 90 days
  age0_con<-subset(toots_con2, age_dat<=0) # initial particles at each location
  
  ## oocysts model
  toots_oocyst2<-subset(toots_oocyst, island_dat == 7) 
  age_df_oocyst<-subset(toots_oocyst2, age_dat==age_par*sec_to_day) # particles that made it to 90 days
  age0_oocyst<-subset(toots_oocyst2, age_dat<=0) # initial particles at each location
  
} else if (island==8){
  
  ## Hawaii/Big Island
  ## constant
  toots_con2<-subset(toots_con, island_dat == 8) 
  age_df_con<-subset(toots_con2, age_dat==age_par*sec_to_day) # particles that made it to 90 days
  age0_con<-subset(toots_con2, age_dat<=0) # initial particles at each location
  
  ## oocysts model
  toots_oocyst2<-subset(toots_oocyst, island_dat == 8) 
  age_df_oocyst<-subset(toots_oocyst2, age_dat==age_par*sec_to_day) # particles that made it to 90 days
  age0_oocyst<-subset(toots_oocyst2, age_dat<=0) # initial particles at each location
  
} else if (island=='MN'){
  
  ## only Maui Nui islands
  toots_con2<-subset(toots_con, island_dat %in% c(4,5,6,7)) 
  age_df_con<-subset(toots_con2, age_dat==age_par*sec_to_day) # particles that made it to 90 days
  age0_con<-subset(toots_con2, age_dat<=0) # initial particles at each location
  
  ## oocysts model
  toots_oocyst2<-subset(toots_oocyst, island_dat %in% c(4,5,6,7)) 
  age_df_oocyst<-subset(toots_oocyst2, age_dat==age_par*sec_to_day) # particles that made it to 90 days
  age0_oocyst<-subset(toots_oocyst2, age_dat<=0) # initial particles at each location
  
} else if (island==1){
  
  ## Niihau
  
  ## constant
  toots_con2<-subset(toots_con, island_dat == 1) 
  age_df_con<-subset(toots_con2, age_dat==age_par*sec_to_day) # particles that made it to 90 days
  age0_con<-subset(toots_con2, age_dat<=0) # initial particles at each location
  
  ## oocysts model
  toots_oocyst2<-subset(toots_oocyst, island_dat == 1) 
  age_df_oocyst<-subset(toots_oocyst2, age_dat==age_par*sec_to_day) # particles that made it to 90 days
  age0_oocyst<-subset(toots_oocyst2, age_dat<=0) # initial particles at each location
  
}

## remove toots_con and toots_oocyst df so it saves space

rm(toots_con, toots_oocyst)

# plot(age_df$lon_dat, age_df$lat_dat, pch=18, cex=0.5, col=rgb(red=1, green=0.6, blue=0.2, alpha=0.2))
# points(age0$lon_dat, age0$lat_dat, pch=18, cex=0.5, col=rgb(red=0, green=0.6, blue=0.2, alpha=0.1))


##########################################################
## AGE0: Calculate number of particles released inside each buffer to use to calculate the end connectivity matrix proportions


##########################################################
## Constant particle release: Create vector of total number of particles released within island buffers

# create empty data frame to store information
tmp_con<-data.frame(
  NI= rep(NA, length.out=nrow(age0_con)),
  KA = rep(NA, length.out=nrow(age0_con)), 
  OA_nwoahu = rep(NA, length.out=nrow(age0_con)), 
  OAdpts_swoahu= rep(NA, length.out=nrow(age0_con)), 
  OA_oahuss = rep(NA, length.out=nrow(age0_con)), 
  OA_kane = rep(NA, length.out=nrow(age0_con)), 
  MN=rep(NA, length.out=nrow(age0_con)), 
  BI=rep(NA, length.out=nrow(age0_con)))


tmp_con[,1]<-1*in.chull(age0_con$lon,age0_con$lat,NIdpts$long,NIdpts$lat) 
sum(tmp_con[,1]) 

tmp_con[,2]<-1*in.chull(age0_con$lon,age0_con$lat,KAdpts$long,KAdpts$lat) 
sum(tmp_con[,2])

tmp_con[,3]<-1*in.chull(age0_con$lon,age0_con$lat,OAdpts_nwoahu$long,OAdpts_nwoahu$lat) 
sum(tmp_con[,3]) 

tmp_con[,4]<-1*in.chull(age0_con$lon,age0_con$lat,OAdpts_swoahu$long,OAdpts_swoahu$lat) 
sum(tmp_con[,4])  

tmp_con[,5]<-1*in.chull(age0_con$lon,age0_con$lat,OAdpts_oahuss$long,OAdpts_oahuss$lat) 
sum(tmp_con[,5])  

tmp_con[,6]<-1*in.chull(age0_con$lon,age0_con$lat,OAdpts_kane$long,OAdpts_kane$lat)
sum(tmp_con[,6])  

tmp_con[,7]<-1*in.chull(age0_con$lon,age0_con$lat,MNdpts$long,MNdpts$lat) 
sum(tmp_con[,7])

tmp_con[,8]<-1*in.chull(age0_con$lon,age0_con$lat,BIdpts$long,BIdpts$lat) 
sum(tmp_con[,8])

# vector of total number of particles released within each buffer locaiton
init_sum_vec_con<-apply(tmp_con, 2, sum)
# sum(init_sum_vec)

##########################################################
## Oocysts model particle release: Create vector of total number of particles released within island buffers

tmp_oocy<-data.frame(
  NI= rep(NA, length.out=nrow(age0_oocyst)),
  KA = rep(NA, length.out=nrow(age0_oocyst)), 
  OA_nwoahu = rep(NA, length.out=nrow(age0_oocyst)), 
  OAdpts_swoahu= rep(NA, length.out=nrow(age0_oocyst)), 
  OA_oahuss = rep(NA, length.out=nrow(age0_oocyst)), 
  OA_kane = rep(NA, length.out=nrow(age0_oocyst)), 
  MN=rep(NA, length.out=nrow(age0_oocyst)), 
  BI=rep(NA, length.out=nrow(age0_oocyst)))


tmp_oocy[,1]<-1*in.chull(age0_oocyst$lon,age0_oocyst$lat,NIdpts$long,NIdpts$lat) 
sum(tmp_oocy[,1]) 

tmp_oocy[,2]<-1*in.chull(age0_oocyst$lon,age0_oocyst$lat,KAdpts$long,KAdpts$lat) 
sum(tmp_oocy[,2])

tmp_oocy[,3]<-1*in.chull(age0_oocyst$lon,age0_oocyst$lat,OAdpts_nwoahu$long,OAdpts_nwoahu$lat) 
sum(tmp_oocy[,3]) 

tmp_oocy[,4]<-1*in.chull(age0_oocyst$lon,age0_oocyst$lat,OAdpts_swoahu$long,OAdpts_swoahu$lat) 
sum(tmp_oocy[,4])  

tmp_oocy[,5]<-1*in.chull(age0_oocyst$lon,age0_oocyst$lat,OAdpts_oahuss$long,OAdpts_oahuss$lat) 
sum(tmp_oocy[,5])  

tmp_oocy[,6]<-1*in.chull(age0_oocyst$lon,age0_oocyst$lat,OAdpts_kane$long,OAdpts_kane$lat)
sum(tmp_oocy[,6])  

tmp_oocy[,7]<-1*in.chull(age0_oocyst$lon,age0_oocyst$lat,MNdpts$long,MNdpts$lat) 
sum(tmp_oocy[,7])

tmp_oocy[,8]<-1*in.chull(age0_oocyst$lon,age0_oocyst$lat,BIdpts$long,BIdpts$lat) 
sum(tmp_oocy[,8])

# vector of total number of particles released within each buffer locaiton
init_sum_vec_oocy<-apply(tmp_oocy, 2, sum)
# sum(init_sum_vec)


##########################################################
## CONSTANT AGE 0: Find points within each buffers polygon of x and y points
# Taking the sum shows the number of times there is true BUT need to divide by the multiplier in front of in.chull to get correct sum


# create empty data frame to store information
tmp_con_init<-data.frame(
  NI= rep(NA, length.out=nrow(age0_con)), 
  KA = rep(NA, length.out=nrow(age0_con)), 
  OA_nwoahu = rep(NA, length.out=nrow(age0_con)), 
  OAdpts_swoahu= rep(NA, length.out=nrow(age0_con)), 
  OA_oahuss = rep(NA, length.out=nrow(age0_con)), 
  OA_kane = rep(NA, length.out=nrow(age0_con)), 
  MN = rep(NA, length.out=nrow(age0_con)), 
  BI = rep(NA, length.out=nrow(age0_con)))


a<-seq(1,8,1) # vector of numbers to represent release buffer


tmp_con_init[,1]<-a[1]*in.chull(age0_con$lon,age0_con$lat,NIdpts$long,NIdpts$lat) 
sum(tmp_con_init[,1]) 

tmp_con_init[,2]<-a[2]*in.chull(age0_con$lon,age0_con$lat,KAdpts$long,KAdpts$lat) 
sum(tmp_con_init[,2])/2 

tmp_con_init[,3]<-a[3]*in.chull(age0_con$lon,age0_con$lat,OAdpts_nwoahu$long,OAdpts_nwoahu$lat) 
sum(tmp_con_init[,3])/3 

tmp_con_init[,4]<-a[4]*in.chull(age0_con$lon,age0_con$lat,OAdpts_swoahu$long,OAdpts_swoahu$lat) 
sum(tmp_con_init[,4])/4  

tmp_con_init[,5]<-a[5]*in.chull(age0_con$lon,age0_con$lat,OAdpts_oahuss$long,OAdpts_oahuss$lat) 
sum(tmp_con_init[,5])/5  

tmp_con_init[,6]<-a[6]*in.chull(age0_con$lon,age0_con$lat,OAdpts_kane$long,OAdpts_kane$lat)
sum(tmp_con_init[,6])/6  

tmp_con_init[,7]<-a[7]*in.chull(age0_con$lon,age0_con$lat,MNdpts$long,MNdpts$lat) 
sum(tmp_con_init[,7])/7

tmp_con_init[,8]<-a[8]*in.chull(age0_con$lon,age0_con$lat,BIdpts$long,BIdpts$lat) 
sum(tmp_con_init[,8])/8


### convert and combine df so that it has the same columns needed to do the connectivity matrix
tmp_con_init[tmp_con_init==0]<-NA
# head(tmp_con_init)

tmp_con_init2<-cbind(age0_con$traj_id, tmp_con_init)
colnames(tmp_con_init2)<-c('traj_id', "NI", "KA", "OA_nwoahu", "OAdpts_swoahu", "OA_oahuss", "OA_kane", "MN", "BI")

test_init_con<-tmp_con_init2 %>% mutate(init_vec = coalesce(NI, KA, OA_nwoahu, OAdpts_swoahu, OA_oahuss, OA_kane, MN, BI))

print(sort(unique(test_init_con$init_vec))) # check to make sure have all buffer zones

## subset out only the ID number and buffer location columns
init_part_con<-test_init_con[,c(1,10)] 
head(init_part_con)

unique(init_part_con$init_vec)
## if there is 
# init_part[which(init_part$init_vec==7),]<-NA


##########################################################
## OOCYSTS MODEL AGE 0: Find points within each buffers polygon of x and y points
# Taking the sum shows the number of times there is true BUT need to divide by the multiplier in front of in.chull to get correct sum


# create empty data frame to store information
tmp_oocy_init<-data.frame(
  NI= rep(NA, length.out=nrow(age0_oocyst)), 
  KA = rep(NA, length.out=nrow(age0_oocyst)), 
  OA_nwoahu = rep(NA, length.out=nrow(age0_oocyst)), 
  OAdpts_swoahu= rep(NA, length.out=nrow(age0_oocyst)), 
  OA_oahuss = rep(NA, length.out=nrow(age0_oocyst)), 
  OA_kane = rep(NA, length.out=nrow(age0_oocyst)), 
  MN = rep(NA, length.out=nrow(age0_oocyst)), 
  BI = rep(NA, length.out=nrow(age0_oocyst)))


a<-seq(1,8,1) # vector of numbers to represent release buffer


tmp_oocy_init[,1]<-a[1]*in.chull(age0_oocyst$lon,age0_oocyst$lat,NIdpts$long,NIdpts$lat) 
sum(tmp_oocy_init[,1]) 

tmp_oocy_init[,2]<-a[2]*in.chull(age0_oocyst$lon,age0_oocyst$lat,KAdpts$long,KAdpts$lat) 
sum(tmp_oocy_init[,2])/2 

tmp_oocy_init[,3]<-a[3]*in.chull(age0_oocyst$lon,age0_oocyst$lat,OAdpts_nwoahu$long,OAdpts_nwoahu$lat) 
sum(tmp_oocy_init[,3])/3 

tmp_oocy_init[,4]<-a[4]*in.chull(age0_oocyst$lon,age0_oocyst$lat,OAdpts_swoahu$long,OAdpts_swoahu$lat) 
sum(tmp_oocy_init[,4])/4  

tmp_oocy_init[,5]<-a[5]*in.chull(age0_oocyst$lon,age0_oocyst$lat,OAdpts_oahuss$long,OAdpts_oahuss$lat) 
sum(tmp_oocy_init[,5])/5  

tmp_oocy_init[,6]<-a[6]*in.chull(age0_oocyst$lon,age0_oocyst$lat,OAdpts_kane$long,OAdpts_kane$lat)
sum(tmp_oocy_init[,6])/6  

tmp_oocy_init[,7]<-a[7]*in.chull(age0_oocyst$lon,age0_oocyst$lat,MNdpts$long,MNdpts$lat) 
sum(tmp_oocy_init[,7])/7

tmp_oocy_init[,8]<-a[8]*in.chull(age0_oocyst$lon,age0_oocyst$lat,BIdpts$long,BIdpts$lat) 
sum(tmp_oocy_init[,8])/8


### convert and combine df so that it has the same columns needed to do the connectivity matrix
tmp_oocy_init[tmp_oocy_init==0]<-NA
# head(tmp_con_init)

tmp_oocy_init2<-cbind(age0_oocyst$traj_id, tmp_oocy_init)
colnames(tmp_oocy_init2)<-c('traj_id', "NI", "KA", "OA_nwoahu", "OAdpts_swoahu", "OA_oahuss", "OA_kane", "MN", "BI")

test_init_oocy<-tmp_oocy_init2 %>% mutate(init_vec = coalesce(NI, KA, OA_nwoahu, OAdpts_swoahu, OA_oahuss, OA_kane, MN, BI))

print(sort(unique(test_init_oocy$init_vec))) # check to make sure have all buffer zones

## subset out only the ID number and buffer location columns
init_part_oocy<-test_init_oocy[,c(1,10)] 
head(init_part_oocy)

unique(init_part_oocy$init_vec)
## if there is 
# init_part[which(init_part$init_vec==7),]<-NA


##########################################################
## CONSTANT AGE N: Find points within each buffers polygon of x and y points
# Taking the sum shows the number of times there is true BUT need to divide by the multiplier in front of in.chull to get correct sum

# create empty data frame to store information
tmp_end_con<-data.frame(
  NI= rep(NA, length.out=nrow(age_df_con)),
  KA = rep(NA, length.out=nrow(age_df_con)), 
  OA_nwoahu = rep(NA, length.out=nrow(age_df_con)), 
  OAdpts_swoahu= rep(NA, length.out=nrow(age_df_con)), 
  OA_oahuss = rep(NA, length.out=nrow(age_df_con)), 
  OA_kane = rep(NA, length.out=nrow(age_df_con)), 
  MN=rep(NA, length.out=nrow(age_df_con)), 
  BI=rep(NA, length.out=nrow(age_df_con)))


a<-seq(1,8,1) # vector of numbers to represent release buffer


tmp_end_con[,1]<-a[1]*in.chull(age_df_con$lon,age_df_con$lat,NIdpts$long,NIdpts$lat) 
sum(tmp_end_con[,1]) 

tmp_end_con[,2]<-a[2]*in.chull(age_df_con$lon,age_df_con$lat,KAdpts$long,KAdpts$lat) 
sum(tmp_end_con[,2])/2 

tmp_end_con[,3]<-a[3]*in.chull(age_df_con$lon,age_df_con$lat,OAdpts_nwoahu$long,OAdpts_nwoahu$lat) 
sum(tmp_end_con[,3])/3 

tmp_end_con[,4]<-a[4]*in.chull(age_df_con$lon,age_df_con$lat,OAdpts_swoahu$long,OAdpts_swoahu$lat) 
sum(tmp_end_con[,4])/4  

tmp_end_con[,5]<-a[5]*in.chull(age_df_con$lon,age_df_con$lat,OAdpts_oahuss$long,OAdpts_oahuss$lat) 
sum(tmp_end_con[,5])/5  

tmp_end_con[,6]<-a[6]*in.chull(age_df_con$lon,age_df_con$lat,OAdpts_kane$long,OAdpts_kane$lat)
sum(tmp_end_con[,6])/6  

tmp_end_con[,7]<-a[7]*in.chull(age_df_con$lon,age_df_con$lat,MNdpts$long,MNdpts$lat) 
sum(tmp_end_con[,7])/7

tmp_end_con[,8]<-a[8]*in.chull(age_df_con$lon,age_df_con$lat,BIdpts$long,BIdpts$lat) 
sum(tmp_end_con[,8])/8


### convert and combine df so that it has the same columns needed to do the connectivity matrix
tmp_end_con[tmp_end_con==0]<-NA
head(tmp_end_con)

tmp_end_con2<-cbind(age_df_con$traj_id, tmp_end_con)
colnames(tmp_end_con2)<-c('traj_id', "NI", "KA", "OA_nwoahu", "OAdpts_swoahu", "OA_oahuss", "OA_kane", "MN", "BI")
head(tmp_end_con2)

test_end_con<-tmp_end_con2 %>% mutate(end_vec = coalesce(NI, KA, OA_nwoahu, OAdpts_swoahu, OA_oahuss, OA_kane, MN, BI)) #  %>% select(a, vec)

head(test_end_con)
print(sort(unique(test_end_con$end_vec))) # check to make sure have all buffer zones

## subset out only the ID number and buffer location columns
end_part_con<-test_end_con[,c(1,10)]
head(end_part_con)


##########################################################
## AGE N: Find points within each buffers polygon of x and y points
# Taking the sum shows the number of times there is true BUT need to divide by the multiplier in front of in.chull to get correct sum


# create empty data frame to store information
tmp_end_oocy<-data.frame(
  NI= rep(NA, length.out=nrow(age_df_oocyst)),
  KA = rep(NA, length.out=nrow(age_df_oocyst)), 
  OA_nwoahu = rep(NA, length.out=nrow(age_df_oocyst)), 
  OAdpts_swoahu= rep(NA, length.out=nrow(age_df_oocyst)), 
  OA_oahuss = rep(NA, length.out=nrow(age_df_oocyst)), 
  OA_kane = rep(NA, length.out=nrow(age_df_oocyst)), 
  MN=rep(NA, length.out=nrow(age_df_oocyst)), 
  BI=rep(NA, length.out=nrow(age_df_oocyst)))


a<-seq(1,8,1) # vector of numbers to represent release buffer


tmp_end_oocy[,1]<-a[1]*in.chull(age_df_oocyst$lon,age_df_oocyst$lat,NIdpts$long,NIdpts$lat) 
sum(tmp_end_oocy[,1]) 

tmp_end_oocy[,2]<-a[2]*in.chull(age_df_oocyst$lon,age_df_oocyst$lat,KAdpts$long,KAdpts$lat) 
sum(tmp_end_oocy[,2])/2 

tmp_end_oocy[,3]<-a[3]*in.chull(age_df_oocyst$lon,age_df_oocyst$lat,OAdpts_nwoahu$long,OAdpts_nwoahu$lat) 
sum(tmp_end_oocy[,3])/3 

tmp_end_oocy[,4]<-a[4]*in.chull(age_df_oocyst$lon,age_df_oocyst$lat,OAdpts_swoahu$long,OAdpts_swoahu$lat) 
sum(tmp_end_oocy[,4])/4  

tmp_end_oocy[,5]<-a[5]*in.chull(age_df_oocyst$lon,age_df_oocyst$lat,OAdpts_oahuss$long,OAdpts_oahuss$lat) 
sum(tmp_end_oocy[,5])/5  

tmp_end_oocy[,6]<-a[6]*in.chull(age_df_oocyst$lon,age_df_oocyst$lat,OAdpts_kane$long,OAdpts_kane$lat)
sum(tmp_end_oocy[,6])/6  

tmp_end_oocy[,7]<-a[7]*in.chull(age_df_oocyst$lon,age_df_oocyst$lat,MNdpts$long,MNdpts$lat) 
sum(tmp_end_oocy[,7])/7

tmp_end_oocy[,8]<-a[8]*in.chull(age_df_oocyst$lon,age_df_oocyst$lat,BIdpts$long,BIdpts$lat) 
sum(tmp_end_oocy[,8])/8


### convert and combine df so that it has the same columns needed to do the connectivity matrix
tmp_end_oocy[tmp_end_oocy==0]<-NA
head(tmp_end_oocy)

tmp_end_oocy2<-cbind(age_df_oocyst$traj_id, tmp_end_oocy)
colnames(tmp_end_oocy2)<-c('traj_id', "NI", "KA", "OA_nwoahu", "OAdpts_swoahu", "OA_oahuss", "OA_kane", "MN", "BI")
head(tmp_end_oocy2)

tmp_end_oocy<-tmp_end_oocy2 %>% mutate(end_vec = coalesce(NI, KA, OA_nwoahu, OAdpts_swoahu, OA_oahuss, OA_kane, MN, BI)) #  %>% select(a, vec)

head(tmp_end_oocy)
print(sort(unique(tmp_end_oocy$end_vec))) # check to make sure have all buffer zones

## subset out only the ID number and buffer location columns
end_part_oocy<-tmp_end_oocy[,c(1,10)]
head(end_part_oocy)



##########################################################
## Merge age0 and age_df together to prep for connectivity matrix analysis for constant and oocysts simulations

all_df_con<-merge(x=init_part_con, y=end_part_con, by=c("traj_id")) 

all_df_oocy<-merge(x=init_part_oocy, y=end_part_oocy, by=c("traj_id")) 

# hist(all_df$end_vec)
# hist(all_df$init_vec, plot=F)

## remove the ID number column to prep for conversion to a table
all_df_con2<-all_df_con[,-1]
# head(all_df_con2)

all_df_oocy2<-all_df_oocy[,-1]
# head(all_df_oocy2)


## convert to a table
tab_all_df_con2<-table(all_df_con2)

tab_all_df_oocy2<-table(all_df_oocy2)

# View(tab_all_df2)
# head(tab_all_df2)
dim(tab_all_df_con2)
dim(tab_all_df_oocy2)



## convert table to a matrix
mat1_con<-as.matrix(tab_all_df_con2) 
dim(mat1_con)

mat1_oocy<-as.matrix(tab_all_df_oocy2) 
dim(mat1_oocy)



##########################################################
# Adding in MISSING rows and columns for for connectivity matrices

## create missing column of 0's
# c8_con<-rep(0, nrow(mat1_con))

# c8_oocy<-rep(0, nrow(mat1_oocy))

## bind missing column to mat1
# mat1_con<-cbind(mat1_con,c8_con)

# mat1_oocy<-cbind(mat1_oocy,c8_oocy)
# 
# ## create missing rows
# r1_con<-rep(0, ncol(mat1_con))
# r2_con<-rep(0, ncol(mat1_con))
# r7_con<-rep(0, ncol(mat1_con))
# r8_con<-rep(0, ncol(mat1_con))
# 
# r1_oocy<-rep(0, ncol(mat1_con))
# r2_oocy<-rep(0, ncol(mat1_con))
# r7_oocy<-rep(0, ncol(mat1_con))
# r8_oocy<-rep(0, ncol(mat1_con))

# 2020 data bind missing rows to mat1
# mat1_con<-rbind(r1_con, r2_con, mat1_con[1:4,],r7_con,r8_con)
# mat1_oocy<-rbind(r1_oocy, r2_oocy, mat1_oocy[1:4,],r7_oocy,r8_oocy)

# for 10 day CM data 
# mat1<-rbind(r1, r2, mat1,r7,r8)
# 
# # rename rownames 
# rownames(mat1_con)<-seq(1,8,1)
# colnames(mat1_con)<-seq(1,8,1)
# 
# rownames(mat1_oocy)<-seq(1,8,1)
# colnames(mat1_oocy)<-seq(1,8,1)

## divide each source location by the total number of particles released by the source grid cell (have to do this above because it is in correct order and I begin to manipulate the matrix below)
# mat_prop<-sweep(mat1,1,init_sum_vec, FUN='/') 
# head(mat_prop)
# dim(mat_prop)

# vec_test<-c(2,2)
# mat_test_sweep<-sweep(mat_test, 2,vec_test, FUN='/')

## reassign any values that are 'Inf' to NA
# mat_prop[mat_prop == 'Inf']<-0

# mat_prop[is.na(mat_prop)]<-0

# vector of release order
grid_order<-seq(1,8,1) 
length(grid_order)

# transpose the matrix
# mat2<-t(mat_prop)
# dev.new()
# image.plot(mat2, main= 'mat2')

# reverse the order of the columns so that image.plot, plots them correctly
# mat3<-mat2[,rev(grid_order)]
# mat3<-mat_prop[,rev(grid_order)]

## normalize connectivity by dividing it by total number of particles released inside of the buffer area
mat_con_sweep<-sweep(mat1_con,1,init_sum_vec_con, FUN='/')

mat_oocy_sweep<-sweep(mat1_oocy,1,init_sum_vec_oocy, FUN='/')


## reorder grid cell so it is plot in the correct order
###  matrix data as a proportions
mat2_con<-mat_con_sweep[,rev(grid_order)]
mat2_oocy<-mat_oocy_sweep[,rev(grid_order)]

###  matrix data as a number
# mat2_con_num<-mat1_con[,rev(grid_order)]
# mat2_oocy_num<-mat1_oocy[,rev(grid_order)]


diff_mat_prop<-(mat2_oocy - mat2_con)*100

# diff_mat_num<-mat2_oocy_num - mat2_con_num



##***NOTE***  
## red means mat2_oocy conn value is greater
## blue means mat2_con conn value is greater


# dev.new()
# image.plot(mat3, main='main3')


# mat1_2<-t(mat1)
# mat1_3<-mat1[,rev(grid_order)]


##########################################################
### Plot connectivity matrix

# zmax<-max(mat3, na.rm=T) # max value of matrix

cols<-rev(c("red", "white" , "blue"))
mypal<-colorRampPalette(cols)(100)

## plot and save connectivity matrix
png(filename = paste0('cmatrix', '_', island_name ,'_', age_par,'days','_', file_con, '.png'),height=8,width=8,units='in', res=200)

# quartz()
brk_dat<-max(abs(diff_mat_prop), na.rm=T)

image.plot(diff_mat_prop, col=mypal, xlab = "Source Locations", ylab= "Sink Locations", xaxt='n', yaxt="n", breaks = seq(-brk_dat, brk_dat, length.out = length(mypal)+1)) #  zlim=c(min(diff_mat_prop, na.rm=T), max(diff_mat_prop, na.rm=T))


dev.off()


##########################################################
### Pau :)


## stop timer
proc.time() - ptm 




