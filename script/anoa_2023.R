## DYNAMIC OCCUPANCY 
## Species : Anoa (Bubalus depressicornis)
##==============================AUTHORS==================## 
# WCS : Alfons Patandung, Arief Rahman
## BTNBNW : Dini Rahmanita, Hendra D Purnama, Haydin R Faizin, Amal Hidayat, Ronny Mokoginta


#============LOAD LIBRARY
#Aktifkan package yang akan digunakan 
library(dplyr) # for data organization
library(magrittr) # for piping
library(ggplot2) # for plotting
library(unmarked) # for occupancy models
library(AICcmodavg) # For GOF tests
library(raster) # for mapping 

#===========SET WD(Directory Kerja)
setwd("D:\\analisis\\dynocc_anoa\\dynocc_anoa")


## DATA PREPARATION 
# Import CSV file 
anoa_23 = read.csv("anoa_23.csv", header = TRUE)

# Separate detection history matrix as y
y =  as.matrix(anoa_23[, 4:23]) 

#Separate detection covariate and observation covariate 
effort = as.matrix(anoa_23[, 24:43])    # camera effort as obs covariate 
presipitasi = as.matrix(anoa_23[, 55:59])
ancaman = as.matrix(anoa_23[, 50:54])

# Covariate for psi, epsilon and gamma(colonization) 
elevation= anoa_23$elevation
slope = anoa_23$slope
jar_sungai = anoa_23$jar_sungai
jar_mukim = anoa_23$jar_mukim
jar_jalan = anoa_23$jar_jalan
jar_hut = anoa_23$jar_hut

# Years as covariate
years = as.character(2019:2023)
years = matrix(years, nrow(y), 5, byrow=TRUE) ## 5 years

# list detection covariate matrix 
det_list = list(year = as.matrix(years),
                curah_hujan = as.matrix(presipitasi), 
                ancaman=as.matrix(ancaman))

#Make data frame for dynamic/multi season occu using "unmarkedMultFrame"
umf_anoa = unmarkedMultFrame(y=y,
                             siteCovs = data.frame(elevation, jar_jalan, slope, jar_sungai, jar_mukim),
                             obsCovs = list(effort=effort),
                             yearlySiteCovs = det_list,
                             numPrimary=5)
