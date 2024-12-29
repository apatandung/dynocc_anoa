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

## FITTING MODEL 
# Constant Model (all parameters ~1) 
fm_null = colext(psiformula= ~ 1, gammaformula = ~ 1, epsilonformula = ~ 1,
                 pformula = ~ 1, umf_anoa)
summary(fm_null)

# Investigate best covariate for p 
# Years as covariate 
fm_1 = colext(psiformula= ~ 1, gammaformula = ~ 1, epsilonformula = ~ 1,
              pformula = ~ year, umf_anoa) 
summary(fm_1)

# distance from threats as covariate 
fm_2 = colext(psiformula= ~ 1, gammaformula = ~ 1, epsilonformula = ~ 1,
              pformula = ~ ancaman, umf_anoa) 
summary(fm_2)

# Effort camera as p covariate 
fm_3 = colext(psiformula= ~ 1, gammaformula = ~ 1, epsilonformula = ~ 1,
              pformula = ~ effort, umf_anoa) 
summary(fm_3)

#Presepitation as covariate p 
fm_4 = colext(psiformula= ~ 1, gammaformula = ~ 1, epsilonformula = ~ 1,
              pformula = ~ presipitasi, umf_anoa)
summary(fm_4)

# Effort+Threats as covariate 
fm_5 = colext(psiformula= ~ 1, gammaformula = ~ 1, epsilonformula = ~ 1,
              pformula = ~ effort+ancaman, umf_anoa) ## Best for p
summary(fm_5)


## ===========Occupancy modelling, insert effort and threats as cov detection ========##
# Elevasi as covariate 
fm_6 = colext(psiformula= ~ elevation, gammaformula = ~ 1, epsilonformula = ~ 1,
              pformula = ~ effort+ancaman, umf_anoa)  ## Best for psi (occupancy)
summary(fm_6)

#Slope as covariate 
fm_7 = colext(psiformula= ~ slope, gammaformula = ~ 1, epsilonformula = ~ 1,
              pformula = ~ effort+ancaman, umf_anoa)
summary(fm_7)

#Sungai sebagai Kovariat 
fm_8 = colext(psiformula= ~ jar_sungai, gammaformula = ~ 1, epsilonformula = ~ 1,
              pformula = ~ effort+ancaman, umf_anoa)
summary(fm_8)

#Jalan sebagai Kovariat  
fm_9= colext(psiformula= ~ jar_jalan, gammaformula = ~ 1, epsilonformula = ~ 1,
             pformula = ~ effort+ancaman, umf_anoa)
summary(fm_9)

# Elevasi+ Jalan sebagai Kovariat 
fm_10 = colext(psiformula= ~ elevation+jar_jalan, gammaformula = ~ 1, epsilonformula = ~ 1,
               pformula = ~ effort+ancaman, umf_anoa) 
summary(fm_10)

# Elevasi+ Sungai sebagai Kovariat 
fm_11 = colext(psiformula= ~ elevation+jar_sungai, gammaformula = ~ 1, epsilonformula = ~ 1,
               pformula = ~ effort+ancaman, umf_anoa) 
summary(fm_11)

#Full 
fm_12 <- colext(psiformula= ~ jar_sungai+slope+jar_jalan+elevation, gammaformula = ~ 1, epsilonformula = ~ 1,
                pformula = ~ effort+ancaman, umf_anoa) 
summary(fm_12)


##=======================PERMODELAN EPSILON(KEPUNAHAN LOKAL) dan GAMMA (KOLONISASI) =================##
## ======Jangan lupa memasukkan "effort+ancaman" dalam pformula dan "elevation" dalam psiformula ====##

# Year as covariate ; kenapa dikurangi(-) 1, karena parameter kolonisasi dan kepunahan lokal diantara tahun pertama dan tahun setelahnya 
# atau dengan kata lain nilainya 2 tahun sekali 
fm_13 = colext(psiformula= ~ elevation, gammaformula = ~ year-1, epsilonformula = ~ year-1,
               pformula = ~ effort+ancaman+year, umf_anoa) 
summary(fm_13)

## Elevasi+Year untuk Kovariat 
fm_14 = colext(psiformula= ~ elevation, gammaformula = ~ elevation+year-1,
               epsilonformula = ~ elevation+year-1,
               pformula = ~ effort+ancaman+year, umf_anoa) # Best 
summary(fm_14)


# Seleksi model 
cand_model = fitList(
  'psi(.)gam(.)eps(.)p(.)' = fm_null,
  'psi(jar_sungai+slope+jar_jalan+elevation)gam(.)eps(.)p(effort+ancaman)' = fm_12, 
  'psi(elevation)gam(year)eps(year)p(effort+ancaman+year)'=fm_13, 
  'psi(elevation)gam(elevation+year)eps(elevation+year)p(effort+ancaman+year)'=fm_14)

model_seleksi = modSel(cand_model)
print(model_seleksi)

# =======3. MODEL ASSESMENT 
# Kurva respon psi(occupancy) terhadap ketinggian 

#Buat sequence/gradian nilai kovariat elevasi 
grad.ele = seq(min(anoa_23$elevation), max(anoa_23$elevation), length=100)

# Buat prediksi hunian awal - prediksi dilakukan secara terpisah untuk setiap gradien
dummyData = data.frame(elevation=grad.ele, year=0)
pred.psi.ele = predict(fm_14, type="psi", newdata=dummyData, appendData=TRUE)

# Kurva respon elevation 
plot(pred.psi.ele[[1]] ~ grad.ele, type = "n", ylim = c(0,1), ylab = "Occupancy probability", xlab = "Elevation [m]") 
polygon(c(grad.ele,rev(grad.ele)), c(pred.psi.ele[,3],rev(pred.psi.ele[,4])), col='grey', border=NA)
lines(pred.psi.ele[[1]] ~ grad.ele, lwd=3, col='blue')

#perhatikan nilai chat, paling baik <=1, opsi kedua <4, jika lebih >4 maka data set tidak bagus 
?mb.gof.test
(gof.results = mb.gof.test(fm_14, nsim=100))
