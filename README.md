### Kriging interpolation using geoR package -----------------------------------------------

#clearing memory space
# Clear the objects from the workspace
rm(list = ls())
gc()
graphics.off()

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load packages
library(geoR)
library(sp)
library(calibrate) # see the number of pairs on semivariogram
library(psych) # descriptive stat
library(Metrics) # calculate RMSE (cross validation)
library(fields)
library(maps)
library(raster)
library(hydroGOF)
library(rgdal)

# ---------------------------------------------------------------------------
# -------------------------- Geostatistical Analysis in geoR -----------------
# ---------------------------------------------------------------------------
# In general, geoR offers much richer possibilities for variogram modeling 
# than gstat.

# Import a table
dados <- read.csv("6 simulaçoes (positivo/pontos_extraidos.csv", sep = ";")

# change values coords 
minX <- 50.5
minY <- 50.5

dados$X <- dados$X - minX
dados$Y <- dados$Y - minY 

## Define soil attribute

soil.attr = "sim6"  # change this argument based in dataframe colnames

## -------------- Convert data frame to geodata object ---------------------

dadosgeo <- as.geodata(data, coords.col = c(1,2), # definir x e y
                     data.col = soil.attr) 


#------------------------ Experimental Variogram -------------------------

var_exp = variog(dadosgeo, uvec=seq(0,1800,l=10)) #lags dependendo da distancia media entre os pontos
plot(var_exp)

### sill =~ variance of the data
var(dadosgeo$data)
ini.raw <- c(0.35,4000) # partial sill and phi (range parameter) estimaded for each soil attribute
nug = 0.2

#Adjusting a model fit with REML
# spherical
var_mod_sph <- likfit(geodata = dadosgeo, ini.cov.pars = ini.raw,
                    nugget = nug, 
                 cov.model = "sph",
                 lik.method = "REML", fix.nugget = F,
                 trend = "1st")
lines(var_mod_sph)
# exponential
var_mod_exp <- likfit(geodata = dadosgeo, ini.cov.pars = ini.raw,
                 #nugget = nug,
                 cov.model = "exp",
                 lik.method = "REML", fix.nugget = F)
lines(var_mod_exp)

# gaussian
var_mod_gau <- likfit(geodata = dadosgeo, ini.cov.pars = ini.raw,
                 #nugget = nug, 
                 cov.model = "gau",
                 lik.method = "REML", fix.nugget = F)
lines(var_mod_gau)

#------------------------------------------------------------------------------------------------------------
#---------------------------- CROSS VALIDATION --------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------

# sph
valid_sph = xvalid(dadosgeo, model = var_mod_sph)
par(mfcol = c(5,2), mar=c(3,3,.5,.5), mgp=c(1.5,.7,0))
plot(valid_sph)
summary(valid_sph)
summary(valid_sph$predicted)

# exp
valid_exp = xvalid(dadosgeo, model = var_mod_exp)
par(mfcol = c(5,2), mar=c(3,3,.5,.5), mgp=c(1.5,.7,0))
plot(valid_exp)
summary(valid_exp)
summary(valid_exp$predicted)

# gau
valid_gau = xvalid(dadosgeo, model = var_mod_gau)
par(mfcol = c(5,2), mar=c(3,3,.5,.5), mgp=c(1.5,.7,0))
plot(valid_gau)
summary(valid_gau)
summary(valid_gau$predicted)

#---------- Cross Validation Stats - Raw Data --------------------------

###RMSE and R2 ----------------------------------------------------------------------
#RMSE 
# Spheric
rmse_sph <- Metrics::rmse(valid_sph$data, valid_sph$predicted)
r2_sph = summary(lm(valid_sph$data ~ valid_sph$predicted))$r.squared  

# EXPONENTIAL
rmse_exp <- Metrics::rmse(valid_exp$data, valid_exp$predicted)
r2_exp = summary(lm(valid_exp$data ~ valid_exp$predicted))$r.squared  

# GAUSSIAN
rmse_gau <- Metrics::rmse(valid_gau$data, valid_gau$predicted)
r2_gau = summary(lm(valid_gau$data ~ valid_gau$predicted))$r.squared  

#####################################    KRIGING      ####################################################

# read boundary
poly <- shapefile("../../boundary/cotorno.shp")
plot(poly)

# convert shp to raster
r = raster(poly, res = 10)

# Input values inside r within boundary dim
rp = rasterize(poly, r, 0)
plot(rp)

# convert raster to DataFrame
grid <- as.data.frame(x = rp, xy = T, na.rm = T)

summary(grid)                 

# change coords values for the same of the dt.geo
# change values coords
grid$x <- grid$x - min(data$x)
grid$y <- grid$y - min(data$y) 

summary(grid)

## collect the prediction grid into a geodata-objec
grid.geodata <- as.geodata(grid, data.col = 3,# ID
                           coords.col = c(1:2)) # coords (x and Y)


### Ordinary kriging

ord.kriging <- krige.conv(geodata = dadosgeo, locations = grid.geodata$coords,
                          krige = krige.control(type.krige = "OK",
                          obj.model = var_mod_exp,
                          trend.d = ~ x + y,
                          trend.l = ~ x + y)) # choose by resume analysis

kriging <- ord.kriging

# input coord, predictions values into k.df
k.df = data.frame(coords = grid[,1:2], predict = kriging$predict) 

colnames(k.df)[3] <- paste(attr, sep = "_") 
# change coords values again

k.df$coords.x <- k.df$coords.x + min(data$x)
k.df$coords.y <- k.df$coords.y + min(data$y)

# rasterize
krig.map <- rasterFromXYZ (k.df) 


## set parameters
crs(krig.map) <- crs(poly)

## plot raster
plot(krig.map, main = paste(attr, "- kriging predictions", sep = " "))

# Export braster as Geotiff

writeRaster(krig.map, "6 simulaçoes (positivo/sim6_geoR_REML",format="GTiff", # SEPARAR POR _
            overwrite = F, NAflag=-9999)
