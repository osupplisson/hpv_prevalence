#Cleaning R objects ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rm(list = ls())
gc()
# Package -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(rmapshaper)
library(dbscan)
library(tidybayes)
library(sf)
library(sfheaders)
library(raster)
library(geodata)
library(janitor)
library(tidyverse)
library(readr)
library(spdep)
library(readxl)
library(tidylog)
library(fmesher)
library(ggthemes)
library(INLA)
library(inlabru)


sd_to_prec <- function(sigma) {
  tibble("sd" = sigma,
         "var" = sigma ^ 2,
         "prec" = 1 / sigma ^ 2)
}
prec_to_var <- function(prec) {
  tibble("prec" = prec,
         "var" = 1 / prec,
         "std" = sqrt(1 / prec))
}
set.seed(123)



# CRS ---------------------------------------------------------------------
#https://epsg.io/2154#google_vignette
#https://epsg.io/2154
FRproj <- CRS(SRS_string = "EPSG:2154")
st_crs("EPSG:2154")$proj4string
kmproj <-  CRS(
  "+proj=lcc +lat_0=46.5 +lon_0=3 +lat_1=49 +lat_2=44 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs"
)
