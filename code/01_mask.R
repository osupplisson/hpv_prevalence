##This code creates the mask for France
##This mask is further used for creating maps and make sure that all datapoint are within the spatial domain
# Packages ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source("hpv/code/00_package.R", echo = T)


#First: Keeping only areas more than 100km² to avoid issue during the fit
#Eras: ile de ré et Belle-ile en mer
carto_region <-
  read_sf(
    "global_raw_data/data_for_cartography/data2024/1_DONNEES_LIVRAISON_2024-02-00156/ADE_3-2_SHP_WGS84G_FRA-ED2024-02-15/REGION.shp"
  )  %>%
  clean_names() %>%
  filter(!nom %in% c("Guadeloupe",
                     "Martinique",
                     "Mayotte",
                     "La Réunion",
                     "Guyane"))

##Removing small islands
fun <- function(x){
  x <- x %>% 
    st_transform(crs = kmproj) 
  x <- st_cast(x$geometry, "POLYGON")
  
  carto_drop <- st_as_sf(x) %>% 
    filter(st_area(x)  < units::set_units(200, "km^2"))
  #ggplot()+gg(carto_drop)
  carto_keep <- st_as_sf(x) %>% 
    filter(st_area(x)  >= units::set_units(200, "km^2"))
  
  #Adding oleron and ré
  ile <- carto_drop %>% 
    filter(st_area(carto_drop) > units::set_units(80, "km^2") &
             st_area(carto_drop) < units::set_units(200, "km^2"))
  
  
  id_ile <- ile %>% 
    st_centroid() %>%
    st_coordinates() %>% 
    as_tibble() %>% 
    mutate(id = row_number()) %>% 
    filter(X < 370) %>% 
    filter(Y < 6600) %>% 
    .$id
  
  marseil <- carto_drop %>% 
    filter(st_area(carto_drop) > units::set_units(3, "km^2"))%>% 
    st_centroid() %>% 
    st_coordinates() %>% 
    as_tibble() %>% 
    mutate(id = row_number()) %>% 
    filter(X > 900) %>% 
    .$id
  
  
  
  
  carto_keep <- carto_keep %>% 
    rbind(ile %>% 
            filter(row_number() %in% id_ile))
  
  carto_keep <- st_union(as(carto_keep, "sf"))
  carto_keep <- sfheaders::sf_remove_holes(carto_keep)
  carto_keep <- st_as_sf(carto_keep)
  carto_keep
}

# Import coastline --------------------------------------------------------
coastlines <-
  read_sf(
    "global_raw_data/Europe_coastline_shapefile/Europe_coastline_poly.shp"
  )  %>%
  clean_names() 

fr_mask <- fun(carto_region)
coastlines <- fun(coastlines)
land <- st_union(coastlines, fr_mask)

save.image("hpv/clean_data/masking_france.rda")

