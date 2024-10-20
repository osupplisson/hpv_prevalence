#This code imports various spatial data

# Packages ----------------------------------------------------------------
source("hpv/code/00_package.R", echo = T)

# Importing regions -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
carto_region <-
  read_sf(
    "global_raw_data/data_for_cartography/data2024/1_DONNEES_LIVRAISON_2024-02-00156/ADE_3-2_SHP_WGS84G_FRA-ED2024-02-15/REGION.shp"
  )  %>%
  clean_names() %>%
  st_as_sf

# Importing districts -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
carto_district <-
  read_sf(
    "global_raw_data/data_for_cartography/data2024/1_DONNEES_LIVRAISON_2024-02-00156/ADE_3-2_SHP_WGS84G_FRA-ED2024-02-15/DEPARTEMENT.shp"
  )  %>%
  clean_names() %>%
  st_as_sf


#Code postal/code insee -------------------------------------------------------
laposte <-
  read_sf("global_raw_data/laposte_hexasmal/laposte_hexasmal.shp") %>%
  clean_names() %>%
  st_as_sf %>%
  rename(insee = code_commun) %>%
  distinct(insee, code_postal)

# Iris code insee -------------------------------------------------------
table_iris <- list.files(
  "global_raw_data/CONTOURS-IRIS_3-0__SHP__FRA_2023-01-01/",
  pattern = ".shp",
  full.names = T,
  all.files = T,
  recursive = T
) %>%
  imap(
    ~ read_sf(.x) %>%
      clean_names() %>%
      st_as_sf() %>%
      rename(insee = insee_com) %>%
      st_transform(crs(carto_district))
  ) %>%
  do.call(rbind, .) %>%
  filter(!st_is_empty(.))

#Code postal/code insee -------------------------------------------------------
data_cp <-
  read_sf("global_raw_data/codes_postaux_V5/codes_postaux_region.shp") %>%
  clean_names() %>%
  st_as_sf %>%
  select(code_postal = id,
         city = lib, 
         geometry)%>%
  mutate(insee_dep = str_sub(code_postal, 1, 2)) %>% 
  st_as_sf %>% 
  st_transform(crs(carto_district))


corsica <- data_cp %>% 
  filter(insee_dep=="20")

corsica_new <- corsica %>% 
  st_centroid() %>% 
  st_join(carto_district %>% select(geometry, insee_dep, nom_dep = nom, insee_reg)) %>% 
  select(-insee_dep.x) %>% 
  rename(insee_dep = insee_dep.y)

corsica <- corsica %>% 
  select(code_postal, geometry) %>% 
  left_join(corsica_new %>% st_drop_geometry(), by = "code_postal") %>% 
  select(code_postal, 
         insee_dep,
         city, geometry)

corsica %>% filter(city == "Calvi")
corsica <- corsica %>% 
  mutate(insee_dep = ifelse(city == "Calvi", "2B", insee_dep))

rm("corsica_new")

data_cp <- data_cp  %>% 
  filter(insee_dep!="20") %>% 
  rbind(corsica) %>% 
  left_join(carto_district %>% st_drop_geometry() %>% select(insee_dep, nom_dep = nom, insee_reg), by = "insee_dep") %>% 
  left_join(carto_region %>% st_drop_geometry() %>% select(insee_reg, nom_region = nom), by = "insee_reg")


corsica %>% 
  filter(is.na(insee_dep))

save.image("hpv/clean_data/spatial_data.rda")
