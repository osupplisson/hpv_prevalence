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


# 
# 
# X <- st_intersects(data_cp)
# 
# A <- data_cp %>% 
#   st_simplify(., preserveTopology = T) %>% 
#   select(id = code_postal, geometry)
# 
# 
# ggplot() + gg(A, fill = 'red', alpha = 0.2)
# 

# plot(st_geometry(idcentroids$geometry), border="blue")
# plot(delaunay_nbobject,  st_centroid(st_geometry(idcentroids$geometry), of_largest_polygon=TRUE), add=TRUE)
#ggplot() + gg(st_as_sf(carto_city$geometry[3])) + gg(st_as_sf(st_centroid(carto_city$geometry[3])))


#
# # Neihborhood matrix ------------------------------------------------------
# id_spatial <- dataset_cp_final %>%
#   distinct(id_intersection, geometry) %>%
#   arrange(id_intersection) %>%
#   mutate(idspace = row_number()) %>%
#   st_as_sf
#
# id_spatial %>%
#   group_by(id_intersection) %>%
#   count() %>%
#   filter(n>1)
#
#
# # Delaunay triangulation neighbours
# delaunay <- tri2nb(st_coordinates(st_centroid(id_spatial$geometry)))
# soi <- graph2nb(soi.graph(delaunay, st_coordinates(st_centroid(id_spatial$geometry))))
# 
# 
# 
# laposte <- laposte %>%
#   filter(insee %in% table_iris$insee)
# 
# #Joining the two dataset
# table_iris_init <- table_iris
# 
# 
# 
# ###List of postcode related to a single insee code-----------------------
# id_within_unique_cp <- laposte %>%
#   st_drop_geometry() %>%
#   group_by(insee) %>%
#   count() %>%
#   ungroup() %>%
#   filter(n == 1) %>%
#   left_join(laposte %>% select(insee, code_postal), by = "insee") %>%
#   select(code_postal, insee)
# 
# cp_contains_unique_or_multiple_id <- laposte %>%
#   st_drop_geometry() %>%
#   group_by(code_postal) %>%
#   count() %>%
#   ungroup() %>%
#   left_join(laposte %>% select(insee, code_postal), by = "code_postal") %>%
#   select(code_postal, insee)
# 
# one_to_one_or_many <- cp_contains_unique_or_multiple_id %>%
#   filter(code_postal %in% id_within_unique_cp$code_postal)
# 
# iris_one_to_one_or_many <- one_to_one_or_many %>%
#   left_join(table_iris %>% select(insee, geometry), by = "insee") %>%
#   st_as_sf %>%
#   st_make_valid() %>%
#   group_by(code_postal) %>%
#   summarise() %>%
#   mutate(id_dataset_cp = row_number())
# 
# # Liste of postcode containing multiple ID and this ID intersect with postcode-------------------------------------------------------------------------
# id_intersecting_multiple_cp <- laposte %>%
#   st_drop_geometry() %>%
#   filter(!insee %in% one_to_one_or_many$insee) %>%
#   filter(!code_postal %in% one_to_one_or_many$code_postal) %>%
#   group_by(insee) %>%
#   count() %>%
#   ungroup() %>%
#   left_join(laposte %>% select(insee, code_postal), by = "insee") %>%
#   select(code_postal, insee)
# 
# 
# iris_many_to_many <- id_intersecting_multiple_cp %>%
#   left_join(table_iris %>% select(insee, geometry), by = "insee") %>%
#   st_as_sf %>%
#   st_make_valid()
# 
# iris_many_to_many_init <- iris_many_to_many
# iris_many_to_many <- iris_many_to_many %>%
#   group_by(code_postal) %>%
#   summarise() %>%
#   mutate(id_carto_city = row_number())
# 
# eval <- 1
# i <- 0
# while (eval != 0) {
#   i <- i + 1
#   print(i)
#   ##Flagging all intersecting geometries
#   intersection <- st_intersects(iris_many_to_many, sparse = T) %>%
#     imap( ~ paste0(.x, collapse = '_')) %>%
#     do.call(rbind, .) %>%
#     tibble(x = .)
#   intersection$id_intersection <- intersection$x[, 1]
#   intersection$x <- NULL
#   
#   #Union of all intersecting geometries
#   union_iris <- iris_many_to_many %>%
#     cbind(intersection) %>%
#     group_by(id_intersection) %>%
#     summarise(geometry = st_union(geometry))
#   
#   
#   #Associate to each id the unioned geometry
#   union_iris <- union_iris %>%
#     separate_wider_delim(
#       id_intersection,
#       delim = "_",
#       names_sep = "_",
#       too_few = "align_start"
#     ) %>%
#     pivot_longer(
#       cols = contains("id_intersection"),
#       names_to = "id",
#       values_to = "values"
#     ) %>%
#     filter(!is.na(values)) %>%
#     mutate(id_carto_city = as.numeric(values)) %>%
#     arrange(id_carto_city) %>%
#     select(id_carto_city, geometry) %>%
#     distinct(id_carto_city, .keep_all = T)
#   
#   # Assigning to each city ID the merged geometry -------------------------
#   iris_many_to_many <- iris_many_to_many %>%
#     st_drop_geometry() %>%
#     select(-one_of(c("id_intersection"))) %>%
#     left_join(union_iris, by = "id_carto_city") %>%
#     add_column(intersection) %>%
#     select(code_postal, id_intersection, geometry)  %>%
#     mutate(id_carto_city = row_number())
#   
#   iris_many_to_many <- iris_many_to_many %>%
#     st_as_sf()
#   
#   iris_many_to_many %>%
#     distinct(id_intersection, geometry) %>%
#     nrow %>%
#     print
#   
#   
#   eval <- st_intersects(iris_many_to_many %>% distinct(id_intersection, geometry),
#                         sparse = T) %>%
#     imap( ~ length(.x)) %>%
#     do.call(rbind, .) %>%
#     tibble(x = .) %>%
#     filter(x > 1) %>%
#     nrow()
#   print(eval)
# }
# 
# 
# ggplot() + gg(iris_many_to_many %>%
#                 filter(as.numeric(str_sub(code_postal, 1, 2)) < 77))
# 
# # Binding -----------------------------------------------------------------
# data_cp <- plyr::rbind.fill(
#   iris_many_to_many,
#   iris_one_to_one_or_many %>%
#     mutate(
#       id_intersection = nrow(iris_many_to_many) + row_number(),
#       id_carto_city = id_intersection
#     ) %>%
#     select(-id_dataset_cp)
# )
# 
# 
# # Adding districts/region -------------------------------------------------
# data_cp <- data_cp %>%
#   mutate(
#     insee_dep = str_sub(code_postal, 1, 2),
#     insee_dep = ifelse(
#       str_sub(code_postal, 1, 2) == "97",
#       str_sub(code_postal, 1, 3),
#       insee_dep
#     )
#   ) %>%
#   left_join(carto_district %>% st_drop_geometry() %>% select(insee_dep, nom, insee_reg),
#             by = "insee_dep") %>%
#   left_join(carto_region %>% st_drop_geometry() %>% select("nom_region" = nom, insee_reg),
#             by = "insee_reg") %>%
#   st_as_sf
# 
# 
# #Handling corsica
# 
# corsica <- data_cp %>%
#   filter(str_sub(code_postal, 1, 2) == 20) %>%
#   select(code_postal, id_intersection, id_carto_city, geometry)
# 
# 
# dpt_cosica <- corsica %>%
#   st_centroid() %>%
#   st_join(carto_district %>% select(insee_dep, nom, geometry, insee_reg)) %>%
#   left_join(carto_region %>% st_drop_geometry() %>% select("nom_region" = nom, insee_reg),
#             by = "insee_reg") %>%
#   st_as_sf
# 
# corsica <- corsica %>%
#   left_join(
#     dpt_cosica %>%
#       select(code_postal, insee_dep, nom, geometry, insee_reg, insee_reg, nom_region) %>% 
#       st_drop_geometry(),
#     by = "code_postal"
#   )
# 
# rm("dpt_cosica")
# data_cp <- data_cp %>%
#   filter(str_sub(code_postal, 1, 2) != 20) %>%
#   rbind(corsica)
# 
# 
# data_cp <- data_cp %>%
#   select(
#     code_postal,
#     "id" = id_carto_city,
#     "id_merged" = id_intersection,
#     nom_region,
#     insee_reg,
#     nom,
#     insee_dep,
#     geometry
#   ) %>%
#   mutate(
#     dpt_numeric = str_sub(code_postal, start = 1, end = 2),
#     dpt_numeric = as.numeric(dpt_numeric),
#     area = ifelse(dpt_numeric < 97, "Metropolitan", "Overseas")
#   )
# 
# data_cp %>% 
#   st_drop_geometry() %>% 
#   group_by(code_postal) %>% 
#   count() %>% 
#   filter(n>1)
# 
# 
# # Adding city name --------------------------------------------------------
# laposte_tmp <-
#   read_sf("global_raw_data/laposte_hexasmal/laposte_hexasmal.shp") %>%
#   clean_names() %>%
#   st_as_sf %>%
#   rename(insee = code_commun) %>% 
#   rename(city = nom_de_la_c) %>% 
#   mutate(city = str_to_sentence(city)) %>% 
#   st_drop_geometry() %>% 
#   distinct(code_postal, .keep_all = T) %>% 
#   select(code_postal, city)
# 
# nrow(data_cp)
# data_cp <- data_cp %>% 
#   left_join(laposte_tmp, by = 'code_postal')
# nrow(data_cp)
# 
# rm("laposte_tmp")