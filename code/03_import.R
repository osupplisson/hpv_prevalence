# Packages ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source("hpv/code/00_package.R", echo = T)

# Import data -------------------------------------------------------------
load("hpv/clean_data/masking_france.rda")
load("hpv/clean_data/spatial_data.rda")

# Importing data ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
raw <- read_excel("hpv/raw_data/Extraction_retravaillée.xlsx") %>%
  clean_names()

table(raw$date_enregistrement)
raw %>%
  filter(str_detect(date_enregistrement, "\\-"))

raw %>%
  filter(str_detect(dateprelev, "\\-"))
# Date for each sampling---------------------------------------------------------------------------------------------------------------------------------------------------------------
raw <- raw %>%
  mutate(
    date_enregistrement = ifelse(
      str_detect(date_enregistrement, "\\-"),
      date_enregistrement,
      paste0(
        str_sub(date_enregistrement, start = 1, end = 4),
        "-",
        str_sub(date_enregistrement, start = 5, end = 6),
        "-",
        str_sub(date_enregistrement, start = 7, end = 8)
      )
    ),
    dateprelev = ifelse(
      str_detect(dateprelev, "\\-"),
      dateprelev,
      paste0(
        str_sub(dateprelev, start = 1, end = 4),
        "-",
        str_sub(dateprelev, start = 5, end = 6),
        "-",
        str_sub(dateprelev, start = 7, end = 8)
      )
    )
  )

raw$dateprelev <- as.Date(raw$dateprelev)
raw$date_enregistrement <- as.Date(raw$date_enregistrement)

# Assigning date ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
raw <- raw %>%
  mutate(
    date =  dplyr::if_else(
      as.character(format(date_enregistrement, format = "%Y")) %in% c("2020", "2021", "2022", "2023"),
      date_enregistrement,
      NA
    ),
    type_date = ifelse(is.na(date), "enregistrement", NA),
    date =  dplyr::if_else(
      is.na(date) &
        as.character(format(dateprelev, format = "%Y")) %in% c("2020", "2021", "2022", "2023"),
      dateprelev,
      date
    ),
    type_date = ifelse(!is.na(date) &
                         is.na(type_date), "enregistrement", type_date),
  )


# Restricting to year -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
raw <- raw  %>%
  filter(!is.na(date)) %>%
  mutate(
    year = as.numeric(format(date, format = "%Y")),
    month = as.numeric(format(date, format = "%m")),
    day = as.numeric(format(date, format = "%d"))
  ) %>%
  filter(year %in% c(2020, 2021, 2022, 2023))


table(raw$type_date, raw$year, useNA = "always")

raw %>%
  group_by(date) %>%
  count() %>%
  ggplot() +
  geom_point(aes(x = date, y = n))

table(raw$type_date, useNA = "always")
table(raw$month, useNA = "always")
raw %>% filter(is.na(date))


raw %>%
  group_by(month, year) %>%
  count() %>%
  view

raw %>% group_by(n_pat) %>% count %>% filter(n > 1) %>% .$n %>% unique

# Rm missing among pri/pro ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
table(raw$year)
table(raw$hppri_45)
colnames(raw$year)
raw <- raw %>%
  dplyr::select(
    pays,
    n_pat,
    age,
    sexe,
    dateprelev,
    date,
    year,
    month,
    day,
    cp_pat,
    hppro_resultat,
    hppri_resultat,
    hppri_bis = hppri_resultats_2020,
    hppro_bis = hpono_resultat
  ) %>%
  filter(!(
    is.na(hppro_resultat) &
      is.na(hppri_resultat) &
      is.na(hppri_bis) &
      is.na(hppro_bis)
  ))

table(raw$year)

table(raw$hppro_resultat)
table(raw$hppri_resultat)
table(raw$hppri_bis)
table(raw$hppro_bis)

raw %>%
  group_by(month, year) %>%
  count() %>% view

# Filter ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
raw <- raw %>%
  filter(
    !str_detect(hppro_resultat, "hpvinh") |
      !str_detect(hppri_resultat, "hpvinh") |
      !str_detect(hppri_bis, "hpvinh") |
      !str_detect(hppro_bis, "hpvinh")
  )

raw <- raw %>%
  filter(sexe != "MASCULIN") %>%
  filter(!is.na(sexe))


# Initial dataset ---------------------------------------------------------
nrow(raw)
raw <- raw %>%
  filter(!is.na(n_pat))  %>%
  filter(!is.na(age)) %>%
  filter(age < 80) %>%
  filter(age >= 15)


nrow(raw)
raw <- raw %>%
  filter(!is.na(cp_pat)) %>%
  mutate(cp_pat = str_remove_all(cp_pat, " ")) %>%
  filter(pays == "FR") %>%
  mutate(dpt = str_sub(cp_pat, start = 1, end = 2)) %>%
  filter(!is.na(as.numeric(dpt))) %>%
  filter(str_length(cp_pat) == 5) %>%
  filter(as.numeric(dpt) %in% seq(1, 97, 1)) %>% 
  filter(!is.na(dpt)) 
nrow(raw)
table(raw$hppro_bis, useNA = "always")
table(raw$hppri_resultat, useNA = "always")
table(raw$year)
table(raw$hppro_resultat, raw$age)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
raw <- raw %>%
  mutate(
    hppri_resultat = case_when(
      is.na(hppri_resultat) &
        !is.na(hppri_bis) ~ hppri_bis,!is.na(hppri_resultat) &
        !is.na(hppri_bis) ~ hppri_resultat,
      TRUE ~ hppri_resultat
    ),
    hppro_resultat = case_when(
      is.na(hppro_resultat) &
        !is.na(hppro_bis) ~ hppro_bis,!is.na(hppro_resultat) &
        !is.na(hppro_bis) ~ hppro_resultat,
      TRUE ~ hppro_resultat
    )
  ) %>%
  select(-hppro_bis, -hppri_bis)


raw %>% filter(is.na(hppri_resultat) & is.na(hppro_resultat))

nrow(raw)
raw <- raw %>%
  filter(!(!is.na(hppro_resultat) & !is.na(hppri_resultat)))
nrow(raw)

# Binding -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
df <- rbind(
  raw %>%
    filter(!is.na(hppro_resultat)) %>%
    select(-starts_with("hppri_")) %>%
    rename_all(~ str_replace(., "hppro_", "hpv_")) %>%
    mutate(type  = "orga"),
  raw %>%
    filter(!is.na(hppri_resultat)) %>%
    rename_all(~ str_replace(., "hppri_", "hpv_")) %>%
    select(-starts_with("hppro_")) %>%
    mutate(type  = "ind")
) %>%
  dplyr::select(n_pat, type, age, cp_pat, date, day, month, year, hpv_resultat)

table(df$type, df$year)
nrow(df)

# Removing organised tests for people aged < 25 ---------------------------------------------------------------------------------------------------------------------------------------------------
df <- df %>%
  filter(!(type == "orga" & age < 25)) %>%
  filter(!(type == "orga" & age > 66))
table(df$age, df$type)



# Unique id per patients --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
df <- df %>%
  left_join(df %>% distinct(n_pat) %>% mutate(id = row_number()), by = 'n_pat') %>%
  dplyr::select(-n_pat)

df %>%
  group_by(id) %>%
  count %>%
  filter(n > 1) %>%
  .$n %>%
  unique


df <- df %>%
  mutate(
    hpv16 = ifelse(
      str_detect(hpv_resultat, "16") |
        str_detect(hpv_resultat, "all"),
      1,
      0
    ),
    hpv18 = ifelse(
      str_detect(hpv_resultat, "18") |
        str_detect(hpv_resultat, "all"),
      1,
      0
    ),
    hpvother = ifelse(
      str_detect(hpv_resultat, "other") |
        str_detect(hpv_resultat, "oth") |
        str_detect(hpv_resultat, "all"),
      1,
      0
    ),
    hpvneg = ifelse(hpv_resultat == "hpvne", 1, 0),
    hpvpos = ifelse(hpv16 == 1 |
                      hpv18 == 1 |
                      hpvother == 1, 1, 0),
    N = 1
  )  %>%
  filter(!(hpvpos == 1 & hpvneg == 1))


table(df$type, df$age)


# Keeping the last test for each individual -------------------------------------------------------------------------------------------------------------------------------------------------------
df <- df %>%
  group_by(id) %>%
  arrange(day, month, year) %>%
  filter(row_number() == max(row_number())) %>%
  ungroup()


# Starting to deal with geometries ----------------------------------------
#Keeping mapping
data_cp_obs <- data_cp

data_cp_obs <- data_cp_obs  %>%
  #Keeping only observed CP or CP in overseas
  filter(code_postal %in% df$cp_pat) %>%
  st_centroid() %>%
  select(code_postal, insee_dep, nom = nom_dep, insee_reg, region_name = nom_region, geometry)


#Removing points for which we could not identify postcode
cp_dom_tom <- laposte %>% filter(str_sub(code_postal, 1,2) == "97") %>% .$code_postal
df <- df %>%
  #Removing all datapoints to which matching with CP failed
  filter(cp_pat %in% data_cp_obs$code_postal | cp_pat %in% cp_dom_tom)




# Splitting DOM/Mainland --------------------------------------------------
df_dom_tom <- df %>%
  rename(code_postal = cp_pat) %>% 
  filter(as.numeric(str_sub(code_postal, 1, 2)) %in% c(97)) %>% 
  mutate(dpt = str_sub(code_postal, 1, 3)) %>% 
  left_join(carto_district %>% st_drop_geometry() %>% select(dpt = insee_dep, nom_dep = nom, insee_reg), by = "dpt") %>%  
  left_join(carto_region %>% st_drop_geometry() %>% select(insee_reg, nom_region = nom), by = "insee_reg") %>% 
  st_drop_geometry()


df <- df %>%
  filter(!as.numeric(str_sub(cp_pat, 1, 2)) %in% c(97))

#Adding information on districts
df <- df %>%
  rename(code_postal = cp_pat) %>%
  left_join(data_cp_obs, by = "code_postal")




carto_district <- carto_district %>%
  filter(str_length(insee_dep) == 2)
#Saving crs
crs_init <- st_crs(df$geometry)



#Changing the projection to meters using the French grid
df <- df %>%
  st_as_sf %>% 
  st_transform(FRproj)

# Removing all datapoints located outside France
#https://rpubs.com/jafet089/880416 cc <- country_codes()

#Removing points outside of the mask
fr_mask <- as(fr_mask, "sf") %>%
  st_transform(FRproj)

coastlines <- as(coastlines, "sf") %>%
  st_transform(FRproj)

land <- as(land, "sf") %>%
  st_transform(FRproj)

crs(fr_mask) == crs(df)
crs(coastlines) == crs(df)


df <- df %>%
  mutate(idBIS = row_number())


# Assigning centroids  ----------------------------------------------------
df <- df %>% 
  st_centroid()

dropped_point <- df[lengths(st_intersects(df, land)) == 0, ]
sum(dropped_point$N)
df <- df %>%
  filter(!idBIS %in% dropped_point$idBIS) %>%
  select(-idBIS)

nrow(df)

###Transforming to coordinates
#longitude and latitude : not a valid coordinate system for spatial modeling, as distances in longitude and distances in latitude are in different scales.
#Modeling should be done in a CRS (Coordinate Reference System) where units along each axis is measured in e.g. kilometers.
#We do not want to use meters as this would result in very large values along the axes, and could cause unstable numerical results.
# Define km projection
fr_mask <- fr_mask %>%
  st_transform(kmproj)

coastlines <- coastlines %>%
  st_transform(kmproj)

land <- land %>%
  st_transform(kmproj)
#Compute boundaries for inlabru
france.bdry <- inla.sp2segment(fr_mask)
df <- df %>%  st_transform(kmproj)


# Adding case -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
case <- df %>%
  st_drop_geometry() %>%
  mutate(case = paste0(hpv16, hpv18, hpvother), idrow = row_number())  %>%
  select(idrow, case) %>%
  mutate(value = 1) %>%
  pivot_wider(
    id_cols = idrow,
    names_from = "case",
    values_from = value,
    values_fill = 0
  ) %>%
  rename_all(function(.x)
    ifelse(.x != "idrow", paste0("case", .x), .x))

df <- df %>%
  mutate(idrow = row_number()) %>%
  left_join(case, by = "idrow")


# Aggregating by strata ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mindate <- min(df$date)
maxdate <- max(df$date)
dfdate <- tibble(date = seq(mindate, maxdate, "day")) %>%
  mutate(week = format(date, '%Y-%U')) %>%
  distinct(week) %>%
  mutate(idtime = row_number())


df <- df %>%
  mutate(week = format(date, '%Y-%U')) %>%
  dplyr::select(-day, -month, -id, -hpv_resultat, -idrow) %>%
  left_join(dfdate, by = "week")


# Aggregation -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
data_bru <- df %>%
  group_by_at(all_of(
    colnames(df) %>% tibble(x = .) %>% filter(!str_detect(x, "hpv") &
                                                x != "N" &
                                                !str_detect(x, "case")) %>% .$x
  )) %>%
  summarise_at(vars(
    colnames(df) %>% tibble(x = .) %>% filter(str_detect(x, "hpv") |
                                                x == "N" |
                                                str_detect(x, "case")) %>% .$x
  ), ~ sum(.)) %>%
  ungroup()

#Changing age class --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
data_bru <- data_bru %>%
  mutate(age_class = cut(
    age,
    c(min(age), 24, seq(29, 60, 10), 66, max(age)),
    ordered = T,
    include.lowest = T
  ))


# Mask city ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
urban_unit <-
  read_sf("global_raw_data/fonds_uu2020_2024/uu2020_2024/uu2020_2024.shp") %>%
  clean_names()  %>%
  st_transform(crs_init)


urban_unit <-
  st_intersection(urban_unit, st_transform(as(fr_mask, "sf"), crs_init))

urban_unit <- urban_unit %>%
  st_transform(kmproj)


rm("case")
rm("centroids")
rm('raw')
rm("tmp")
rm("df")
rm("exporting_results_model")
rm("pp_check_function")
sum(data_bru$N)


# Saving ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
save.image("hpv/clean_data/aggregated_data.rda")
load("hpv/clean_data/aggregated_data.rda")

#Unique ID of centroids
data_cp_metropolitan <- data_cp %>%
  st_as_sf() %>%
  st_transform(kmproj)

dropped_point <- data_cp_metropolitan[lengths(st_intersects(data_cp_metropolitan, land)) == 0, ]
nrow(dropped_point)

data_cp_final <- data_cp_metropolitan %>%
  filter(!code_postal %in% dropped_point$code_postal) %>% 
  mutate(idspace = row_number())

data_bru <- data_bru %>%
  left_join(data_cp_final %>%
              st_drop_geometry() %>%
              select(code_postal, idspace),
            by = 'code_postal')

# # Delaunay triangulation neighbours
delaunay <- tri2nb(st_coordinates(st_centroid(data_cp_final$geometry)))
soi <- graph2nb(soi.graph(delaunay, st_coordinates(st_centroid(data_cp_final$geometry))))
queenfirstorder <- poly2nb(as(data_cp_final$geometry, "Spatial"),queen = TRUE)
queensecondorder <- spdep::nblag(neighbours = queenfirstorder, maxlag = 2)
queensecondorder <- nblag_cumul(queensecondorder)
# Saving ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
save.image("hpv/clean_data/aggregated_data.rda")

#ggplot() + gg(st_transform(st_as_sf(carto_region), "lambert_globe"))
# 
# 
# ggplot() + 
#   gg(idcentroidspred %>% 
#        mutate(obs = ifelse(code_postal %in% data_bru$code_postal, "A", "B")), 
#      aes(fill = obs), alpha = 0.2)+
#   scale_fill_manual(values = c("red", "blue"))
# 
# 
# A <- idcentroidspred %>% 
#   mutate(obs = ifelse(code_postal %in% data_bru$code_postal, "A", "B")) %>% 
#   distinct(id_merged,.keep_all = T)
# 
# A <- data_cp %>% 
#   st_simplify(., preserveTopology = T, dTolerance = 1)
# 
# 
# B <- A%>% 
#   filter(as.numeric(insee_dep)==49) %>% 
#   distinct(id_merged, .keep_all = T) %>% 
#   filter(row_number()%in% c(4,14))
# 
# ggplot() + 
#   gg(B, fill = 'red', alpha = 0.2)
# 
# 
# 
# A <- A %>% 
#   group_by(city) %>% 
#   summarise()
# 
# st_overlaps(A, A)
