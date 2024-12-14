# Packages ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source("hpv/code/00_package.R", echo = T)
# Read aggregated data ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load("hpv/clean_data/aggregated_data.rda")

# Read functions ----------------------------------------------------------
source("hpv/code/99_post_fit_function.R", echo = T)

# Delaunay matrix ---------------------------------------------------------
delaunay_w <- delaunay
soi_w <- soi
queen1_w <- queenfirstorder
soi <- as(nb2mat(soi, style = 'B', zero.policy = TRUE), 'Matrix')
delaunay <- as(nb2mat(delaunay, style = 'B', zero.policy = TRUE), 'Matrix')
queenfirstorder <- as(nb2mat(queenfirstorder, style = 'B', zero.policy = TRUE),
                      'Matrix')
queensecondorder <- as(nb2mat(queensecondorder, style = 'B', zero.policy = TRUE),
                       'Matrix')



# plot(st_geometry(fr_mask), border="grey")
# plot(delaunay_w, st_coordinates(st_centroid(idcentroids_unique$geometry)), add=TRUE)
# plot(st_geometry(fr_mask), border="grey")
# plot(soi_w, st_coordinates(st_centroid(idcentroids_unique$geometry)), add=TRUE)

# utils -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sd_to_prec <- function(sigma) {
  tibble(
    "sd" = sigma,
    "var" = sigma ^ 2,
    "prec" = 1 / sigma ^ 2,
    "0025quantile" = qnorm(0.025, 0, sigma),
    "0975quantile" = qnorm(0.975, 0, sigma)
  )
}


# Filtering out <30 >66 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sum(data_bru$N)
data_bru %>%
  filter(!(age >= 30 & age < 67)) %>%
  .$N %>%
  sum

data_bru <- data_bru %>%
  filter(age >= 30 & age < 67)

sum(data_bru$N)

data_bru <- data_bru %>%
  mutate(hpv1618 = N - case000 - case001) %>%
  st_as_sf() %>%
  st_transform(kmproj)
print(sum(data_bru$hpv1618) / sum(data_bru$N))

dfidage <-
  tibble(age = seq(min(data_bru$age), max(data_bru$age), 1)) %>%
  arrange(age) %>%
  mutate(idage = row_number())

# Building mesh -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# https://rpubs.com/jafet089/886687
# https://onlinelibrary.wiley.com/doi/epdf/10.1002/ece3.4789
# https://github.com/inlabru-org/inlabru/issues/57
# https://rpubs.com/jafet089/880416
if (!file.exists("hpv/clean_data/output/mesh.RDS")) {
  locx1 <- data_bru %>%
    st_centroid %>%
    dplyr::select(code_postal, geometry)
  
  locx2 <- data_cp_final %>%
    st_centroid %>%
    filter(!code_postal %in% locx1$code_postal)
  
  locx <- plyr::rbind.fill(locx1, locx2) %>%
    distinct(code_postal, .keep_all = T)
  
  mesh_input <- fmesher::fm_mesh_2d_inla(
    boundary = france.bdry,
    max.edge = c(5, 10),
    cutoff = 5,
    offset = c(5, 10),
    loc = locx$geometry,
    crs = fmesher::fm_crs(kmproj),
    min.angle = 28
  )
  print(mesh_input$n)
  saveRDS(mesh_input, "hpv/clean_data/output/mesh.RDS")
} else{
  mesh_input <- readRDS("hpv/clean_data/output/mesh.RDS")
}
range_space <- round(0.5 * diff(range(mesh_input$loc[, 2])))
print(range_space)
print(mesh_input$n)



# print(coarser_mesh$n)
# coarser_mesh$n*max(dfdate$idtime)
# https://becarioprecario.bitbucket.io/spde-gitbook/ch-nonstationarity.html#ch:barrier

#ggplot() + gg(coastlines)
#fr_mask
water.tri <- inla.over_sp_mesh(land, y = mesh_input, type = "centroid")
num.tri <- length(mesh_input$graph$tv[, 1])
barrier.tri <- setdiff(1:num.tri, water.tri)
poly.barrier <- inla.barrier.polygon(mesh_input, barrier.triangles = barrier.tri)
poly.barrier <- as(poly.barrier, "sf") %>% st_set_crs(kmproj)


# Mesh --------------------------------------------------------------------
plot_mesh.barrier <- ggplot() +
  gg(mesh_input, fill = "black")

if (!file.exists("hpv/clean_data/output/post_mesh.RDS")) {
  plot_mesh.barrier <- ggplot() +
    gg(coastlines, fill = "grey", alpha = 0.25) +
    gg(mesh_input, fill = "black") +
    gg(poly.barrier, fill = "red", alpha = 0.25) +
    gg(
      data = data_bru  %>% distinct(code_postal, .keep_all = T),
      col = 'darkgreen',
      alpha = 0.25
    ) +
    labs(fill = element_blank(), caption = "Blue lines indicate France borders. \n Red shaded area indicates Oceans and Seas, acting as natural barrier.") +
    theme_map() +
    lims(
      x = layer_scales(plot_mesh.barrier)$x$range$range,
      y = layer_scales(plot_mesh.barrier)$y$range$range
    )
  
  saveRDS(plot_mesh.barrier, file = "hpv/clean_data/output/post_mesh.RDS")
}

# Pivoting ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
data_bru <- data_bru %>%
  mutate(
    intercept = 1,
    opportunistic = ifelse(type == "ind", 1, 0),
    testgroup = opportunistic + 1
  )


dfbinomial <- data_bru %>%
  pivot_longer(c(hpv1618, hpvother),
               names_to = "virus",
               values_to = "result") %>%
  mutate(
    id = row_number(),
    hpv1618 = ifelse(virus == "hpv1618", 1, 0),
    hpvother = ifelse(virus == "hpvother", 1, 0),
    virusgroup = hpvother + 1,
    intercept = 1
  ) %>%
  mutate(
    intercept_level = case_when(
      virusgroup == 0 & opportunistic == 0 ~ 1,
      virusgroup == 1 & opportunistic == 0 ~ 2,
      virusgroup == 0 & opportunistic == 1 ~ 3,
      virusgroup == 1 & opportunistic == 1 ~ 4
    )
  )

# Save --------------------------------------------------------------------
saveRDS(dfbinomial, file = "hpv/clean_data/output/df_for_fit.RDS")

# Components ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# https://inla.r-inla-download.org/r-inla.org/doc/latent/copy.pdf
function_bru_elements <- function(prec_fixed,
                                  sd_pc_input,
                                  prec_corr,
                                  range.fraction.input,
                                  typespace) {
  fx <- function(x) {
    exp(x) / (1 + exp(x))
  }
  ggplot(tibble(x = fx(rnorm(
    10000, 0, prec_to_var(0.4)$std
  ))), aes(x)) +
    geom_density(colour = "blue")
  
  # Mesh for space ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  if (str_detect(typespace, "barrier")) {
    spdespace <- INLAspacetime::barrierModel.define(
      mesh = mesh_input,
      barrier.triangles = barrier.tri,
      prior.sigma = c(sd_pc_input, 0.01),
      prior.range = c(range_space, 0.99),
      range.fraction = range.fraction.input,
      constr = FALSE
    )
  } else if (str_detect(typespace, "stationary")) {
    spdespace <- inla.spde2.pcmatern(
      mesh_input,
      prior.sigma = c(sd_pc_input, 0.01),
      prior.range = c(range_space, 0.99),
      constr = F
    )
  } else if (str_detect(typespace, "anisotropicrspde")) {
    #Very long and difficult to fit, not possible given our computing capabilities
    spdespace <- rSPDE::rspde.anistropic2d(mesh = mesh_input)
  }
  
  # https://stats.stackexchange.com/questions/454647/relation-between-gaussian-processes-and-gaussian-markov-random-fields
  fx <- function(x, n) {
    (exp(x) - 1) / (exp(x) + n - 1)
  }
  ggplot(tibble(x = fx(rnorm(
    10000, 0, prec_to_var(0.2)$std
  ), 2)), aes(x)) +
    geom_density(colour = "blue")
  
  
  ggplot(tibble(x = rnorm(10000, 0, prec_to_var(2)$std))) +
    geom_function(
      fun = function(x) {
        (exp(x) - 1) / (exp(x) + 4 - 1)
      },
      colour = "red"
    )
  
  
  fxpos <- function(x, n) {
    exp(x) / (1 + exp(x))
  }
  ggplot(tibble(x = plogis(inla.pc.rcor0(1000000, 0.5, 0.5))), aes(x)) +
    geom_density(colour = "blue")
  
  
  # Components ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Baseline component ------------------------------------------------------
  cmp <- ~ 0 +
    intercept(
      intercept,
      model = "linear",
      mean.linear = 0,
      prec.linear = prec_fixed
    ) +
    intercepthpvother(virusgroup, model = "factor_contrast", hyper = list(prec = list(
      prior = 'pc.prec', param = c(u = 1, a = 0.5)
    ))) +
    interceptopportunistic(opportunistic,
                           model = "factor_contrast",
                           hyper = list(prec = list(
                             prior = 'pc.prec', param = c(u = 1, a = 0.5)
                           ))) +
    intercepthpvotheropportunistic(opportunistic * (virusgroup - 1),
                                   model = "factor_contrast",
                                   hyper = list(prec = list(
                                     prior = 'pc.prec', param = c(u = 1, a = 0.5)
                                   )))
  
  
  # Age ---------------------------------------------------------------------
  cmp <- update(
    cmp,
    ~ . + fage(
      (age - 29),
      group = virusgroup,
      model = "rw2",
      ngroup = 2,
      scale.model = TRUE,
      constr = TRUE,
      hyper = list(prec = list(
        prior = 'pc.prec', param = c(u = sd_pc_input, a = 0.01)
      )),
      control.group = list(model = "exchangeable", hyper = list(rho = list(
        prior = "normal", param = c(0, prec_corr)
      )))
    ) + fageopportunistic(
      (age - 29),
      group = virusgroup,
      weights = opportunistic,
      model = "rw2",
      ngroup = 2,
      scale.model = TRUE,
      constr = TRUE,
      hyper = list(prec = list(
        prior = 'pc.prec', param = c(u = sd_pc_input, a = 0.01)
      )),
      control.group = list(model = "exchangeable", hyper = list(rho = list(
        prior = "normal", param = c(0, prec_corr)
      )))
    )
  )
  
  
  
  # Space and time ----------------------------------------------------------
  if (str_detect(typespace, "spacetime")) {
    #Does not work, INLA fails
    cmp <- update(
      cmp,
      ~ . +
        fspacetime(
          idspace,
          group = idtime,
          replicate = virusgroup,
          model = "bym2",
          graph = W,
          scale.model = TRUE,
          adjust.for.con.comp = TRUE,
          constr = TRUE,
          hyper = list(
            prec = list(
              prior = 'pc.prec',
              param = c(u = sd_pc_input, a = 0.01)
            ),
            phi = list(prior = 'pc', param =  c(u = 0.5, a = 0.5))
          ),
          control.group = list(model = "ar1", hyper = list(rho = list(
            prior = 'pc.cor0', param = c(0.5, 0.5)
          )))
        ) +
        fspacetimeopportunistic(
          idspace,
          weights = opportunistic,
          group = idtime,
          replicate = virusgroup,
          model = "bym2",
          graph = W,
          scale.model = TRUE,
          adjust.for.con.comp = TRUE,
          constr = TRUE,
          hyper = list(
            prec = list(
              prior = 'pc.prec',
              param = c(u = sd_pc_input, a = 0.01)
            ),
            phi = list(prior = 'pc', param =  c(u = 0.5, a = 0.5))
          ),
          control.group = list(model = "ar1", hyper = list(rho = list(
            prior = 'pc.cor0', param = c(0.5, 0.5)
          )))
        )
    )
  }
  else{
    ###Space
    if (str_detect(typespace, "bym2")) {
      #GMRF component, model 3 to 6
      cmp <- update(
        cmp,
        ~ . +
          fspace(
            idspace,
            group = virusgroup,
            ngroup = 2,
            model = "bym2",
            graph = W,
            scale.model = TRUE,
            adjust.for.con.comp = TRUE,
            constr = TRUE,
            hyper = list(
              prec = list(
                prior = 'pc.prec',
                param = c(u = sd_pc_input, a = 0.01)
              ),
              phi = list(prior = 'pc', param =  c(u = 0.5, a = 0.5))
            ),
            control.group = list(model = "exchangeable", hyper = list(rho = list(
              prior = "normal", param = c(0, prec_corr)
            )))
          ) +
          fspaceopportunistic(
            idspace,
            weights = opportunistic,
            group = virusgroup,
            ngroup = 2,
            model = "bym2",
            graph = W,
            scale.model = TRUE,
            adjust.for.con.comp = TRUE,
            constr = TRUE,
            hyper = list(
              prec = list(
                prior = 'pc.prec',
                param = c(u = sd_pc_input, a = 0.01)
              ),
              phi = list(prior = 'pc', param =  c(u = 0.5, a = 0.5))
            ),
            control.group = list(model = "exchangeable", hyper = list(rho = list(
              prior = "normal", param = c(0, prec_corr)
            )))
          )
      )
    }
    else{
      #GRF component, model 1 and 2
      cmp <- update(
        cmp,
        ~ . +
          fspace(
            geometry,
            model = spdespace,
            group = virusgroup,
            ngroup = 2,
            control.group = list(model = "exchangeable", hyper = list(rho = list(
              prior = "normal", param = c(0, prec_corr)
            )))
          ) +
          fspaceopportunistic(
            geometry,
            weights = opportunistic,
            model = spdespace,
            group = virusgroup,
            ngroup = 2,
            control.group = list(model = "exchangeable", hyper = list(rho = list(
              prior = "normal", param = c(0, prec_corr)
            )))
          )
      )
    }
    ##Time
    cmp <- update(
      cmp,
      ~ . + ftime(
        idtime,
        group = virusgroup,
        model = "rw2",
        ngroup = 2,
        scale.model = TRUE,
        constr = TRUE,
        hyper = list(prec = list(
          prior = 'pc.prec',
          param = c(u = sd_pc_input, a = 0.01)
        )),
        control.group = list(model = "exchangeable", hyper = list(rho = list(
          prior = "normal", param = c(0, prec_corr)
        )))
      ) + ftimeopportunistic(
        idtime,
        group = virusgroup,
        weights = opportunistic,
        model = "rw2",
        ngroup = 2,
        scale.model = TRUE,
        constr = TRUE,
        hyper = list(prec = list(
          prior = 'pc.prec',
          param = c(u = sd_pc_input, a = 0.01)
        )),
        control.group = list(model = "exchangeable", hyper = list(rho = list(
          prior = "normal", param = c(0, prec_corr)
        )))
      )
    )
  }
  
  
  
  # Formula -----------------------------------------------------------------
  fml <- result ~ .
  
  # Likelihood ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  lik1 <- like(
    data = dfbinomial,
    family = "binomial",
    Ntrials = N,
    formula = fml
  )
  
  
  
  print(cmp)
  print(lik1$expr)
  
  # Out ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  out <- list("cmp" = cmp, "lik1" = lik1)
  return(out)
}



# Importing mapping for regions -------------------------------------------------------------------------------------------------------------------------------------------------------------------
carto_region <-
  read_sf(
    "global_raw_data/data_for_cartography/data2024/1_DONNEES_LIVRAISON_2024-02-00156/ADE_3-2_SHP_WGS84G_FRA-ED2024-02-15/REGION.shp"
  )  %>%
  clean_names()


# Building pixels -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
dim_grid <- c(500, 500)
#dim_grid <- c(10, 10)

# https://inlabru-org.github.io/inlabru/reference/predict.bru.html
### Using the kmproj CRS take forever, should not be used
### Instead: transform both mesh and mask with new crs
### Then transform the object again

mask_for_sampling <- data_cp_final %>%
  filter(!is.na(idspace)) %>%
  distinct(idspace, .keep_all = T) %>%
  st_transform(crs_init)


suppressWarnings({
  pixel_pred <- fm_pixels(
    fm_transform(mesh_input, crs_init),
    format = "sf",
    mask = mask_for_sampling,
    minimal = TRUE,
    dims = dim_grid
  ) %>%
    st_as_sf() %>%
    st_transform(kmproj)
})
nrow(pixel_pred)

pixel_pred <- pixel_pred %>%
  st_join(mask_for_sampling %>%
            st_as_sf() %>%
            st_transform(kmproj))

#ggplot() + gg(st_transform(as(fr_mask, "sf")), crs_init)+gg(mask_city, color = 'black')
# Pixel pred for Paris region ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
dim_grid <- c(300, 300)
#dim_grid <- c(10, 10)

suppressWarnings({
  pixel_pred_paris <- fm_pixels(
    fm_transform(mesh_input, crs_init),
    format = "sf",
    mask = mask_for_sampling,
    minimal = TRUE,
    dims = dim_grid,
    xlim = st_bbox(
      carto_region %>%
        filter(str_detect(nom, "le-de-France")) %>%
        st_transform(crs_init) %>%
        .$geometry
    )[c(1, 3)],
    ylim = st_bbox(
      carto_region %>%
        filter(str_detect(nom, "le-de-France")) %>%
        st_transform(crs_init) %>%
        .$geometry
    )[c(2, 4)]
  ) %>%
    st_as_sf() %>%
    st_transform(kmproj)
})

nrow(pixel_pred_paris)
pixel_pred_paris <- pixel_pred_paris  %>%
  st_join(mask_for_sampling %>%
            st_as_sf() %>%
            st_transform(kmproj))

pixel_pred_paris <- pixel_pred_paris %>%
  tidylog::filter(str_detect(nom_region, "le-de-France"))

pixel_pred <- pixel_pred %>%
  rbind(pixel_pred_paris) %>%
  distinct(geometry, .keep_all = T) %>%
  mutate(idgeometry = row_number())
nrow(pixel_pred)


# Centroids of major French cities ----------------------------------------------------------------------------------------------------------------------------------------------------------------
centroid_major_city <- data_cp_final %>%
  filter(
    code_postal %in% c(
      "75001",
      "69001",
      "13001",
      "59000",
      "33000",
      "06000",
      "35000",
      "31000",
      "34000",
      "44000",
      "67000"
    )
  ) %>%
  st_transform(kmproj) %>%
  dplyr::select(
    code_postal,
    geometry,
    idspace,
    insee_dep,
    insee_reg,
    "district" = nom_dep,
    region = nom_region,
    city
  ) %>%
  st_centroid() %>%
  mutate(city = ifelse(city == "Marseille-I", "Marseille", city))


table(centroid_major_city$city)

ggplot() + gg(fr_mask) + gg(centroid_major_city)


saveRDS(centroid_major_city, file = "hpv/clean_data/output/centroid_major_city.RDS")


# Modifying dataset of centroids to include additional spatial var --------
idcentroidspred <- data_cp_final %>%
  filter(!is.na(idspace)) %>%
  distinct(idspace, .keep_all = T)

idcentroidspred$geometry <- ms_simplify(idcentroidspred$geometry,
                                        keep = 0.05,
                                        keep_shapes = T)
idcentroidspred <- idcentroidspred %>%
  st_as_sf

idcentroidspred <- idcentroidspred %>%
  rename(district = nom_dep, region = nom_region) %>%
  select(code_postal,
         city,
         idspace,
         insee_dep,
         district,
         insee_reg,
         region,
         geometry) %>%
  mutate(idgeometry = row_number())

pixel_pred <- pixel_pred %>%
  rename(district = nom_dep, region = nom_region) %>%
  select(code_postal,
         city,
         idspace,
         insee_dep,
         district,
         insee_reg,
         region,
         geometry) %>%
  mutate(idgeometry = row_number())

centroid_major_city <- centroid_major_city %>%
  select(code_postal,
         city,
         idspace,
         insee_dep,
         district,
         insee_reg,
         region,
         geometry) %>%
  mutate(idgeometry = row_number())

print(colnames(centroid_major_city))
print(colnames(idcentroidspred))
print(colnames(pixel_pred))
# Two first models ---------------------------------------------------------

# BRU options -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
input_bru_options <- bru_options(
  bru_verbose = 4,
  verbose = T,
  bru_initial = NULL,
  control.compute = list(waic = TRUE, dic = TRUE, config = TRUE),
  control.inla = list(
    strategy = "auto",
    int.strategy = "auto",
    optimise.strategy = "smart",
    fast = TRUE,
    compute.initial.values = T,
    stencil = 9,
    dz = 0.1,
    diff.logdens = 0.1,
    numint.maxfeval = 800000000,
    stupid.search = TRUE,
    restart = 10
  ),
  num.threads = "30:1",
  safe = T
)

pass <- F

setname <- c(
  "bym2queenfirstorder",
  "bym2delaunay",
  "bym2queensecondorder",
  "bym2soi",
  "stationary",
  "barrier"
)


if (pass != T) {
  for (name in setname) {
    pathfit <- paste0("hpv/clean_data/fit/fit_", name, ".RDS")
    pathpostfit <- paste0("hpv/clean_data/output/post_fit_", name, ".RDS")
    pathpp <- paste0("hpv/clean_data/output/pp_", name, ".RDS")
    pathgcpo  <- paste0("hpv/clean_data/output/gcpo_", name, ".RDS")
    print(pathfit)
    print(pathpostfit)
    print(pathpp)
    print(pathgcpo)
    
    if (str_detect(name, "delaunay")) {
      W <- delaunay
    } else if (str_detect(name, "soi")) {
      W <- soi
    } else if (str_detect(name, "queenfirstorder")) {
      W <- queenfirstorder
    } else if (str_detect(name, "queensecondorder")) {
      W <- queensecondorder
    }
    # Fit ---------------------------------------------------------------------
    if (!file.exists(pathfit)) {
      x <- function_bru_elements(
        prec_fixed = 1 / 100,
        sd_pc_input = 10,
        prec_corr = 0.2,
        range.fraction.input = 0.1,
        typespace = name
      )
      print(x$cmp)
      #Fit
      fit <- bru(components = x$cmp, x$lik1, options = input_bru_options)
      #Saving
      saveRDS(fit, pathfit)
      print(summary(fit))
    }
    rm("fit")
    gc()
    if (!file.exists(pathpostfit) |
        !file.exists(pathpp) |
        !file.exists(pathgcpo)) {
      fit <- readRDS(pathfit)
      print(summary(fit))
    }
    gc()
    # GCPO ----------------------------------------------------------------------
    if (!file.exists(pathgcpo)) {
      print("GCPO")
      gcpo <- inla.group.cv(
        fit,
        group.cv = NULL,
        num.level.sets = 32,
        size.max = 32,
        strategy = "posterior"
      )
      saveRDS(gcpo, file = pathgcpo)
      rm("gcpo")
    }
    
    fml <- ~ tibble(
      p = plogis(
        intercept +
          intercepthpvother +
          interceptopportunistic +
          intercepthpvotheropportunistic +
          fspace +
          fspaceopportunistic +
          fage +
          fageopportunistic +
          ftime +
          ftimeopportunistic
      )
    )

    # Postfit -----------------------------------------------------------------
    if (!file.exists(pathpostfit)) {
      print("Postfit")
      post_fit <- exporting_results_model(
        fit_input = fit,
        Nsample = 1500,
        Name = 1500,
        namemodel = name
      )
      saveRDS(post_fit, file = pathpostfit)
      rm("post_fit")
    }
    gc()
    
    # PP ----------------------------------------------------------------------
    if (!file.exists(pathpp)) {
      print("PP")
      pp <- pp_check_function(fit_input = fit, Nsample = 2000)
      saveRDS(pp, file = pathpp)
      rm("pp")
    }
    gc()
    rm("fit")
    rm("fml")
  }
  gc()
}

if (!file.exists("hpv/clean_data/output/gcpo.RDS")) {
  barrier <- readRDS("hpv/clean_data/output/gcpo_barrier.RDS")
  stationary <- readRDS("hpv/clean_data/output/gcpo_stationary.RDS")
  delaunay <- readRDS("hpv/clean_data/output/gcpo_bym2delaunay.RDS")
  soi <- readRDS("hpv/clean_data/output/gcpo_bym2soi.RDS")
  queen1 <- readRDS("hpv/clean_data/output/gcpo_bym2queenfirstorder.RDS")
  queen2 <- readRDS("hpv/clean_data/output/gcpo_bym2queensecondorder.RDS")
  
  res <- tibble(
    model = "gcpo_barrier.RDS",
    ls = c(-mean(log(barrier$cv))),
    space = "grf",
    type = "barrier",
    n = "Model 1"
  ) %>%
    add_row(
      tibble(
        model = "gcpo_stationary.RDS",
        ls = c(-mean(log(stationary$cv))),
        space = "grf",
        type = "stationary",
        n = "Model 2"
      )
    ) %>%
    add_row(
      tibble(
        model = "gcpo_bym2delaunay.RDS",
        ls = c(-mean(log(delaunay$cv))),
        space = "bym2",
        type = "delaunay",
        n = "Model 3"
      )
    ) %>%
    add_row(tibble(
      model = "gcpo_bym2soi.RDS",
      ls = c(-mean(log(soi$cv))),
      space = "bym2",
      type = "soi",
      n = "Model 4"
    )) %>%
    add_row(
      tibble(
        model = "gcpo_bym2queenfirstorder.RDS",
        ls = c(-mean(log(queen1$cv))),
        space = "bym2",
        type = "queenfirstorder",
        n = "Model 5"
      )
    ) %>%
    add_row(
      tibble(
        model = "gcpo_bym2queensecondtorder.RDS",
        ls = c(-mean(log(queen2$cv))),
        space = "bym2",
        type = "queensecondorder",
        n = "Model 6"
      )
    )
  saveRDS(res, "hpv/clean_data/output/gcpo.RDS")
} else{
  res <- readRDS("hpv/clean_data/output/gcpo.RDS")
}

# Model selection ---------------------------------------------------------
#NB: WAIC: Smaller values are better, lower out of sample deviance.
#NB: log-cpo is a positively oriented score, i.e. a large value is “good”.
#logarithmic score = - 1/N sum_{i to N}log(gcpo_i)
#We want to maximise 1/N sum_{i to N}log(gcpo_i) or minimise - 1/N sum_{i to N}log(gcpo_i)

res <- res %>%
  filter(ls == min(ls))
path_fit_init <- paste0("hpv/clean_data/fit/", str_replace(res$model, "gcpo", "fit"))
name <- str_remove(str_replace(res$model, "gcpo", "fit"), ".RDS")

print(res)

if (str_detect(res$type, "delaunay")) {
  W <- delaunay
} else if (str_detect(res$type, "soi")) {
  W <- soi
} else if (str_detect(name, "queenfirstorder")) {
  W <- queenfirstorder
} else if (str_detect(name, "queensecondorder")) {
  W <- queensecondorder
}


fml <- ~ tibble(
  p = plogis(
    intercept +
      intercepthpvother +
      interceptopportunistic +
      intercepthpvotheropportunistic +
      fspace +
      fspaceopportunistic +
      fage +
      fageopportunistic +
      ftime +
      ftimeopportunistic
  )
)

gc()

# Sensitivitiy analyses: changing priors for correlation parameters -------
for (sensi in c("sensi1", "sensi2")) {
  print(name)
  pathfit <- paste0("hpv/clean_data/fit/fit_", sensi, ".RDS")
  pathpostfit <- paste0("hpv/clean_data/output/post_fit_", sensi, ".RDS")
  pathpp <- paste0("hpv/clean_data/output/pp_", sensi, ".RDS")
  print(pathfit)
  print(pathpostfit)
  print(pathpp)
  
  
  # Fit ---------------------------------------------------------------------
  if (!file.exists(pathfit)) {
    #If sensi1 or sensi2: increasing the precision of the gaussian prior for the latent correlation parameters
    prec_corr_input <- case_when(sensi == "sensi1" ~ 0.4, sensi == "sensi2" ~ 0.8, TRUE ~ 0.2)
    print(prec_corr_input)
    #Importing baseline fit, used as starting values
    fit <- readRDS(path_fit_init)
    input_bru_options[["control.mode"]] <- list(theta = fit$mode$theta, x = fit$mode$x)
    rm("fit")
    gc()
    
    #Formula call
    x <- function_bru_elements(
      prec_fixed = 1 / 100,
      sd_pc_input = 10,
      prec_corr = prec_corr_input,
      range.fraction.input = 0.1,
      typespace =  paste0(res$space, "_", res$type)
    )
    #Fit
    fit <- bru(components = x$cmp, x$lik1, options = input_bru_options)
    #Saving
    saveRDS(fit, pathfit)
    print(summary(fit))
    rm("fit")
  }
  gc()
  if (!file.exists(pathpostfit) |
      !file.exists(pathpp)) {
    fit <- readRDS(pathfit)
  }
  gc()
  
  
  # Postfit -----------------------------------------------------------------
  if (!file.exists(pathpostfit)) {
    print("Postfit")
    post_fit <- exporting_results_model(
      fit_input = fit,
      Nsample = 1500,
      Name = 1500,
      namemodel = res$space
    )
    saveRDS(post_fit, file = pathpostfit)
    rm("post_fit")
  }
  gc()
}
rm("fit")
gc()