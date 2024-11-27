# Function ----------------------------------------------------------------
# PP ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pp_check_function <- function(fit_input, Nsample) {
  dfbinomial <- dfbinomial %>%
    mutate(age_class = base::cut(
      age,
      breaks = c(30, 39, 49, 59, max(dfbinomial$age)),
      include.lowest = TRUE,
      ordered_result = TRUE
    ))
  # Posterior predictive ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  draw_pp <- generate(
    object = fit_input,
    newdata = dfbinomial,
    n.samples = Nsample,
    seed = 123,
    formula = fml
  )
  
  
  
  library(tidytable)
  
  df <- draw_pp  %>%
    imap( ~ .x  %>%
            data.table() %>%
            mutate(id = row_number(), draws = .y)) %>%
    do.call(rbind, .)  %>%
    left_join(data.table(dfbinomial), by = "id") %>%
    rowwise() %>%
    mutate(ysim = rbinom(n = 1, size = N, prob = p)) %>%
    ungroup() %>%
    mutate(yobs = result) %>%
    ungroup()
  
  print(colnames(df))
  
  summary_function <- function(vector) {
    df %>%
      group_by(c(vector, "draws")) %>%
      summarise(value = sum(ysim) / sum(N)) %>%
      group_by(vector) %>%
      summary_draw() %>%
      tibble() %>%
      left_join(
        dfbinomial  %>%
          group_by(vector) %>%
          summarise(obs = sum(result) / sum(N), N = sum(N)),
        by = c(vector)
      )
  }
  
  out <- list()
  x <- c("virusgroup", "year", "age_class", "region_name", "type")
  list_x <- do.call(c, lapply(seq_along(x), combn, x = x, simplify = FALSE))
  print(length(list_x))
  for (i in 1:length(list_x)) {
    print(i)
    var <- list_x[[i]]
    loop <- summary_function(var)
    print(loop)
    out[[i]] <- loop
    gc()
  }
  
  detach("package:tidytable", unload = TRUE)
  
  return(out)
}

exporting_results_model <- function(fit_input, Nsample, Name, namemodel) {
  
  summary.fixed <- fit_input$summary.fixed
  summary.random <- fit_input$summary.random
  summary.hyperpar <- fit_input$summary.hyperpar
  
  # Hyperparameters ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  correlation <- fit_input$marginals.hyperpar %>%
    imap(function(.x, .y) {
      .x
    }) %>%
    purrr::discard(is_null)
  
  print("AME")
  
  # AME for tests---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  #Keeping only one line by stratum, aggregating the corresponding number of tests
  dfbinomialtmp <- dfbinomial %>%
    group_by(
      geometry,
      idspace,
      age,
      idtime,
      year,
      hpvother,
      virusgroup,
      virus,
      insee_dep,
      region_name
    ) %>%
    summarise(N = sum(N)) %>%
    ungroup() %>%
    mutate(id = row_number())
  
  colnames(dfbinomialtmp)
  
  #Duplicating the dataframe: one block for organised, one block for opportunistic
  datatmp <- rbind(
    dfbinomialtmp %>%
      mutate(opportunistic = 0, test = "Organised"),
    dfbinomialtmp %>%
      mutate(opportunistic = 1, test = "Opportunistic")
  ) %>%
    mutate(testgroup = opportunistic + 1, id_join = row_number()) %>%
    mutate(
      age_class = base::cut(
        age,
        breaks = c(30, 39, 49, 59, max(dfbinomialtmp$age)),
        include.lowest = TRUE,
        ordered_result = TRUE
      ),
      intercept = 1
    )
  
  # Total number of tests ---------------------------------------------------------
  #We want to have the number of tests performed stratified by dimension
  #Becaus eah test yield a result for HPV16/19 and other genotypes, the number of tests is duplicated
  #Filtering in only HPV1618
  nall <- dfbinomial %>%
    st_drop_geometry() %>%
    filter(virus == "hpv1618") %>%
    summarise(N = sum(N)) %>%
    ungroup
  
  nyear <- dfbinomial %>%
    st_drop_geometry() %>%
    filter(virus == "hpv1618") %>%
    group_by(year) %>%
    summarise(N = sum(N))  %>%
    ungroup
  
  nage <- dfbinomial %>%
    st_drop_geometry() %>%
    filter(virus == "hpv1618") %>%
    mutate(age_class = base::cut(
      age,
      breaks = c(30, 39, 49, 59, max(dfbinomialtmp$age)),
      include.lowest = TRUE,
      ordered_result = TRUE
    )) %>%
    group_by(age_class) %>%
    summarise(N = sum(N)) %>%
    ungroup
  
  nageyear <- dfbinomial %>%
    st_drop_geometry() %>%
    mutate(age_class = base::cut(
      age,
      breaks = c(30, 39, 49, 59, max(dfbinomialtmp$age)),
      include.lowest = TRUE,
      ordered_result = TRUE
    )) %>%
    filter(virus == "hpv1618") %>%
    group_by(age_class, year) %>%
    summarise(N = sum(N)) %>%
    ungroup
  
  
  print(Name)
  #Draws from joint
  tmp <- generate(
    object = fit_input,
    newdata = datatmp,
    n.samples = Name,
    seed = 123,
    formula = fml
  )
  #  saveRDS(tmp, "hpv/clean_data/test.RDS")
  print(nrow(tmp[[1]]))
  print(colnames(tmp[[1]]))
  
  
  #Adding variables to draws
  #Computing expecting number
  library(tidytable)
  tmp_bis <- tmp %>%
    imap( ~ .x %>%
            mutate(id_join = row_number(), draws = .y)) %>%
    do.call(rbind, .)  %>%
    left_join(datatmp %>% st_drop_geometry() %>% data.table, by = "id_join") %>%
    mutate(ENpos = p * N)
  
  #Creating df for the ame
  data_for_ame_test <- tmp_bis %>%
    ungroup() %>%
    select(id,
           draws,
           age,
           age_class,
           year,
           virus,
           insee_dep,
           region_name,
           test,
           ENpos) %>%
    pivot_wider(names_from = test, values_from = ENpos) %>%
    mutate(value = Opportunistic - Organised)
  
  print(colnames(data_for_ame_test))
  
  
  #Expected prevalence for each pathway
  print("Summary")
  data_for_expected_prevalence <- data_for_ame_test %>%
    mutate(test = "Opportunistic", value = Opportunistic) %>%
    select(-Organised, -Opportunistic) %>%
    rbind(
      data_for_ame_test %>%
        mutate(test = "Organised", value = Organised) %>%
        select(-Organised, -Opportunistic)
    )
  
  print("Prevalence")
  print("Total")
  expected_prevalence <- data_for_expected_prevalence %>%
    group_by(draws, virus, test) %>%
    summarise(value = sum(value)) %>%
    mutate(value = value / nall$N) %>%
    group_by(virus, test) %>%
    summary_draw() %>%
    tibble()
  
  print(expected_prevalence)
  
  print("Year")
  expected_prevalence_year <- data_for_expected_prevalence %>%
    group_by(draws, virus, test, year) %>%
    summarise(value = sum(value)) %>%
    left_join(nyear, by = c('year')) %>%
    mutate(value = value / N) %>%
    group_by(virus, test, year) %>%
    summary_draw() %>%
    tibble()
  print(expected_prevalence_year)
  
  print("Age")
  expected_prevalence_age <- data_for_expected_prevalence %>%
    group_by(draws, virus, test, age_class) %>%
    summarise(value = sum(value)) %>%
    left_join(nage, by = c('age_class')) %>%
    mutate(value = value / N) %>%
    group_by(virus, age_class) %>%
    summary_draw() %>%
    tibble()
  print(expected_prevalence_age)
  
  
  print("Year-Age")
  expected_prevalence_year_age <- data_for_expected_prevalence %>%
    group_by(draws, virus, test, year, age_class) %>%
    summarise(value = sum(value)) %>%
    left_join(nageyear, by = c("year", "age_class")) %>%
    mutate(value = value / N) %>%
    group_by(virus, test, year, age_class) %>%
    summary_draw() %>%
    tibble()
  
  print(expected_prevalence_year_age)
  
  #AME-Rather marginal expected difference b/w opportunistic and organised
  ame_test <- data_for_ame_test %>%
    group_by(draws, virus) %>%
    summarise(value = sum(value) / nall$N) %>%
    group_by(virus) %>%
    summary_draw() %>%
    tibble()
  print(ame_test)
  
  ame_test_year <- data_for_ame_test %>%
    group_by(draws, virus, year) %>%
    summarise(value = sum(value)) %>%
    left_join(nyear, by = c("year")) %>%
    group_by(virus, year) %>%
    mutate(value = value / N) %>%
    summary_draw() %>%
    tibble()
  print(ame_test_year)
  
  ame_test_age <- data_for_ame_test  %>%
    group_by(draws, virus, age_class) %>%
    summarise(value = sum(value)) %>%
    left_join(nage, by = c("age_class")) %>%
    mutate(value = value / N) %>%
    group_by(virus, age_class) %>%
    summary_draw() %>%
    tibble()
  print(ame_test_age)
  
  ame_test_age_year <- data_for_ame_test  %>%
    group_by(draws, virus, age_class, year) %>%
    summarise(value = sum(value)) %>%
    left_join(nageyear, by = c("year", "age_class")) %>%
    group_by(virus, age_class, year) %>%
    mutate(value = value / N) %>%
    summary_draw() %>%
    tibble()
  print(ame_test_age_year)
  
  ame_test_after2020 <- data_for_ame_test %>%
    filter(year > 2020)  %>%
    group_by(draws, virus) %>%
    summarise(value = sum(value) / nall$N) %>%
    group_by(virus) %>%
    summary_draw() %>%
    tibble()
  print(ame_test_after2020)
  
  ame_test_after2021 <- data_for_ame_test %>%
    filter(year > 2021)  %>%
    group_by(draws, virus) %>%
    summarise(value = sum(value) / nall$N) %>%
    group_by(virus) %>%
    summary_draw() %>%
    tibble()
  print(ame_test_after2021)
  
  detach("package:tidytable", unload = TRUE)
  rm("tmp")
  rm("tmp_bis")
  rm("datatmp")
  rm("data_for_ame_test")
  rm("dfbinomialtmp")
  
  gc()
  
  print("Risk level")
  ##Space -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  generate_risk_level <- function(df_input, age_input, idtime_input) {
    x <- crossing(
      df_input,
      expand.grid(age = age_input),
      expand.grid(idtime = idtime_input),
      expand.grid(tibble(
        hpv1618 = 1, hpvother = 0
      )),
      expand.grid(opportunistic = c(0, 1))
    ) %>%
      st_as_sf() %>%
      st_transform(crs(kmproj)) %>%
      select(-one_of(c("id")))
    
    x <- rbind(x, x %>% mutate(hpv1618 = 0, hpvother = 1))
    
    #Adding lastly required covariates
    datatmp <- x %>%
      mutate(
        intercept = 1,
        testgroup = opportunistic + 1,
        virusgroup = ifelse(hpv1618 == 1, 1, 2),
        id = row_number()
      )
    
    print(colnames(datatmp))
    print(nrow(datatmp))
    
    state <- evaluate_state(
      model = fit_input$bru_info$model,
      result = fit_input,
      data = datatmp,
      property = 'sample',
      n = Nsample,
      seed = 123
    )
    
    
    tmp <- inlabru::evaluate_model(
      model = fit_input$bru_info$model,
      state = state,
      data = datatmp,
      predictor =  fml
    )
    
    
    library(tidytable)
    tmp_new <- tmp %>%
      imap(~ data.table(.x) %>%
             mutate(id = row_number(), draws = .y) %>% 
             rename(value = p)) %>%
      do.call(rbind, .) %>%
      left_join(datatmp %>% st_drop_geometry() %>% data.table(), by = "id")
    
    virusspecific <- tmp_new %>%
      group_by(id) %>%
      summary_draw() %>%
      tibble()   %>%
      left_join(datatmp, by = "id") %>%
      st_as_sf()
    
    print(virusspecific)
    
    difference_virus <- tmp_new %>%
      select(-id, -hpv1618, -hpvother, -opportunistic) %>%
      pivot_wider(names_from = virusgroup, values_from = value) %>%
      rename(hpv1618 = `1`, hpvother = `2`) %>%
      mutate(value = hpvother -  hpv1618) %>%
      group_by(idgeometry, testgroup, age) %>%
      summary_draw() %>%
      tibble()   %>%
      left_join(df_input %>%
                  select(one_of(colnames(pixel_pred))), by = "idgeometry") %>%
      st_as_sf()
    
    print(difference_virus)
    
    difference_test <- tmp_new %>%
      select(-id, -hpv1618, -hpvother, -opportunistic) %>%
      pivot_wider(names_from = testgroup, values_from = value) %>%
      rename(organised = `1`, opportunistic = `2`) %>%
      mutate(value = opportunistic -  organised) %>%
      group_by(idgeometry, virusgroup, age) %>%
      summary_draw() %>%
      tibble()   %>%
      left_join(df_input %>%
                  select(one_of(colnames(pixel_pred))), by = "idgeometry") %>%
      st_as_sf()
    
    print(difference_test)
    detach("package:tidytable", unload = TRUE)
    
    out <- list(
      "virusspecific" = virusspecific,
      "difference_virus" = difference_virus,
      "difference_test" = difference_test
    )
    return(out)
  }
  print(namemodel)
  if (str_detect(namemodel, "bym2")) {
    spatial_domain <- idcentroidspred
  }
  else{
    spatial_domain <- pixel_pred
  }
  
  print(colnames(spatial_domain))
  
  
  baseline_space_30 <-   generate_risk_level(
    df_input = spatial_domain,
    age_input = 30,
    idtime_input = c(max(dfdate$idtime))
  )
  gc()
  baseline_space_paris_30 <- baseline_space_30 %>%
    imap( ~ .x %>% filter(insee_reg == 11))
  gc()
  
  baseline_space_66 <-   generate_risk_level(
    df_input = spatial_domain,
    age_input = 66,
    idtime_input = c(max(dfdate$idtime))
  )
  gc()
  baseline_space_paris_66 <- baseline_space_66 %>%
    imap( ~ .x %>% filter(insee_reg == 11))
  gc()
  
  
  print("French city")
  # Age/Time/sampling in major french city ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  datatmp <- crossing(
    centroid_major_city,
    expand.grid(
      age = seq(min(dfbinomial$age), max(dfbinomial$age)),
      opportunistic = c(0, 1),
      idtime = c(min(data_bru$idtime), max(dfdate$idtime)),
      virusgroup = c(1, 2)
    )
  ) %>%
    st_as_sf() %>%
    st_transform(crs(kmproj)) %>%
    mutate(
      intercept = 1,
      testgroup = opportunistic + 1,
      hpv1618 = ifelse(virusgroup == 1, 1, 0),
      hpvother = 1 - hpv1618,
      id = row_number()
    ) 
  
  #Generate predictions
  effect_main_cities <- generate(
    fit_input,
    datatmp,
    n.samples = Nsample,
    probs = vector_quantiles,
    seed = 123,
    formula = fml
  )
  
  
  library(tidytable)
  
  effect_main_cities <- effect_main_cities %>%
    imap(~ data.table(.x) %>%
           mutate(id = row_number(), draws = .y)) %>%
    do.call(rbind, .) %>%
    left_join(data.table(datatmp %>% st_drop_geometry()), by = "id")
  
  
  summary_city <- effect_main_cities %>%
    rename(value = p) %>%
    group_by(city, idspace, age, idtime, virusgroup, testgroup) %>%
    summary_draw() %>%
    tibble()
  
  
  difference_virus_summary_city <- effect_main_cities %>%
    select(-id, -hpv1618, -hpvother) %>%
    rename(value = p) %>%
    mutate(virusgroup = ifelse(virusgroup == 1, "hpv1618", "hpvother")) %>%
    data.table() %>%
    pivot_wider(names_from = virusgroup, values_from = value) %>%
    mutate(value = hpvother - hpv1618) %>%
    group_by(city, idspace, age, idtime, testgroup) %>%
    summary_draw() %>%
    tibble()
  
  
  difference_virus_test_summary_city <- effect_main_cities %>%
    select(-id, -hpv1618, -hpvother, -intercept, -opportunistic) %>%
    rename(value = p) %>%
    mutate(virusgroup = ifelse(virusgroup == 1, "hpv1618", "hpvother")) %>%
    data.table() %>%
    pivot_wider(names_from = virusgroup, values_from = value) %>%
    mutate(value = hpvother - hpv1618) %>%
    dplyr::select(-hpvother, -hpv1618) %>%
    group_by(city, idspace, age, idtime, draws) %>%
    arrange(testgroup) %>%
    mutate(value = value - value[1]) %>%
    filter(testgroup == 2) %>%
    ungroup() %>%
    group_by(city, idspace, age, idtime) %>%
    summary_draw() %>%
    mutate(label = "Difference between test (opportunistic-organised) of the difference between virus (other - 1618)") %>%
    tibble
  
  
  difference_test_summary_city <- effect_main_cities %>%
    select(-id, -hpv1618, -hpvother, -opportunistic) %>%
    rename(value = p) %>%
    mutate(testgroup = ifelse(testgroup == 1, "organised", "opportunistic")) %>%
    data.table() %>%
    pivot_wider(names_from = testgroup, values_from = value) %>%
    mutate(value = opportunistic - organised) %>%
    group_by(city, idspace, age, idtime, virusgroup) %>%
    summary_draw() %>%
    tibble()
  
  tmp <- effect_main_cities %>%
    rename(value = p) %>%
    arrange(age) %>%
    group_by(draws, city, idspace, idtime, virusgroup, testgroup) %>%
    mutate(value_ref = value[1]) %>%
    ungroup()
  
  #Age
  summary_city_age <- rbind(
    tmp %>%
      data.table() %>%
      mutate(value = value - value_ref) %>%
      group_by(city, idspace, age, idtime, virusgroup, testgroup) %>%
      summary_draw() %>%
      mutate(metric = "rd"),
    tmp %>%
      data.table() %>%
      mutate(value = value / value_ref) %>%
      group_by(city, idspace, age, idtime, virusgroup, testgroup) %>%
      summary_draw() %>%
      mutate(metric = "rr"),
    tmp %>%
      data.table() %>%
      mutate(value = (value / (1 - value)) / (value_ref / (1 - value_ref))) %>%
      group_by(city, idspace, age, idtime, virusgroup, testgroup) %>%
      summary_draw() %>%
      mutate(metric = "or")
  ) %>%
    tibble()
  
  
  #Year
  tmp <- effect_main_cities %>%
    rename(value = p) %>%
    filter(idtime %in% c(min(idtime), max(idtime))) %>%
    arrange(idtime) %>%
    group_by(draws, city, idspace, age, virusgroup, testgroup) %>%
    mutate(value_ref = value[1]) %>%
    ungroup()
  
  summary_city_year <- rbind(
    tmp %>%
      mutate(value = value - value_ref) %>%
      group_by(city, idspace, age, idtime, virusgroup, testgroup) %>%
      summary_draw() %>%
      mutate(metric = "rd"),
    tmp %>%
      mutate(value = value / value_ref) %>%
      group_by(city, idspace, age, idtime, virusgroup, testgroup) %>%
      summary_draw() %>%
      mutate(metric = "rr"),
    tmp %>%
      mutate(value = (value / (1 - value)) / (value_ref / (1 - value_ref))) %>%
      group_by(city, idspace, age, idtime, virusgroup, testgroup) %>%
      summary_draw() %>%
      mutate(metric = "or")
  ) %>%
    tibble()
  rm("tmp")
  rm("datatmp")
  gc()
  detach("package:tidytable", unload = TRUE)
  
  rm("effect_main_cities")
  
  
  # Out ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  out <- list(
    "expected_prevalence" = expected_prevalence,
    "expected_prevalence_year" = expected_prevalence_year,
    "expected_prevalence_age" = expected_prevalence_age,
    "expected_prevalence_year_age" = expected_prevalence_year_age,
    "ame_test" = ame_test,
    "ame_test_year" = ame_test_year,
    "ame_test_age" = ame_test_age,
    "ame_test_age_year" = ame_test_age_year,
    "ame_test_after2020" = ame_test_after2020,
    "ame_test_after2021" = ame_test_after2021,
    "summary_city" = summary_city,
    "difference_virus_summary_city" = difference_virus_summary_city,
    "difference_virus_test_summary_city" = difference_virus_test_summary_city,
    "difference_test_summary_city" = difference_test_summary_city,
    "summary_city_age" = summary_city_age,
    "summary_city_year" = summary_city_year,
    "baseline_space_30" = baseline_space_30,
    "baseline_space_paris_30" = baseline_space_paris_30,
    "baseline_space_66" = baseline_space_66,
    "baseline_space_paris_66" = baseline_space_paris_66,
    "correlation" = correlation,
    "summary.fixed" = summary.fixed,
    "summary.random" = summary.random,
    "summary.hyperpar" = summary.hyperpar
  )
  
  
  
  return(out)
  
}



#ggplot()+gg(pixel_pred, aes(color = type))
# Useful function ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
summary_draw <- function(df) {
  df %>%
    summarise(
      mean = mean(value, na.rm = T),
      median = median(value, na.rm = T),
      qi_lb = tidybayes::qi(value, .width = 0.95, na.rm = T)[1],
      qi_ub = tidybayes::qi(value, .width = 0.95, na.rm = T)[2],
      qi_lb_09 = tidybayes::qi(value, .width = 0.9, na.rm = T)[1],
      qi_ub_09 = tidybayes::qi(value, .width = 0.9, na.rm = T)[2],
      qi_lb_08 = tidybayes::qi(value, .width = 0.8, na.rm = T)[1],
      qi_ub_08 = tidybayes::qi(value, .width = 0.8, na.rm = T)[2],
      qi_lb_07 = tidybayes::qi(value, .width = 0.7, na.rm = T)[1],
      qi_ub_07 = tidybayes::qi(value, .width = 0.7, na.rm = T)[2],
      qi_lb_06 = tidybayes::qi(value, .width = 0.6, na.rm = T)[1],
      qi_ub_06 = tidybayes::qi(value, .width = 0.6, na.rm = T)[2],
      qi_lb_05 = tidybayes::qi(value, .width = 0.5, na.rm = T)[1],
      qi_ub_05 = tidybayes::qi(value, .width = 0.5, na.rm = T)[2],
      min = min(value, na.rm = T),
      max = max(value, na.rm = T),
      N_higher_0 = sum(value > 0, na.rm = T),
      N_higher_1 = sum(value > 1, na.rm = T),
      higher_0 = mean(value > 0, na.rm = T),
      higher_05 = mean(value > 0.5, na.rm = T),
      higher_1 = mean(value > 1, na.rm = T),
      lower_0 = mean(value < 0, na.rm = T),
      lower_05 = mean(value < 0.5, na.rm = T),
      lower_1 = mean(value < 1, na.rm = T),
      higher_2 = mean(value > 2, na.rm = T)
    ) %>%
    ungroup()
}
