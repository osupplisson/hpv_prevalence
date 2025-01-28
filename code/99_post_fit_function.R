fun_format <- function(input) {
  if ("testgroup" %in% colnames(input)) {
    input <- input %>%
      mutate(test = ifelse(testgroup == 1, "Organised", "Opportunistic"))
  }
  
  if ("virusgroup" %in% colnames(input)) {
    input <- input %>%
      mutate(genotypes = ifelse(virusgroup == 1, "HPV16/18", "Other genotypes"))
  }
  
  if ("virus" %in% colnames(input)) {
    input <- input %>%
      mutate(genotypes = ifelse(virus == "hpv1618", "HPV16/18", "Other genotypes"))
  }
  
  if ("idyear" %in% colnames(input)) {
    input <- input %>%
      mutate(year = case_when(
        idyear == 1 ~ 2020,
        idyear == 2 ~ 2021,
        idyear == 3 ~ 2022,
        idyear == 4 ~ 2023,
        TRUE ~ NA
      ))
  }
  
  
  
  input
}

# Function ----------------------------------------------------------------
summary_function <- function(input, vector, type) {
  if (type == "pp") {
    tmp <- input %>%
      mutate(val = ysim)
  } else{
    tmp <- input %>%
      mutate(val = yexp)
  }
  tmp %>%
    group_by(c(vector, "draws")) %>%
    summarise(value = sum(val) / sum(N)) %>%
    group_by(vector) %>%
    summary_draw() %>%
    tibble() %>%
    left_join(dfbinomial  %>%
                group_by(vector) %>%
                summarise(obs = sum(result) / sum(N), N = sum(N)),
              by = c(vector)) %>%
    fun_format
}


# PP ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pp_check_function <- function(fit_input, stateinput, type = "pp") {
  dfbinomial <- dfbinomial %>%
    mutate(age_class = base::cut(
      age,
      breaks = c(30, 39, 49, 59, max(dfbinomial$age)),
      include.lowest = TRUE,
      ordered_result = TRUE
    ))
  
  # Posterior predictive ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  print("Evaluate model")
  draw_pp <- inlabru::evaluate_model(
    model = fit_input$bru_info$model,
    state = stateinput,
    data = dfbinomial,
    predictor = fmlpred
  )
  
  
  library(tidytable)
  setDTthreads(10)
  df <- draw_pp  %>%
    imap( ~ .x  %>%
            data.table() %>%
            mutate(id = row_number(), draws = .y)) %>%
    do.call(rbind, .)  %>%
    left_join(data.table(dfbinomial), by = "id") %>%
    rowwise() %>%
    mutate(ysim = rbinom(n = 1, size = N, prob = p), yexp = p * N) %>%
    ungroup() %>%
    mutate(yobs = result) %>%
    ungroup()
  
  print(colnames(df))
  
  
  out <- list()
  x <- c("virusgroup", "year", "age_class", "region_name", "type")
  list_x <- do.call(c, lapply(seq_along(x), combn, x = x, simplify = FALSE))
  print(length(list_x))
  for (i in 1:length(list_x)) {
    print(i)
    var <- list_x[[i]]
    loop <- summary_function(df, var, type)
    print(loop)
    out[[i]] <- loop
    gc()
  }
  
  detach("package:tidytable", unload = TRUE)
  
  return(out)
}

exporting_results_model <- function(fit_input, stateinput, namemodel, runmap = T) {
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
  
  
  #Evaluate from joint
  print("Evaluate model")
  tmp <- inlabru::evaluate_model(
    model = fit_input$bru_info$model,
    state = stateinput[seq(1, min(length(stateinput), 1550))],
    data = datatmp,
    predictor = fmlpred
  )
  
  
  #Adding variables to draws
  #Computing expecting number
  library(tidytable)
  setDTthreads(10)
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
  prevalence <- list()
  prevalence[["expected_prevalence"]] <- data_for_expected_prevalence %>%
    group_by(draws, virus, test) %>%
    summarise(value = sum(value)) %>%
    mutate(value = value / nall$N) %>%
    group_by(virus, test) %>%
    summary_draw() %>%
    tibble() %>%
    fun_format
  
  print(prevalence[["expected_prevalence"]])
  
  print("Year")
  prevalence[["expected_prevalence_year"]] <- data_for_expected_prevalence %>%
    group_by(draws, virus, test, year) %>%
    summarise(value = sum(value)) %>%
    left_join(nyear, by = c('year')) %>%
    mutate(value = value / N) %>%
    group_by(virus, test, year) %>%
    summary_draw() %>%
    tibble() %>%
    fun_format
  
  print(prevalence[["expected_prevalence_year"]])
  
  print("Age")
  prevalence[["expected_prevalence_age"]] <- data_for_expected_prevalence %>%
    group_by(draws, virus, test, age_class) %>%
    summarise(value = sum(value)) %>%
    left_join(nage, by = c('age_class')) %>%
    mutate(value = value / N) %>%
    group_by(virus, age_class) %>%
    summary_draw() %>%
    tibble() %>%
    fun_format
  
  print(prevalence[["expected_prevalence_age"]])
  
  
  print("Year-Age")
  prevalence[["expected_prevalence_year_age"]] <- data_for_expected_prevalence %>%
    group_by(draws, virus, test, year, age_class) %>%
    summarise(value = sum(value)) %>%
    left_join(nageyear, by = c("year", "age_class")) %>%
    mutate(value = value / N) %>%
    group_by(virus, test, year, age_class) %>%
    summary_draw() %>%
    tibble() %>%
    fun_format
  
  print(prevalence[["expected_prevalence_year_age"]])
  
  #AME-Rather marginal expected difference b/w opportunistic and organised
  ame <- list()
  ame[["ame_test"]] <- data_for_ame_test %>%
    group_by(draws, virus) %>%
    summarise(value = sum(value) / nall$N) %>%
    group_by(virus) %>%
    summary_draw() %>%
    tibble() %>%
    fun_format
  
  print(ame[["ame_test"]])
  
  ame[["ame_test_year"]] <- data_for_ame_test %>%
    group_by(draws, virus, year) %>%
    summarise(value = sum(value)) %>%
    left_join(nyear, by = c("year")) %>%
    group_by(virus, year) %>%
    mutate(value = value / N) %>%
    summary_draw() %>%
    tibble() %>%
    fun_format
  
  print(ame[["ame_test_year"]])
  
  ame[["ame_test_age"]] <- data_for_ame_test  %>%
    group_by(draws, virus, age_class) %>%
    summarise(value = sum(value)) %>%
    left_join(nage, by = c("age_class")) %>%
    mutate(value = value / N) %>%
    group_by(virus, age_class) %>%
    summary_draw() %>%
    tibble() %>%
    fun_format
  
  print(ame[["ame_test_age"]])
  
  ame[["ame_test_age_year"]] <- data_for_ame_test  %>%
    group_by(draws, virus, age_class, year) %>%
    summarise(value = sum(value)) %>%
    left_join(nageyear, by = c("year", "age_class")) %>%
    group_by(virus, age_class, year) %>%
    mutate(value = value / N) %>%
    summary_draw() %>%
    tibble() %>%
    fun_format
  
  print(ame[["ame_test_age_year"]])
  
  
  detach("package:tidytable", unload = TRUE)
  rm("tmp")
  rm("tmp_bis")
  rm("datatmp")
  rm("data_for_ame_test")
  rm("dfbinomialtmp")
  
  gc()
  
  print("Risk level")
  ##Space -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  generate_risk_level <- function(df_input,
                                  age_input,
                                  idtime_input,
                                  city = F) {
    datatmp <- crossing(
      df_input,
      expand.grid(age = age_input),
      expand.grid(idtime = idtime_input),
      expand.grid(virusgroup = c(1, 2)),
      expand.grid(opportunistic = c(0, 1))
    ) %>%
      st_as_sf() %>%
      st_transform(crs(kmproj)) %>%
      select(-one_of(c("id"))) %>%
      mutate(intercept = 1,
             testgroup = opportunistic + 1,
             id = row_number())
    table(datatmp$opportunistic, datatmp$testgroup)
    table(datatmp$opportunistic, datatmp$testgroup)
    #Should be 4
    datatmp %>% st_drop_geometry %>% group_by(idgeometry) %>% count() %>% .$n %>% unique
    print(colnames(datatmp))
    
    print(nrow(datatmp))
    print("Evaluate model")
    tmp <- inlabru::evaluate_model(
      model = fit_input$bru_info$model,
      state = stateinput,
      data = datatmp,
      predictor =  fmlpred
    )
    
    
    library(tidytable)
    setDTthreads(10)
    tmp_new <- tmp %>%
      imap(~ data.table(.x) %>%
             mutate(id = row_number(), draws = .y) %>%
             rename(value = p)) %>%
      do.call(rbind, .) %>%
      left_join(datatmp %>% st_drop_geometry() %>% data.table(), by = "id")
    
    print(colnames(tmp_new))
    
    ##Expected prevalence by post code
    virusspecific <- tmp_new %>%
      st_drop_geometry() %>%
      as_tidytable() %>%
      group_by(id) %>%
      summary_draw() %>%
      tibble()   %>%
      left_join(datatmp, by = "id") %>%
      st_as_sf()
    
    print(virusspecific)
    
    #Difference between virus within screening type
    difference_virus <- tmp_new %>%
      as_tidytable() %>%
      st_drop_geometry() %>%
      group_by(idgeometry, age, idtime, testgroup, draws) %>%
      arrange(virusgroup) %>%
      mutate(value = value - value[1]) %>%
      group_by(idgeometry, testgroup, virusgroup, age, idtime) %>%
      summary_draw() %>%
      filter(virusgroup == 2) %>%
      tibble()   %>%
      left_join(df_input %>%
                  select(one_of(colnames(pixel_pred))), by = "idgeometry") %>%
      st_as_sf()
    
    print(difference_virus)
    
    #Difference between test within virus
    difference_test <- tmp_new %>%
      st_drop_geometry() %>%
      as_tidytable() %>%
      group_by(idgeometry, age, idtime, virusgroup, draws) %>%
      arrange(testgroup) %>%
      mutate(value = value - value[1]) %>%
      group_by(idgeometry, testgroup, virusgroup, age, idtime) %>%
      summary_draw() %>%
      tibble()   %>%
      filter(testgroup == 2) %>%
      left_join(df_input %>%
                  select(one_of(colnames(pixel_pred))), by = "idgeometry") %>%
      st_as_sf()
    print(difference_test)
    
    #Difference (between virus within screening type) between screening type
    difference_virus_test <- tmp_new %>%
      st_drop_geometry() %>%
      as_tidytable() %>%
      group_by(idgeometry, age, idtime, testgroup, draws) %>%
      arrange(virusgroup) %>%
      mutate(value = value - value[1]) %>%
      group_by(idgeometry, age, idtime, virusgroup, draws) %>%
      arrange(testgroup) %>%
      mutate(value = value - value[1]) %>%
      group_by(idgeometry, age, idtime, testgroup, virusgroup) %>%
      summary_draw() %>%
      filter(testgroup == 2 & virusgroup == 2) %>%
      mutate(metric = "Difference in difference (between screening between genotypes within screening)") %>%
      tibble()   %>%
      left_join(df_input %>%
                  select(one_of(colnames(pixel_pred))), by = "idgeometry") %>%
      st_as_sf()
    
    print(difference_virus_test)
    
    
    
    ###Count
    count_difference_test <- tmp_new %>%
      st_drop_geometry() %>%
      as_tidytable() %>%
      group_by(idgeometry, age, idtime, virusgroup, draws) %>%
      arrange(testgroup) %>%
      mutate(value = value - value[1]) %>%
      mutate(value = ifelse(value > 0, 1, 0)) %>%
      group_by(virusgroup, testgroup, age, idtime, draws) %>%
      summarise(value = sum(value)) %>%
      group_by(virusgroup, testgroup, age, idtime) %>%
      summary_draw() %>%
      filter(testgroup == 2) %>%
      tibble()
    
    
    print(count_difference_test)
    
    #Counting the number of district for which we found a greater expected prevalence for other genotypes than HPV16/18
    #For each age and time and screening pathway
    count_difference_virus <- tmp_new %>%
      st_drop_geometry() %>%
      as_tidytable() %>%
      group_by(idgeometry, testgroup, age, idtime, draws) %>%
      arrange(virusgroup) %>%
      mutate(value = value -  value[1]) %>%
      mutate(value = ifelse(value > 0, 1, 0)) %>%
      group_by(virusgroup, testgroup, age, idtime, draws) %>%
      summarise(value = sum(value)) %>%
      group_by(virusgroup, testgroup, age, idtime) %>%
      summary_draw() %>%
      tibble() %>%
      filter(virusgroup == 2)
    
    print(count_difference_virus)
    
    #Counting the number of district for which we found a greater expected prevalence for other genotypes than HPV16/18
    #For each age and time
    count_difference_virus_test <- tmp_new %>%
      st_drop_geometry() %>%
      as_tidytable() %>%
      group_by(idgeometry, testgroup, age, idtime, draws) %>%
      arrange(virusgroup) %>%
      mutate(value = value -  value[1]) %>%
      mutate(value = ifelse(value > 0, 1, 0)) %>%
      group_by(testgroup, virusgroup, age, idtime, draws) %>%
      summarise(value = sum(value)) %>%
      filter(virusgroup == 2) %>%
      ungroup() %>%
      group_by(age, idtime, draws) %>%
      arrange(testgroup) %>%
      mutate(value = value -  value[1]) %>%
      filter(testgroup == 2) %>%
      group_by(age, idtime) %>%
      summary_draw() %>%
      tibble() %>%
      mutate(var = "Difference in difference (between screening between genotypes within screening)")
    
    
    print(count_difference_virus_test)
    
    detach("package:tidytable", unload = TRUE)
    
    
    out <- list(
      "virusspecific" = virusspecific %>% fun_format,
      "difference_virus" = difference_virus %>% fun_format,
      "difference_test" = difference_test %>% fun_format,
      "difference_virus_test" = difference_virus_test %>% fun_format,
      "count_difference_test" = count_difference_test %>% fun_format,
      "count_difference_virus" = count_difference_virus %>% fun_format,
      "count_difference_virus_test" = count_difference_virus_test %>% fun_format
    )
    
    
    if (city == T) {
      print("Part special city")
      #RD/RR/OR for age in main cities
      summary_city_age_tmp <- tmp_new %>%
        st_drop_geometry() %>%
        group_by(draws, city, idspace, idtime, virusgroup, testgroup) %>%
        arrange(age) %>%
        mutate(
          rd = value - value[1],
          rr = value / value[1],
          or = (value / (1 - value)) / (value[1] / (1 - value[1]))
        ) %>%
        ungroup()
      
      summary_city_age <- plyr::rbind.fill(
        summary_city_age_tmp %>%
          mutate(value = rd) %>%
          group_by(city, idspace, idtime, virusgroup, testgroup, age) %>%
          summary_draw() %>%
          mutate(metric = "rd"),
        summary_city_age_tmp %>%
          mutate(value = rr) %>%
          group_by(city, idspace, idtime, virusgroup, testgroup, age) %>%
          summary_draw() %>%
          mutate(metric = "rr"),
        summary_city_age_tmp %>%
          mutate(value = or) %>%
          group_by(city, idspace, idtime, virusgroup, testgroup, age) %>%
          summary_draw() %>%
          mutate(metric = "or")
      ) %>%
        tibble()
      
      
      #RD/RR/OR for time in main cities
      summary_city_time_tmp <- tmp_new %>%
        st_drop_geometry() %>%
        group_by(draws, city, idspace, age, virusgroup, testgroup) %>%
        arrange(idtime) %>%
        mutate(
          rd = value - value[1],
          rr = value / value[1],
          or = (value / (1 - value)) / (value[1] / (1 - value[1]))
        ) %>%
        ungroup()
      
      
      summary_city_time <- plyr::rbind.fill(
        summary_city_time_tmp %>%
          mutate(value = rd) %>%
          group_by(city, idspace, idtime, virusgroup, testgroup, age) %>%
          summary_draw() %>%
          mutate(metric = "rd"),
        summary_city_time_tmp %>%
          mutate(value = rr) %>%
          group_by(city, idspace, idtime, virusgroup, testgroup, age) %>%
          summary_draw() %>%
          mutate(metric = "rr"),
        summary_city_time_tmp %>%
          mutate(value = or) %>%
          group_by(city, idspace, idtime, virusgroup, testgroup, age) %>%
          summary_draw() %>%
          mutate(metric = "or")
      ) %>%
        tibble()
      
      
      out[["summary_city_age"]] <- summary_city_age %>% fun_format
      out[["summary_city_year"]] <- summary_city_time %>% fun_format
    }
    rm("tmp_new")
    gc()
    return(out)
  }
  
  
  
  print("Restriction to French cities")
  print(namemodel)
  
  if (str_detect(namemodel, "bym2")) {
    spatial_domain <- idcentroidspred
  }
  else{
    spatial_domain <- idcentroidspred %>%
      st_centroid()
  }
  
  
  generate_diff_expected_prevalence_number_cp <- function(df_input, age_input, idtime_input) {
    datatmp <- crossing(
      df_input,
      expand.grid(age = age_input),
      expand.grid(idtime = idtime_input),
      expand.grid(virusgroup = c(1, 2)),
      expand.grid(opportunistic = c(0, 1))
    ) %>%
      st_as_sf() %>%
      st_transform(crs(kmproj)) %>%
      select(-one_of(c("id"))) %>%
      mutate(intercept = 1,
             testgroup = opportunistic + 1,
             id = row_number())
    
    
    print(colnames(datatmp))
    print(nrow(datatmp))
    
    print("Evaluate model")
    tmp <- inlabru::evaluate_model(
      model = fit_input$bru_info$model,
      state = stateinput,
      data = datatmp,
      predictor =  fmlpred
    )
    
    library(tidytable)
    setDTthreads(10)
    tmp_new <- tmp %>%
      imap(~ data.table(.x) %>%
             mutate(id = row_number(), draws = .y) %>%
             rename(value = p)) %>%
      do.call(rbind, .) %>%
      left_join(datatmp %>%
                  st_drop_geometry() %>%
                  data.table(), by = "id")
    
    print(colnames(tmp_new))
    
    #Difference between screening pathway within virus
    count_test <- tmp_new %>%
      st_drop_geometry() %>%
      as_tidytable() %>%
      group_by(idgeometry, virusgroup, age, idtime, draws) %>%
      arrange(testgroup) %>%
      mutate(value = value - value[1]) %>%
      mutate(value = ifelse(value > 0, 1, 0)) %>%
      group_by(testgroup, virusgroup, age, idtime, draws) %>%
      summarise(value = sum(value)) %>%
      filter(testgroup == 2) %>%
      group_by(virusgroup, idtime, draws) %>%
      arrange(age) %>%
      mutate(value = value - value[1]) %>%
      filter(age > 30) %>%
      group_by(virusgroup, age, idtime) %>%
      summary_draw() %>%
      tibble()
    
    print(count_test)
    
    #Difference between genotype within screening pathway
    count_virus <- tmp_new %>%
      st_drop_geometry() %>%
      as_tidytable() %>%
      group_by(idgeometry, testgroup, age, idtime, draws) %>%
      arrange(virusgroup) %>%
      mutate(value = value - value[1]) %>%
      mutate(value = ifelse(value > 0, 1, 0)) %>%
      group_by(testgroup, virusgroup, age, idtime, draws) %>%
      summarise(value = sum(value)) %>%
      filter(virusgroup == 2) %>%
      group_by(testgroup, idtime, draws) %>%
      arrange(age) %>%
      mutate(value = value - value[1]) %>%
      filter(age > 30) %>%
      group_by(testgroup, age, idtime) %>%
      summary_draw() %>%
      tibble()
    
    print(count_virus)
    
    ###"Difference in difference"
    diff_between_screenig_between_prev_within_screening <- tmp_new %>%
      st_drop_geometry() %>%
      as_tidytable() %>%
      group_by(idgeometry, testgroup, age, idtime, draws) %>%
      arrange(virusgroup) %>%
      mutate(value = value - value[1]) %>%
      mutate(value = ifelse(value > 0, 1, 0)) %>%
      group_by(testgroup, virusgroup, age, idtime, draws) %>%
      summarise(value = sum(value)) %>%
      filter(virusgroup == 2) %>%
      group_by(age, idtime, draws) %>%
      arrange(testgroup) %>%
      mutate(value = value -  value[1]) %>%
      filter(testgroup == 2) %>%
      group_by(idtime, draws) %>%
      arrange(age) %>%
      mutate(value = value - value[1]) %>%
      group_by(age, idtime) %>%
      summary_draw() %>%
      tibble() %>%
      filter(age > 30)
    
    
    print(diff_between_screenig_between_prev_within_screening)
    detach("package:tidytable", unload = TRUE)
    
    return(
      list(
        "diff_prev_between_screening" = count_test %>% fun_format,
        "diff_prev_within_screening" = count_virus %>% fun_format,
        "diff_between_screenig_between_prev_within_screening" = diff_between_screenig_between_prev_within_screening %>% fun_format
      )
    )
  }
  
  
  list_main_city <- generate_risk_level(
    df_input = centroid_major_city,
    age_input = seq(min(dfbinomial$age), max(dfbinomial$age)),
    idtime_input = c(min(data_bru$idtime), max(dfdate$idtime)),
    city = T
  )
  
  map <- list()
  map_paris <- list()
  change_inflation_age <- list()
  if (runmap == T) {
    print("Mapping")
    
    for (a in seq(30, 66)) {
      print(a)
      map[[a]] <-   generate_risk_level(
        df_input = spatial_domain,
        age_input = a,
        idtime_input = c(max(dfdate$idtime))
      )
      
      map_paris[[a]] <-  map[[a]] %>%
        imap(function(.x, .y) {
          if (!str_detect(.y, "count")) {
            .x %>% filter(insee_reg == "11")
          }
          else{
            NULL
          }
        })
      
      gc()
    }
    
    print("Difference in count")
    if (runmap == T) {
      for (a in seq(31, 66)) {
        print(a)
        change_inflation_age[[a]] <- generate_diff_expected_prevalence_number_cp(
          df_input = spatial_domain,
          age_input = c(30, a),
          idtime_input = c(max(dfdate$idtime))
        )
        gc()
      }
    }
  }
  
  
  
  
  
  gc()
  
  # Out ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  out <- list(
    "correlation" = correlation,
    "summary.fixed" = summary.fixed,
    "summary.random" = summary.random,
    "summary.hyperpar" = summary.hyperpar,
    "prevalence" = prevalence,
    "ame" = ame,
    "map" = map,
    "map_paris" = map_paris,
    "change_inflation_age" = change_inflation_age,
    "list_main_city" = list_main_city
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
      qi_lb = quantile(value, probs = 1 - 0.975, na.rm = T),
      qi_ub = quantile(value, probs = 0.975, na.rm = T),
      qi_lb_099 = quantile(value, probs = 1 - 0.99, na.rm = T),
      qi_ub_099 = quantile(value, probs = 0.99, na.rm = T),
      qi_lb_095 = quantile(value, probs = 1 - 0.95, na.rm = T),
      qi_ub_095 = quantile(value, probs = 0.95, na.rm = T),
      qi_lb_09 = quantile(value, probs = 1 - 0.9, na.rm = T),
      qi_ub_09 = quantile(value, probs = 0.9, na.rm = T),
      qi_lb_08 = quantile(value, probs = 1 - 0.8, na.rm = T),
      qi_ub_08 = quantile(value, probs = 0.8, na.rm = T),
      qi_lb_07 = quantile(value, probs = 1 - 0.7, na.rm = T),
      qi_ub_07 = quantile(value, probs = 0.7, na.rm = T),
      qi_lb_06 = quantile(value, probs = 1 - 0.6, na.rm = T),
      qi_ub_06 = quantile(value, probs = 0.6, na.rm = T),
      qi_lb_05 = quantile(value, probs = 1 - 0.5, na.rm = T),
      qi_ub_05 = quantile(value, probs = 0.5, na.rm = T),
      min = min(value, na.rm = T),
      max = max(value, na.rm = T),
      N_draws = max(row_number()),
      N_higher_0 = sum(value > 0, na.rm = T),
      N_higher_05 = sum(value > 0.5, na.rm = T),
      N_higher_1 = sum(value > 1, na.rm = T),
      N_lower_0 = sum(value < 0, na.rm = T),
      N_lower_05 = sum(value < 0.5, na.rm = T),
      N_lower_1 = sum(value < 1, na.rm = T)
    ) %>%
    ungroup() %>%
    mutate(
      higher_0 = N_higher_0 / N_draws,
      higher_05 = N_higher_05 / N_draws,
      higher_1 = N_higher_1 / N_draws,
      lower_0 = N_lower_0 / N_draws,
      lower_05 = N_lower_05 / N_draws,
      lower_1 = N_lower_1 / N_draws
    )
}
