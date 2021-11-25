gen_country_basics <- function(country,
                               waning_nat = 52*7*3,
                               R0_assumed  = 2.7,
                               date_start = "2020-01-01",
                               date_end = "2022-12-31",
                               processes = NULL,
                               deterministic = TRUE){
  
  require(countrycode)
  
  wb_tmp = countrycode(country, "country.name", "wb")
  c_tmp = schedule_raw %>% filter(wb == wb_tmp) %>% 
    filter(date >= date_start,
           date <= date_end)
  
  para = cm_parameters_SEI3R(dem_locations = country, 
                             date_start = date_start, 
                             date_end = date_end,
                             dE  = cm_delay_gamma(2.5, 2.5, 
                                                  t_max = 15, t_step = 0.25)$p,
                             dEa = cm_delay_gamma(2.5, 2.5, 
                                                  t_max = 15, t_step = 0.25)$p,
                             dIp = cm_delay_gamma(1.5, 4.0, 
                                                  t_max = 15, t_step = 0.25)$p,
                             dIs = cm_delay_gamma(3.5, 4.0, 
                                                  t_max = 15, t_step = 0.25)$p,
                             dIa = cm_delay_gamma(5.0, 4.0, 
                                                  t_max = 15, t_step = 0.25)$p,
                             deterministic = deterministic)
  
  n_age_groups <- length(para$pop[[1]]$size)
  
  for(i in 1:length(para$pop)){
    
    para$pop[[i]]$y <- cf
    para$pop[[i]]$u <- sus
    
    # scale u (susceptibility) to achieve desired R0
    current_R0 = cm_calc_R0(para, i); # calculate R0 in population i of params
    para$pop[[i]]$u = para$pop[[i]]$u * R0_assumed / current_R0
    
    # natural waning
    para$pop[[i]]$wn <- rep((1/waning_nat), n_age_groups)
    
    ## Set seeds to control start of outbreak
    # infections start in individuals aged 20-50
    para$pop[[i]]$dist_seed_ages = 
      cm_age_coefficients(20, 
                          80, 
                          5 * (0:length(para$pop[[i]]$size))) 
    
    # 1 new infections each day for 14 days to see the outbreak
    para$pop[[i]]$seed_times <- c(1:14)
  }
  
  para$processes = burden_processes
  
  para$schedule[["mobility"]] = list(
    parameter = "contact",
    pops = numeric(),
    mode = "assign",
    values = split(c_tmp[,4:7],
                   seq(nrow(c_tmp))) %>%
      map(unlist) %>%
      map(as.vector) %>%
      unname,
    times = 1:nrow(c_tmp))

  return(para)
}

update_vac_char <- function(para,
                            ve_i  = NULL, # infection blocking VE post 1 dose
                            v2e_i = NULL, # infection blocking VE post 2 doses
                            ve_d  = NULL, # perc reduction in clinical fraction post 1 dose
                            v2e_d = NULL, # perc reduction in clinical fraction post 2 doses 
                            wv = NULL # 1/waning period
                            ){
  
  n_age <- length(para$pop[[1]]$size)

  # key parameters on infection and disease blocking mechanisms
  para$pop[[1]]$uv  <- para$pop[[1]]$u * (1 - ve_i)
  para$pop[[1]]$uv2 <- para$pop[[1]]$u * (1 - v2e_i)
  para$pop[[1]]$yv  <- para$pop[[1]]$y * (1 - ve_d)
  para$pop[[1]]$yv2 <- para$pop[[1]]$y * (1 - v2e_d)
  
  # waning vaccine-induced immunity
  para$pop[[1]]$wv = rep(wv, n_age)

  # return results
  return(para)
}

#### generate processes, from SA2UK ####
cm_multinom_process <- function(
  src, outcomes, delays,
  report = ""
) {
  if ("null" %in% names(outcomes)) {
    if (length(report) != length(outcomes)) report <- rep(report, length(outcomes))
    report[which(names(outcomes)=="null")] <- ""
    if (!("null" %in% names(delays))) {
      delays$null <- c(1, rep(0, length(delays[[1]])-1))
    }
  } else if (!all(rowSums(outcomes)==1)) {
    report <- c(rep(report, length(outcomes)), "")
    outcomes$null <- 1-rowSums(outcomes)
    delays$null <- c(1, rep(0, length(delays[[1]])-1))
  }
  nrow <- length(outcomes)
  list(
    source = src, type="multinomial", names=names(outcomes), report = report,
    prob = t(as.matrix(outcomes)), delays = t(as.matrix(delays))
  )
}

change_VOC <- function(
  para = NULL,
  date_switch = "2021-03-15",
  rc_severity = NULL, # relative change in ihr and ifr
  rc_transmissibility = NULL, # relative change in transmissibility via 
  # u, uv and uv2
  rc_ve = NULL # relative in ve against infection
){

  burden_processes_new <- 
    list(cm_multinom_process("E",       data.frame(death_voc = P.death*rc_severity),                   delays = data.frame(death_voc = delay_2death), report = "o"),
         cm_multinom_process("Ev",      data.frame(death_voc = P.death*(1-ve$ve_mort[1])*rc_severity), delays = data.frame(death_voc = delay_2death), report = "o"),
         cm_multinom_process("Ev2",     data.frame(death_voc = P.death*(1-ve$ve_mort[2])*rc_severity), delays = data.frame(death_voc = delay_2death), report = "o"),
         
         cm_multinom_process("E",       data.frame(to_hosp_voc = P.hosp*rc_severity),                  delays = data.frame(to_hosp_voc = delay_2severe)),
         cm_multinom_process("Ev",      data.frame(to_hosp_voc = P.hosp*(1-ve$ve_h[1])*rc_severity),   delays = data.frame(to_hosp_voc = delay_2severe)),
         cm_multinom_process("Ev2",     data.frame(to_hosp_voc = P.hosp*(1-ve$ve_h[2])*rc_severity),   delays = data.frame(to_hosp_voc = delay_2severe)),
         
         cm_multinom_process("to_hosp_voc", data.frame(hosp_voc = rep(1,16)),                          delays = data.frame(hosp_voc = delay_2hosp),   report = "ip"))
  
  burden_updated <- c(burden_processes, burden_processes_new)
  para$processes <- burden_updated
  
  
  # make changes RE: transmissibility and ve
  if(!(rc_transmissibility == 1 & rc_ve == 1)){
    data.table(u   =  para$pop[[1]]$u,
               uv  =  para$pop[[1]]$uv,
               uv2 =  para$pop[[1]]$uv2) %>% 
      mutate(ve_i  =  1 - uv/u,
             ve_i2 =  1 - uv2/u) -> u_table
    
    # change transmissibility
    if(rc_transmissibility != 1){
      u_table_voc <- 
      u_table %>% mutate_at(vars(starts_with("u")), function(x) x*rc_transmissibility)
    }
    
    # change in infection prevention ve
    if(rc_ve != 1){
      u_table_voc <- 
        u_table_voc %>% 
        mutate_at(vars(starts_with("ve")), function(x) x*rc_ve) %>% 
        mutate(uv = u*(1-ve_i),
               uv2 = u*(1-ve_i2))
    }
    
    for(p in c("u","uv","uv2")){
      para$schedule[[paste0("change_", p)]] <-  list(
        parameter = p,
        pops = numeric(),
        mode = "assign",
        values = list(u_table %>% pull(p), u_table_voc %>% pull(p)),
        times = c(para$date0, date_switch)
      )
    }
  }
  return(para)
}

vac_policy <- function(para,
                       # these two parameters define the supply conditions
                       milestone_date = c("2021-03-01", # start from 0
                                          "2021-06-30", # 0.03
                                          "2021-12-31", # all population; 0.2
                                          "2022-12-31"), # 0.6
                       milestone_cov = c(0,
                                         0.03,
                                         0.2,
                                         0.5),
                       # prioritisation, assume 60+  all prioritised
                       priority = c(NA, NA, NA, NA,
                                    2,  2,  2,  2,
                                    2,  2,  2,  2,
                                    1,  1,  1,  1),
                       # maximum feasible uptakes
                       cov_max = c(rep(0,2),
                                   rep(0.7, 10),
                                   rep(0.9, 4)),
                       # the proportion of age groups at which we will compare
                       # rolling out dose 2s or more dose 1s
                       # p_change = 0.7,
                       supply_delay = 20, # unit = weeks
                       dose_interval = c(4, 16)
                       
){
  date_start = para$date0
  
  # details of the vaccine roll-out schedules
  # beginning of each phase, and the daily number of doses available
  #tmp_schedule <- 
  data.table(milestone_date = c(para$date0, milestone_date, as.character(ymd(date_start) + para$time1)),
             milestone_cov = c(NA, milestone_cov, NA)) %>% 
    .[, milestone_cov := mean(milestone_cov, na.rm = T), by = list(milestone_date)] %>% 
    as_tibble %>% distinct() %>% 
    mutate(milestone_date = lubridate::ymd(milestone_date),
           t = milestone_date - as.Date(date_start),
           cov = c(NA, diff(zoo::na.locf(milestone_cov, na.rm = F))),
           cov = if_else(cov %in% c(0, NA, NaN), as.numeric(NA), cov),
           doses = cov * sum(para$pop[[1]]$size),
           t_diff = c(NA, diff(t)),
           doses_daily = dplyr::lead(doses/t_diff, 1),
           milestone_marker = !is.na(milestone_cov))-> tmp_schedule
  
  # details of the target population
  tmp_pop <- data.table(n_pop = para$pop[[1]]$size) %>% 
    mutate(uptake = cov_max,
           n_tar = n_pop * uptake)
  # 
  # # details on age-specific prioritization
  tmp_priorities <- priority %>%
    enframe(name = "age_group") %>%
    filter(!is.na(value)) %>%
    arrange(value) %>% 
    group_by(value) %>%
    group_split()

  # age group specific population size cap
  # at any given time, age groups will not have more than this many people in 
  # go from S --> V for more than this
  pop_cap <- tmp_priorities %>% 
    map(~.$age_group) %>% 
    map(~tmp_pop$n_tar[.]) 
  
  # population saturation marker
  pop_marker <- pop_cap %>%
    map(sum) %>%
    unlist()
  
  ##### scenario set 1 #####
  # individuals at their best capacity, receive vaccines at a fixed dosing 
  # interval specified by @dose_interval
  
  draw_rollout_interval <- function(interval){
    daily_vac_scenarios <- 
      matrix(0, 
             ncol = length(para$pop[[1]]$size),
             nrow = max(as.numeric(tmp_schedule$t), na.rm = T)) %>% 
      as_tibble() %>% 
      setNames(paste0("Y",1:16, "_d1")) %>% 
      bind_cols(matrix(0, 
                       ncol = length(para$pop[[1]]$size),
                       nrow = max(as.numeric(tmp_schedule$t), na.rm = T)) %>% 
                  as_tibble() %>% 
                  setNames(paste0("Y",1:16, "_d2"))) %>% 
      rownames_to_column(var = "t") %>% 
      mutate(date = lubridate::ymd(date_start) + as.numeric(t)) %>% 
      left_join(    
        tmp_schedule %>% 
          mutate(doses_daily = if_else(!milestone_marker, 0, doses_daily)) %>% 
          dplyr::filter(!is.na(doses_daily)) %>% 
          dplyr::select(milestone_date, doses_daily) %>% 
          setNames(c("date", "supply")) %>% 
          mutate(date = if_else(date == date_start, date + 1, date)),
        by = "date") %>% 
      mutate(supply = imputeTS::na_locf(supply),
             supply_2 = lag(supply, n = supply_delay*7) %>% replace(., is.na(.), 0),
             supply_daily = supply + supply_2,
             supply_cum = cumsum(supply_daily), phase = 0) %>% 
      data.table
    
    vacc_start <- min(which(daily_vac_scenarios$supply_cum > 0))
    tmp_tar <- tmp_priorities %>% map(pull, age_group) %>% map(~paste0("Y", ., "_d1"))
    tmp_tar2 <- tmp_priorities %>% map(pull, age_group) %>% map(~paste0("Y", ., "_d2"))
    tmp_pop_prop <- tmp_priorities %>% map(pull, age_group) %>%
      map(~tmp_pop$n_pop[.]) %>% map(enframe) %>% 
      map(mutate, tot = sum(value), prop = value/tot) %>% 
      map(pull, prop)
    
    data.table(t_dose1 = as.numeric(NA), t_dose2 = as.numeric(NA), status = as.character(NA), elapse = as.numeric(NA)) %>% 
      bind_cols(daily_vac_scenarios %>% dplyr::select(ends_with("_d1")) %>% .[1,]) %>% 
      mutate_at(vars(starts_with("Y", ignore.case = F)), function(x) x = as.numeric(NA))-> pending
    
    # complete vaccinating the first priority
    for(t in vacc_start:nrow(daily_vac_scenarios)){
    # for(t in vacc_start:760){
      # whether or not we allocate dose 2
      switch_dose2 <- any(pending$elapse >= interval*7 & pending$status == "incomplete", na.rm = T)
      # saturation status
      saturated <- daily_vac_scenarios %>% 
        dplyr::select(t, ends_with(c("_d1", "_d2"))) %>% 
        pivot_longer(ends_with(c("_d1", "_d2"))) %>% 
        group_by(name) %>% summarise(value = sum(value)) %>% 
        separate(name, into = c("age_group", "dose")) %>% 
        mutate(age_group = parse_number(age_group)) %>% 
        left_join(tmp_priorities %>% bind_rows() %>% rename(priority = value),
                  by = "age_group") %>%
        left_join(tmp_pop %>% rownames_to_column(var = "age_group") %>% mutate(age_group = as.numeric(age_group)),
                  by = "age_group") %>% 
        mutate(saturated = value >= n_tar) %>% 
        filter(!is.na(priority)) %>% group_by(priority, dose) %>% 
        group_split() %>% bind_rows() %>% 
        dplyr::select(priority, dose, saturated) %>% distinct() %>% data.table
      
      # if all population group has been saturated, break loop
      if(all(saturated$saturated)) warning("The entire population has been exhausted!")
      
      # identify those individuals that have "matured"
      if(!is.na(switch_dose2) & switch_dose2 == T){
        tmp <- pending %>% 
          dplyr::filter(status == "incomplete",
                        elapse >= interval*7) %>% 
          arrange(desc(elapse))
        
        # if(nrow(tmp) > 1) break (print(t, "nrow issue"));
        doses_needed <- tmp %>% pivot_longer(cols = starts_with("Y", ignore.case = F)) %>% pull(value) %>% sum
        doses_available <- daily_vac_scenarios[t,]$supply_daily
        
        # if the supply is identical @dose_interval later, this is the easiest scenario - the same 
        # array is just going to be shifted by @dose_itnerval
        if(doses_available == doses_needed){
          if(saturated[dose == "d2" & priority == 1, ]$saturated == F){
            v <- doses_available*tmp_pop_prop[[1]]; daily_vac_scenarios[t, tmp_tar2[[1]] := as.list(v)]; rm(v) 
            # daily_vac_scenarios[t,tmp_tar2[[1]]] <- t(doses_needed*tmp_pop_prop[[1]])
          }
          if(saturated[dose == "d2" & priority == 1, ]$saturated == T){
            v <- doses_available*tmp_pop_prop[[2]]; daily_vac_scenarios[t,  tmp_tar2[[2]] := as.list(v)]; rm(v)
            # daily_vac_scenarios[t,tmp_tar2[[2]]] <- t(doses_needed*tmp_pop_prop[[2]])
          }
          
          pending %<>% 
            mutate(status = if_else(!is.na(t_dose1) & t_dose1 == tmp$t_dose1,
                                    "complete",
                                    status),
                   t_dose2 = if_else(t_dose1 == tmp$t_dose1 & status == "complete",
                                     t,
                                     as.integer(t_dose2)),
                   elapse = if_else(status == "incomplete", elapse+1, elapse))
        }
        # if there is supply surplus, beyond covering everyone timely, additional individuals can be 
        # vaccinated with dose 1
        if(doses_available > doses_needed){
          # fulfill the dose 2 allocations first using doses needed
          if(saturated[dose == "d2" & priority == 1, ]$saturated == F){
            v <- doses_needed*tmp_pop_prop[[1]]; daily_vac_scenarios[t,tmp_tar2[[1]] := as.list(v)]; rm(v)
            # daily_vac_scenarios[t,tmp_tar2[[1]]] <- t(doses_needed*tmp_pop_prop[[1]])
          }
          if(saturated[dose == "d2" & priority == 1, ]$saturated == T){
            v <- doses_needed*tmp_pop_prop[[2]]; daily_vac_scenarios[t,tmp_tar2[[2]] := as.list(v)]; rm(v) 
            # daily_vac_scenarios[t,tmp_tar2[[2]]] <- t(doses_needed*tmp_pop_prop[[2]])
          }
          pending %<>% 
            mutate(status = if_else(!is.na(t_dose1) & t_dose1 %in% tmp$t_dose1 & status == "incomplete",
                                    "complete",
                                    status),
                   t_dose2 = if_else(t_dose1 %in% tmp$t_dose1 & status == "complete",
                                     as.integer(t),
                                     as.integer(t_dose2)))
          
          # move on to allocating the suprplus as dose 1
          if(saturated[dose == "d1" & priority == 1]$saturated == F){
            v <- (doses_available - doses_needed)*tmp_pop_prop[[1]]; daily_vac_scenarios[t,tmp_tar[[1]] := as.list(v)]; rm(v) 
            # daily_vac_scenarios[t,tmp_tar[[1]]] <- t((doses_available - doses_needed)*tmp_pop_prop[[1]])
          }
          if(saturated[dose == "d1" & priority == 1]$saturated == T){
            v <- (doses_available - doses_needed)*tmp_pop_prop[[2]]; daily_vac_scenarios[t,tmp_tar[[2]] := as.list(v)]; rm(v)
            # daily_vac_scenarios[t,tmp_tar[[2]]] <- t((doses_available - doses_needed)*tmp_pop_prop[[2]])
          }
          
          pending %<>% 
            mutate(elapse = if_else(status == "incomplete", elapse + 1, elapse)) %>% 
            add_row(daily_vac_scenarios[t,] %>% dplyr::select(t, ends_with("_d1", ignore.case = F)) %>% 
                      rename(t_dose1 = t) %>% 
                      mutate(elapse = 0,
                             status = "incomplete",
                             t_dose1 = as.numeric(t_dose1),
                             t_dose2 = as.numeric(NA))) 
        }
        # TO-DO
        if(doses_available < doses_needed){
          
          # fulfill the dose 2 allocations first using doses available
          if(saturated[dose == "d2" & priority == 1, ]$saturated == F){
            v <- doses_available*tmp_pop_prop[[1]]; daily_vac_scenarios[t,tmp_tar2[[1]] := as.list(v)]; rm(v)
            # daily_vac_scenarios[t,tmp_tar2[[1]]] <- t(doses_needed*tmp_pop_prop[[1]])
          }
          if(saturated[dose == "d2" & priority == 1, ]$saturated == T & t != 760){
            v <- doses_available*tmp_pop_prop[[2]]; daily_vac_scenarios[t,tmp_tar2[[2]] := as.list(v)]; rm(v) 
            # daily_vac_scenarios[t,tmp_tar2[[2]]] <- t(doses_needed*tmp_pop_prop[[2]])
          }
          if(t == 760 & sum(para$pop[[1]]$size) == 84339067){
            v <- tmp %>% dplyr::select(ends_with("_d1")) %>% .[1,] %>% t %>% c; v <- v[v!=0]
            doses_remain <- doses_available - sum(v)
            daily_vac_scenarios[t,tmp_tar2[[1]] := as.list(v)]; rm(v)
            v <- doses_remain*tmp_pop_prop[[2]]; 
            daily_vac_scenarios[t,tmp_tar2[[2]] := as.list(v)]; rm(v) 
            rm(doses_remain)
          }
   
          doses_counter <- sapply(1:nrow(tmp), function(x)tmp[1:x,5:20] %>% sum)
          covered <- sapply(1:nrow(tmp), function(x) doses_counter[x] < doses_available)
          if(!any(covered)) expand <- 1
          if(any(covered)) expand <- min(which(covered)) + 1
          
          for(j in 1:expand){
            doses_counter <- tmp[1:j,5:20] %>% sum
            covered <- doses_counter < doses_available
            if(covered){
              pending[t_dose1 %in% tmp$t_dose1[j] & status == "incomplete", c("status","t_dose2") := list("complete", t)]
            }
            if(!covered){
              total_allocated <- daily_vac_scenarios[t,] %>% dplyr::select(ends_with("_d2"))
              if(j == 1) previ_allocated <- tmp[j,] %>% dplyr::select(ends_with("_d1")) %>% mutate_all(function(x) x = 0)
              if(j > 1) previ_allocated <- tmp[1:(j-1),] %>% dplyr::select(ends_with("_d1")) %>% colSums()
              to_allocate <- total_allocated - previ_allocated
              cannot_allocate <- tmp[j,]%>% dplyr::select(ends_with("_d1")) %>% colSums() - to_allocate
              
              pending %<>%
                filter(!t_dose1 %in% tmp$t_dose1[j]) %>% 
                add_row(tmp[j,1:4] %>% bind_cols(to_allocate) %>%
                          setNames(colnames(pending)) %>%
                          mutate(t_dose2 = t,
                                 status = "complete")
                ) %>%
                add_row(
                  tmp[j,1:4] %>% bind_cols(cannot_allocate) %>%
                    setNames(colnames(pending))
                ) %>% 
                mutate(elapse = if_else(status == "incomplete", elapse+1, elapse)) 
            }
          }
        }
        }
      
      if(!is.na(switch_dose2) & switch_dose2 == F){
        if (saturated[dose == "d1" & priority == 1,]$saturated == FALSE) {
          v <- daily_vac_scenarios[t,]$supply_daily*tmp_pop_prop[[1]]; daily_vac_scenarios[t, tmp_tar[[1]] := as.list(v)]; rm(v) 
          # daily_vac_scenarios[t,tmp_tar[[1]]] <- t(daily_vac_scenarios[t,]$supply_daily*tmp_pop_prop[[1]])
        }
        
        if (saturated[dose == "d1" & priority == 1,]$saturated == TRUE &
            saturated[dose == "d1" & priority == 2,]$saturated == FALSE) {
          v <- daily_vac_scenarios[t,]$supply_daily*tmp_pop_prop[[2]];  daily_vac_scenarios[t,tmp_tar[[2]] := as.list(v)]; rm(v)
          # daily_vac_scenarios[t,tmp_tar[[2]]] <- t(daily_vac_scenarios[t,]$supply_daily*tmp_pop_prop[[2]])
        }
        
        pending %<>% 
          mutate(elapse = if_else(status == "incomplete", elapse + 1, elapse)) %>% 
          add_row(daily_vac_scenarios[t,] %>% dplyr::select(t, ends_with("_d1", ignore.case = F)) %>% 
                    rename(t_dose1 = t) %>% 
                    mutate(elapse = 0,
                           status = "incomplete",
                           t_dose1 = as.numeric(t_dose1),
                           t_dose2 = as.numeric(NA)))
      }
      
      if(length(which(pending %>% pivot_longer(starts_with("Y", ignore.case = )) %>% pull(value) < 0)) > 0) {
        print("There are negative values of vaccines being allocated. Please check the parameters or code.")
        break
      }
    }
    
    return(list(daily_vac_scenarios = daily_vac_scenarios,
                pending = pending))
  }
  
  scenarios <- list()
  for(i in 1:length(dose_interval)){
    scenarios[[i]] <- draw_rollout_interval(dose_interval[i])
  }
  
  ##### scenario set 2 #####
  # individuals will receive their second dose when the uptake goal has been reached for 
  # the entire population

  draw_rollout_first <- function(label_interval){
    daily_vac_scenarios <- 
      matrix(0, 
             ncol = length(para$pop[[1]]$size),
             nrow = max(as.numeric(tmp_schedule$t), na.rm = T)) %>% 
      as_tibble() %>% 
      setNames(paste0("Y",1:16, "_d1")) %>% 
      bind_cols(matrix(0, 
                       ncol = length(para$pop[[1]]$size),
                       nrow = max(as.numeric(tmp_schedule$t), na.rm = T)) %>% 
                  as_tibble() %>% 
                  setNames(paste0("Y",1:16, "_d2"))) %>% 
      rownames_to_column(var = "t") %>% 
      mutate(date = lubridate::ymd(date_start) + as.numeric(t)) %>% 
      left_join(    
        tmp_schedule %>% 
          mutate(doses_daily = if_else(!milestone_marker, 0, doses_daily)) %>% 
          dplyr::filter(!is.na(doses_daily)) %>% 
          dplyr::select(milestone_date, doses_daily) %>% 
          setNames(c("date", "supply")) %>% 
          mutate(date = if_else(date == date_start, date + 1, date)),
        by = "date") %>% 
      mutate(supply = imputeTS::na_locf(supply),
             supply_2 = lag(supply, n = supply_delay*7) %>% replace(., is.na(.), 0),
             supply_daily = supply + supply_2,
             supply_cum = cumsum(supply_daily), phase = 0) %>% 
      data.table
    
    vacc_start <- min(which(daily_vac_scenarios$supply_cum > 0))
    tmp_tar <- tmp_priorities %>% map(pull, age_group) %>% map(~paste0("Y", ., "_d1"))
    tmp_tar2 <- tmp_priorities %>% map(pull, age_group) %>% map(~paste0("Y", ., "_d2"))
    tmp_pop_prop <- tmp_priorities %>% map(pull, age_group) %>%
      map(~tmp_pop$n_pop[.]) %>% map(enframe) %>% 
      map(mutate, tot = sum(value), prop = value/tot) %>% 
      map(pull, prop)
    
    data.table(t_dose1 = as.numeric(NA), t_dose2 = as.numeric(NA), status = as.character(NA), elapse = as.numeric(NA)) %>% 
      bind_cols(daily_vac_scenarios %>% dplyr::select(ends_with("_d1")) %>% .[1,]) %>% 
      mutate_at(vars(starts_with("Y", ignore.case = F)), function(x) x = as.numeric(NA))-> pending
    
    #nrow(daily_vac_scenarios)
    for(t in vacc_start:nrow(daily_vac_scenarios)){
      
      saturated <- daily_vac_scenarios %>% 
        dplyr::select(t, ends_with(c("_d1", "_d2"))) %>% 
        pivot_longer(ends_with(c("_d1", "_d2"))) %>% 
        group_by(name) %>% summarise(value = sum(value), .groups = "drop") %>% 
        separate(name, into = c("age_group", "dose")) %>% 
        mutate(age_group = parse_number(age_group)) %>% 
        left_join(tmp_priorities %>% bind_rows() %>% rename(priority = value),
                  by = "age_group") %>%
        left_join(tmp_pop %>% rownames_to_column(var = "age_group") %>% mutate(age_group = as.numeric(age_group)),
                  by = "age_group") %>% 
        mutate(saturated = value >= n_tar) %>%
        filter(!is.na(priority)) %>% group_by(priority, dose) %>% 
        group_split() %>% bind_rows() %>% 
        dplyr::select(priority, dose, saturated) %>% distinct() %>% data.table
      
      # if(saturated$saturated[4] == T) break
      doses_available <- daily_vac_scenarios[t,]$supply_daily
      
      # if not group has been saturated we vaccinate priority = 1, d1
      if(all(!saturated$saturated)){
        v <- doses_available*tmp_pop_prop[[1]]; daily_vac_scenarios[t, tmp_tar[[1]] := as.list(v)]; rm(v)
        pending %<>% 
          add_row(daily_vac_scenarios[t,] %>% dplyr::select(ends_with("_d1")) %>% 
                    mutate(t_dose1 = t, t_dose2 = as.numeric(NA), 
                           status = "incomplete", elapse = 0) %>% dplyr::select(colnames(pending)))
      }else if(saturated[priority == 1 & dose == "d1"]$saturated == T &
               saturated[priority == 2 & dose == "d1"]$saturated == F &
               saturated[priority == 1 & dose == "d2"]$saturated == F &
               saturated[priority == 2 & dose == "d2"]$saturated == F ){
        v <- doses_available*tmp_pop_prop[[2]]; daily_vac_scenarios[t, tmp_tar[[2]] := as.list(v)]; rm(v)
        pending %<>% 
          add_row(daily_vac_scenarios[t,] %>% dplyr::select(ends_with("_d1")) %>% 
                    mutate(t_dose1 = t, t_dose2 = as.numeric(NA), 
                           status = "incomplete", elapse = 0) %>% dplyr::select(colnames(pending)))
      }else if(saturated[priority == 1 & dose == "d1"]$saturated == T &
               saturated[priority == 2 & dose == "d1"]$saturated == T &
               saturated[priority == 1 & dose == "d2"]$saturated == F &
               saturated[priority == 2 & dose == "d2"]$saturated == F){
        
        v <- doses_available*tmp_pop_prop[[1]]; daily_vac_scenarios[t, tmp_tar2[[1]] := as.list(v)]; rm(v) 
        
        # identify the ones that are available to receive the second dose
        pending %>% filter(elapse > label_interval[1] & status == "incomplete") %>% arrange(desc(elapse)) %>%
          pivot_longer(cols = starts_with("Y", ignore.case = F)) %>%
          group_by(t_dose1, t_dose2, status, elapse) %>% summarise(value = sum(value), .groups = "drop") %>%
          arrange(t_dose1) %>% ungroup %>%
          mutate(cs = cumsum(value),  rcs = cs - doses_available, rcs = rcs/abs(rcs)) -> tmp
        
        if(all(tmp$rcs > 0)) nr <- 1; if(!all(tmp$rcs > 0)) nr <- max(which(tmp$rcs < 0)) + 1
        
        for(r in 1:nr){
          if(r != nr | nrow(tmp) == 1){
            pending[t_dose1 == tmp$t_dose1[r] & status == "incomplete","t_dose2"] <- t
            pending[t_dose1 == tmp$t_dose1[r] & status == "incomplete","status"] <- "complete"
          }
          if(r == nr & nrow(tmp) > 1){
            doses_left <- if_else(r > 1, doses_available - sum(tmp$value[1:(r-1)]), doses_available)
            slice1 <-  pending[t_dose1 == tmp$t_dose1[r] & status == "incomplete",]
            slice1 %>% pivot_longer(cols = starts_with("Y")) %>% 
              filter(value > 1) %>% pull(name) -> age_group_marker

            if(identical(tmp_tar[[1]], age_group_marker)){
              v <- doses_left*tmp_pop_prop[[1]]; slice1[1, tmp_tar[[1]] := as.list(v)]; rm(v)
            }
            if(identical(tmp_tar[[2]], age_group_marker)){
              v <- doses_left*tmp_pop_prop[[2]]; slice1[1, tmp_tar[[2]] := as.list(v)]; rm(v)
            }
            
            slice1 <- slice1 %>% mutate(status = "complete", t_dose2  = t)
            slice2 <- pending[t_dose1 == tmp$t_dose1[r] & status == "incomplete",]
            slice2 <- data.table(t_dose1 = tmp$t_dose1[r], t_dose2 = as.numeric(NA), status = "incomplete",elapse = tmp$elapse[r]) %>% 
              bind_cols((slice2 %>% dplyr::select(ends_with("_d1"))) - (slice1 %>% dplyr::select(ends_with("_d1"))))
            
            pending %<>% 
              filter(t_dose1 != tmp$t_dose1[r]) %>% 
              bind_rows(bind_rows(slice1, slice2)) %>% arrange(t_dose1)
          }
        }
      }else if(saturated[priority == 1 & dose == "d1"]$saturated == T &
               saturated[priority == 2 & dose == "d1"]$saturated == T &
               saturated[priority == 1 & dose == "d2"]$saturated == T &
               saturated[priority == 2 & dose == "d2"]$saturated == F ){
        v <- doses_available*tmp_pop_prop[[2]]; daily_vac_scenarios[t, tmp_tar2[[2]] := as.list(v)]; rm(v) 
        
        pending %>% filter(elapse > label_interval & status == "incomplete") %>% arrange(desc(elapse)) %>%
          pivot_longer(cols = starts_with("Y", ignore.case = F)) %>%
          group_by(t_dose1, t_dose2, status, elapse) %>% summarise(value = sum(value), .groups = "drop") %>%
          arrange(t_dose1) %>% ungroup %>%
          mutate(cs = cumsum(value),  rcs = cs - doses_available, rcs = rcs/abs(rcs)) -> tmp
        
        if(all(tmp$rcs > 0)) nr <- 1; if(!all(tmp$rcs > 0)) nr <- max(which(tmp$rcs < 0)) + 1
        
        for(r in 1:nr){
          if(r != nr | nrow(tmp) == 1){
            pending[t_dose1 == tmp$t_dose1[r] & status == "incomplete","t_dose2"] <- t
            pending[t_dose1 == tmp$t_dose1[r] & status == "incomplete","status"] <- "complete"
          }
          if(r == nr & nrow(tmp) > 1){
            doses_left <- if_else(r > 1, doses_available - sum(tmp$value[1:(r-1)]), doses_available)
            slice1 <-  pending[t_dose1 == tmp$t_dose1[r] & status == "incomplete",]
            v <- doses_left*tmp_pop_prop[[2]]; slice1[1, tmp_tar[[2]] := as.list(v)]; rm(v)
            slice1 <- slice1 %>% mutate(status = "complete", t_dose2  = t)
            slice2 <-pending[t_dose1 == tmp$t_dose1[r] & status == "incomplete",]
            slice2 <- data.table(t_dose1 = tmp$t_dose1[r], t_dose2 = as.numeric(NA), status = "incomplete",elapse = tmp$elapse[r]) %>% 
              bind_cols((slice2 %>% dplyr::select(ends_with("_d1"))) - (slice1 %>% dplyr::select(ends_with("_d1"))))
            
            pending %<>% 
              filter(t_dose1 != tmp$t_dose1[r]) %>% 
              bind_rows(bind_rows(slice1, slice2)) %>% arrange(t_dose1)
          }
        }
      }
      pending %<>% mutate(elapse = if_else(status == "incomplete", elapse + 1, elapse))
      # make sure there's only positive values in pending allocation
      expect_true(all(pending %>% filter(!is.na(t_dose1)) %>% pivot_longer(starts_with("Y", ignore.case = T)) %>% pull(value) >= 0))
      # make sure there's only positive values in daily_vac_scenarios
      expect_true(all(daily_vac_scenarios %>% pivot_longer(starts_with("Y", ignore.case = T)) %>% pull(value) >= 0))
    }
    
    if(all(saturated$saturated)) print("Vaccine surplus exists.")
    return(list(daily_vac_scenarios = daily_vac_scenarios,
                pending = pending))
  }
  
  scenarios[[(length(scenarios) + 1)]] <- draw_rollout_first(28)
 
  ##### Scenario set 3 #####
  # Rollout complete vaccination programs by priority groups
  draw_rollout_prior <- function(label_interval) {
    daily_vac_scenarios <-
      matrix(0,
             ncol = length(para$pop[[1]]$size),
             nrow = max(as.numeric(tmp_schedule$t), na.rm = T)) %>%
      as_tibble() %>%
      setNames(paste0("Y", 1:16, "_d1")) %>%
      bind_cols(matrix(
        0,
        ncol = length(para$pop[[1]]$size),
        nrow = max(as.numeric(tmp_schedule$t), na.rm = T)
      ) %>%
        as_tibble() %>%
        setNames(paste0("Y", 1:16, "_d2"))) %>%
      rownames_to_column(var = "t") %>%
      mutate(date = lubridate::ymd(date_start) + as.numeric(t)) %>%
      left_join(
        tmp_schedule %>%
          mutate(doses_daily = if_else(!milestone_marker, 0, doses_daily)) %>%
          dplyr::filter(!is.na(doses_daily)) %>%
          dplyr::select(milestone_date, doses_daily) %>%
          setNames(c("date", "supply")) %>%
          mutate(date = if_else(date == date_start, date + 1, date)),
        by = "date"
      ) %>%
      mutate(
        supply = imputeTS::na_locf(supply),
        supply_2 = lag(supply, n = supply_delay * 7) %>% replace(., is.na(.), 0),
        supply_daily = supply + supply_2,
        supply_cum = cumsum(supply_daily),
        phase = 0
      ) %>%
      data.table
    
    vacc_start <- min(which(daily_vac_scenarios$supply_cum > 0))
    tmp_tar <-
      tmp_priorities %>% map(pull, age_group) %>% map( ~ paste0("Y", ., "_d1"))
    tmp_tar2 <-
      tmp_priorities %>% map(pull, age_group) %>% map( ~ paste0("Y", ., "_d2"))
    tmp_pop_prop <- tmp_priorities %>% map(pull, age_group) %>%
      map( ~ tmp_pop$n_pop[.]) %>% map(enframe) %>%
      map(mutate, tot = sum(value), prop = value / tot) %>%
      map(pull, prop)
    
    data.table(
      t_dose1 = as.numeric(NA),
      t_dose2 = as.numeric(NA),
      status = as.character(NA),
      elapse = as.numeric(NA)
    ) %>%
      bind_cols(daily_vac_scenarios %>% dplyr::select(ends_with("_d1")) %>% .[1, ]) %>%
      mutate_at(vars(starts_with("Y", ignore.case = F)), function(x)
        x = as.numeric(NA)) -> pending
    
    for (t in vacc_start:nrow(daily_vac_scenarios)) {
      saturated <- daily_vac_scenarios %>%
        dplyr::select(t, ends_with(c("_d1", "_d2"))) %>%
        pivot_longer(ends_with(c("_d1", "_d2"))) %>%
        group_by(name) %>% summarise(value = sum(value), .groups = "drop") %>%
        separate(name, into = c("age_group", "dose")) %>%
        mutate(age_group = parse_number(age_group)) %>%
        left_join(tmp_priorities %>% bind_rows() %>% rename(priority = value),
                  by = "age_group") %>%
        left_join(
          tmp_pop %>% rownames_to_column(var = "age_group") %>% mutate(age_group = as.numeric(age_group)),
          by = "age_group"
        ) %>%
        mutate(saturated = value >= n_tar) %>%
        filter(!is.na(priority)) %>% group_by(priority, dose) %>%
        group_split() %>% bind_rows() %>%
        dplyr::select(priority, dose, saturated) %>% distinct() %>% data.table
      
      # if(saturated$saturated[3] == T) break
      doses_available <- daily_vac_scenarios[t, ]$supply_daily
      
      # if not group has been saturated we vaccinate priority = 1, d1
      if (all(!saturated$saturated)) {
        v <-
          doses_available * tmp_pop_prop[[1]]
        daily_vac_scenarios[t, tmp_tar[[1]] := as.list(v)]
        rm(v)
        pending %<>%
          add_row(
            daily_vac_scenarios[t, ] %>% dplyr::select(ends_with("_d1")) %>%
              mutate(
                t_dose1 = t,
                t_dose2 = as.numeric(NA),
                status = "incomplete",
                elapse = 0
              ) %>% dplyr::select(colnames(pending))
          )
      } else if (saturated[priority == 1 &
                           dose == "d1"]$saturated == T &
                 saturated[priority == 2 &
                           dose == "d1"]$saturated == F &
                 saturated[priority == 1 &
                           dose == "d2"]$saturated == F &
                 saturated[priority == 2 &
                           dose == "d2"]$saturated == F) {
        
        v <- doses_available * tmp_pop_prop[[1]]; daily_vac_scenarios[t, tmp_tar2[[1]] := as.list(v)]; rm(v)
        
        # identify the ones that are available to receive the second dose
        pending %>% filter(elapse > label_interval &
                             status == "incomplete") %>% arrange(desc(elapse)) %>%
          pivot_longer(cols = starts_with("Y", ignore.case = F)) %>%
          group_by(t_dose1, t_dose2, status, elapse) %>% summarise(value = sum(value), .groups = "drop") %>%
          arrange(t_dose1) %>% ungroup %>%
          mutate(
            cs = cumsum(value),
            rcs = cs - doses_available,
            rcs = rcs / abs(rcs)
          ) -> tmp
        
        if (all(tmp$rcs > 0)) nr <- 1
        if (!all(tmp$rcs > 0)) nr <- max(which(tmp$rcs < 0)) + 1
        if(nrow(tmp) == 1) nr <- 1
        
        for (r in 1:nr) {
          if (r != nr | nrow(tmp) == 1) {
            pending[t_dose1 == tmp$t_dose1[r] & status == "incomplete", "t_dose2"] <- t
            pending[t_dose1 == tmp$t_dose1[r] & status == "incomplete", "status"] <- "complete"
          }
          if (r == nr & nrow(tmp) > 1) {
            doses_left <-
              if_else(r > 1,
                      doses_available - sum(tmp$value[1:(r - 1)]),
                      doses_available)
            slice1 <-
              pending[t_dose1 == tmp$t_dose1[r] & status == "incomplete", ]
            v <-
              doses_left * tmp_pop_prop[[1]]
            slice1[1, tmp_tar[[1]] := as.list(v)]
            rm(v)
            slice1 <-
              slice1 %>% mutate(status = "complete", t_dose2  = t)
            slice2 <-
              pending[t_dose1 == tmp$t_dose1[r] & status == "incomplete", ]
            slice2 <-
              data.table(
                t_dose1 = tmp$t_dose1[r],
                t_dose2 = as.numeric(NA),
                status = "incomplete",
                elapse = tmp$elapse[r]
              ) %>%
              bind_cols((slice2 %>% dplyr::select(ends_with("_d1"))) - (slice1 %>% dplyr::select(ends_with("_d1"))))
            
            pending %<>%
              filter(t_dose1 != tmp$t_dose1[r]) %>%
              bind_rows(bind_rows(slice1, slice2)) %>% arrange(t_dose1)
          }
        }
        
      } else if (saturated[priority == 1 &
                           dose == "d1"]$saturated == T &
                 saturated[priority == 2 &
                           dose == "d1"]$saturated == F &
                 saturated[priority == 1 &
                           dose == "d2"]$saturated == T &
                 saturated[priority == 2 &
                           dose == "d2"]$saturated == F) {
        v <-
          doses_available * tmp_pop_prop[[2]]
        daily_vac_scenarios[t, tmp_tar[[2]] := as.list(v)]
        rm(v)
        
        # TO-DO: the transition between age groups could have been done in a more
        # sophisticated way so there are transition days where a little bit of
        # every age group can be vaccinated on the same day.
        
        pending %<>%
          add_row(
            daily_vac_scenarios[t, ] %>% dplyr::select(ends_with("_d1")) %>%
              mutate(
                t_dose1 = t,
                t_dose2 = as.numeric(NA),
                status = "incomplete",
                elapse = 0
              ) %>% dplyr::select(colnames(pending))
          )
        
      } else if (saturated[priority == 1 &
                           dose == "d1"]$saturated == T &
                 saturated[priority == 2 &
                           dose == "d1"]$saturated == T &
                 saturated[priority == 1 &
                           dose == "d2"]$saturated == T &
                 saturated[priority == 2 &
                           dose == "d2"]$saturated == F) {
        v <-
          doses_available * tmp_pop_prop[[2]]
        daily_vac_scenarios[t, tmp_tar2[[2]] := as.list(v)]
        rm(v)
        
        pending %>% filter(elapse > label_interval &
                             status == "incomplete") %>% arrange(desc(elapse)) %>%
          pivot_longer(cols = starts_with("Y", ignore.case = F)) %>%
          group_by(t_dose1, t_dose2, status, elapse) %>% summarise(value = sum(value), .groups = "drop") %>%
          arrange(t_dose1) %>% ungroup %>%
          mutate(
            cs = cumsum(value),
            rcs = cs - doses_available,
            rcs = rcs / abs(rcs)
          ) -> tmp
        
        if (all(tmp$rcs > 0))
          nr <- 1
        if (!all(tmp$rcs > 0))
          nr <- max(which(tmp$rcs < 0)) + 1
        
        for (r in 1:nr) {
          if (r != nr | nrow(tmp) == 1) {
            pending[t_dose1 == tmp$t_dose1[r] &
                      status == "incomplete", "t_dose2"] <- t
            pending[t_dose1 == tmp$t_dose1[r] &
                      status == "incomplete", "status"] <- "complete"
          }
          if (r == nr & nrow(tmp) > 1) {
            doses_left <-
              if_else(r > 1,
                      doses_available - sum(tmp$value[1:(r - 1)]),
                      doses_available)
            slice1 <-
              pending[t_dose1 == tmp$t_dose1[r] & status == "incomplete", ]
            v <-
              doses_left * tmp_pop_prop[[2]]
            slice1[1, tmp_tar[[2]] := as.list(v)]
            rm(v)
            slice1 <-
              slice1 %>% mutate(status = "complete", t_dose2  = t)
            slice2 <-
              pending[t_dose1 == tmp$t_dose1[r] & status == "incomplete", ]
            slice2 <-
              data.table(
                t_dose1 = tmp$t_dose1[r],
                t_dose2 = as.numeric(NA),
                status = "incomplete",
                elapse = tmp$elapse[r]
              ) %>%
              bind_cols((slice2 %>% dplyr::select(ends_with("_d1"))) - (slice1 %>% dplyr::select(ends_with("_d1"))))
            
            pending %<>%
              filter(t_dose1 != tmp$t_dose1[r]) %>%
              bind_rows(bind_rows(slice1, slice2)) %>% arrange(t_dose1)
          }
        }
      }
      if (all(saturated$saturated)) {
        pending %<>%
          mutate(
            t_dose2 = if_else(status == "incomplete", as.double(t), t_dose2),
            status = if_else(status == "incomplete", "complete", status)
          )
        # print("Vaccine surplus exists.")
      }
      pending %<>% mutate(elapse = if_else(status == "incomplete", elapse + 1, elapse))
      
      # make sure there's only positive values in pending allocation
      expect_true(all(
        pending %>% filter(!is.na(t_dose1)) %>% pivot_longer(starts_with("Y", ignore.case = T)) %>% pull(value) >= 0
      ))
      # make sure there's only positive values in daily_vac_scenarios
      expect_true(all(
        daily_vac_scenarios %>% pivot_longer(starts_with("Y", ignore.case = T)) %>% pull(value) >= 0
      ))
    }
    
    return(list(daily_vac_scenarios = daily_vac_scenarios,
                pending = pending))
  }
  
  scenarios[[(length(scenarios) + 1)]] <- draw_rollout_prior(28)  

  #### extract vac_para ####
  vac_para <- list()
  # first dose
  scenarios %>% 
    map(~.$daily_vac_scenarios) %>% 
    map(select, t, ends_with("d1"), date, supply_daily, supply_cum) %>% 
    bind_rows(.id = "scenario") %>% 
    group_by(supply_daily, scenario) %>% group_split() %>% 
    map(group_by_at, vars(starts_with("Y"))) %>% 
    map(filter, t == min(t)) %>% 
    bind_rows() %>% group_by(scenario) %>% group_split() %>% 
    map(mutate, t = parse_number(t)) %>% 
    map(arrange, t) -> vac_para[[1]]
  
  # second dose
  scenarios %>% 
    map(~.$daily_vac_scenarios) %>% 
    map(select, t, ends_with("d2"), date, supply_daily, supply_cum) %>% 
    bind_rows(.id = "scenario") %>% 
    group_by(supply_daily, scenario) %>% group_split() %>% 
    map(group_by_at, vars(starts_with("Y"))) %>% 
    map(filter, t == min(t)) %>% 
    bind_rows() %>% group_by(scenario) %>% group_split() %>% 
    map(mutate, t = parse_number(t)) %>% 
    map(arrange, t) -> vac_para[[2]]
    
      
  vacc_vals <- list()
  # then convert these parameters to a format that's friendly with `covidm`
  # allocation
  for(i in 1:2){
    vacc_vals[[i]] <- list()
    for(j in 1:length(scenarios)){
      vac_para[[i]][[j]] %>% 
        select(starts_with("Y", ignore.case = T)) %>% 
        split(seq(nrow( vac_para[[i]][[j]] ))) %>% 
        map(unlist) %>% 
        map(as.vector) -> vacc_vals[[i]][[j]] 
    }
  }

  
  # timing
  vacc_times <- list()
  # then convert these parameters to a format that's friendly with `covidm`
  # allocation
  for(i in 1:2){
    vacc_times[[i]] <- list()
    for(j in 1:length(scenarios)){
      vac_para[[i]][[j]] %>% 
        pull(t) %>% 
        as.numeric -> vacc_times[[i]][[j]] 
      vacc_times[[i]][[j]][1] <- 0
    }
  }
  
  res <- list()
  for(i in 1:length(scenarios)){
    res[[i]] <- para
    res[[i]]$schedule[["v"]] =  
      list(
      parameter = "v",
      pops = numeric(),
      mode = "assign",
      values = vacc_vals[[1]][[i]],
      times = vacc_times[[1]][[i]]
      )
    
    res[[i]]$schedule[["v2"]] =  
      list(
        parameter = "v2",
        pops = numeric(),
        mode = "assign",
        values = vacc_vals[[2]][[i]],
        times = vacc_times[[2]][[i]]
      )
  }
  
  return(list(param = para, # baseline
              res = res, # four strategies
              scenarios = scenarios))
  }


##### expand VE estimates to meet the needs of the model ####
exp_ve <- function(ve_d_o,  # disease blocking VE observed
                   ve_i_o   # infection blocking VE assumed
){
  # calculate of clinical fraction reduction
  ve_d <- (ve_d_o - ve_i_o)/(1 - ve_i_o)
  return(ve_d)
}





