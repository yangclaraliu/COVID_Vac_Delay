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
                          50, 
                          5 * (0:length(para$pop[[i]]$size))) 
    
    # 1 new infections each day for 14 days to see the outbreak
    para$pop[[i]]$seed_times <- c(1:14)
  }
  
  para$processes = processes
  
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
  date_swtich = "2021-03-15",
  rc_severity = 1.5, # relative change in ihr and ifr
  rc_transmissibility = 1.5, # relative change in transmissibility via 
  # u, uv and uv2
  rc_ve = 0.4 # relative in ve against infection
){

  if(rc_severity != 1){
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
  }
  
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
        times = c(para$date0, date_swtich)
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
                       p_change = 0.7,
                       supply_delay = 20 # unit = weeks
                       
){
  date_start = para$date0
  
  # details of the vaccine roll-out schedules
  # beginning of each phase, and the daily number of doses available
  #tmp_schedule <- 
  data.frame(milestone_date = c(para$date0, 
                                milestone_date, 
                                as.character(ymd(date_start) + para$time1)),
             milestone_cov = c(NA, 
                               milestone_cov, 
                               NA)) %>% 
    group_by(milestone_date) %>% 
    summarise(milestone_cov = mean(milestone_cov, 
                                   na.rm = T),
              .groups = "drop") %>% 
    ungroup %>% 
    mutate(milestone_date = lubridate::ymd(milestone_date),
           t = milestone_date - as.Date(para$date0),
           cov = c(NA, diff(zoo::na.locf(milestone_cov, 
                                         na.rm = F))),
           cov = if_else(cov %in% c(0, NA, NaN), 
                         as.numeric(NA), 
                         cov),
           # assume the goal is to fully cover x% of the population given a two
           # dose schedule
           doses = cov * sum(para$pop[[1]]$size),
           t_diff = c(NA, diff(t)) ,
           doses_daily = dplyr::lead(doses/t_diff, 1),
           milestone_marker = !is.na(milestone_cov)) -> tmp_schedule
  
  # details of the target population
  tmp_pop <- data.frame(n_pop = para$pop[[1]]$size) %>% 
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
  
  # create empty grid to score vaccination plans
  # these "YXX_dX" are *cumulative*
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
           supply_cum = cumsum(supply_daily), phase = 0) -> daily_vac
  
  daily_vac %>% 
    dplyr::select(date, supply, supply_2, supply_daily, supply_cum) %>%
    pivot_longer(starts_with("supply")) %>% 
    mutate(metric = if_else(name %in% c("supply", "supply_2", "supply_daily"), "daily", "cumulative"),
           name = case_when(name == "supply_daily" ~ "Total Daily Supply",
                            name == "supply_2" ~ "Daily Supply of Dose 2",
                            name == "supply" ~ "Daily Suply of Dose 1",
                            name == "supply_cum" ~ "Total Cumulative Supply")) %>% 
    filter(metric == "daily") %>% 
    ggplot(., aes(x = date, y = value, group = name, color = name)) +
    geom_line(size = 2) +
    geom_vline(xintercept = ymd(milestone_date), linetype = 2) +
    facet_wrap(~name, ncol = 1) +
    theme_bw() +
    labs(x = "", y="", color = "") +
    theme(legend.position = "top") -> p_supply
  
  ##### scenario 1 #####
  # (1) complete vaccinating 60+ with dose 1, reaching the uptake goal 
  # (2) complete vaccinating 60+ with dose 2
  # (3) vaccinate other adults with dose 1
  
  daily_vac_scenarios <- list()
  daily_vac_scenarios[[1]] <- daily_vac
  tmp_tar <- c(paste0("Y", tmp_priorities[[1]]$age_group,"_d1"))
  tmp_tar2 <- c(paste0("Y", tmp_priorities[[1]]$age_group,"_d2"))
  tmp_pop_prop <- tmp_pop$n_pop[tmp_priorities[[1]]$age_group]/sum(tmp_pop$n_pop[tmp_priorities[[1]]$age_group])
  t_marker <- which(daily_vac_scenarios[[1]]$supply_cum < pop_marker[[1]])
  for(i in 1:length(tmp_tar)){
    daily_vac_scenarios[[1]][t_marker,tmp_tar[i]] <-
      daily_vac_scenarios[[1]][t_marker, "supply_daily"]*tmp_pop_prop[i]
  }
  # remainder of dose 1
  daily_vac_scenarios[[1]][max(t_marker)+1,tmp_tar] <- 
    as.list(pop_cap[[1]]  - 
              (daily_vac_scenarios[[1]][,tmp_tar] %>% 
                 colSums()))
  # move on to dose 2
  daily_vac_scenarios[[1]][max(t_marker)+1,tmp_tar2] <- 
    as.list(tmp_pop_prop*unlist(daily_vac_scenarios[[1]][max(t_marker)+1,"supply_daily"] - 
                                  sum(daily_vac_scenarios[[1]][max(t_marker)+1,tmp_tar])))
  # allocate dose 2 to 60+
  t_marker2 <- which(daily_vac_scenarios[[1]]$supply_cum < pop_marker[[1]]*2) %>% 
    .[!.%in%t_marker] %>% 
    .[.!=(max(t_marker) + 1)]
  for(i in 1:length(tmp_tar2)){
    daily_vac_scenarios[[1]][t_marker2,tmp_tar2[i]] <-
      daily_vac_scenarios[[1]][t_marker2, "supply_daily"]*tmp_pop_prop[i]
  }
  daily_vac_scenarios[[1]][max(t_marker2)+1,tmp_tar2] <- 
    as.list(pop_cap[[1]]  - 
              (daily_vac_scenarios[[1]][,tmp_tar2] %>% 
                 colSums()))
  # move on to dose 1 for other adults
  tmp2_tar <- c(paste0("Y", tmp_priorities[[2]]$age_group,"_d1"))
  tmp2_tar2 <- c(paste0("Y", tmp_priorities[[2]]$age_group,"_d2"))
  tmp2_pop_prop <- tmp_pop$n_pop[tmp_priorities[[2]]$age_group]/sum(tmp_pop$n_pop[tmp_priorities[[2]]$age_group])
  daily_vac_scenarios[[1]][max(t_marker2)+1,tmp2_tar] <- 
    as.list(tmp2_pop_prop*unlist(daily_vac_scenarios[[1]][max(t_marker2)+1,"supply_daily"] - 
                                   sum(daily_vac_scenarios[[1]][max(t_marker2)+1,tmp_tar2])))
  # allocate dose 1 to other adults
  t_marker3 <- which(daily_vac_scenarios[[1]]$supply_cum < pop_marker[[1]]*2 +
                       pop_marker[[2]]) %>% 
    .[!.%in%t_marker2] %>% 
    .[.!=(max(t_marker2) + 1)] %>% 
    .[!.%in%t_marker] %>% 
    .[.!=(max(t_marker) + 1)]
  for(i in 1:length(tmp2_tar)){
    daily_vac_scenarios[[1]][t_marker3,tmp2_tar[i]] <-
      daily_vac_scenarios[[1]][t_marker3, "supply_daily"]*tmp2_pop_prop[i]
  }
  
  # should we end the allocation sequence now?
  end_now <- !max(t_marker3) < max(daily_vac_scenarios[[1]]$t)
  if(!end_now){
    t_marker4 <- which(daily_vac_scenarios[[1]]$supply_cum <= pop_marker[[1]]*2 + 
                         pop_marker[[2]]*2) %>% .[. > max(t_marker3)]
    for(i in 1:length(tmp2_tar2)){
      daily_vac_scenarios[[1]][t_marker4,tmp2_tar2[i]] <-
        daily_vac_scenarios[[1]][t_marker4, "supply_daily"]*tmp2_pop_prop[i]
    }
  }

  ##### scenario 2 #####
  # (1) complete vaccinating 60+ with dose 1, reaching the uptake goal 
  # (2) vaccinate other adults with dose 1
  # (3) complete vaccinating 60+ with dose 2
  
  daily_vac_scenarios[[2]] <- daily_vac
  t_marker <- which(daily_vac_scenarios[[2]]$supply_cum < pop_marker[[1]])
  for(i in 1:length(tmp_tar)){
    daily_vac_scenarios[[2]][t_marker,tmp_tar[i]] <-
      daily_vac_scenarios[[2]][t_marker, "supply_daily"]*tmp_pop_prop[i]
  }
  # remainder of dose 1
  daily_vac_scenarios[[2]][max(t_marker)+1,tmp_tar] <- 
    as.list(pop_cap[[1]]  - 
              (daily_vac_scenarios[[2]][,tmp_tar] %>% 
                 colSums()))
  
  # move on to dose 1 for other adults
  daily_vac_scenarios[[2]][max(t_marker)+1,tmp2_tar] <- 
    as.list(tmp2_pop_prop*unlist(daily_vac_scenarios[[2]][max(t_marker)+1,"supply_daily"] - 
                                   sum(daily_vac_scenarios[[2]][max(t_marker)+1,tmp_tar])))
  
  # allocate dose 1 to other adults
  t_marker2 <- which(daily_vac_scenarios[[2]]$supply_cum < pop_marker[[1]] +
                       pop_marker[[2]]) %>% 
    .[!.%in%t_marker] %>% 
    .[.!=(max(t_marker) + 1)]
  for(i in 1:length(tmp2_tar)){
    daily_vac_scenarios[[2]][t_marker2,tmp2_tar[i]] <-
      daily_vac_scenarios[[2]][t_marker2, "supply_daily"]*tmp2_pop_prop[i]
  }
  
  # move on to dose 2 for 60+
  daily_vac_scenarios[[2]][max(t_marker2)+1,tmp_tar2] <- 
    as.list(tmp_pop_prop*unlist(daily_vac_scenarios[[2]][max(t_marker2)+1,"supply_daily"] - 
                                  sum(daily_vac_scenarios[[2]][max(t_marker2)+1,tmp2_tar])))
  
  # allocate dose 2 to 60+
  t_marker3 <- which(daily_vac_scenarios[[2]]$supply_cum < pop_marker[[1]]*2 +
                       pop_marker[[2]]) %>% 
    .[!.%in%t_marker2] %>% 
    .[.!=(max(t_marker2) + 1)] %>% 
    .[!.%in%t_marker] %>% 
    .[.!=(max(t_marker) + 1)]
  for(i in 1:length(tmp_tar2)){
    daily_vac_scenarios[[2]][t_marker3,tmp_tar2[i]] <-
      daily_vac_scenarios[[2]][t_marker3, "supply_daily"]*tmp_pop_prop[i]
  }
  
  
  # should we end the allocation sequence now?
  end_now <- !max(t_marker3) < max(daily_vac_scenarios[[2]]$t)
  if(!end_now){
    t_marker4 <- which(daily_vac_scenarios[[2]]$supply_cum <= pop_marker[[1]]*2 + 
                         pop_marker[[2]]*2) %>% .[. > max(t_marker3)]
    for(i in 1:length(tmp2_tar2)){
      daily_vac_scenarios[[2]][t_marker4,tmp2_tar2[i]] <-
        daily_vac_scenarios[[2]][t_marker4, "supply_daily"]*tmp2_pop_prop[i]
    }
  }
  
  ##### scenario 3 #####
  # (1) start vaccinating the first phase older adults for dose 2
  # (2) catch up on the remaining older adults 60+, for both dose 1 and 2
  #     need to make sure the duration between dose 1 and 2 are more than 4 weeks
  # (3) complete vaccinating other adults for dose 1
  
  # finish vaccinating the first p_change% of 60+
  daily_vac_scenarios[[3]] <- daily_vac
  t_marker <- which(daily_vac_scenarios[[3]]$supply_cum < pop_marker[[1]] * p_change)
  for(i in 1:length(tmp_tar)){
    daily_vac_scenarios[[3]][t_marker,tmp_tar[i]] <-
      daily_vac_scenarios[[3]][t_marker, "supply_daily"]*tmp_pop_prop[i]
  }
  
  # start vaccinating some older adults for dose 2 before reaching the uptake goal
  t_marker2 <- which(daily_vac_scenarios[[3]]$supply_cum < pop_marker[[1]] * p_change * 2) %>% 
    .[!.%in%t_marker] 
  for(i in 1:length(tmp_tar2)){
    daily_vac_scenarios[[3]][t_marker2,tmp_tar2[i]] <-
      daily_vac_scenarios[[3]][t_marker2, "supply_daily"]*tmp_pop_prop[i]
  }
  
  # catch up the adults who has not yet been vaccinated for dose 1
  t_marker3 <- which(daily_vac_scenarios[[3]]$supply_cum < (pop_marker[[1]] * p_change *2 +
                                                              pop_marker[[1]] * (1-p_change))) %>% 
    .[!.%in%t_marker] %>% 
    .[!.%in%t_marker2]
  
  if(length(t_marker3)<28) print("Not enough time between doses!")
  for(i in 1:length(tmp_tar)){
    daily_vac_scenarios[[3]][t_marker3,tmp_tar[i]] <-
      daily_vac_scenarios[[3]][t_marker3, "supply_daily"]*tmp_pop_prop[i]
  }
  
  # allocate the remaining dose 1 to 60+
  daily_vac_scenarios[[3]][max(t_marker3)+1,tmp_tar] <- 
    as.list(pop_cap[[1]] - daily_vac_scenarios[[3]][,tmp_tar] %>% colSums())
  
  daily_vac_scenarios[[3]][max(t_marker3)+1,tmp_tar2] <-
    (daily_vac_scenarios[[3]][max(t_marker3+1), "supply_daily"] - 
       sum(daily_vac_scenarios[[3]][max(t_marker3+1),tmp_tar]))
  
  # catch up the older adults who has not yet been vaccinated for dose 2
  t_marker4 <- which(daily_vac_scenarios[[3]]$supply_cum < (pop_marker[[1]] *2)) %>% 
    .[!.%in%t_marker] %>% 
    .[!.%in%t_marker2] %>% 
    .[!.%in%t_marker3] 
  
  for(i in 1:length(tmp_tar)){
    daily_vac_scenarios[[3]][t_marker4,tmp_tar2[i]] <-
      daily_vac_scenarios[[3]][t_marker4, "supply_daily"]*tmp_pop_prop[i]
  }
  
  # start vaccinating all adults
  t_marker5 <- which(daily_vac_scenarios[[3]]$supply_cum < (pop_marker[[1]]*2 +
                                                              pop_marker[[2]])) %>% 
    .[!.%in%t_marker] %>% 
    .[!.%in%t_marker2] %>% 
    .[!.%in%t_marker3] %>% 
    .[!.%in%t_marker4] 
  
  for(i in 1:length(tmp2_tar)){
    daily_vac_scenarios[[3]][t_marker5,tmp2_tar[i]] <-
      daily_vac_scenarios[[3]][t_marker5, "supply_daily"]*tmp2_pop_prop[i]
  }
  
  # should we end the allocation sequence now?
  end_now <- !max(t_marker5) < max(daily_vac_scenarios[[3]]$t)
  if(!end_now){
    t_marker6 <- which(daily_vac_scenarios[[3]]$supply_cum <= pop_marker[[1]]*2 + 
                         pop_marker[[2]]*2) %>% .[. > max(t_marker5)]
    for(i in 1:length(tmp2_tar2)){
      daily_vac_scenarios[[3]][t_marker6,tmp2_tar2[i]] <-
        daily_vac_scenarios[[3]][t_marker6, "supply_daily"]*tmp2_pop_prop[i]
    }
  }
  
  
  # daily_vac_scenarios[[3]] %>% 
  #   dplyr::select(t, starts_with("Y", ignore.case = T)) %>% 
  #   pivot_longer(starts_with("Y", ignore.case = T)) %>% 
  #   separate(name,into = c("ag", "dose"), sep = "_") %>% 
  #   mutate(t = as.numeric(t),
  #          ag = parse_number(ag)) %>% 
  #   ggplot(., aes(x = t, y = value, color = dose)) +
  #   geom_line() + facet_wrap(~ag)
  

  #### extract vac_para ####
  vac_para <- list()
  # first dose
  daily_vac_scenarios %>% 
    map(select, t, ends_with("d1"), date, supply_daily, supply_cum) %>% 
    bind_rows(.id = "scenario") %>% 
    group_by(supply_daily, scenario) %>% group_split() %>% 
    map(group_by_at, vars(starts_with("Y"))) %>% 
    map(filter, t == min(t)) %>% 
    bind_rows() %>% group_by(scenario) %>% group_split() %>% 
    map(mutate, t = parse_number(t)) %>% 
    map(arrange, t) -> vac_para[[1]]
  
  # second dose
  daily_vac_scenarios %>% 
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
    for(j in 1:3){
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
    for(j in 1:3){
      vac_para[[i]][[j]] %>% 
        pull(t) %>% 
        as.numeric -> vacc_times[[i]][[j]] 
      vacc_times[[i]][[j]][1] <- 0
    }
  }
  
  res <- list()
  for(i in 1:3){
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
              res = res, # three strategies
              p_supply = p_supply,
              daily_vac_scenarios = daily_vac_scenarios))
  }


##### expand VE estimates to meet the needs of the model ####
exp_ve <- function(ve_d_o,  # disease blocking VE observed
                   ve_i_o   # infection blocking VE assumed
){
  # calculate of clinical fraction reduction
  ve_d <- (ve_d_o - ve_i_o)/(1 - ve_i_o)
  return(ve_d)
}





