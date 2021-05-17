gen_country_basics <- function(country,
                               waning_nat = 45*7,
                               R0_assumed  = 2.7,
                               date_start = "2020-01-01",
                               date_end = "2022-12-31",
                               # s_A = 0.1,
                               deterministic = TRUE){
  
  require(countrycode)
  
  wb_tmp = countrycode(country, "country.name", "wb")
  # c_tmp = schedule_raw %>% 
  #   filter(wb == wb_tmp,
  #          date >= lubridate::ymd(date_start),
  #          date <= lubridate::ymd(date_end))
  # 
  # t_run <- nrow(c_tmp)
  
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
  
  para$processes = burden_processes
  
  # para$schedule[["mobility"]] = list(
  #   parameter = "contact",
  #   pops = numeric(),
  #   mode = "assign",
  #   values = split(c_tmp[,4:7],
  #                  seq(nrow(c_tmp))) %>%
  #     map(unlist) %>%
  #     map(as.vector) %>%
  #     unname,
  #   times = 1:t_run)
  
  return(para)
}

update_vac_char <- function(para,
                            waning_vac = 105,
                            waning_vac2 = 365*3,
                            ve_i = NULL,
                            ve_i2 = NULL,
                            ve_d = NULL,
                            ve_d2 = NULL){
  
  n_age <- length(para$pop[[1]]$size)
  
  # currently closing two-dose setup
  # para$pop[[1]]$ev2 <- rep(0, n_age)
  
  # 
  # para$pop[[1]]$pd_ri <- rep(1, n_age)
  # para$pop[[1]]$pi_r <- rep(1, n_age)
  
  # key parameters on infection and disease blocking mechanisms
  para$pop[[1]]$ei_v <- rep(ve_i, n_age)
  para$pop[[1]]$ed_vi <- rep(ve_d, n_age)
  para$pop[[1]]$ei_v2 <- rep(ve_i2, n_age)
  para$pop[[1]]$ed_vi2 <- rep(ve_d2, n_age)
  
  # waning vaccine-induced immunity
  para$pop[[1]]$wv = rep((1/waning_vac), n_age)
  para$pop[[1]]$wv2 = rep((1/waning_vac2), n_age)
  
  
  # return results
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
                                         0.4),
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
                                para$time1),
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
  
  # scenario 1
  # (1) complete vaccinating 60+ with dose 1, reaching the uptake goal 
  # (2) complete vaccinating 60+ with dose 2
  # (3) vaccinate other adults with dose 1
 source("code/prioritisation/scenario1.R")
  
  # scenario 2
  # (1) complete vaccinating 60+ with dose 1, reaching the uptake goal 
  # (2) vaccinate other adults with dose 1
  # (3) complete vaccinating 60+ with dose 2
  source("code/prioritisation/scenario2.R")
  
  daily_vac_scenarios[[2]] %>%
    pivot_longer(cols = starts_with("Y", ignore.case = F)) %>%
    separate(name, into = c("ag", "dose"), sep = "_") %>%
    mutate(ag = parse_number(ag)) %>%
    ggplot(., aes(x = date, y = value)) +
    # geom_bar(stat = "identity")+
    geom_line()+
    facet_grid(ag~dose)
   
  # tmp_tar <- c(paste0("Y", tmp_priorities[[1]]$age_group,"_d1"))
  # for(i in 1:length(tmp_tar)){
  #   daily_vac[daily_vac$supply_cum < pop_marker * p_change,tmp_tar[i]] <-
  #     daily_vac[daily_vac$supply_cum < pop_marker * p_change, "supply_daily"]*
  #     tmp_pop$n_pop[tmp_priorities[[1]]$age_group[i]]/sum(tmp_pop$n_pop[tmp_priorities[[1]]$age_group])
  # }
  # 
  # 
  # has uptake target been met yet?
  # uptake_tar <- p_change >= 1 # TRUE: met, FALSE: not met
  # if(!uptake_tart){
  #   
  # }
  # 
  # daily_vac_scenarios <- list()
  # # delay scenarios among top priority groups
  # tmp_tar <- c(paste0("Y", tmp_priorities[[1]]$age_group,"_d1"),
  #              paste0("Y", tmp_priorities[[1]]$age_group,"_d2"))
  # ag_tmp <- tmp_tar[i] %>% strsplit(split = "_") %>% unlist %>% .[1]
  # # i - total number of age groups in this priority level
  # # j - number of delays we'd like to explore
  # 
  # # 
  # daily_vac_scenarios[[1]] <- list()
  # for(j in 1:length(dose_interval)){
  #   daily_vac_scenarios[[1]][[j]] <- daily_vac
  #   for(i in 1:((length(tmp_tar))/2)){
  #     daily_vac_scenarios[[1]][[j]][,tmp_tar[i+4]] <- 
  #       lag(daily_vac[,tmp_tar[i]], dose_interval[j]*7) %>% replace(.,is.na(.), 0)
  #     daily_vac_scenarios[[1]][[j]] %<>% 
  #       mutate(supply_tmp = (!!as.name(tmp_tar[i])),
  #              !!tmp_tar[i]  := supply_tmp - !!sym(tmp_tar[i+4])) %>% 
  #       dplyr::select(-supply_tmp)
  #   }
  # }
  # daily_vac_scenarios[[1]] %<>% bind_rows(.id = "dose_interval") %>% 
  #   left_join(dose_interval %>% enframe(name = "dose_interval") %>% 
  #               mutate(dose_interval = as.character(dose_interval)),
  #             by = "dose_interval") %>% 
  #   dplyr::select(-dose_interval) %>% 
  #   rename(dose_interval = value)
# 
#   t_change <- min(which(daily_vac_scenarios[[1]]$supply_cum > pop_marker*p_change))
#     
# 
#   
#   pivot_longer(starts_with("Y", ignore.case = F)) %>%
#     separate(name, into = c("ag", "dose"), sep = "_") %>% 
#     filter(ag %in% paste0("Y",13:16)) %>% 
#     ggplot(., aes(x = date, y = value, group = ag, color = ag)) +
#     geom_line()+
#     facet_grid(dose~dose_interval)
  
  

  # daily_vac_scenarios[[1]] %>% 
  #   ggplot(., aes(x = date, y = supply)) +
  #   geom_line()
  
  # t_current min(which(daily_vac$supply_cum >= pop_marker*p_change))
  # if(is.infinite(tail(date_marker, 1))) {exhaust = F} else {exhaust = T}
  # if(all(is.infinite(date_marker)) & length(date_marker) == 2) {np = T} else {np = F}
  
  date_marker <- c(as.numeric(tmp_schedule$t[2]),
                   date_marker,
                   nrow(daily_vac)) %>% unique %>% sort %>% .[!is.infinite(.)]
  
  for(i in 1:(length(date_marker)-1)) {
    if(i == 1) {
      daily_vac$phase[date_marker[i]] <- i
    }
    daily_vac$phase[(date_marker[i]+1):(date_marker[i+1])]<- i
  }
  
  if(np) {end <- 1} else {end <- length(tmp_priorities)}
  
  for(i in 1:end){
    tmp_ag <- paste0("Y",tmp_priorities[[i]]$age_group)
    group_tot <- para$pop[[1]]$size[tmp_priorities[[i]]$age_group]
    group_prop_ag <- group_tot/sum(group_tot)
    
    if(length(tmp_priorities) > i){
      tmp_ag_next <-  paste0("Y",tmp_priorities[[i+1]]$age_group)
      group_tot_next <- para$pop[[1]]$size[tmp_priorities[[i+1]]$age_group]
      group_prop_next <- group_tot_next/sum(group_tot_next)
    }
    
    if(i == 1){
      daily_vac[((date_marker[i]):(date_marker[i+1]-1)),tmp_ag] <- 
        sapply(1:length(group_prop_ag), 
               function(x) group_prop_ag[x]*
                 daily_vac[((date_marker[i]):(date_marker[i+1]-1)),
                           "supply"]) %>% 
        bind_cols() %>% 
        setNames(tmp_ag)
    } 
    
    if(i > 1 & date_marker[i] < nrow(daily_vac)){
      daily_vac[((date_marker[i]+1):(date_marker[i+1]-1)),tmp_ag] <- 
        sapply(1:length(group_prop_ag), 
               function(x) group_prop_ag[x]*
                 daily_vac[((date_marker[i]+1):(date_marker[i+1]-1)),
                           "supply"]) %>% 
        bind_cols() %>% 
        setNames(tmp_ag)
    }
    
    if(date_marker[i] < nrow(daily_vac)){
      remainder <- sum(pop_cap[[i]]) - daily_vac[,tmp_ag] %>% 
        colSums(., na.rm = T) %>% sum
      supply <- daily_vac[date_marker[i+1],"supply"]
      daily_vac[date_marker[i+1], tmp_ag] <- as.list(remainder*group_prop_ag)
      
      if(length(tmp_priorities) > i){
        daily_vac[date_marker[i+1], tmp_ag_next] <- 
          as.list(group_prop_next*unlist((supply - remainder)))
        
      }
    }
    
  }
  
  # Highlight the exception
  # If the vaccine supply has yet to saturated the population, we will replace 
  # the last row with the row before, because the last row cannot mop up what's 
  # happening before hand
  
  if(!exhaust){
    daily_vac[nrow(daily_vac), paste0("Y",1:16)] <- 
      daily_vac[nrow(daily_vac)-1, paste0("Y",1:16)] 
  }
  
  # putting all vaccine policy related raw parameters back together
  daily_vac %>% 
    mutate(t = as.numeric(t)) %>% 
    group_by(phase) %>% group_split() %>% 
    map(group_by_at, vars(starts_with("Y"))) %>% 
    map(filter, t == min(t)) %>% 
    bind_rows() -> vac_para
  
  # then convert these parameters to a format that's friendly with `covidm`
  # allocation
  vac_para %>% 
    dplyr::select(starts_with("Y")) %>% 
    split(seq(nrow(vac_para))) %>% 
    map(unlist) %>% 
    map(as.vector) -> vacc_vals
  
  # timing
  vacc_times <- vac_para$t %>% as.numeric; vacc_times[1] <- 0
  
  para$schedule[["vaccination"]] = list(
    parameter = "v",
    pops = numeric(),
    mode = "assign",
    values = vacc_vals,
    times = vacc_times)
  
  daily_vac %<>% 
    rename(doses_daily = supply,
           supply = supply_cum) %>% 
    mutate(date = as.Date(date_start) + as.numeric(t)) %>% 
    dplyr::select(-phase)
  
  return(list(param = para, 
              supply = tmp_schedule,
              vac_para = vac_para,
              daily_vac = daily_vac))
}
