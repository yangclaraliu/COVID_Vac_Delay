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
  
  # scenario 3
  # (1) start vaccinating the first phase older adults for dose 2
  # (2) catch up on the remaining older adults 60+, for both dose 1 and 2
  #     need to make sure the duration between dose 1 and 2 are more than 4 weeks
  # (3) complete vaccinating other adults for dose 1
  source("code/prioritisation/scenario3.R")
  
  # daily_vac_scenarios[[3]] %>%
  #   pivot_longer(cols = starts_with("Y", ignore.case = F)) %>%
  #   separate(name, into = c("ag", "dose"), sep = "_") %>%
  #   mutate(ag = parse_number(ag)) %>%
  #   filter(ag > 4) %>% 
  #   ggplot(., aes(x = date, y = value)) +
  #   # geom_bar(stat = "identity")+
  #   geom_line()+
  #   # geom_point() +
  #   facet_grid(ag~dose)
  
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
