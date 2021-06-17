gen_country_basics <- function(country,
                               waning_nat = 52*7*3,
                               R0_assumed  = 2.7,
                               date_start = "2020-01-01",
                               date_end = "2022-12-31",
                               deterministic = TRUE){
  
  require(countrycode)
  
  wb_tmp = countrycode(country, "country.name", "wb")
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
                            ve_i  = NULL, # infection blocking VE post 1 dose
                            v2e_i = NULL, # infection blocking VE post 2 doses
                            ve_d  = NULL, # perc reduction in clinical fraction post 1 dose
                            v2e_d = NULL, # perc reduction in clinical fraction post 2 doses 
                            wv = NULL # 1/waning period
                            ){
  
  n_age <- length(para$pop[[1]]$size)

  # key parameters on infection and disease blocking mechanisms
  para$pop[[1]]$uv  <- rep(para$pop[[1]]$u * (1 - ve_i),  n_age)
  para$pop[[1]]$uv2 <- rep(para$pop[[1]]$u * (1 - v2e_i), n_age)
  para$pop[[1]]$yv  <- rep(para$pop[[1]]$y * (1 - ve_d),  n_age)
  para$pop[[1]]$yv2 <- rep(para$pop[[1]]$y * (1 - v2e_d), n_age)
  
  # waning vaccine-induced immunity
  para$pop[[1]]$wv = rep(wv, n_age)

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
  

  daily_vac_scenarios[[1]] %>% 
    dplyr::select(t, starts_with("Y", ignore.case = T)) %>% 
    pivot_longer(starts_with("Y", ignore.case = T)) %>% 
    separate(name,into = c("ag", "dose"), sep = "_") %>% 
    mutate(t = as.numeric(t),
           ag = parse_number(ag)) %>% 
    ggplot(., aes(x = t, y = value, color = dose)) +
    geom_line() + facet_wrap(~ag)
  
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
  
  # create empty rows to record things
  matrix(0, 
         ncol = length(para$pop[[1]]$size),
         nrow = max(as.numeric(tmp_schedule$t), na.rm = T)) %>% 
    as_tibble() %>% 
    setNames(paste0("Y",1:16, "_d1_S")) %>%
    # record the one that effectively still remains in V
    bind_cols(matrix(0, 
                     ncol = length(para$pop[[1]]$size),
                     nrow = max(as.numeric(tmp_schedule$t), na.rm = T)) %>% 
                as_tibble() %>% 
                setNames(paste0("Y",1:16, "_d1_V"))) %>%
    bind_cols(matrix(0, 
                     ncol = length(para$pop[[1]]$size),
                     nrow = max(as.numeric(tmp_schedule$t), na.rm = T)) %>% 
                as_tibble() %>% 
                setNames(paste0("Y",1:16, "_d1_SV2"))) %>%
    bind_cols(matrix(0, 
                     ncol = length(para$pop[[1]]$size),
                     nrow = max(as.numeric(tmp_schedule$t), na.rm = T)) %>% 
                as_tibble() %>% 
                setNames(paste0("Y",1:16, "_d1_VV2"))) %>%
    bind_cols(., daily_vac_scenarios[[2]]) -> daily_vac_scenarios[[2]]
  
  # allocate vaccines
  for(i in 5:length(para$pop[[1]]$group_names)){
    tag1 <- paste0("Y",i,"_d1")
    tag2 <- paste0("Y",i,"_d2")
    tag3 <- paste0("Y",i,"_d1_S")
    tag4 <- paste0("Y",i,"_d1_V")
    tag5 <- paste0("Y",i,"_d1_SV2")
    tag6 <- paste0("Y",i,"_d1_VV2")
    p_wane <- para$pop[[1]]$wv[1]#1-exp(-para$pop[[1]]$wv[1])
    
    # phase 1, only vaccinating 60+ 
    for(j in 2:max(t_marker2)){
      # new V = newly vaccinated + those who didn't wane
      daily_vac_scenarios[[2]][j, tag4] <- 
        daily_vac_scenarios[[2]][j, tag1] + daily_vac_scenarios[[2]][j-1, tag4]*(1-p_wane)
      # new S waned back from V1
      daily_vac_scenarios[[2]][j, tag3] <- 
        daily_vac_scenarios[[2]][j-1, tag3] +
        daily_vac_scenarios[[2]][j-1, tag4]*p_wane
    }
    
    if(i >= 13){
      k = max(t_marker2) + 1
      repeat{
        daily_vac_scenarios[[2]][k, tag4] <-
          (daily_vac_scenarios[[2]][k-1, tag4])*(1-p_wane)
        
        daily_vac_scenarios[[2]][k, tag3] <-
          daily_vac_scenarios[[2]][k-1, tag3] +
          (daily_vac_scenarios[[2]][k, tag4])*(p_wane) -
          daily_vac_scenarios[[2]][k, tag2]
        
        daily_vac_scenarios[[2]][k, tag5] <-
          daily_vac_scenarios[[2]][k, tag2]
        
        if(daily_vac_scenarios[[2]][k, tag3] < 0 |
           is.na(daily_vac_scenarios[[2]][k, tag3])|
           k >= nrow(daily_vac_scenarios[[2]])){
          
          daily_vac_scenarios[[2]][k, tag3] <-
            daily_vac_scenarios[[2]][k-1, tag3] +
            (daily_vac_scenarios[[2]][k-1, tag4])*(p_wane)
          
          daily_vac_scenarios[[2]][k, tag4] <-
            (daily_vac_scenarios[[2]][k-1, tag4])*(1-p_wane)
          
          
          daily_vac_scenarios[[2]][k, tag5] <-
            daily_vac_scenarios[[2]][k, tag2]
          
          daily_vac_scenarios[[2]][k, tag6] <- 0
          
          
          break
        }
        k = k+1
      }
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
  
  matrix(0, 
         ncol = length(para$pop[[1]]$size),
         nrow = max(as.numeric(tmp_schedule$t), na.rm = T)) %>% 
    as_tibble() %>% 
    setNames(paste0("Y",1:16, "_d1_S")) %>%
    # record the one that effectively still remains in V
    bind_cols(matrix(0, 
                     ncol = length(para$pop[[1]]$size),
                     nrow = max(as.numeric(tmp_schedule$t), na.rm = T)) %>% 
                as_tibble() %>% 
                setNames(paste0("Y",1:16, "_d1_V"))) %>%
    bind_cols(matrix(0, 
                     ncol = length(para$pop[[1]]$size),
                     nrow = max(as.numeric(tmp_schedule$t), na.rm = T)) %>% 
                as_tibble() %>% 
                setNames(paste0("Y",1:16, "_d1_SV2"))) %>%
    bind_cols(matrix(0, 
                     ncol = length(para$pop[[1]]$size),
                     nrow = max(as.numeric(tmp_schedule$t), na.rm = T)) %>% 
                as_tibble() %>% 
                setNames(paste0("Y",1:16, "_d1_VV2"))) %>%
    bind_cols(., daily_vac_scenarios[[3]]) -> daily_vac_scenarios[[3]]
  
  for(i in 5:length(para$pop[[1]]$group_names)){
    tag1 <- paste0("Y",i,"_d1")
    tag2 <- paste0("Y",i,"_d2")
    tag3 <- paste0("Y",i,"_d1_S")
    tag4 <- paste0("Y",i,"_d1_V")
    tag5 <- paste0("Y",i,"_d1_SV2")
    tag6 <- paste0("Y",i,"_d1_VV2")
    p_wane <- para$pop[[1]]$wv[1]#1-exp(-para$pop[[1]]$wv[1])
    
    # phase 1, only vaccinating 60+ 
    for(j in 2:max(t_marker)){
      # new V = newly vaccinated + those who didn't wane
      daily_vac_scenarios[[3]][j, tag4] <- 
        daily_vac_scenarios[[3]][j, tag1] + daily_vac_scenarios[[3]][j-1, tag4]*(1-p_wane)
      # new S waned back from V1
      daily_vac_scenarios[[3]][j, tag3] <- 
        daily_vac_scenarios[[3]][j-1, tag3]+
        daily_vac_scenarios[[3]][j-1, tag4]*p_wane
    }
    
    # convert those waned from V1 to S from daily to cumulative measure
    # daily_vac_scenarios[[3]][, tag3]  <- 
    #   cumsum(daily_vac_scenarios[[3]][, tag3])
    
    if(i >= 13){
      # phase 2 start vaccinating the vaccinated individuals who has lost their 
      # immunity
      k = max(t_marker) + 1
      repeat{
        daily_vac_scenarios[[3]][k, tag4] <-
          (daily_vac_scenarios[[3]][k-1, tag4])*(1-p_wane)
        
        daily_vac_scenarios[[3]][k, tag3] <-
          daily_vac_scenarios[[3]][k-1, tag3] +
          (daily_vac_scenarios[[3]][k, tag4])*(p_wane) -
          daily_vac_scenarios[[3]][k, tag2]
        
        daily_vac_scenarios[[3]][k, tag5] <-
          daily_vac_scenarios[[3]][k, tag2]
        
        if(daily_vac_scenarios[[3]][k, tag3] < 0 |
           is.na(daily_vac_scenarios[[3]][k, tag3])){
          daily_vac_scenarios[[3]][k, tag3] <- 0
          
          daily_vac_scenarios[[3]][k, tag5] <-
            daily_vac_scenarios[[3]][k-1, tag3] +
            (daily_vac_scenarios[[3]][k, tag4])*(p_wane)
          
          daily_vac_scenarios[[3]][k, tag6] <-
            daily_vac_scenarios[[3]][k, tag2] -
            (daily_vac_scenarios[[3]][k, tag5])
          
          daily_vac_scenarios[[3]][k, tag4] <-
            (daily_vac_scenarios[[3]][k-1, tag4])*(1-p_wane) -
            daily_vac_scenarios[[3]][k, tag6]
          k = k+1
          break
        }
        k = k+1
      }
      
      daily_vac_scenarios[[3]][k:nrow(daily_vac_scenarios[[3]]), tag3] <- 0
      
      repeat{
        # flux changes to V, waning out and then vaccinate out
        daily_vac_scenarios[[3]][k, tag4] <-
          (daily_vac_scenarios[[3]][k-1, tag4])*(1-p_wane) -
          (daily_vac_scenarios[[3]][k, tag2] -
             daily_vac_scenarios[[3]][k-1, tag4]*(p_wane))
        
        daily_vac_scenarios[[3]][k, tag3] <-
          (daily_vac_scenarios[[3]][k-1, tag4])*(p_wane) -
          (daily_vac_scenarios[[3]][k-1, tag4])*(p_wane)
        
        # flux change to SV2 path
        daily_vac_scenarios[[3]][k, tag5] <-
          (daily_vac_scenarios[[3]][k-1, tag4])*(p_wane)
        
        # flux change to VV2 path
        daily_vac_scenarios[[3]][k, tag6] <-
          (daily_vac_scenarios[[3]][k, tag2]) -
          (daily_vac_scenarios[[3]][k-1, tag4])*(p_wane)
        
        if(daily_vac_scenarios[[3]][k, tag4] < 0 |
           is.na(daily_vac_scenarios[[3]][k, tag4])){
          k = k + 1
          daily_vac_scenarios[[3]][k, tag6] <-
            (daily_vac_scenarios[[3]][k, tag2])
          break
        }
        k = k + 1
      }
      
      for(j in t_marker3){
        # new V = newly vaccinated + those who didn't wane
        daily_vac_scenarios[[3]][j, tag4] <-
          daily_vac_scenarios[[3]][j, tag1] + daily_vac_scenarios[[3]][j-1, tag4]*(1-p_wane)
        # new S waned back from V1
        daily_vac_scenarios[[3]][j, tag3] <-
          daily_vac_scenarios[[3]][j-1, tag3] +
          daily_vac_scenarios[[3]][j-1, tag4]*p_wane
      }
      
      k = max(t_marker3) + 1
      repeat{
        daily_vac_scenarios[[3]][k, tag4] <-
          (daily_vac_scenarios[[3]][k-1, tag4])*(1-p_wane)
        
        daily_vac_scenarios[[3]][k, tag3] <-
          daily_vac_scenarios[[3]][k-1, tag3] +
          (daily_vac_scenarios[[3]][k, tag4])*(p_wane) -
          daily_vac_scenarios[[3]][k, tag2]
        
        daily_vac_scenarios[[3]][k, tag5] <-
          daily_vac_scenarios[[3]][k, tag2]
        
        if(daily_vac_scenarios[[3]][k, tag3] < 0 |
           is.na(daily_vac_scenarios[[3]][k, tag3])){
          daily_vac_scenarios[[3]][k, tag3] <- 0
          
          daily_vac_scenarios[[3]][k, tag5] <-
            daily_vac_scenarios[[3]][k-1, tag3] +
            (daily_vac_scenarios[[3]][k, tag4])*(p_wane)
          
          daily_vac_scenarios[[3]][k, tag6] <-
            daily_vac_scenarios[[3]][k, tag2] -
            (daily_vac_scenarios[[3]][k, tag5])
          
          daily_vac_scenarios[[3]][k, tag4] <-
            (daily_vac_scenarios[[3]][k-1, tag4])*(1-p_wane) -
            daily_vac_scenarios[[3]][k, tag6]
          k = k+1
          break
        }
        k = k+1
      }
      
      repeat{
        # flux changes to V, waning out and then vaccinate out
        daily_vac_scenarios[[3]][k, tag4] <-
          (daily_vac_scenarios[[3]][k-1, tag4])*(1-p_wane) -
          (daily_vac_scenarios[[3]][k, tag2] - 
             daily_vac_scenarios[[3]][k-1, tag4]*(p_wane))
        
        daily_vac_scenarios[[3]][k, tag3] <-
          (daily_vac_scenarios[[3]][k-1, tag4])*(p_wane) -
          (daily_vac_scenarios[[3]][k-1, tag4])*(p_wane)
        
        # flux change to SV2 path
        daily_vac_scenarios[[3]][k, tag5] <-
          (daily_vac_scenarios[[3]][k-1, tag4])*(p_wane)
        
        # flux change to VV2 path 
        daily_vac_scenarios[[3]][k, tag6] <-
          (daily_vac_scenarios[[3]][k, tag2]) - 
          (daily_vac_scenarios[[3]][k-1, tag4])*(p_wane)
        
        if(daily_vac_scenarios[[3]][k, tag4] < 0 |
           is.na(daily_vac_scenarios[[3]][k, tag4])){
          k = k + 1
          daily_vac_scenarios[[3]][k, tag6] <-
            (daily_vac_scenarios[[3]][k, tag2])
          break
        }
        k = k + 1
      }
    }
  }

  #### extract vac_para ####
  vac_para <- list()
  # S --> V1
  daily_vac_scenarios %>% 
    map(select, t, ends_with("d1"), date, supply_daily, supply_cum) %>% 
    bind_rows(.id = "scenario") %>% 
    group_by(supply_daily, scenario) %>% group_split() %>% 
    map(group_by_at, vars(starts_with("Y"))) %>% 
    map(filter, t == min(t)) %>% 
    bind_rows() %>% group_by(scenario) %>% group_split() -> vac_para[[1]]
  
  # S --> V2
  daily_vac_scenarios %>% 
    map(select, t, ends_with("SV2"), date, supply_daily, supply_cum) %>% 
    bind_rows(.id = "scenario") %>% 
    group_by(supply_daily, scenario) %>% group_split() %>% 
    map(group_by_at, vars(starts_with("Y"))) %>% 
    map(filter, t == min(t)) %>% 
    bind_rows() %>% group_by(scenario) %>% group_split() -> vac_para[[2]]
  
  # V --> V2
  daily_vac_scenarios %>% 
    map(select, t, ends_with("VV2"), date, supply_daily, supply_cum) %>% 
    bind_rows(.id = "scenario") %>% 
    group_by(supply_daily, scenario) %>% group_split() %>% 
    map(group_by_at, vars(starts_with("Y"))) %>% 
    map(filter, t == min(t)) %>% 
    bind_rows() %>% group_by(scenario) %>% group_split() -> vac_para[[3]]

  vacc_vals <- list()
  # then convert these parameters to a format that's friendly with `covidm`
  # allocation
  for(i in 1:3){
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
  for(i in 1:3){
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
    res[[i]]$schedule[["SV1"]] =  
      list(
      parameter = "v",
      pops = numeric(),
      mode = "assign",
      values = vacc_vals[[1]][[i]],
      times = vacc_times[[1]][[i]]
      )
    
    res[[i]]$schedule[["SV2"]] =  
      list(
        parameter = "v2",
        pops = numeric(),
        mode = "assign",
        values = vacc_vals[[2]][[i]],
        times = vacc_times[[2]][[i]]
      )
    
    res[[i]]$schedule[["VV2"]] =  
      list(
        parameter = "v12",
        pops = numeric(),
        mode = "assign",
        values = vacc_vals[[3]][[i]],
        times = vacc_times[[3]][[i]]
      )
  }
  
  return(list(param = para, 
              res = res,
              p_supply = p_supply,
              daily_vac_scenarios = daily_vac_scenarios))
}

change_ve <- function(
  para = NULL,
  transition_start = "2021-03-15",
  transition_end = "2021-06-15",
  ve_reduction = NULL
){
  require(lubridate)
  require(magrittr)
  require(tidyverse)
  
  tmp_start <- para$res[[1]]$date0
  tmp_end <- para$res[[1]]$time1
  
  
  data.table(date = seq(ymd(tmp_start), ymd(tmp_end), 1)) %>% 
    mutate(value = case_when(date >= tmp_start & date <= ymd(transition_start) ~ 1,
                             date >= ymd(transition_end) ~ 1-ve_reduction)) -> tmp
  
  tmp %<>% 
    left_join(seq(ymd(transition_start), 
                  ymd(transition_end), 
                  1) %>% enframe %>% 
                mutate(max = name - round(max(name)/2),
                       sig = sigmoid::sigmoid(max, SoftMax = T),
                       sig = round(sig, 2)) %>% 
                rename(date = value),
              by = "date"
    ) %>%   
    rowid_to_column(var = "t") %>% 
    mutate(value = if_else(is.na(value), 1-sig*0.2, value),
           t = as.numeric(t) - 1) %>% 
    dplyr::select(-name, -max, -sig) %>% 
    group_by(value) %>% group_split() %>% 
    map(~.[1,]) %>% bind_rows() %>% arrange(date)
  
  n_para <- length(para$res)
  for(i in 1:n_para){
    para$res[[i]]$schedule[["ei_v_change"]] <- list(
      parameter = "ei_v",
      pops = numeric(),
      mode = "multiply",
      values = tmp$value %>% map(rep, 16),
      times = tmp$t
    )
    
    para$res[[i]]$schedule[["ei_v2_change"]] <- list(
      parameter = "ei_v2",
      pops = numeric(),
      mode = "multiply",
      values = tmp$value %>% map(rep, 16),
      times = tmp$t
    )
  }
  return(para)
}




