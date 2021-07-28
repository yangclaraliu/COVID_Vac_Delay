params_2$daily_vac_scenarios %>%
  bind_rows(.id = "ROS") %>% 
  dplyr::select(-c(supply_cum, supply_daily, supply, supply_2, phase)) %>%   
  pivot_longer(cols = starts_with("Y", ignore.case = T)) %>% 
  separate(name, into = c("ag", "dose")) %>% 
  filter(value != 0) %>% 
  mutate(ag = factor(ag, 
                     levels = paste0("Y", 1:16),
                     labels = params$pop[[1]]$group_names),
         ROS = factor(ROS,
                      levels = 1:3,
                      labels = paste0("Roll out strategy: ",1:3)),
         dose = factor(dose,
                       levels = c("d1", "d2"),
                       labels = c("1st dose",
                                  "2nd dose"))
         ) %>% 
  ggplot(., aes(x = date, y = ag, color = dose)) +
  geom_point() +
  labs(x = "", y = "Age Group", color = "") +
  facet_wrap(~ROS, ncol = 1) +
  ggsci::scale_color_lancet() +
  theme_bw() +
  theme(strip.text = element_text(size = 12)) -> p

ggsave("figs/ROS.png", p, height = 9, width = 16)

tmp %>% 
  pivot_longer(cols = starts_with("uv"),
               values_to = "uv") %>% 
  separate(name, into = c("var", "ag")) %>% 
  mutate(ag = factor(ag,
                     levels = paste0("Y",1:16),
                     labels = params$pop[[1]]$group_names
                     )) %>% 
    left_join(data.table(u = params_2$res[[1]]$pop[[1]]$u,
                       ag = params$pop[[1]]$group_names)) %>%
  mutate(ag = factor(ag,
                     levels = params$pop[[1]]$group_names)) %>% 
  ggplot(., aes(x = date, y = uv, group = ag)) +
  geom_line() +
  geom_hline(aes(yintercept = u)) + facet_wrap(~ag, ncol = 4) +
  theme_bw() +
  labs(x = "", y = "Susceptibility") +
  theme(strip.text = element_text(size = 12))-> p

ggsave("figs/sus_change.png", p, height = 9, width = 16)

# Albania
DEoptim2_selected_debug <- readRDS("~/GitHub/covid19_vac_prioritisation_europe/data/intermediate/DEoptim2_selected_debug.rds")
DEoptim3_selected_debug <- readRDS("~/GitHub/covid19_vac_prioritisation_europe/data/intermediate/DEoptim3_selected_debug.rds")



predict_outbreak <- function(delay = 20,
                             under_report = F,
                             wane_vac = 120){
  
  if(under_report == F){
    i <- which(DEoptim2_selected_debug$country_name == "Georgia")
    tmp <- gen_country_basics(country = DEoptim3_selected_debug$country_name[i],
                              waning_nat = 52*7*3,
                              R0_assumed  = DEoptim2_selected_debug$r[i],
                              date_start = DEoptim2_selected_debug$date[i],
                              date_end = "2022-12-31",
                              processes = burden_processes,
                              deterministic = T) 
  }
  
  if(under_report == T){
    i <- which(DEoptim3_selected_debug$country_name == "Georgia")
    tmp <- gen_country_basics(country = DEoptim3_selected_debug$country_name[i],
                              waning_nat = 52*7*3,
                              R0_assumed  = DEoptim3_selected_debug$r[i],
                              date_start = DEoptim3_selected_debug$date[i],
                              date_end = "2022-12-31",
                              processes = burden_processes,
                              deterministic = T) 
    
      tmp$processes[[1]]$prob[1,] <- 
        tmp$processes[[1]]$prob[1,] * DEoptim3_selected_debug$par3[i]
      tmp$processes[[2]]$prob[1,] <- 
        tmp$processes[[2]]$prob[1,] * DEoptim3_selected_debug$par3[i]
      tmp$processes[[3]]$prob[1,] <- 
        tmp$processes[[3]]$prob[1,] * DEoptim3_selected_debug$par3[i]
  }

  require(magrittr)
  tmp %<>% 
    update_vac_char(.,
                    ve_i   = ve$ve_i_o[1],  # infection blocking VE post 1 dose
                    v2e_i  = ve$ve_i_o[2],  # infection blocking VE post 2 doses
                    ve_d   = ve$ve_d[1],    # clinical fraction among breakthrough post 1 dose
                    v2e_d  = ve$ve_d[2],    # clinical fraction among breakthrough post 2 doses
                    wv = 1/wane_vac) %>% 
    vac_policy(.,
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
               cov_max = c(rep(0,4),
                           rep(0.7, 8),
                           rep(0.9, 4)),
               # the proportion of age groups at which we will compare
               # rolling out dose 2s or more dose 1s
               p_change = 0.7,
               supply_delay = delay # unit = weeks
    )
  return(tmp)
}


16:40 %>% 
  map(~predict_outbreak(delay = .)) -> all_params

lapply(1:25, function(x){
  all_params[[x]]$res %>% 
    map(cm_simulate) %>% 
    lapply(., "[[", "dynamics") %>% 
    bind_rows(.id = "ROS") %>% 
    filter(compartment == "death_o") %>% 
    group_by(ROS) %>% summarise(value = sum(value))
}) -> vary_delay

vary_delay %>% 
  map(mutate, rk = rank(value)) %>% 
  bind_rows(., .id = "delay") %>% 
  mutate(delay = as.numeric(delay) + 15,
         rk = factor(rk)) %>% 
  ggplot(., aes(x = delay, 
                y = value, 
                group = ROS,
                color = ROS)) +
  geom_line(size = 2) +
  theme_bw() +
  ggsci::scale_color_lancet() +
  labs(x = "Supply Delay", y = "COVID-19 Mortality") -> p

ggsave("figs/vary_delay_compare.png", height = 9, width = 8)


seq(90, 52*7*3,30) %>% 
  map(~predict_outbreak(delay = 20,
                        wane_vac = .)) -> all_params


lapply(1:length(all_params), function(x){
  all_params[[x]]$res %>% 
    map(cm_simulate) %>% 
    lapply(., "[[", "dynamics") %>% 
    bind_rows(.id = "ROS") %>% 
    filter(compartment == "death_o") %>% 
    group_by(ROS) %>% summarise(value = sum(value))
}) -> vary_wv

vary_wv %>% 
  map(mutate, rk = rank(value)) %>% 
  bind_rows(., .id = "wv_index") %>% 
  left_join(data.frame(wv_value = seq(90, 52*7*3,30)) %>% 
              rownames_to_column(var = "wv_index")) %>% 
  ggplot(., aes(x = wv_value, y = value, group = ROS, color = ROS)) +
  geom_line(size = 2) +
  theme_bw() +
  ggsci::scale_color_lancet() +
  labs(x = "1st Dose Waning", y = "COVID-19 Mortality") -> p

ggsave("figs/vary_wv.png", p, height = 9, width = 16)

# incorporate change ve
seq(0.1, 0.9, 0.1) %>% 
  map(~change_ve(params, ve_reduction = .)) -> all_params

lapply(1:9, function(x){
  all_params[[x]]$res %>% 
    map(cm_simulate) %>% 
    lapply(., "[[", "dynamics") %>% 
    bind_rows(.id = "ROS") %>% 
    filter(compartment == "death_o") %>% 
    group_by(ROS) %>% summarise(value = sum(value))
}) -> vary_ve

vary_ve %>% 
  map(mutate, rk = rank(value)) %>% 
  bind_rows(., .id = "ve_index") %>% 
  left_join(data.frame(ve_value = seq(0.1, 0.9,0.1)) %>% 
              rownames_to_column(var = "ve_index")) %>% 
  ggplot(., aes(x = ve_value, y = value, group = ROS, color = ROS)) +
  geom_line(size = 2) +
  theme_bw() +
  ggsci::scale_color_lancet() +
  labs(x = "VE reduction", y = "COVID-19 Mortality") +
  scale_y_log10()-> p

ggsave("figs/vary_ve.png", p, height = 9, width = 16)

# delay + high prevalence
16:40 %>% 
  map(~predict_outbreak(delay = .,
                        under_report = T)) -> all_params

lapply(1:25, function(x){
  all_params[[x]]$res %>% 
    map(cm_simulate) %>% 
    lapply(., "[[", "dynamics") %>% 
    bind_rows(.id = "ROS") %>% 
    filter(compartment == "death_o") %>% 
    group_by(ROS) %>% summarise(value = sum(value)) %>% 
    mutate(tot = value/DEoptim3_selected_debug$par3[i])
}) -> vary_delay_ur

vary_delay_ur %>% 
  map(mutate, rk = rank(value)) %>% 
  bind_rows(., .id = "delay") %>% 
  mutate(delay = as.numeric(delay) + 15,
         rk = factor(rk)) %>% 
  ggplot(., aes(x = delay, 
                y = tot, 
                group = ROS,
                color = ROS)) +
  geom_line(size = 2) +
  theme_bw() +
  ggsci::scale_color_lancet() +
  labs(x = "Supply Delay", y = "COVID-19 Mortality")+  
  coord_cartesian(ylim = c(6600, 12200)) -> p

ggsave("figs/vary_delay_ur.png", height = 9, width = 8)


gen_country_basics(country = DEoptim3_selected_debug$country_name[i],
                   waning_nat = 52*7*3,
                   R0_assumed  = DEoptim2_selected_debug$r[i],
                   date_start = DEoptim2_selected_debug$date[i],
                   date_end = "2022-12-31",
                   processes = burden_processes,
                   deterministic = T) -> params
res_d <- cm_simulate(params)
gen_country_basics(country = DEoptim3_selected_debug$country_name[i],
                   waning_nat = 52*7*3,
                   R0_assumed  = DEoptim2_selected_debug$r[i],
                   date_start = DEoptim2_selected_debug$date[i],
                   date_end = "2022-12-31",
                   processes = burden_processes,
                   deterministic = F) -> params
res <- cm_simulate(params, n = 100)

res$dynamics %>%
  filter(compartment == "death_o") %>%
  bind_rows(res_d$dynamics %>%
              filter(compartment == "death_o") %>%
              mutate(run = 0)) %>%
  mutate(date = ymd(DEoptim2_selected_debug$date[i]) + t) %>%
  group_by(run, t, population, date) %>%
  summarise(death_predicted = sum(value)) %>%
  left_join(owid_epi %>%
              filter(loc == "Georgia") %>%
              dplyr::select(-wb, -cases),
            by = c("date", "population" = "loc")) %>%
  mutate(flag = if_else(run == 0, "deterministic", "stochastic")) %>%
  pivot_wider(names_from = flag, values_from = death_predicted) %>%
  filter(!(is.na(deaths) & is.na(deterministic) & is.na(stochastic))) -> tmp

tmp %>% 
  filter(date <= ymd("2021-06-30")) %>% 
  ggplot(., aes(x = date)) +
  geom_line(aes(y = stochastic, group = run), alpha = 0.2, color = "grey70") +
  geom_line(aes(y = deterministic, group = run)) +
  geom_point(aes(y = deaths), color = "red") +
  theme_bw() +
  labs(x = "", y = "Daily COVID-19 Mortality") -> p

ggsave("figs/fit.png", height = 9, width = 16)
