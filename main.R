library(tidyverse); library(data.table)

##### load covidm #####
cm_path <- "code/covidm_for_fitting/"
cm_force_rebuild <- F
cm_build_verbose <- T
cm_version <- 2
source(paste0(cm_path, "/R/covidm.R"))

source("code/util_functions.R")
source("code/0_LoadData.R")

params <- gen_country_basics(country = "Belarus",
                             waning_nat = 52*7*3,
                             R0_assumed  = 2.7,
                             date_start = "2020-01-01",
                             date_end = "2022-12-31",
                             processes = burden_processes,
                             deterministic = TRUE)

params_1 <-   update_vac_char(params,
                              ve_i   = ve$ve_i_o[1],  # infection blocking VE post 1 dose
                              v2e_i  = ve$ve_i_o[2],  # infection blocking VE post 2 doses
                              ve_d   = ve$ve_d[1],    # clinical fraction among breakthrough post 1 dose
                              v2e_d  = ve$ve_d[2],    # clinical fraction among breakthrough post 2 doses
                              wv = 1/120) # 1/ waning duration
params_2 <- vac_policy(params_1,
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
                       supply_delay = 20 # unit = weeks
                       
)

res <- cm_simulate(params_2$res[[1]])
res$dynamics %>% filter(compartment == "severe_i") %>% mutate(tot = sum(value))

res$dynamics %>% 
  filter(compartment %in% c("cases","death_o")) %>% 
  ggplot(., aes(x = t, y = value, group = group)) +
  geom_line() +
  facet_wrap(~compartment, scales = "free") 

res$dynamics %>% 
  filter(compartment %in% c("S", "Sw", "Sv", "Sv2",
                            "R", "Rv", "Rv2")) %>% 
  mutate(broad = substr(compartment, 1, 1) %>% 
           factor(., levels = c("S","R"))) %>% 
  ggplot(., aes(x = t, 
                y = value,
                color = compartment,
                group = compartment)) +
  geom_line() +
  facet_grid(broad~group)

res$dynamics %>% 
  filter(compartment %in% c("Sw", "Sv", "Sv2",
                            "E", "Ev", "Ev2")) %>% 
  # group_by(t, compartment) %>% summarise(value = sum(value)) %>% 
  # mutate(compartment = factor(compartment, levels = c("Sv2", "Sv", "Sw"))) %>% 
  ggplot(., aes(x = t, y = value, group = compartment, color = compartment)) +
  geom_line() + facet_wrap(~group)
  
geom_bar(stat = "identity")

res$dynamics %>% 
  filter(compartment %in% c("S","Sw", "Sv", "Sv2",
                            "Ip","Ia","Is",
                            "E", "Ev", "Ev2",
                            "R", "Rv","Rv2"
                            )) %>%  
  group_by(t, compartment) %>% summarise(value = sum(value)) %>% 
  mutate(broad = substr(compartment, 1, 1),
         broad = if_else(compartment %in% c("Sw", "Sv","Sv2"), 
                         "V", 
                         as.character(broad)),
         broad = factor(broad,
                        levels = c("S","V","E","I","R"))) %>% 
  ggplot(., aes(x = t, 
                y = value,
                color = compartment,
                group = compartment)) +
  geom_bar(stat = "identity") +
  facet_wrap(~broad, ncol = 1, scale = "free")


# params <- gen_country_basics("Thailand") 
# 
# cm_simulate(params)
# 
# 
# 

#                     ) 
# 
# params <- vac_policy(para = params,
#                      # these two parameters define the supply conditions
#                      milestone_date = c("2021-01-01", # start from 0
#                                         "2021-06-30", # 0.03
#                                         "2021-12-31", # all population; 0.2
#                                         "2022-12-31"), # 0.6
#                      milestone_cov = c(0,
#                                        0.03,
#                                        0.2,
#                                        0.5),
#                      # prioritisation, assume 60+  all prioritised
#                      priority = c(NA, NA, NA, NA,
#                                   2,  2,  2,  2,
#                                   2,  2,  2,  2,
#                                   1,  1,  1,  1),
#                      # maximum feasible uptakes
#                      cov_max = c(rep(0,2),
#                                  rep(0.7, 10),
#                                  rep(0.9, 4)),
#                      # the proportion of age groups at which we will compare
#                      # rolling out dose 2s or more dose 1s. This value needs to
#                      # be smaller than the maximum uptake of the corresponding #
#                      # age group.
#                      p_change = 0.7,
#                      supply_delay = 20 # unit = weeks
# )
# 
# params2 <-   change_ve(params,
#                        transition_start = "2021-03-15",
#                        transition_end = "2021-06-15",
#                        ve_reduction = 0.5) 
# 
# params3 <-   change_ve(params,
#                       transition_start = "2021-03-15",
#                       transition_end = "2021-06-15",
#                       ve_reduction = 0.9) 
# 
# run1 = cm_simulate(params2$res[[1]])
# run2 = cm_simulate(params2$res[[2]])
# run3 = cm_simulate(params2$res[[3]])
# 
# run4 = cm_simulate(params3$res[[1]])
# run5 = cm_simulate(params3$res[[2]])
# run6 = cm_simulate(params3$res[[3]])
# 
# run1$dynamics %>% mutate(run = 1) %>% 
#   bind_rows(run2$dynamics %>% mutate(run = 2)) %>% 
#   bind_rows(run3$dynamics %>% mutate(run = 3)) %>% 
#   bind_rows(run4$dynamics %>% mutate(run = 4)) %>% 
#   bind_rows(run5$dynamics %>% mutate(run = 5)) %>% 
#   bind_rows(run6$dynamics %>% mutate(run = 6)) %>% 
#   separate(group, into = c("LL","UL"), remove = F) %>% 
#   mutate(LL = as.numeric(LL),
#          date = ymd(params2$param$date0) + t,
#          reduction = if_else(run %in% 4:6, 
#                              "Large Reduction", 
#                              "Small Reduction") %>% 
#            factor,
#          scenario = case_when(run %in% c(1,4) ~ 1,
#                               run %in% c(2,5) ~ 2,
#                               run %in% c(3,6) ~ 3) %>% factor,
#          run = factor(run)) %>% 
#   filter(compartment == "death_o",
#          date > "2021-01-01",
#          LL > 50) %>% 
#   ggplot(., aes(x = date, y = value, group = reduction, color = reduction)) +
#   geom_line(size = 1.5) +
#   facet_grid(scenario~group, scales = "free")
# 
