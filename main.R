library(tidyverse); library(data.table)

##### load covidm #####
cm_path <- "code/covidm_for_fitting/"
cm_force_rebuild <- F
cm_build_verbose <- T
cm_version <- 2
source(paste0(cm_path, "/R/covidm.R"))
source("code/util_functions.R")
source("code/0_LoadData.R")
load("data/intermediate/params_3_vp_18.rdata")
# load(paste0(path_dropbox,"params_vp_all.rdata"))

# 
# #
# res_by_country_2 %>% 
#   map(filter, grepl("death|hosp|cases", compartment)) %>% 
#   map(group_by, scenario, t, compartment) %>% 
#   map(summarise, value = sum(value)) %>% 
#   map(pivot_wider, names_from = compartment, values_from = value) -> res_2
# 
# res_by_country_3 %>% 
#   map(filter, grepl("death|hosp|cases", compartment)) %>% 
#   map(group_by, scenario, t, compartment) %>% 
#   map(summarise, value = sum(value)) %>% 
#   map(pivot_wider, names_from = compartment, values_from = value) -> res_3
# 
# res_by_country_2_sw %>% 
#   map(filter, grepl("death|hosp|cases", compartment)) %>% 
#   map(group_by, scenario, t, compartment) %>% 
#   map(summarise, value = sum(value)) %>% 
#   map(pivot_wider, names_from = compartment, values_from = value) -> res_2_sw
# 
# res_by_country_3_sw %>% 
#   map(filter, grepl("death|hosp|cases", compartment)) %>% 
#   map(group_by, scenario, t, compartment) %>% 
#   map(summarise, value = sum(value)) %>% 
#   map(pivot_wider, names_from = compartment, values_from = value) -> res_3_sw
# 
# for(i in 12:nrow(model_selected)){
#   res_2[[i]] %<>% 
#     mutate(date = model_selected$t[i] + ymd("2019-12-01") + t)
#   
#   res_3[[i]] %<>% 
#     mutate(date = model_selected$t[i] + ymd("2019-12-01") + t)
# }
# 
# for(i in 1:11){
#   res_2_sw[[i]] %<>% 
#     mutate(date = model_selected$t[i] + ymd("2019-12-01") + t)
#   
#   res_3_sw[[i]] %<>% 
#     mutate(date = model_selected$t[i] + ymd("2019-12-01") + t)
# }
# 
# save(res_2, res_3, res_2_sw, res_3_sw,
#      file = paste0(path_dropbox,"res_vp.rdata"))
# 
# 
# res_2 %<>% 
#   map(mutate, VOC = if_else(date <= "2021-04-15", F, T)) %>% 
#   map(mutate, death_w_voc = if_else(VOC == T, death_voc_o, death_o)) %>% 
#   map(mutate, hosp_w_voc = if_else(VOC == T, hosp_voc_i, hosp_i)) %>% 
#   map(filter, date >= "2021-03-01")
# 
# res_3 %<>% 
#   map(mutate, VOC = if_else(date <= "2021-04-15", F, T)) %>% 
#   map(mutate, death_w_voc = if_else(VOC == T, death_voc_o, death_o)) %>% 
#   map(mutate, hosp_w_voc = if_else(VOC == T, hosp_voc_i, hosp_i)) %>% 
#   map(filter, date >= "2021-03-01")
# 
# 
# res_2_sw %<>% 
#   map(mutate, VOC = if_else(date <= "2021-04-15", F, T)) %>% 
#   map(mutate, death_w_voc = if_else(VOC == T, death_voc_o, death_o)) %>% 
#   map(mutate, hosp_w_voc = if_else(VOC == T, hosp_voc_i, hosp_i)) %>% 
#   map(filter, date >= "2021-03-01")
# 
# res_3_sw %<>% 
#   map(mutate, VOC = if_else(date <= "2021-04-15", F, T)) %>% 
#   map(mutate, death_w_voc = if_else(VOC == T, death_voc_o, death_o)) %>% 
#   map(mutate, hosp_w_voc = if_else(VOC == T, hosp_voc_i, hosp_i)) %>% 
#   map(filter, date >= "2021-03-01")
#   
# save(res_2, res_3, res_2_sw, res_3_sw,
#      file = paste0(path_dropbox,"res_vp.rdata"))
# 
# 
# res_3_sw %>% 
#   map(dplyr::select, -death_voc_o, -hosp_voc_i, -death_o, -hosp_i) %>% 
#   map(group_by, scenario) %>% 
#   map(summarise, 
#       cases = sum(cases),
#       death = sum(death_w_voc),
#       hosp = sum(hosp_w_voc)) %>% 
#   bind_rows(.id = "country_index") %>% 
#   left_join(model_selected %>% 
#               mutate(country_index = 1:13) %>% 
#               dplyr::select(country_index, country_name) %>% 
#               mutate(country_index = as.character(country_index)),
#             by = "country_index") %>% 
#   group_by(country_index) %>% group_split() %>% 
#   map(mutate, 
#       cases_p = cases/min(cases),
#       death_p = death/min(death),
#       hosp_p = hosp/min(hosp)) %>% 
#   bind_rows() %>% 
#   pivot_longer(cols = ends_with("_p")) %>% 
#   mutate(scenario = factor(scenario, levels = 1:7,
#                            labels = c(paste0(seq(4,20,4),"w"),
#                                       "Prior_Cov",
#                                       "Prior_Comp")),
#          name = factor(name,  levels = c("cases_p", "hosp_p", "death_p"),
#                        labels = c("Infections", "Severe Cases", "Mortality"))) %>% 
#   ggplot(., aes(x = scenario, y = value, color = name)) +
#   geom_boxplot() +
#   labs(y = "Prevention Protential Lost\n(Target Outcome)/(Best Achievable Outcome)") +
#   theme_cowplot() +
#   scale_color_lancet() +
#   labs(color = "", x = "Dosing Interval", title = "First Dose Waning = 120 Days") +
#   theme(legend.position = "top",
#         legend.text = element_text(size = 18),
#         title = element_text(size = 18),
#         axis.text = element_text(size = 14))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # %>% 
#   dplyr::select(-death_voc_o, -hosp_voc_i, -death_o, -hosp_i) %>%
#   # group_by(scenario) %>% group_split() %>%
#   # map(~mutate_at(., vars(c(cases, death_w_voc,hosp_w_voc)), cumsum)) %>% bind_rows() %>%
#   pivot_longer(cols = c(starts_with("death"), starts_with("hosp"), "cases")) %>% 
#   mutate(VOC = if_else(date >= "2021-04-15", F, T),
#          name = factor(name,
#                        levels = c("cases", "hosp_w_voc", "death_w_voc"),
#                        labels = c("Cases", "Severe Cases", "Deaths"))) %>% ungroup %>% 
#   # separate(name, into = c("metric", "stat", "voc"), remove = F) %>% 
#   ggplot(., aes(x = date, y = value, group = scenario, color = scenario)) +
#   geom_line(size = 1.2) +
#   facet_wrap(~name, ncol = 1, scales = "free") + theme_bw() +
#   labs(x = "", y = "") +
#   theme(strip.text = element_text(size = 16))
# 
# res_by_country[[1]] %>% 
#   bind_rows(.id = "scenario") %>% 
#   filter(grepl("death|hosp|cases", compartment)) %>% 
#   mutate(date = model_selected$t[6] + ymd("2019-12-01") + t) %>% 
#   mutate(VOC = if_else(date <= "2021-04-15", F, T)) %>% 
#   pivot_wider(names_from = compartment, values_from = value) %>% 
#   group_by(scenario, date, VOC) %>% 
#   summarise(death_o = sum(death_o),
#             death_voc_o = sum(death_voc_o),
#             hosp_i = sum(hosp_i),
#             hosp_voc_i = sum(hosp_voc_i),
#             cases = sum(cases)) %>% 
#   mutate(death_w_voc = if_else(VOC == T, death_voc_o, death_o),
#          hosp_w_voc = if_else(VOC == T, hosp_voc_i, hosp_i)) %>%
#   filter(date >= "2021-03-01") %>% 
#   dplyr::select(-death_voc_o, -hosp_voc_i, -death_o, -hosp_i) %>%
#   # group_by(scenario) %>% group_split() %>%
#   # map(~mutate_at(., vars(c(cases, death_w_voc,hosp_w_voc)), cumsum)) %>% bind_rows() %>%
#   pivot_longer(cols = c(starts_with("death"), starts_with("hosp"), "cases")) %>% 
#   group_by(scenario, name) %>% 
#   summarise(value = sum(value)) %>% 
#   filter(scenario %in% 1:5) -> test
# 
# test %>% 
#   mutate(scenario = factor(scenario, levels = c(1:5),
#                            labels = seq(4,20,4))) %>% 
#   ggplot(., aes(x = scenario, y = value, group = name)) +
#   geom_line() +
#   facet_wrap(~name, scales = "free")


# res %>% 
#   map(~.$dynamics) %>% 
#   bind_rows(.id = "scenario") %>% 
#   filter(compartment == "death_o") %>% 
#   group_by(scenario, t) %>% 
#   summarise(value = sum(value)) %>% 
#   ggplot(., aes(x = t, y = value, group = scenario, color = scenario)) +
#   geom_line()


# res$dynamics %>% 
#   filter(compartment %in% c("cases","death_o")) %>% 
#   ggplot(., aes(x = t, y = value, group = group)) +
#   geom_line() +
#   facet_wrap(~compartment, scales = "free") 
# 
# res$dynamics %>% 
#   filter(compartment %in% c("S", "Sw", "Sv", "Sv2",
#                             "R", "Rv", "Rv2")) %>% 
#   mutate(broad = substr(compartment, 1, 1) %>% 
#            factor(., levels = c("S","R"))) %>% 

#   ggplot(., aes(x = t, 
#                 y = value,
#                 color = compartment,
#                 group = compartment)) +
#   geom_line() +
#   facet_grid(broad~group)

# res$dynamics %>% 
#   filter(compartment %in% c("Sw", "Sv", "Sv2",
#                             "E", "Ev", "Ev2")) %>% 
#   # group_by(t, compartment) %>% summarise(value = sum(value)) %>% 
#   # mutate(compartment = factor(compartment, levels = c("Sv2", "Sv", "Sw"))) %>% 
#   ggplot(., aes(x = t, y = value, group = compartment, color = compartment)) +
#   geom_line() + facet_wrap(~group)
#   
# geom_bar(stat = "identity")
# 
# res$dynamics %>% 
#   filter(compartment %in% c("S","Sw", "Sv", "Sv2",
#                             "Ip","Ia","Is",
#                             "E", "Ev", "Ev2",
#                             "R", "Rv","Rv2"
#                             )) %>%  
#   group_by(t, compartment) %>% summarise(value = sum(value)) %>% 
#   mutate(broad = substr(compartment, 1, 1),
#          broad = if_else(compartment %in% c("Sw", "Sv","Sv2"), 
#                          "V", 
#                          as.character(broad)),
#          broad = factor(broad,
#                         levels = c("S","V","E","I","R"))) %>% 
#   ggplot(., aes(x = t, 
#                 y = value,
#                 color = compartment,
#                 group = compartment)) +
#   geom_bar(stat = "identity") +
#   facet_wrap(~broad, ncol = 1, scale = "free")


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
