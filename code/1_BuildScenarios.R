params_3 <- params_3_VOC <- list()
for(i in seq_len(nrow(model_selected_ur))){
  params_3[[i]] <-   params_3_VOC[[i]] <- gen_country_basics(country = model_selected_ur$country_name[i],
                                      waning_nat = 52*7*3,
                                      R0_assumed  = model_selected_ur$r[i],
                                      date_start = as.character(ymd("2019-12-01") + model_selected_ur$t[i]),
                                      date_end = "2022-12-31",
                                      processes = burden_processes,
                                      deterministic = TRUE)

  add_vac <- function(params){
    res <- params %>%
      update_vac_char(.,
                      ve_i   = ve$ve_i_o[1],  # infection blocking VE post 1 dose
                      v2e_i  = ve$ve_i_o[2],  # infection blocking VE post 2 doses
                      ve_d   = ve$ve_d[1],    # clinical fraction among breakthrough post 1 dose
                      v2e_d  = ve$ve_d[2],    # clinical fraction among breakthrough post 2 doses
                      wv = 1/360)
    return(res);rm(res)
  }

  add_vac_VOC <- function(params){
    res <- params %>%
      update_vac_char(.,
                      ve_i   = ve$ve_i_o[1],  # infection blocking VE post 1 dose
                      v2e_i  = ve$ve_i_o[2],  # infection blocking VE post 2 doses
                      ve_d   = ve$ve_d[1],    # clinical fraction among breakthrough post 1 dose
                      v2e_d  = ve$ve_d[2],    # clinical fraction among breakthrough post 2 doses
                      wv = 1/360) %>% # 1/ waning duration
      change_VOC(.,
                 date_switch = "2021-04-15",
                 rc_severity = 1.5,
                 rc_transmissibility = 1.5,
                 rc_ve = 0.4)
    return(res);rm(res)
  }

  params_3[[i]] %<>% add_vac(.)
  params_3_VOC[[i]] %<>% add_vac_VOC(.)
}


# # load(paste0("data/intermediate/params_3_vp_18.rdata"))
# n_scenario <- length(params_3_vp[[1]]$res)
# 
# for(i in 1:nrow(model_selected_ur)) {
#   params_3_VOC_vp[[i]] <- list()
#   params_3_VOC_vp[[i]][["param"]] <- params_3_VOC[[i]]
#   params_3_VOC_vp[[i]][["res"]] <- list()
#   for(j in 1:n_scenario){
#     # main body is from @params_2_VOC created above, which doesn't include vaccine roll-out schedule
#     params_3_VOC_vp[[i]][["res"]][[j]] <-  params_3_VOC[[i]]
#     # the vaccine roll-out schedule follows the corresponding countries in params_2_vp
#     params_3_VOC_vp[[i]][["res"]][[j]]$schedule[["v"]] <- params_3_vp[[i]]$res[[j]]$schedule[["v"]]
#     params_3_VOC_vp[[i]][["res"]][[j]]$schedule[["v2"]] <- params_3_vp[[i]]$res[[j]]$schedule[["v2"]]
#   }
#   params_3_VOC_vp[[i]][["scenarios"]] <- params_3_vp[[i]]$scenarios
# }

# save(params_3_vp, params_3_VOC_vp, file = paste0("data/intermediate/params_vp_18.rdata"))

add_vp <- function(params,
                   ms_cov = c(0,0.03,0.2,0.5),
                   delay = 24){
  res <- vac_policy(params,
                    # these two parameters define the supply conditions
                    milestone_date = c("2021-03-01", # start from 0
                                       "2021-06-30", # 0.03
                                       "2021-12-31", # all population; 0.2
                                       "2022-12-31"), # 0.6
                    milestone_cov = ms_cov,
                    # prioritisation, assume 60+  all prioritised
                    priority = c(NA, NA, NA, NA,
                                 2,  2,  2,  2,
                                 2,  2,  2,  2,
                                 1,  1,  1,  1),
                    # maximum feasible uptakes
                    cov_max = c(rep(0,4),
                                rep(0.7, 8),
                                rep(0.9, 4)),
                    supply_delay = as.numeric(delay), # unit = weeks
                    dose_interval = c(4, 8, 12, 16, 20))
  return(res);rm(res)
}

params_3_vp <-  params_3_VOC_vp <- list()
for(i in 1:length(euro_inuse)){
  which(model_selected_ur$iso3c ==   euro_inuse[i]) -> j
  params_3_vp[[i]] <- add_vp(params_3[[j]])
  params_3_VOC_vp[[i]] <- add_vp(params_3_VOC[[j]])
  print(round(i*100/length(euro_inuse),0))
  save(params_3_vp, params_3_VOC_vp,
       file = "data/intermediate/params_3_vp_delay.rdata") #this means immune delay
}
# 
# for(i in 1:13){
#   lapply(params_3_VOC_vp[[i]]$scenarios, "[[", "daily_vac_scenarios") %>%
#     map(arrange, date) %>%
#     map(mutate_at, vars(starts_with("Y", ignore.case = T)), cumsum) %>%
#     bind_rows(.id = "scenario") %>%
#     select(-starts_with("supply")) %>%
#     pivot_longer(starts_with("Y")) %>%
#     separate(name, into = c("ag", "dose")) %>%
#     group_by(scenario, date, dose) %>%
#     summarise(value = sum(value), .groups = "drop") -> x
# 
# 
#   lapply(params_3_VOC_vp[[i]]$scenarios, "[[", "daily_vac_scenarios") %>%
#     map(arrange, date) %>%
#     map(mutate_at, vars(starts_with("supply", ignore.case = T)), cumsum) %>%
#     bind_rows(.id = "scenario") %>%
#     dplyr::select(date, starts_with("supply"), -supply_cum) -> y
# 
#   y %>% filter(supply > 10000) %>% .[1,] %>% pull(date) -> t1
#   y %>% filter(supply_2 > 10000) %>% .[1,] %>% pull(date) -> t2
# 
#   print(t2-t1)
# 
#   x %>%
#     left_join(y, by = "date") %>%
#     filter(date >= "2021-01-01") %>%
#     mutate(pop =  vac_denom %>% filter(iso3c == euro_inuse[i]) %>% pull(tot),
#            value = value/pop,
#            supply = supply/pop,
#            supply_2 = supply_2/pop) %>%
#     ggplot(., aes(x = date, y = value, group = dose, color = dose)) +
#     geom_line() +
#     geom_line(aes(y = supply), color = "black", linetype = 2) +
#     geom_line(aes(y = supply_2), color = "black", linetype = 3) +
#     facet_wrap(~scenario, nrow = 1) +
#     labs(title = euro_inuse[i]) -> p
# 
#   ggsave(paste0("figs/supplemental/Delay_24wk_Sanity/",i,
#                 "_",euro_inuse[i],".png"),
#          p, width = 20, height = 6)
# }



# sanity check
# lapply(1:18, function(y){
# lapply(1:13, function(y){
#   lapply(1:7, function(x){
#     params_3_vp[[y]]$scenarios[[x]]$daily_vac_scenarios %>%
#       dplyr::select(t, starts_with("Y", ignore.case = F)) %>%
#       pivot_longer(starts_with("Y", ignore.case = F)) %>%
#       filter(value < 0)
#   })
# })


# for(i in 1:13){
#   lapply(params_3_vp[[i]]$res,"[[","schedule") %>% lapply(., "[[", "v") %>% 
#     lapply(., "[[","values") %>% unlist %>% .[.<0] %>% print
#   lapply(params_3_vp[[i]]$res,"[[","schedule") %>% lapply(., "[[", "v2") %>% 
#     lapply(., "[[","values") %>% unlist %>% .[.<0] %>% print
# }
# 
# for(i in 1:13){
#   lapply(params_3_vp[[i]]$scenarios, "[[", "daily_vac_scenarios") %>%
#     map(mutate_at, vars(starts_with("Y", ignore.case = F)), cumsum) %>%
#     bind_rows(.id = "scenarios") %>%
#     dplyr::select(-starts_with("supply"), -phase, -t) %>%
#     pivot_longer(cols = starts_with("Y")) %>%
#     separate(name, into = c("ag","dose")) %>%
#     group_by(date, scenarios, dose)  -> tmp
# 
#   tmp %>%
#     # filter(date >= "2021-03-01") %>%
#     ggplot(., aes(x = date, y = value, group = interaction(ag, dose), color = dose)) +
#     geom_line() +
#     facet_wrap(~scenarios) +
#     theme_bw() +
#     labs(title = model_selected_ur$country_name[i]) -> tmp_p
# 
#   ggsave(paste0("figs/intermediate/ROS_ag/",
#                 model_selected_ur$country_name[i],
#                 ".png"),
#          plot = tmp_p,
#          width = 15, height = 10)
# 
#   tmp %>%
#     summarise(value = sum(value))%>%
#     # filter(date >= "2021-03-01") %>%
#     ggplot(., aes(x = date, y = value, group = interaction(dose), color = dose)) +
#     geom_line() +
#     facet_wrap(~scenarios) +
#     theme_bw() +
#     labs(title = model_selected_ur$country_name[i]) -> tmp_p
# 
#   ggsave(paste0("figs/intermediate/ROS/",
#                 model_selected_ur$country_name[i],
#                 ".png"),
#          plot = tmp_p,
#          width = 15, height = 10)
# 
# }

  