# load(paste0(path_dropbox, "params_vp_all.rdata"))
# load("data/intermediate/params_vp_18.rdata")

n_scenario <- params_3_VOC_vp[[1]]$res %>% length

CJ(ve_i = seq(0.35,0.95,0.1), v2e_i = seq(0.35,0.95,0.1),
   ve_d_o = seq(0.35,0.95,0.1), v2e_d_o = seq(0.35,0.95,0.1)) %>% 
  filter(ve_i <= v2e_i,
         ve_d_o <= v2e_d_o,
         ve_i <= ve_d_o,
         v2e_i <= v2e_d_o) %>% 
  mutate(ve_d =  (ve_d_o - ve_i)/(1 - ve_i),
         v2e_d =  (v2e_d_o - v2e_i)/(1 - v2e_i)) %>% 
  rownames_to_column(var = "ve_set") -> ve_new

n_vac = nrow(ve_new)
n_country <- which(model_selected_ur$iso3c %in% euro_inuse)

# require(parall)
# no_cores <- detectCores() - 1
# cluster <- makeCluster(no_cores)
# plan(multisession)
pb <- progress_bar$new(total = n_vac*length(euro_inuse))
res_3_ve <- res_3_VOC_ve <- list()

for(i in n_country){
# for(i in 1){
  res_3_ve[[i]] <- res_3_VOC_ve[[i]] <- list()
  for(j in 1:n_vac){
  # for(j in 1){
    params_3_vp[[i]]$res %>% 
      map(~update_vac_char(.,
                           ve_i   = ve_new$ve_i[j],  
                           v2e_i  = ve_new$v2e_i[j], 
                           ve_d   = ve_new$ve_d[j], 
                           v2e_d  = ve_new$v2e_d[j],
                           wv = 1/360)) %>% 
      map(~  change_VOC(.,
                        date_switch = "2021-04-15",
                        rc_severity = 1.5,
                        rc_transmissibility = 1.5,
                        rc_ve = 0.4)) %>% 
      map(cm_simulate) %>% map(~.[["dynamics"]]) %>% 
      bind_rows(.id = "scenarios") %>% 
      filter(grepl("cases_reported|hosp_i|hosp_voc_i|death_", compartment)) %>% 
      mutate(date = ymd("2019-12-01") + model_selected_ur$t[i] + t,
             VOC = if_else(date <= "2021-04-15",F,T)) %>% 
      filter(date >= "2021-03-01") %>% 
      .[,sum(value), by = list(scenarios, population, compartment, VOC)] -> res_3_ve[[i]][[j]]
    
    params_3_VOC_vp[[i]]$res %>% 
      map(~update_vac_char(.,
                           ve_i   = ve_new$ve_i[j],  
                           v2e_i  = ve_new$v2e_i[j], 
                           ve_d   = ve_new$ve_d[j], 
                           v2e_d  = ve_new$v2e_d[j],
                           wv = 1/360)) %>% 
      map(~  change_VOC(.,
                        date_switch = "2021-04-15",
                        rc_severity = 1.5,
                        rc_transmissibility = 1.5,
                        rc_ve = 0.4)) %>% 
      map(cm_simulate) %>% map(~.[["dynamics"]]) %>% 
      bind_rows(.id = "scenarios") %>% 
      filter(grepl("cases_reported|hosp_i|hosp_voc_i|death_", compartment)) %>% 
      mutate(date = ymd("2019-12-01") + model_selected_ur$t[i] + t,
             VOC = if_else(date <= "2021-04-15",F,T)) %>% 
      filter(date >= "2021-03-01") %>% 
      .[,sum(value), by = list(scenarios, population, compartment, VOC)] -> res_3_VOC_ve[[i]][[j]]
    pb$tick()
  }
  res_3_ve[[i]] %<>% bind_rows(., .id = "ve_set")
  res_3_VOC_ve[[i]] %<>% bind_rows(., .id = "ve_set")
}
# write_rds(res_3_ve, paste0(path_dropbox,"intermediate/SA_VE_120.rds"))
# write_rds(res_3_VOC_ve, paste0(path_dropbox,"intermediate/SA_VE_120_VOC.rds"))
write_rds(res_3_ve, paste0("data/intermediate/SA_VE_360_v4.rds"))
write_rds(res_3_VOC_ve, paste0("data/intermediate/SA_VE_360_VOC_v4.rds"))

pb <- progress_bar$new(total = n_vac*nrow(model_selected_ur))
res_3_ve <- res_3_VOC_ve <- list()
for(i in n_country){
  res_3_ve[[i]] <- res_3_VOC_ve[[i]] <- list()
  for(j in 1:n_vac){
    params_3_vp[[i]]$res %>% 
      map(~update_vac_char(.,
                           ve_i   = ve_new$ve_i[j],  
                           v2e_i  = ve_new$v2e_i[j], 
                           ve_d   = ve_new$ve_d[j], 
                           v2e_d  = ve_new$v2e_d[j],
                           wv = 1/120)) %>% 
      map(~  change_VOC(.,
                        date_switch = "2021-04-15",
                        rc_severity = 1.5,
                        rc_transmissibility = 1.5,
                        rc_ve = 0.4)) %>% 
      map(cm_simulate) %>% map(~.[["dynamics"]]) %>% 
      bind_rows(.id = "scenarios") %>% 
      filter(grepl("cases_reported|hosp_i|hosp_voc_i|death_", compartment)) %>% 
      mutate(date = ymd("2019-12-01") + model_selected_ur$t[i] + t,
             VOC = if_else(date <= "2021-04-15",F,T)) %>% 
      filter(date >= "2021-03-01") %>% 
      .[,sum(value), by = list(scenarios, population, compartment, VOC)] -> res_3_ve[[i]][[j]]
    
    params_3_VOC_vp[[i]]$res %>% 
      map(~update_vac_char(.,
                           ve_i   = ve_new$ve_i[j],  
                           v2e_i  = ve_new$v2e_i[j], 
                           ve_d   = ve_new$ve_d[j], 
                           v2e_d  = ve_new$v2e_d[j],
                           wv = 1/120)) %>% 
      map(~  change_VOC(.,
                        date_switch = "2021-04-15",
                        rc_severity = 1.5,
                        rc_transmissibility = 1.5,
                        rc_ve = 0.4)) %>% 
      map(cm_simulate) %>% map(~.[["dynamics"]]) %>% 
      bind_rows(.id = "scenarios") %>% 
      filter(grepl("cases_reported|hosp_i|hosp_voc_i|death_", compartment)) %>% 
      mutate(date = ymd("2019-12-01") + model_selected_ur$t[i] + t,
             VOC = if_else(date <= "2021-04-15",F,T)) %>% 
      filter(date >= "2021-03-01") %>% 
      .[,sum(value), by = list(scenarios, population, compartment, VOC)] -> res_3_VOC_ve[[i]][[j]]
    pb$tick()
  }
  res_3_ve[[i]] %<>% bind_rows(., .id = "ve_set")
  res_3_VOC_ve[[i]] %<>% bind_rows(., .id = "ve_set")
}
# write_rds(res_3_ve, paste0(path_dropbox,"intermediate/SA_VE_120.rds"))
# write_rds(res_3_VOC_ve, paste0(path_dropbox,"intermediate/SA_VE_120_VOC.rds"))
write_rds(res_3_ve, paste0("data/intermediate/SA_VE_120_v4.rds"))
write_rds(res_3_VOC_ve, paste0("data/intermediate/SA_VE_120_VOC_v4.rds"))













# 
# pb <- progress_bar$new(total = n_vac*nrow(model_selected))
# lapply(1:nrow(model_selected), function(y){
#   lapply(1:n_vac, function(x){
#     pb$tick()
#     params_3_VOC_ve[[y]][[x]] %>% map(cm_simulate) %>% map(~.[["dynamics"]]) %>% 
#       bind_rows(.id = "scenarios") %>% data.table %>% 
#       filter(grepl("cases_reported|hosp_i|hosp_voc_i|death_", compartment)) %>% 
#       mutate(date = ymd("2019-12-01") + model_selected$t[y] + t,
#              VOC = if_else(date <= "2021-04-15",F,T)) %>% 
#       filter(date >= "2021-03-01") %>% 
#       group_by(scenarios, population, compartment, VOC) %>% 
#       summarise(value = sum(value), .groups = "drop")
#   }) %>% bind_rows(.id = "ve_set")
# }) -> test
# write_rds(test, paste0(path_dropbox,"intermediate/SA_VE_360_VOC.rds"))

















# 
# params_3_VOC <- params_3_ve <- list()
# for(i in nrow(model_selected)){
#   params_2_ve[[i]] <- params_3_ve[[i]] <- list()
#   for(j in 1:nrow(ve_new)){
#     params_3_vp[[i]]$res %>% 
#       map(~update_vac_char(.,
#                            ve_i   = ve_new$ve_i[j],  
#                            v2e_i  = ve_new$v2e_i[j], 
#                            ve_d   = ve_new$ve_d[j], 
#                            v2e_d  = ve_new$v2e_d[j],
#                            wv = 1/120))%>% 
#       map(~  change_VOC(.,
#                         date_swtich = "2021-04-15",
#                         rc_severity = 1.5,
#                         rc_transmissibility = 1.5,
#                         rc_ve = 0.4)) -> params_3_ve[[i]][[j]]
#     
#     params_3_VOC_vp[[i]]$res %>% 
#       map(~update_vac_char(.,
#                            ve_i   = ve_new$ve_i[j],  
#                            v2e_i  = ve_new$v2e_i[j], 
#                            ve_d   = ve_new$ve_d[j], 
#                            v2e_d  = ve_new$v2e_d[j],
#                            wv = 1/120)) %>% 
#       map(~  change_VOC(.,
#                         date_swtich = "2021-04-15",
#                         rc_severity = 1.5,
#                         rc_transmissibility = 1.5,
#                         rc_ve = 0.4)) -> params_2_ve[[i]][[j]]
#   }
# }
# 
# pb <- progress_bar$new(total = n_vac*nrow(model_selected))
# lapply(10, function(y){
#   lapply(1:n_vac, function(x){
#     pb$tick()
#     params_3_ve[[y]][[x]] %>% map(cm_simulate) %>% lapply(., "[[", "dynamics") %>% 
#       bind_rows(.id = "scenarios") %>% filter(grepl("cases_reported|hosp_i|hosp_voc_i|death_", compartment)) %>% 
#       mutate(date = ymd("2019-12-01") + model_selected$t[y] + t,
#              VOC = if_else(date <= "2021-04-15",F,T)) %>% 
#       filter(date >= "2021-03-01") %>% 
#       group_by(scenarios, population, compartment, VOC) %>% 
#       summarise(value = sum(value), .groups = "drop")
#   }) %>% bind_rows(.id = "ve_set")
# }) -> test
# write_rds(test, paste0(path_dropbox,"intermediate/SA_VE_120.rds"))
# 
# #### plotting results ####
# SA_VE_120 <- read_rds(paste0(path_dropbox,"intermediate/SA_VE_120_Serbia.rds"))
# SA_VE_360 <- read_rds(paste0(path_dropbox,"intermediate/SA_VE_360_Serbia.rds"))
# 
# SA_VE_120 %>% 
#   map(separate, col = "compartment", into = c("metric","voc","seg2")) %>% 
#   bind_rows() %>% 
#   filter(!(metric == "death" & voc == "o" & VOC == T)) %>% 
#   filter(!(metric == "death" & voc == "voc" & VOC == F))%>% 
#   filter(!(metric == "hosp" & voc == "o" & VOC == T)) %>% 
#   filter(!(metric == "hosp" & voc == "voc" & VOC == F)) %>% 
#   dplyr::select(-voc, -seg2) %>% 
#   group_by(ve_set, scenarios, population, metric) %>% 
#   summarise(value = sum(value), .groups = "drop") %>% 
#   left_join(ve_new, by = "ve_set") %>% 
#   mutate(ve_set = as.numeric(ve_set)) -> tmp
# 
# tmp %>% 
#   filter(metric == "death") %>% 
#   ggplot(., aes(x = scenarios, y = value, group = ve_set)) +
#   geom_point() +
#   geom_line(alpha = 0.2) +
#   facet_wrap(~population, scales = "free") 
# 
# tmp %>% 
#   filter(metric == "death" & population == "Serbia") %>% 
#   ggplot(., aes(x = ve_i, y = v2e_i, fill = value, base = 10)) +
#   geom_tile() +
#   scale_fill_viridis() +
#   facet_grid(~scenarios)
# 
# tmp
# 
# 
# ve_new %>% 
#   mutate(ve_set = as.numeric(ve_set)) %>% 
#   arrange(ve_set)
# 
# 
# i = 5
# 
# params_3_vp[[i]]$res %>% 
#   map(~update_vac_char(.,
#                        ve_i   = ve_new$ve_i[334],  
#                        v2e_i  = ve_new$v2e_i[334], 
#                        ve_d   = ve_new$ve_d[334], 
#                        v2e_d  = ve_new$v2e_d[334],
#                        wv = 1/360)) %>% 
#   map(cm_simulate) -> test_high
# 
# params_3_vp[[i]]$res %>% 
#   map(~update_vac_char(.,
#                        ve_i   = ve_new$ve_i[4],  
#                        v2e_i  = ve_new$v2e_i[4], 
#                        ve_d   = ve_new$ve_d[4], 
#                        v2e_d  = ve_new$v2e_d[4],
#                        wv = 1/360)) %>% 
#   map(cm_simulate) -> test_low
# 
# test_high[[1]]$dynamics %>% filter(compartment == "death_o") %>% 
#   group_by(population) %>% summarise(value = sum(value))
# 
# test_low[[1]]$dynamics %>% filter(compartment == "death_o") %>% 
#   group_by(population) %>% summarise(value = sum(value))
# 
# test_high[[1]]$dynamics %>% mutate(test = "high") %>% 
#   bind_rows(test_low[[1]]$dynamics %>% mutate(test = "low")) %>%
#   # pivot_wider(names_from = test, values_from = value) %>% 
#   filter(compartment == "death_o") %>%
#   mutate(date = ymd("2019-12-01") + t + model_selected$t[10]) %>% 
#   group_by(t, test, compartment, date) %>% 
#   summarise(value = sum(value)) %>%
#   ggplot(., aes(x = date, y = value, group = test, color = test)) +
#   geom_line() 
# 
