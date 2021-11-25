res_3 <- res_3_VOC <- list()
n_scenario <- params_3_vp[[1]]$res %>% length

pb <- progress::progress_bar$new(total = nrow(model_selected_ur))
for(i in 1:nrow(model_selected_ur)){
  pb$tick()
  # res_2[[i]] <- lapply(1:n_scenario, function(x) cm_simulate(params_2_vp[[i]]$res[[x]]))
  # res_2[[i]]  <- lapply(res_2[[i]], "[[" ,"dynamics") %>% 
  #   bind_rows(.id = "scenario")
  
  res_3[[i]] <- lapply(1:n_scenario, function(x) cm_simulate(params_3_vp[[i]]$res[[x]]))
  res_3[[i]]  <- lapply(res_3[[i]], "[[" ,"dynamics") %>% 
    bind_rows(.id = "scenario")
  
  # res_2_VOC[[i]] <- lapply(1:n_scenario, function(x) cm_simulate(params_2_VOC_vp[[i]]$res[[x]]))
  # res_2_VOC[[i]]  <- lapply(res_2_VOC[[i]], "[[" ,"dynamics") %>% 
  #   bind_rows(.id = "scenario")
  
  res_3_VOC[[i]] <- lapply(1:n_scenario, function(x) cm_simulate(params_3_VOC_vp[[i]]$res[[x]]))
  res_3_VOC[[i]]  <- lapply(res_3_VOC[[i]], "[[" ,"dynamics") %>% 
    bind_rows(.id = "scenario")
}


#
sum_endpoints <- function(res_list,
                          # date of VOC emergence
                          ds = NULL){
  rl_name <- deparse(substitute(res_list))
  VOC <- grepl("VOC", rl_name)
  
  if(VOC == F){
    # w/o VOC, get cases, hosp_i, death_o
    res_list %>% bind_rows() %>% 
      filter(compartment %in% c("cases", "hosp_i", "death_o")) %>% 
      group_by(compartment) %>% group_split() -> tmp
    
    tmp[[1]] %>% group_by(scenario, population, compartment) %>% 
      summarise(value = sum(value), .groups = "drop") %>% 
      mutate(group = factor(NA)) -> seg1
    
    tmp[[2]] %>% group_by(scenario, population, compartment) %>% 
      summarise(value = sum(value), .groups = "drop") %>% 
      mutate(group = factor(NA)) -> seg2
    
    tmp[[3]] %>% group_by(scenario, group, population, compartment) %>% 
      summarise(value = sum(value), .groups = "drop") %>% 
      dplyr::select(colnames(seg2))  -> seg3
    
    bind_rows(seg1, seg2, seg3) -> tmp
    return(tmp)
  }
  
  if(VOC == T){
    # w/o VOC, get cases, hosp_i, hosp_voc_i, death_o, death_voc_o
    res_list %>% bind_rows() %>% 
      filter(compartment %in% c("cases", "hosp_i", "death_o",
                                "hosp_voc_i", "death_voc_o")) %>% 
      left_join(model_selected_ur[,c("country_name", "t")] %>% 
                  rename(t_start = t), 
                by = c("population" = "country_name")) %>% 
      mutate(metric = substr(compartment, 1,4),
             date = t_start + t + ymd("2019-12-01"),
             metric = factor(metric,
                             levels = c("case", "hosp", "death"))) %>% 
      group_by(metric) %>% group_split() -> tmp
    
    tmp[[1]] %>% group_by(scenario, t_start, population, compartment) %>% 
      summarise(value = sum(value), .groups = "drop") %>% 
      mutate(group = factor(NA)) -> seg1
    
    tmp[[2]] %>% 
      mutate(VOC = date >= ds) %>% 
      filter(!(compartment == "hosp_voc_i" & VOC == F)) %>% 
      filter(!(compartment == "hosp_i" & VOC == T)) %>% 
      group_by(scenario, t_start, population) %>% 
      summarise(value = sum(value), .groups = "drop") %>% 
      mutate(group = factor(NA), compartment = "hosp") %>% 
      dplyr::select(colnames(seg1)) -> seg2
    
    tmp[[3]] %>% 
      mutate(VOC = date >= ds) %>% 
      filter(!(compartment == "death_voc_o" & VOC == F)) %>% 
      filter(!(compartment == "death_o" & VOC == T)) %>% 
      group_by(scenario, t_start, group, population) %>% 
      summarise(value = sum(value), .groups = "drop") %>%
      mutate(compartment = "death") %>% 
      dplyr::select(colnames(seg2))  -> seg3
    
    bind_rows(seg1, seg2, seg3) -> tmp
    return(tmp)
  }
}

dyna <- list()
# dyna[["res_2"]] <- sum_endpoints(res_2)
dyna[["res_3"]] <- sum_endpoints(res_3)
# dyna[["res_2_VOC"]] <- sum_endpoints(res_2_VOC, ds = "2021-04-15")
dyna[["res_3_VOC"]] <- sum_endpoints(res_3_VOC, ds = "2021-04-15")

# rm(res_2, res_3, res_2_VOC, res_3_VOC)
write_rds(dyna, "data/intermediate/res_baseline_v3.rds")

#### shorter immunity duration ####
pb <- progress::progress_bar$new(total = nrow(model_selected_ur))

# params_2_vp_sw <- params_2_vp
params_3_vp_sw <- params_3_vp
# params_2_VOC_vp_sw <- params_2_VOC_vp
params_3_VOC_vp_sw <- params_3_VOC_vp

for(j in 1:nrow(model_selected_ur)){
  for(i in 1:n_scenario){
    params_2_vp_sw[[j]]$res[[i]]$pop[[1]]$wv <- rep(1/120, 16)
    params_3_vp_sw[[j]]$res[[i]]$pop[[1]]$wv <- rep(1/120, 16)
    params_2_VOC_vp_sw[[j]]$res[[i]]$pop[[1]]$wv <- rep(1/120, 16)
    params_3_VOC_vp_sw[[j]]$res[[i]]$pop[[1]]$wv <- rep(1/120, 16)
  }
}

res_3_sw <- res_3_VOC_sw <-  list()

pb <- progress::progress_bar$new(total = nrow(model_selected_ur))
for(i in 1:nrow(model_selected_ur)){
  pb$tick()
  
  res_3_sw[[i]] <- lapply(1:n_scenario, function(x) cm_simulate(params_3_vp_sw[[i]]$res[[x]]))
  res_3_sw[[i]]  <- lapply(res_3_sw[[i]], "[[" ,"dynamics") %>% 
    bind_rows(.id = "scenario")
  
  res_3_VOC_sw[[i]] <- lapply(1:n_scenario, function(x) cm_simulate(params_3_VOC_vp_sw[[i]]$res[[x]]))
  res_3_VOC_sw[[i]]  <- lapply(res_3_VOC_sw[[i]], "[[" ,"dynamics") %>% 
    bind_rows(.id = "scenario")
}

dyna_sw <- list()
dyna_sw[["res_3_sw"]] <- sum_endpoints(res_3_sw)
dyna_sw[["res_3_VOC_sw"]] <- sum_endpoints(res_3_VOC_sw, ds = "2021-04-15")

# rm(res_2_sw, res_3_sw, res_2_VOC_sw, res_3_VOC_sw)
write_rds(dyna_sw, paste0("data/intermediate/res_sw_v3.rds"))
