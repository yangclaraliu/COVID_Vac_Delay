# scenario 1
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
  bind_cols(., daily_vac_scenarios[[1]]) -> daily_vac_scenarios[[1]]

for(i in 5:length(para$pop[[1]]$group_names)){
  tag1 <- paste0("Y",i,"_d1")
  tag2 <- paste0("Y",i,"_d2")
  tag3 <- paste0("Y",i,"_d1_S")
  tag4 <- paste0("Y",i,"_d1_V")
  tag5 <- paste0("Y",i,"_d1_SV2")
  tag6 <- paste0("Y",i,"_d1_VV2")
  
  for(j in 2:max(t_marker)){
    # new S
    daily_vac_scenarios[[1]][j, tag3] <- 
      (daily_vac_scenarios[[1]][j-1, tag4] +
         daily_vac_scenarios[[1]][j, tag1])*(para$pop[[1]]$wv[1])
    # new V
    daily_vac_scenarios[[1]][j, tag4] <- 
      (daily_vac_scenarios[[1]][j-1, tag4] +
         daily_vac_scenarios[[1]][j, tag1])*(1-para$pop[[1]]$wv[1])
  }
  
  daily_vac_scenarios[[1]][, tag3]  <- 
    cumsum(daily_vac_scenarios[[1]][, tag3])
  
  daily_vac_scenarios[[1]][max(t_marker) + 1, tag5] <- 
    daily_vac_scenarios[[1]][max(t_marker) + 1, tag2]
  daily_vac_scenarios[[1]][max(t_marker) + 1, tag2] <- 0
  daily_vac_scenarios[[1]][max(t_marker) + 1, tag4] <-
    daily_vac_scenarios[[1]][max(t_marker), tag4] 
  
  if(i >= 13){
    k = max(t_marker) + 2
    repeat{
      daily_vac_scenarios[[1]][k, tag4] <-
        (daily_vac_scenarios[[1]][k-1, tag4])*(1-para$pop[[1]]$wv[1])
      
      daily_vac_scenarios[[1]][k, tag3] <-
        daily_vac_scenarios[[1]][k-1, tag3] +
        (daily_vac_scenarios[[1]][k, tag4])*(para$pop[[1]]$wv[1]) -
        daily_vac_scenarios[[1]][k, tag2]
      
      daily_vac_scenarios[[1]][k, tag5] <-
        daily_vac_scenarios[[1]][k, tag2]
      
      if(daily_vac_scenarios[[1]][k, tag3] < 0 |
         is.na(daily_vac_scenarios[[1]][k, tag3])){
        daily_vac_scenarios[[1]][k, tag3] <- 0
        
        daily_vac_scenarios[[1]][k, tag5] <-
          daily_vac_scenarios[[1]][k-1, tag3] +
          (daily_vac_scenarios[[1]][k, tag4])*(para$pop[[1]]$wv[1])
        
        daily_vac_scenarios[[1]][k, tag6] <-
          daily_vac_scenarios[[1]][k, tag2] -
          (daily_vac_scenarios[[1]][k, tag5])
        
        daily_vac_scenarios[[1]][k, tag4] <-
          (daily_vac_scenarios[[1]][k-1, tag4])*(1-para$pop[[1]]$wv[1]) - 
          daily_vac_scenarios[[1]][k, tag6]
        k = k+1
        break
      }
      k = k+1
    }
    daily_vac_scenarios[[1]][k:nrow(daily_vac_scenarios[[1]]), tag3] <- 0
    
    repeat{
      # flux changes to V, waning out and then vaccinate out
      daily_vac_scenarios[[1]][k, tag4] <-
        (daily_vac_scenarios[[1]][k-1, tag4])*(1-para$pop[[1]]$wv[1]) -
        daily_vac_scenarios[[1]][k, tag2]
      
      # flux change to SV2 path
      daily_vac_scenarios[[1]][k, tag5] <-
        (daily_vac_scenarios[[1]][k, tag4])*(para$pop[[1]]$wv[1])
      
      # flux change to VV2 path 
      daily_vac_scenarios[[1]][k, tag6] <-
        (daily_vac_scenarios[[1]][k, tag2]) - 
        (daily_vac_scenarios[[1]][k, tag4])*(para$pop[[1]]$wv[1])
      
      if(k > max(t_marker2)){
        break
      }
      
      k = k + 1
    }
  }
  
}


daily_vac_scenarios[[1]] %>%
  dplyr::select(ends_with(c("_d1", "_d2", "SV2", "VV2", "V", "S")), date) %>%
  pivot_longer(cols = starts_with(c("Y"), ignore.case = F)) %>%
  separate(name, into = c("ag", "dose", "metric"), sep = "_") %>%
  mutate(ag = parse_number(ag) %>% factor, 
         value = if_else(value < 0, 0, value)) %>%
  filter(ag %in% c(13:16)) %>%
  unite("metric", c(dose, metric))%>% 
  mutate(metric = factor(metric,
                         levels = c("d1_NA", "d2_NA", "d1_S", "d1_V",
                                    "d1_SV2", "d1_VV2"),
                         labels = c("Daily Dose 1 Available",
                                    "Daily Dose 2 Available",
                                    "S from Dose 1 Waning",
                                    "Protected by V1",
                                    "S --> V2",
                                    "V1 --> V2"))) %>% 
  #filter(metric %in% c("d1_eff", "d1_NA")) %>%
  ggplot(., aes(x = date, y = value, color = ag, group = ag)) +
  # geom_bar(stat = "identity")+
  geom_line()+
  # geom_point() +
  facet_wrap(~metric, ncol = 1, scales = "free")

