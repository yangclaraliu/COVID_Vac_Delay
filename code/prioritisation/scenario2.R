# scenario 2
# (1) complete vaccinating 60+ with dose 1, reaching the uptake goal 
# (2) vaccinate other adults with dose 1
# (3) complete vaccinating 60+ with dose 2

daily_vac_scenarios[[2]] <- daily_vac
# tmp_tar <- c(paste0("Y", tmp_priorities[[1]]$age_group,"_d1"))
# tmp_tar2 <- c(paste0("Y", tmp_priorities[[1]]$age_group,"_d2"))
# tmp_pop_prop <- tmp_pop$n_pop[tmp_priorities[[1]]$age_group]/sum(tmp_pop$n_pop[tmp_priorities[[1]]$age_group])
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
# tmp2_tar <- c(paste0("Y", tmp_priorities[[2]]$age_group,"_d1"))
# tmp2_tar2 <- c(paste0("Y", tmp_priorities[[2]]$age_group,"_d2"))
# tmp2_pop_prop <- tmp_pop$n_pop[tmp_priorities[[2]]$age_group]/sum(tmp_pop$n_pop[tmp_priorities[[2]]$age_group])
daily_vac_scenarios[[2]][max(t_marker)+1,tmp2_tar] <- 
  as.list(tmp2_pop_prop*unlist(daily_vac_scenarios[[2]][max(t_marker2)+1,"supply_daily"] - 
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

rm(t_marker, t_marker2, t_marker3)
