# scenario 3
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
