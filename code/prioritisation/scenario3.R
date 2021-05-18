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
# 
# daily_vac_scenarios[[3]][k, ] %>% t
# daily_vac_scenarios[[3]] %>%
#   dplyr::select(ends_with(c("_d1", "_d2", "S","V","SV2", "VV2")), date) %>%
#   pivot_longer(cols = starts_with(c("Y"), ignore.case = F)) %>%
#   separate(name, into = c("ag", "dose", "metric"), sep = "_") %>%
#   filter(ag > 4) %>%
#   mutate(ag = parse_number(ag) %>% factor) %>%
#   unite("metric", c(dose, metric))%>%
#   #filter(metric %in% c("d1_eff", "d1_NA")) %>%
#   ggplot(., aes(x = date, y = value, color = ag, group = ag)) +
#   # geom_bar(stat = "identity")+
#   geom_line()+
#   # geom_point() +
#   facet_wrap(~metric, ncol = 1, scales = "free")
# 
# a <- round((daily_vac_scenarios[[3]]$Y16_d2),2)
# b <- round(daily_vac_scenarios[[3]]$Y16_d1_SV2 + (daily_vac_scenarios[[3]]$Y16_d1_VV2),2)
# sum(a)
# sum(b)
# sum(a)-sum(b)
# 
# round((daily_vac_scenarios[[3]]$Y16_d1_V),2) %>% tail
# 
# where_diff <- which(!a==b)
# a[where_diff]
# b[where_diff]
# 
# 
