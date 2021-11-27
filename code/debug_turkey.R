# Turkey
i = 12
params <- gen_country_basics(country = model_selected$country_name[i],
                                 waning_nat = 52*7*3,
                                 R0_assumed  = model_selected$r[i],
                                 date_start = as.character(ymd("2019-12-01") + model_selected$t[i]),
                                 date_end = "2022-12-31",
                                 processes = burden_processes,
                                 deterministic = TRUE) %>% 
  update_vac_char(.,
                  ve_i   = 0,  # infection blocking VE post 1 dose
                  v2e_i  = 0,  # infection blocking VE post 2 doses
                  ve_d   = 0,    # clinical fraction among breakthrough post 1 dose
                  v2e_d  = 0,    # clinical fraction among breakthrough post 2 doses
                  wv = 1/360)  %>% # 1/ waning duration
  change_VOC(.,
             date_switch = "2021-04-15",
             rc_severity = 1.5,
             rc_transmissibility = 1.5,
             rc_ve = 0.4)

para = params
milestone_date = c("2021-03-01", # start from 0
                   "2021-06-30", # 0.03
                   "2021-12-31", # all population; 0.2
                   "2022-12-31")
milestone_cov = c(0,
                  0.03,
                  0.2,
                  0.5)
priority = c(NA, NA, NA, NA,
             2,  2,  2,  2,
             2,  2,  2,  2,
             1,  1,  1,  1)
cov_max = c(rep(0,4),
            rep(0.7, 8),
            rep(0.9, 4))
supply_delay = 24
dose_interval = 4


daily_vac_scenarios %>% 
  mutate_at(vars(starts_with("Y", ignore.case = F)), cumsum) %>% 
  pivot_longer(cols = starts_with("Y", ignore.case = F)) %>% 
  separate(name, into = c("ag","dose")) %>% 
  group_by(date, dose) %>% summarise(value = sum(value)) %>% 
  ggplot(., aes(x = date, y = value, group = dose)) +
  geom_line()


test <- cm_simulate(params)[["dynamics"]]
test1 <- cm_simulate(params_VOC)[["dynamics"]]


bind_rows(test, test1, .id = "group") %>% 
  group_by(t, compartment, group) %>% summarise(value = sum(value)) %>% 
  filter(grepl("death_o", compartment)) %>% 
  mutate(date = ymd("2019-12-01") + t + model_selected$t[i]) %>% 
  ggplot(., aes(x = date, y = value, color = group)) +
  geom_line() +
  geom_vline(xintercept = ymd("2021-04-15"))



m <- ve_new[,c(2,3,6,7)]
params_VOC <- gen_country_basics(country = model_selected$country_name[i],
                                 waning_nat = 52*7*3,
                                 R0_assumed  = model_selected$r[i],
                                 date_start = as.character(ymd("2019-12-01") + model_selected$t[i]),
                                 date_end = "2022-12-31",
                                 processes = burden_processes,
                                 deterministic = TRUE) 

para <- params_VOC

ve_i  = m[1,1] # infection blocking VE post 1 dose
v2e_i = m[1,2] # infection blocking VE post 2 doses
ve_d  = m[1,3] # perc reduction in clinical fraction post 1 dose
v2e_d = m[1,4] # perc reduction in clinical fraction post 2 doses 

date_switch = "2021-04-15"

para$pop[[1]]$uv  <- para$pop[[1]]$u * (1 - ve_i)
para$pop[[1]]$uv2 <- para$pop[[1]]$u * (1 - v2e_i)
para$pop[[1]]$yv  <- para$pop[[1]]$y * (1 - ve_d)
para$pop[[1]]$yv2 <- para$pop[[1]]$y * (1 - v2e_d)

# make changes RE: transmissibility and ve
if(!(rc_transmissibility == 1 & rc_ve == 1)){
  data.table(u   =  para$pop[[1]]$u,
             uv  =  para$pop[[1]]$uv,
             uv2 =  para$pop[[1]]$uv2) %>% 
    mutate(ve_i  =  1 - uv/u,
           ve_i2 =  1 - uv2/u) -> u_table
  
  # change transmissibility
  if(rc_transmissibility != 1){
    u_table_voc <- 
      u_table %>% mutate_at(vars(starts_with("u")), function(x) x*rc_transmissibility)
  }
  
  # change in infection prevention ve
  if(rc_ve != 1){
    u_table_voc <- 
      u_table_voc %>% 
      mutate_at(vars(starts_with("ve")), function(x) x*rc_ve) %>% 
      mutate(uv = u*(1-ve_i),
             uv2 = u*(1-ve_i2))
  }
  
  for(p in c("u","uv","uv2")){
    para$schedule[[paste0("change_", p)]] <-  list(
      parameter = p,
      pops = numeric(),
      mode = "assign",
      values = list(u_table %>% pull(p), u_table_voc %>% pull(p)),
      times = c(para$date0, date_switch)
    )
  }
}

