#### Response to reviewer, Round 1 #####
# 2022/02/02

#### Low supply condition scenario ####
low_supply <- read_rds("data/R1_low_supply.rds")

low_supply$res %>% 
  map(cm_simulate)

lapply(low_supply$scenarios, "[[", "daily_vac_scenarios") %>% 
  map(mutate_at, vars(starts_with("Y", ignore.case = F)), cumsum) %>% 
  map(pivot_longer, cols = starts_with("Y", ignore.case = F)) %>% 
  map(separate, name, into = c("ag", "dose")) %>% 
  map(dplyr::select, -starts_with("supply"), -phase) %>% 
  bind_rows(.id = "strategy") %>% 
  group_by(strategy, date, dose) %>% 
  summarise(value = sum(value), .groups = "drop") %>% 
  mutate(iso3c = "ALB") %>% 
  left_join(vac_denom, by = "iso3c") %>% 
  mutate(cov = value/tot,
         strategy = factor(strategy,
                           levels = c(1:5, 7, 6),
                           labels = strategy_labels),
         dose = factor(dose, levels = c("d1", "d2"),
                       labels = c("First", "Second"))) -> tmp

ref <- data.table(date = ymd(c("2021-03-01", # start from 0
                           "2021-06-30", # 0.03
                           "2021-12-31", # all population; 0.2
                           "2022-12-31")),
                  val = c(0, 0.01, 0.05, 0.1), dose = "d1")

ref %<>% 
  mutate(date = date + 24*7,
         dose = "d2") %>% 
  bind_rows(ref) %>% 
  mutate(dose = factor(dose,
                       levels = c("d1", "d2"),
                       labels = c("First", "Second")))
  
tmp %>% 
  ggplot(., aes(x = date, y = cov)) +
  geom_line(aes(color = dose)) +
  geom_line(data = ref,
            aes(x = date, y = val, group = dose, linetype = dose)) +
  scale_linetype_manual(values = c(2,3)) +
  facet_wrap(~strategy, nrow = 2) +
  theme_bw() +
  theme(legend.position = "top") +
  labs(x = "Date", y = "Vaccine Coverage", color = "Allocated", linetype = "Supply") +
  ggsci::scale_color_lancet() -> p

ggsave("figs/R1_low_supply.png", p, width = 12, height = 6)

#### sensitivity analyses around wv #####
# Loading dyna_sw_R1 and dyna_VOC_sw_R1
load("data/intermediate/res_sw_R1.rdata")
# Loading 
res_sw_v4 <- read_rds("data/intermediate/res_sw_v4.rds")
res_v4 <- read_rds("data/intermediate/res_baseline_v4.rds")

res_v4$res_3 %>% 
  mutate(wv = 360) %>% 
  bind_rows(res_sw_v4$res_3_sw %>% 
              mutate(wv = 120)) %>% 
  bind_rows(dyna_sw_R1 %>% 
              setNames(c(60,90)) %>% 
              bind_rows(.id = "wv") %>% 
              mutate(wv = as.numeric(wv))) -> tmp

res_v4$res_3_VOC %>% 
  mutate(wv = 360) %>% 
  bind_rows(res_sw_v4$res_3_VOC_sw %>% 
              mutate(wv = 120)) %>% 
  bind_rows(dyna_VOC_sw_R1 %>% 
              setNames(c(60,90)) %>% 
              bind_rows(.id = "wv") %>% 
              mutate(wv = as.numeric(wv))) -> tmp_VOC

tmp_VOC %>% 
  filter(compartment == "death") %>% 
  group_by(scenario, population, wv) %>% 
  summarise(value = sum(value), .groups = "drop") %>% 
  mutate(scenario = factor(scenario, levels = c(1:5,7,6), labels = strategy_labels),
         wv = factor(wv)) %>% 
  mutate(di = case_when(scenario == "A1" ~ 4,
                        scenario == "A2" ~ 8,
                        scenario == "A3" ~ 12,
                        scenario == "A4" ~ 16,
                        scenario == "A5" ~ 20,
                        scenario == "B1" ~ 25,
                        scenario == "B2" ~ 50)) %>% 
  group_by(population, wv) %>% mutate(rk = rank(value)) -> tab

tab %>% 
  filter(rk == 1) %>% 
  dplyr::select(scenario, population, di, wv) %>% 
  rename(optimal = di) %>% 
  right_join(tab, by = c("population", "wv")) %>% 
  select(-rk) %>% 
  ggplot(., aes(color = wv)) +
  geom_line(aes(x = di, y = value, group = wv)) +
  geom_point(aes(x = di, y = value)) +
  geom_vline(aes(xintercept = optimal, color = wv, group = wv), show.legend = F) +
  facet_wrap(~population, scales = "free", nrow = 5) +
  ggsci::scale_color_lancet() +
  scale_x_continuous(breaks = c(4,8,12,16,20,25,50),
                     labels = c("A1\n4w", "A2\n8w", "A3\n12w", "A4\n16w",
                                "A5\n20w","B1\n22-34w","B2\n43-56w")) +
  theme_bw() + 
  theme(legend.position = "top") +
  labs(x = "Dosing Interval Strategy", color = "First Dose Waning Duration",
       y = "Cumulative Mortality by end of '22") -> p


ggsave("figs/R1_shorter_wv.png", p, width = 15, height = 10)

#### additional delay functions ####
additional_delay <- c(12,52)

# params_3_vp <-  params_3_VOC_vp <- list()
for(m in 1:length(additional_delay)){
  params_3_vp[[m]] <- params_3_VOC_vp[[m]] <- list()
  # delay = 12, i = 3, 6, 12
  # delay = 52, i = 5, 11, 13
  for(i in 12:length(euro_inuse)){
    which(model_selected_ur$iso3c ==  euro_inuse[i]) -> j
    params_3_vp[[m]][[i]] <- add_vp(params = params_3[[j]], 
                                    delay = additional_delay[m])
    params_3_VOC_vp[[m]][[i]] <- add_vp(params_3_VOC[[j]], delay = additional_delay[m])
    print(round(i*100/length(euro_inuse),0))
    save(params_3_vp, params_3_VOC_vp,
         file = "data/intermediate/params_3_vp_additional_add_supply_delay.rdata")
  }
}

w1 <- 1:13; w1 <- w1[!w1 %in% c(3,6,12)]
w2 <- 1:13; w2 <- w2[!w2 %in% c(5,11,13)]

for(i in w1){
  lapply(params_3_vp[[1]][[i]]$scenarios, "[[", "daily_vac_scenarios") %>% 
    map(arrange, date) %>% 
    map(mutate_at, vars(starts_with("Y", ignore.case = T)), cumsum) %>% 
    bind_rows(.id = "scenario") %>%
    select(-starts_with("supply")) %>% 
    pivot_longer(starts_with("Y")) %>% 
    separate(name, into = c("ag", "dose")) %>% 
    group_by(scenario, date, dose) %>% 
    summarise(value = sum(value), .groups = "drop") -> x
 
  
  lapply(params_3_vp[[1]][[i]]$scenarios, "[[", "daily_vac_scenarios") %>% 
    map(arrange, date) %>% 
    map(mutate_at, vars(starts_with("supply", ignore.case = T)), cumsum) %>% 
    bind_rows(.id = "scenario") %>% 
    dplyr::select(date, starts_with("supply"), -supply_cum) -> y
  
  y %>% filter(supply > 10000) %>% .[1,] %>% pull(date) -> t1
  y %>% filter(supply_2 > 10000) %>% .[1,] %>% pull(date) -> t2
  
  print(t2-t1)
  
  x %>% 
    left_join(y, by = "date") %>% 
    filter(date >= "2021-01-01") %>%
    mutate(pop =  vac_denom %>% filter(iso3c == euro_inuse[i]) %>% pull(tot),
           value = value/pop,
           supply = supply/pop,
           supply_2 = supply_2/pop) %>% 
    ggplot(., aes(x = date, y = value, group = dose, color = dose)) +
    geom_line() +
    geom_line(aes(y = supply), color = "black", linetype = 2) +
    geom_line(aes(y = supply_2), color = "black", linetype = 3) +
    facet_wrap(~scenario, nrow = 1) +
    labs(title = euro_inuse[i]) -> p
  
  ggsave(paste0("figs/supplemental/Delay_12wk_Sanity/",i,
                "_",euro_inuse[i],".png"),
         p, width = 20, height = 6)
}

for(i in w2){
  lapply(params_3_vp[[2]][[i]]$scenarios, "[[", "daily_vac_scenarios") %>% 
    map(arrange, date) %>% 
    map(mutate_at, vars(starts_with("Y", ignore.case = T)), cumsum) %>% 
    bind_rows(.id = "scenario") %>%
    select(-starts_with("supply")) %>% 
    pivot_longer(starts_with("Y")) %>% 
    separate(name, into = c("ag", "dose")) %>% 
    group_by(scenario, date, dose) %>% 
    summarise(value = sum(value), .groups = "drop") -> x
  
  
  lapply(params_3_vp[[2]][[i]]$scenarios, "[[", "daily_vac_scenarios") %>% 
    map(arrange, date) %>% 
    map(mutate_at, vars(starts_with("supply", ignore.case = T)), cumsum) %>% 
    bind_rows(.id = "scenario") %>% 
    dplyr::select(date, starts_with("supply"), -supply_cum) -> y
  
  y %>% filter(supply > 10000) %>% .[1,] %>% pull(date) -> t1
  y %>% filter(supply_2 > 10000) %>% .[1,] %>% pull(date) -> t2
  
  print(t2-t1)
  
  x %>% 
    left_join(y, by = "date") %>% 
    filter(date >= "2021-01-01") %>%
    mutate(pop =  vac_denom %>% filter(iso3c == euro_inuse[i]) %>% pull(tot),
           value = value/pop,
           supply = supply/pop,
           supply_2 = supply_2/pop) %>% 
    ggplot(., aes(x = date, y = value, group = dose, color = dose)) +
    geom_line() +
    geom_line(aes(y = supply), color = "black", linetype = 2) +
    geom_line(aes(y = supply_2), color = "black", linetype = 3) +
    facet_wrap(~scenario, nrow = 1) +
    labs(title = euro_inuse[i]) -> p
  
  ggsave(paste0("figs/supplemental/Delay_52wk_Sanity/",i,
                "_",euro_inuse[i],".png"),
         p, width = 20, height = 6)
}
 