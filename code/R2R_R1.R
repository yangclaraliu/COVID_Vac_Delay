#### Response to reviewer, Round 1 #####
# 2022/02/02

# 
low_supply <- read_rds("data/R1_low_supply.rds")

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
