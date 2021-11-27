res_3 %>% 
  bind_rows() %>% 
  filter(compartment == "S") %>% 
  group_by(population, t) %>% 
  summarise(value = sum(value)) %>% 
  left_join(model_selected_ur, by = c("population" = "country_name")) -> tmp

tmp %>% 
  mutate(date = t.x + t.y + ymd("2019-12-01")) %>% 
  filter(t.x == 0 | date == "2021-03-01") %>% 
  group_by(iso3c) %>% group_split() %>% map(add_column, cat = c("start","end")) %>% 
  bind_rows() %>% 
  dplyr::select(population, iso3c, value, cat) %>% 
  pivot_wider(names_from = cat, values_from = value) %>% 
  mutate(diff = (start-end)/start) %>% 
  filter(iso3c %in% euro_inuse) %>% 
  pull(diff) %>% range

#### inline statistics ####
owid_vac %>% 
  filter(wb %in% euro_inuse) %>% 
  left_join(vac_denom, by = "wb") %>% 
  mutate(p = people_fully_vaccinated/tot) %>% 
  ggplot(., aes(x = date, y = p, group = wb)) +
  geom_line() +
  facet_wrap(~location) +
  geom_hline(yintercept = c(0.2, 0.3))

# 
l <- c("death","death_o")
res_baseline[["res_3_VOC"]] %>% 
  filter(compartment %in% l) %>%
  group_by(scenario, population, compartment) %>% 
  summarise(value = sum(value)) %>% 
  pivot_wider(names_from = scenario, values_from = value) %>% 
  mutate(B1 = `7`,
         A1 = `1`) %>% 
  pivot_longer(cols = as.character(c(1:7)),
               names_to = "scenario",
               values_to = "value") %>% 
  mutate(relative = value/B1 - 1,
         scenario = factor(scenario, 
                           levels = c(1:5,7,6), 
                           labels = strategy_labels),
         di = case_when(scenario == "A1" ~ 4,
                        scenario == "A2" ~ 8,
                        scenario == "A3" ~ 12,
                        scenario == "A4" ~ 16,
                        scenario == "A5" ~ 20,
                        scenario == "B1" ~ 25,
                        scenario == "B2" ~ 50),
         wb = countrycode(population, "country.name", "wb")) -> tmp 
  group_by(wb) %>% 
  mutate(rk = rank(relative)) %>% 
  dplyr::select(wb, relative, rk) %>% 
  pivot_wider(names_from = rk, values_from = relative) %>% 
  mutate(diff = abs(`1`-`2`)) %>% pull(diff) %>% range

tmp  %>%
  filter(wb %in% euro_inuse) %>% 
  group_by(population) %>% 
  filter(scenario == "B2")
  mutate(rk = rank(relative)) %>% 
  filter(rk == 1)
  
res_sw$res_3_sw %>% 
    mutate(iso3c = countrycode(population, "country.name", "iso3c")) %>% 
    filter(compartment %in% l,
           iso3c %in% euro_inuse) %>%
    group_by(scenario, population, compartment) %>% 
    summarise(value = sum(value)) %>% 
    pivot_wider(names_from = scenario, values_from = value) %>% 
    mutate(B1 = `7`,
           A1 = `1`) %>% 
    pivot_longer(cols = as.character(c(1:7)),
                 names_to = "scenario",
                 values_to = "value") %>% 
    mutate(relative = value/B1 - 1,
           scenario = factor(scenario, 
                             levels = c(1:5,7,6), 
                             labels = strategy_labels),
           di = case_when(scenario == "A1" ~ 4,
                          scenario == "A2" ~ 8,
                          scenario == "A3" ~ 12,
                          scenario == "A4" ~ 16,
                          scenario == "A5" ~ 20,
                          scenario == "B1" ~ 25,
                          scenario == "B2" ~ 50),
           wb = countrycode(population, "country.name", "wb")) %>% 
  group_by(population) %>% mutate(rk = rank(relative)) %>% filter(rk == 1) %>% pull(scenario) %>% table


cm_populations %>% 
  filter(location_type == 4) %>% 
  mutate(wb = countrycode(name, "country.name", "wb")) -> pop_structure

pop_structure %>% 
  filter(wb %in% euro_inuse) %>% 
  separate(age, into = c("LL","UL")) %>% 
  mutate(ag = case_when(LL < 20 ~ "children",
                        LL >= 20 & LL < 60 ~ "adults",
                        LL >= 60 ~ "older adults")) %>% 
  group_by(wb, ag) %>%
  mutate(tot = sum(m + f)) %>% 
  summarise(tot = sum(tot)) %>% 
  pivot_wider(names_from = ag,
              values_from = tot) %>% 
  mutate(tot = adults + children + `older adults`,
         p_a = adults/tot) %>% 
  arrange(p_a)

non_S_3_debug %>% filter(wb %in% euro_inuse) %>%  arrange(p)


# Russia
SA_VE_VOC_360 %>% 
  bind_rows() %>% 
  filter(grepl("death", compartment)) %>% 
  filter(!(compartment == "hosp_voc_i" & VOC == F)) %>% 
  filter(!(compartment == "hosp_i" & VOC == T)) %>% 
  filter(!(compartment == "death_voc_o" & VOC == F)) %>% 
  filter(!(compartment == "death_o" & VOC == T)) %>% 
  group_by(ve_set, scenarios, population) %>% 
  summarise(value = sum(V1)) %>% 
  mutate(status = "w/ VOC") %>% #) %>% 
  pivot_wider(names_from = scenarios, values_from = value) %>% 
  mutate(B1 = `7`) %>% 
  pivot_longer(cols = as.character(c(1:7)),
               names_to = "scenarios",
               values_to = "value") %>%
  # filter(scenarios %in% c(6,7)) %>% 
  mutate(relative = value/B1,
         scenarios = factor(scenarios,
                            levels = c(1:5,7,6),
                            labels = strategy_labels)) %>% 
  filter(population == "Russian Federation",
         scenarios %in% c("B1","B2")) %>% 
  dplyr::select(-value, -B1) %>% 
  pivot_wider(names_from = "scenarios",
              values_from = relative) %>% 
  left_join(ve_new, by = "ve_set") %>% 
  mutate(res = if_else(B1 <= B2, "B1", "B2")) -> tmp

tmp %>% 
  pivot_longer(cols = c("ve_i", "v2e_i","ve_d_o","v2e_d_o")) %>% 
  ggplot(., aes(x = name, y = value, group = ve_set)) +
  geom_line() +
  facet_wrap(~res)


tmp %>% 
  mutate(diff = B1 - B2) %>% 
  # mutate(res = factor(res)) %>% 
  ggplot(., aes(x = ve_i, y = v2e_i,  fill = B1 - B2, color = res)) +
  geom_tile() +
  facet_grid(ve_d_o ~ v2e_d_o) +
  scale_fill_distiller(palette = "Spectral")

all <- list()
for(i in which(model_selected_ur$iso3c %in% euro_inuse)){
  n <- model_selected_ur$country_name[i]
  cm_matrices[[n]][[1]] + cm_matrices[[n]][[2]] + cm_matrices[[n]][[3]] + cm_matrices[[n]][[4]]  -> all[[n]]
}


p_list <- list()
for(i in 1:length(all)){
  (all$`Bosnia and Herzegovina` - all[[i]]) %>%
    set_colnames(ag_labels) %>% 
    set_rownames(ag_labels) %>% 
    as.data.frame %>% 
    rownames_to_column(var = "V1") %>% 
    pivot_longer(cols = ag_labels) %>% 
    mutate(V1 = factor(V1, levels = ag_labels, labels = 1:16) %>% as.numeric,
           name = factor(name, levels = ag_labels, labels = 1:16) %>% as.numeric) %>% 
    ggplot(., aes(x = V1, y = name, fill = value)) +
    geom_tile() -> p_list[[i]]
}
do.call("grid.arrange", c(p_list, ncol=4))

