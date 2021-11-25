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
         wb = countrycode(population, "country.name", "wb")) %>% 
  group_by(wb) %>% 
  mutate(rk = rank(relative)) %>% 
  dplyr::select(wb, relative, rk) %>% 
  pivot_wider(names_from = rk, values_from = relative) %>% 
  mutate(diff = abs(`1`-`2`)) %>% pull(diff) %>% range

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

