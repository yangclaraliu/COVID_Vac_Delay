strategy_labels <- c(paste0("A",1:5),paste0("B", c(1,2)))

# 
# owid_vac %>% 
#   # filter(iso_code %in% euro_inuse) %>% 
#   filter(iso_code %in% euro_lmic,
#          !(iso_code %in% "TKM")) %>% 
#   left_join(vac_denom, by = c("iso_code" = "iso3c")) %>% 
#   mutate(cov = people_vaccinated/tot) %>% 
#   filter(!is.na(cov)) %>% 
#   ggplot(., aes(x = date, y = cov, group = iso_code)) +
#   geom_line() +
#   geom_point() +
#   geom_line(data = data.frame(date = ymd(c("2021-03-01", 
#                                            "2021-06-30",
#                                            "2021-12-31")),
#                               cov = c(0, 0.03, 0.2),
#                               iso_code = "mean"),
#             aes(x = date, y = cov),
#             color = "red") +
#   theme_cowplot() +
#   labs(x = "", y = "% of popullation vaccinated") +
#   # facet_wrap(~iso_code) +
#   geom_hline(yintercept = 0.2, linetype = 2) -> p
#   # filter(cov <= 0.1) %>% 
#   # group_by(iso_code) %>% mutate(rk_date = rank(desc(date))) %>% 
#   # filter(rk_date == 1) %>% 
#   # pull(date) %>% sort
#   # ggplot(., aes(x = date)) +
#   # geom_histogram()
# 
# ggsave("figs/supplemental/observed_coverage.png")

#### Figures ####
lapply(1:13, function(y) {
  lapply(c(1:7), function(x) params_3_vp[[y]]$scenarios[[x]]$daily_vac_scenarios) %>% 
    setNames(c(1:7)) %>% bind_rows(.id = "strategy") %>% 
    mutate(iso3c = model_selected_ur$iso3c[[y]],
           pop = sum(params_3_vp[[y]]$param$pop[[1]]$size))
}) %>% 
  bind_rows() %>% 
  dplyr::select(date, strategy,starts_with("Y"), iso3c, pop) %>% 
  group_by(strategy, iso3c,pop) %>% 
  group_split() %>% map(arrange, date) %>% 
  map(mutate_at, vars(starts_with("Y")), cumsum) %>% bind_rows() %>% 
  pivot_longer(cols = starts_with("Y")) %>%
  mutate(dose = if_else(grepl("d1",name), "Dose 1", "Dose 2")) %>% 
  #separate(name, into = c("ag","dose")) %>% 
  group_by(date, strategy, dose, iso3c, pop) %>% 
  summarise(value = sum(value)) %>% 
  mutate(cov = value/pop,
         strategy = factor(strategy,
                           levels = c(1:5, 7, 6),
                           labels = strategy_labels)) -> tmp

lapply(1, function(y) {
  lapply(c(1), function(x) params_3_vp[[y]]$scenarios[[x]]$daily_vac_scenarios) %>% 
    setNames(c(1)) %>% bind_rows(.id = "strategy") %>% 
    mutate(wb = model_selected_ur$iso3c[[y]],
           pop = sum(params_3_vp[[y]]$param$pop[[1]]$size))
}) %>% 
  .[[1]] %>% 
  select( starts_with("supply"), pop, date) %>% 
  mutate_at(vars(c("supply", "supply_2")), cumsum) %>% 
  mutate(supply = supply/pop,
         supply_2 = supply_2/pop) -> tmp_supply

# tmp %>% 
#   filter(strategy %in% strategy_labels[c(1,5:7)]) %>% 
#   left_join(tmp_supply, by = "date") %>%
#   filter(date >= "2021-02-01",
#          iso3c != "MNE",
#          iso3c != "TJK",
#          strategy == "A1") %>% 
#   ggplot(., aes(x = date)) +
#   geom_line(size = 1.2, alpha = 0.8, aes(y = cov, group = interaction(iso3c, dose), color = dose)) +
#   geom_line(aes(y = supply, linetype = "Dose 1"), size = 1.5) +
#   geom_line(aes(y = supply_2, linetype = "Dose 2"), size = 1.5) +
#   scale_linetype_manual(values = c(2,3)) +
#   facet_wrap(~iso3c) +
#   ggsci::scale_color_lancet() +
#   cowplot::theme_cowplot() +
#   labs(color = "Allocated:", y = "Vaccine Coverage by Dose", 
#        x = "Year", linetype = "Supply:") +
#   theme(strip.background = element_rect(fill = NA, color = "black"),
#         strip.text = element_text(size = 16),
#         legend.text = element_text(size =16),
#         legend.position = "top",
#         legend.key.width = unit(2,"cm")) 

tmp %>% 
  filter(strategy %in% strategy_labels[c(1,5:7)],
         iso3c %in% euro_inuse
         ) %>% 
  left_join(tmp_supply, by = "date") %>%
  filter(date >= "2021-02-01") %>% 
  ggplot(., aes(x = date)) +
  geom_line(size = 1.2, alpha = 0.8, aes(y = cov, group = interaction(iso3c, dose), color = dose)) +
  geom_line(aes(y = supply, linetype = "Dose 1"), size = 1.5) +
  geom_line(aes(y = supply_2, linetype = "Dose 2"), size = 1.5) +
  scale_linetype_manual(values = c(2,3)) +
  scale_x_continuous(breaks = ymd(c("2021-01-01", "2022-01-01")),
                     labels = c(2021, 2022)) +
  facet_wrap(~strategy) +
  ggsci::scale_color_lancet() +
  cowplot::theme_cowplot() +
  labs(color = "Allocated:", y = "Vaccine Coverage by Dose", 
       x = "Year", linetype = "Supply:") +
  theme(strip.background = element_rect(fill = NA, color = "black"),
        strip.text = element_text(size = 16),
        legend.text = element_text(size =16),
        legend.position = "top",
        legend.key.width = unit(2,"cm")) -> p1

ve %>% 
  dplyr::select(-ve_d) %>% 
  add_column(ve_onward = c(0.38, 0.5)) %>% 
  rownames_to_column(var = "dose") %>% 
  pivot_longer(cols = starts_with("ve"),
               names_to = "ve_type") %>%
  mutate(dose = factor(dose, levels = c(1,2), 
                       labels = c("1st","2nd")),
         ve_type = factor(ve_type,
                          levels = c("ve_onward",
                                     "ve_i_o",
                                     "ve_d_o",
                                     "ve_h",
                                     "ve_mort"),
                          labels = c("Onward\nTransmission",
                                     "Infection",
                                     "Disease",
                                     "Severe Case",
                                     "Mortality"))) %>% 
  ggplot(., aes(x = dose, y = value, group = ve_type, pch = ve_type)) +
  geom_line() + 
  geom_hline(yintercept = seq(0.35, 0.95, 0.1), linetype = 2) +
  scale_y_continuous(breaks = seq(0.35, 0.95, 0.1)) +
  geom_point(size = 3) + 
  labs(x = "Dose", y = "Vaccine Efficacy\nby Target", pch = "") +
  theme_cowplot() +
  theme(legend.position = "right") -> p2

data.frame(strategy = strategy_labels, 
           d1 = 0.65,
           d2 = c(0.65, 0.75, 0.75,
                  0.85, 0.85,
                  0.95, 0.95)) %>% 
  pivot_longer(cols = starts_with("d"),
               names_to = "dose") %>% 
  mutate(dose = factor(dose, levels = c("d1","d2"), 
                       labels = c("1st","2nd"))) %>% 
  ggplot(., aes(x = dose,  y = value, group = strategy)) +
  geom_line() +
  # geom_point() + 
  geom_hline(yintercept = seq(0.65, 0.95, 0.1), linetype = 2) +
  scale_y_continuous(breaks = seq(0.65, 0.95, 0.1)) +
  geom_label(x = "2nd", y = 0.65, label = "A1", size = 6) +
  geom_label(x = "2nd", y = 0.75, label = "A2/3", size = 6) +
  geom_label(x = "2nd", y = 0.85, label = "A4/5", size = 6) +
  geom_label(x = "2nd", y = 0.95, label = "B1/2", size = 6) +
  labs(x = "Dose", 
       y = "Infection- and Disease-reducing\nVaccine Efficacy\nby Dosing Interval") +
  theme_cowplot() -> p3


# owid_vac %>% 
#   rename(d1 = people_vaccinated,
#          d2 = people_fully_vaccinated) %>% 
#   dplyr::select(wb, date, d1, d2) %>% 
#   pivot_longer(cols = c("d1","d2"),
#                names_to = "dose") %>% 
#   filter(!is.na(value)) %>% 
#   right_join(tmp %>% ungroup %>% dplyr::select(wb,pop) %>% distinct(),
#              by = "wb") %>% 
#   mutate(p = value/pop) %>% 
#   ggplot(.,aes(x = date, y = p, group = interaction(wb, dose), color = dose)) +
#   geom_line(size = 1.2) +
#   scale_color_lancet(labels = c(("Dose1"),"Dose2")) +
#   theme_cowplot() +
#   theme(legend.position = "top") +
#   labs(y = "Vaccine Coverage", x = "2021", color = "") -> p7
# ggsave("figs/supplemental_fig3.png",
#        p7)

font_sizes <- theme(axis.text = element_text(size = 18),
                    axis.title = element_text(size = 20),
                    legend.text = element_text(size = 18),
                    legend.title = element_text(size = 18),
                    strip.text = element_text(size = 18))

p4 <- plot_grid(p2 + 
                  guides(pch = guide_legend(nrow = 6)) +
                  font_sizes,
                p3 +
                  font_sizes, 
                ncol = 1,
                labels = c("(a)",
                           "(b)"),
                align = "v",
                axis = "l",
                label_size = 20)

p5 <- plot_grid(p1 + 
                  guides(color = guide_legend(nrow = 2),
                         linetype = guide_legend(nrow = 2)) +
                  font_sizes,
                labels = c("(c)"),
                label_size = 20)

p6 <- plot_grid(p4,
                p5, ncol = 2, 
                rel_widths = c(1,1.5))

ggsave("figs/fig2.png", p6, height = 10, width = 18)
 
#### vaccine supply curves ####
# params_2_vp[[1]]$scenarios[[i]]$daily_vac_scenarios %>% 
#   mutate(t = as.numeric(t)) %>% 
#   dplyr::select(t, starts_with("supply")) %>% 
#   mutate_at(vars(starts_with("supply")), cumsum) %>% 
#   dplyr::select(-supply_cum) %>% 
#   pivot_longer(cols = starts_with("supply")) %>% 
#   mutate(name = factor(name, labels = c("Dose 1", "Dose 2", "Total"))) %>% 
#   ggplot(., aes(x = t, y = value, group = name, linetype = name)) +
#   geom_line(size = 1.2) +
#   theme_bw() +
#   labs(x = "", y = "", title = paste0("Total Vaccine Supply"), linetype = "") +
#   theme(title = element_text(size = 14)) 

#### effective dosing interval ####
lapply(1:13, function(x){
  lapply(params_3_vp[[x]]$scenarios, "[[", "pending") %>% 
    bind_rows(.id = "scenario") %>%
    pivot_longer(cols = starts_with("Y")) %>% 
    filter(elapse >= 28 & !is.na(t_dose2)) %>% 
    group_by(status, scenario, elapse) %>% summarise(value = sum(value)) %>% 
    mutate(scenario = factor(scenario)) %>% group_by(scenario) %>% group_split()
}) -> tmp

tmp_unroll <- list()
for(j in 1:length(params_3_vp)){
  sapply(1:length(tmp[[j]]), function(y){
    sapply(1:nrow(tmp[[j]][[y]]), 
           function(x) rep(tmp[[j]][[y]]$elapse[x], tmp[[j]][[y]]$value[x]/1000)) 
  }) %>% 
    map(c) %>% 
    map(unlist) -> tmp_unroll[[j]]
}

tmp_unroll %>% map(map, data.table) %>% 
  map(bind_rows, .id = "scenario") %>% bind_rows(.id = "country_index") %>% 
  filter(scenario %in%  c(6,7)) %>% 
  left_join(model_selected_ur %>%
              mutate(country_index = 1:18) %>%
              dplyr::select(country_index, country_name) %>%
              mutate(country_index = as.character(country_index)),
            by = "country_index") -> tmp

rm(tmp_unroll)

require(ggridges)

tmp %>% 
  mutate(scenario = factor(scenario, levels = c(7,6),
                           labels = c("B1","B2")),
         iso3c = countrycode(country_name, "country.name", "iso3c")) %>% 
  filter(iso3c %in% euro_inuse) %>% 
  group_by(scenario, country_name) %>% 
  summarise(v1 = mean(V1)/7,
            v2 = median(V1)/7) %>% 
  group_by(scenario) %>% 
  summarise(mean_LL = min(v1),
            mean_UL = max(v1),
            md_LL = min(v2),
            md_UL = max(v2))

tmp %>% 
  mutate(scenario = factor(scenario, levels = c(7,6),
                           labels = c("B1","B2")),
         iso3c = countrycode(country_name, "country.name", "iso3c")) %>% 
  filter(iso3c %in% euro_inuse) %>% 
  ggplot(., aes(x = V1, y = country_name, fill = scenario)) +
  geom_density_ridges(alpha = 0.75) +
  labs(y = "", color = "", x = "Dosing Interval (Days)") +
  theme_cowplot() +
  theme(legend.position = "top",
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)) +
  ggsci::scale_fill_lancet() +
  labs(fill = "") +
  scale_y_discrete(limits = rev) +
  geom_vline(xintercept = c(4*7, 20*7)) +
  scale_x_continuous(breaks = c(7*4, 7*20, 7*30, 40*7,50*7,60*7),
                     labels = c("A1", "A5", "30wk", "40wk", "50wk","60wk")) -> p 


# lapply(1:nrow(model_selected), function(x) {
#   tmp_unroll[[x]] %>% map(summary) %>% bind_rows(.id = "scenario")
# }) %>% bind_rows(.id = "country_index") -> tmp_unroll
# 
# model_selected %>% 
#   mutate(country_index = 1:13) %>% 
#   dplyr::select(country_index, country_name) %>% 
#   mutate(country_index = as.character(country_index)) %>% 
#   right_join(tmp_unroll, by = "country_index") %>% 
#   filter(scenario %in% c(6,7)) %>% 
#   mutate(scenario = factor(scenario, levels = c(6,7),
#                            labels = c("B1",
#                                       "B2"))) -> p_tmp
# dplyr::select(scenario, Mean) %>% group_by(scenario) %>% group_split() %>% map(pull, Mean) %>% map(range)
  
# ggplot(p_tmp, aes(x = country_name, y = Mean, color = scenario)) +
#   geom_point(size = 5, pch = 15) +
#   geom_segment(aes(x = country_name, xend = country_name, 
#                    y = `1st Qu.`, yend = `3rd Qu.`)) +
#   coord_flip() +
  # labs(x = "", color = "", y = "Dosing Interval (Days)") +
  # theme_cowplot() +
  # theme(legend.position = "none",
  #       axis.text = element_text(size = 14),
  #       axis.title = element_text(size = 14)) +
  # ggsci::scale_color_lancet() +
  # scale_x_discrete(limits = rev) +
  # geom_hline(yintercept = c(4*7, 20*7)) +
  # geom_label(x ="Georgia", y = 20*7, label = "A5", color = "black", size = 8) +
  # geom_label(x ="Georgia", y = 4*7, label = "A1", color = "black", size = 8) +
  # geom_label(x ="Georgia", y = 250, label = "B1", color = ggsci::pal_lancet()(2)[2], size = 8) +
  # geom_label(x ="Georgia", y = 300, label = "B2", color = ggsci::pal_lancet()(2)[1], size = 8) -> p
#     wow, this is really good code please use it


ggsave("figs/fig3_new.png", p, height = 8, width = 12)

#### baseline results ####
# update to `res_baseline_v4.rds` on 2022/02/03 upon R1
res_baseline <- read_rds(paste0("data/intermediate/res_baseline_v4.rds"))
outcome_by_strategy <- function(res_list, tab_name, outcome){
  if(outcome == "death") {l <- c("death", "death_o"); st = ""}
  if(outcome == "hosp") {l <- c("hosp", "hosp_i"); st = ""}
  if(outcome == "cases") {l <- c("cases"); st = ""}
  
  rl_name <- deparse(substitute(res_list))
  VOC <- grepl("VOC", tab_name)
  if(VOC) {
    st <- paste0(st, "w/ VOC")
    } else {
    st <- paste0(st, "w/o VOC")    
  }
  
  
  res_list[[tab_name]] %>% 
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
    ggplot(., aes(x = di, y = relative)) +
    # facet_wrap(~population, ncol = 5, scales = "free") +
    # geom_boxplot(outlier.shape = NA) +
    geom_line(aes(group = population), alpha = 0.5) +
    # geom_jitter() +
    geom_point(size = 2, alpha = 0.5) +
    scale_x_continuous(breaks = c(4,8,12,16,20,25,50),
                       labels = c("A1\n4w", "A2\n8w", "A3\n12w", "A4\n16w", 
                                  "A5\n20w","B1\n22-34w","B2\n43-56w")) +
    scale_y_continuous(breaks = c(0, 0.15, 0.3), labels = c("0%","15%","30%"),
                       limits = c(-0.1, 0.5)
                       ) +
    # geom_dl(aes(label = wb), 
    #         method = list(dl.combine("first.points", "last.points")))
    geom_vline(xintercept = 25, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2 ) +
    theme_cowplot() +
    theme(axis.title = element_text(size = 16),
          # title = element_text(size = 18),
          plot.subtitle = element_text(size = 16)) +
    labs(x = "Dosing Interval Strategy",
         subtitle = st,
        #  title = "Relative Differences in Cumulative Outcome Compared to B1",
         y = "%Differences Compared to B1") -> p

  return(p)  
}

p1 <- outcome_by_strategy(res_baseline, "res_3", "death")
p2 <- outcome_by_strategy(res_baseline, "res_3_VOC", "death")
title <- ggdraw() + draw_label("Relative Differences in Cumulative COVID-19 Mortality Compared to B1\nwv = 360 days",
                               fontface = "bold", x = 0, hjust = 0, size = 18)
p3 <- plot_grid(title, p2,p1,ncol = 1, rel_heights = c(0.1, 1,1), labels = c("","(A)","(B)"),
               label_fontface = "bold", label_size = 20)

res_baseline[["res_3"]] %>% 
  filter(compartment %in% c("death", "death_o")) %>%
  group_by(scenario, population, compartment) %>% 
  summarise(value = sum(value)) %>% 
  pivot_wider(names_from = scenario, values_from = value) %>% 
  mutate(B1 = `7`) %>% 
  pivot_longer(cols = as.character(c(1:7)),
               names_to = "scenario",
               values_to = "value") %>% 
  mutate(relative = value/B1) %>% 
  # filter(scenario == 1) %>% pull(relative) %>% summary
  filter(scenario %in% c(3:5)) %>% group_by(scenario) %>% group_split() %>% map(pull, relative) %>% map(summary)

# ggsave("figs/res_3_baseline_R1.png", p, width = 10, height = 14)


#### shorter waning ####
# update to `res_sw_v4.rds` on 2022/02/03 upon R1
res_sw <- read_rds("data/intermediate/res_sw_v4.rds")
p4 <- outcome_by_strategy(res_sw, "res_3_sw", "death")
p5 <- outcome_by_strategy(res_sw, "res_3_VOC_sw", "death")
title <- ggdraw() + draw_label("Relative Differences in Cumulative COVID-19 Mortality Compared to B1\nwv = 120 days",
                               fontface = "bold", x = 0, hjust = 0, size = 18)
p6 <- plot_grid(title, p5,p4,ncol = 1, rel_heights = c(0.1, 1,1), labels = c("","(A)","(B)"),
               label_fontface = "bold", label_size = 20)

plot_grid(p2 + labs(subtitle = "w/ VOC\n1st dose wane at 360 days"), 
          p1 + labs(subtitle = "w/o VOC\n1st dose wane at 360 days"),  
          p5 + labs(subtitle = "w/ VOC\n1st dose wane at 120 days"), 
          p4 + labs(subtitle = "w/o VOC\n1st dose wane at 120 days"),
          labels = c("(a)","(b)","(c)","(d)"),
          nrow = 2, byrow =T) -> p7

# update to `fig4_v4.png` on 2022/02/03 upon R1
ggsave("figs/fig4_v4.png", p7, width = 12, height = 9)

res_sw[["res_3_sw"]] %>% 
  filter(compartment %in% c("death", "death_o")) %>%
  group_by(scenario, population, compartment) %>% 
  summarise(value = sum(value)) %>% 
  pivot_wider(names_from = scenario, values_from = value) %>% 
  mutate(B1 = `7`) %>% 
  pivot_longer(cols = as.character(c(1:7)),
               names_to = "scenario",
               values_to = "value") %>% 
  mutate(relative = value/B1) %>% 
  group_by(population) %>% group_split() %>% map(mutate, rk = rank(relative)) %>% 
  bind_rows() %>% group_by(scenario, rk) %>% tally %>% filter(rk == 1) %>% group_by(scenario) %>% group_split()
  

##### longer interval ~ higher efficacy
ve_new <- read_rds((paste0(path_dropbox,"intermediate/ve_new.rds")))
ve_new[ve_i == 0.65 & ve_d == 0 & v2e_d == 0] %>% 
  slice(rep(1:n(), each = 2)) %>% 
  select(ve_set) %>% mutate(scenarios = c(NA, 1:5,7,6)) %>% 
  filter(!is.na(scenarios)) -> profiles

# SA_VE_120 <- read_rds(paste0(path_dropbox,"intermediate/SA_VE_120.rds"))
# SA_VE_VOC_120 <- read_rds(paste0(path_dropbox,"intermediate/SA_VE_120_VOC.rds"))
# 
# SA_VE_360 <- read_rds(paste0(path_dropbox,"intermediate/SA_VE_360.rds"))
# SA_VE_VOC_360 <- read_rds(paste0(path_dropbox,"intermediate/SA_VE_360_VOC.rds"))

SA_VE_360 %>% 
  bind_rows() %>% 
  filter(!grepl("voc", compartment),
         ve_set %in% profiles$ve_set,
         compartment == "death_o"
         ) %>% 
  group_by(ve_set, scenarios, population, compartment) %>% 
  summarise(value = sum(V1)) %>% 
  mutate(scenarios = factor(scenarios, levels = c(1:5,7,6))) %>% 
  filter((ve_set == profiles$ve_set[1] & scenarios == profiles$scenarios[1]) |
           (ve_set == profiles$ve_set[2] & scenarios == profiles$scenarios[2]) |
           (ve_set == profiles$ve_set[3] & scenarios == profiles$scenarios[3]) |
           (ve_set == profiles$ve_set[4] & scenarios == profiles$scenarios[4]) |
           (ve_set == profiles$ve_set[5] & scenarios == profiles$scenarios[5]) |
           (ve_set == profiles$ve_set[6] & scenarios == profiles$scenarios[6]) |
           (ve_set == profiles$ve_set[7] & scenarios == profiles$scenarios[7])) %>% 
  ungroup() %>% dplyr::select(-ve_set) %>% 
  pivot_wider(names_from = scenarios, values_from = value) %>% 
  mutate(B1 = `7`) %>% 
  pivot_longer(cols = as.character(c(1:7)),
               names_to = "scenarios",
               values_to = "value") %>% 
  mutate(relative = value/B1,
         scenarios = factor(scenarios,
                            levels = c(1:5,7,6),
                            labels = strategy_labels)) -> tmp

res_baseline_v3[["res_3"]] %>% 
  filter(compartment %in% c("death", "death_o")) %>%
  group_by(scenario, population, compartment) %>% 
  summarise(value = sum(value)) %>% 
  pivot_wider(names_from = scenario, values_from = value) %>% 
  mutate(B1 = `7`) %>% 
  pivot_longer(cols = as.character(c(1:7)),
               names_to = "scenarios",
               values_to = "value") %>% 
  mutate(relative = value/B1,
         scenarios = factor(scenarios, levels = c(1:5,7,6), labels = strategy_labels)) %>% 
  dplyr::select(colnames(tmp)) %>% 
  mutate(status = "original") %>% 
  bind_rows(tmp %>% mutate(status = "updated")) %>% 
  mutate(relative = relative - 1) %>% 
  ggplot(., aes(x = scenarios, 
                y = relative,
                # group = population,
                color = status)) +
  geom_boxplot(outlier.shape = NA)+
  # geom_line(alpha = 0.5, size = 2) +
  # geom_jitter() +
  # geom_point(size = 2, alpha = 0.5) +
  ggsci::scale_color_lancet(labels = c("Baseline VE",
                                       "Dosing-interval Sensitive VE")) +
  geom_vline(xintercept = 6, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2 ) +
  scale_y_continuous(breaks = c(0,0.15,0.3),
                     limits = c(0,0.5),
                     labels = c("0%","15%","30%")) +
  theme_cowplot() +
  theme(axis.title = element_text(size = 16),
        legend.position = "none",
        # title = element_text(size = 18),
        plot.subtitle = element_text(size = 16)) +
  labs(x = "Dosing Interval Strategy",
       subtitle = "w/o VOC\n1st Dose wane at 360 days",
       color = "",
       #  title = "Relative Differences in Cumulative Outcome Compared to B1",
       y = "%Differences Compared to B1") -> p2

SA_VE_360_VOC %>% 
  bind_rows() %>% 
  filter(grepl("death", compartment)) %>% 
  filter(!(compartment == "hosp_voc_i" & VOC == F)) %>% 
  filter(!(compartment == "hosp_i" & VOC == T)) %>% 
  filter(!(compartment == "death_voc_o" & VOC == F)) %>% 
  filter(!(compartment == "death_o" & VOC == T)) %>% 
  group_by(ve_set, scenarios, population) %>% 
  summarise(value = sum(V1)) %>% 
  mutate(scenarios = factor(scenarios, levels = c(1:5,7,6))) %>% 
  filter((ve_set == profiles$ve_set[1] & scenarios == profiles$scenarios[1]) |
           (ve_set == profiles$ve_set[2] & scenarios == profiles$scenarios[2]) |
           (ve_set == profiles$ve_set[3] & scenarios == profiles$scenarios[3]) |
           (ve_set == profiles$ve_set[4] & scenarios == profiles$scenarios[4]) |
           (ve_set == profiles$ve_set[5] & scenarios == profiles$scenarios[5]) |
           (ve_set == profiles$ve_set[6] & scenarios == profiles$scenarios[6]) |
           (ve_set == profiles$ve_set[7] & scenarios == profiles$scenarios[7])) %>% 
  ungroup() %>% dplyr::select(-ve_set) %>% 
  pivot_wider(names_from = scenarios, values_from = value) %>% 
  mutate(B1 = `7`) %>% 
  pivot_longer(cols = as.character(c(1:7)),
               names_to = "scenarios",
               values_to = "value") %>% 
  mutate(relative = value/B1,
         scenarios = factor(scenarios,
                            levels = c(1:5,7,6),
                            labels = strategy_labels)) -> tmp

res_baseline_v3[["res_3_VOC"]] %>% 
  filter(compartment %in% c("death", "death_o"),
         population %in% model_selected_ur$country_name) %>%
  group_by(scenario, population, compartment) %>% 
  summarise(value = sum(value)) %>% 
  pivot_wider(names_from = scenario, values_from = value) %>% 
  mutate(B1 = `7`) %>% 
  pivot_longer(cols = as.character(c(1:7)),
               names_to = "scenarios",
               values_to = "value") %>% 
  mutate(relative = value/B1,
         scenarios = factor(scenarios, levels = c(1:5,7,6), labels = strategy_labels)) %>% 
  dplyr::select(colnames(tmp)) %>% 
  mutate(status = "original") %>% 
  # group_by(population, status) %>% group_split() %>% 
  # map(mutate, rk = rank(relative)) %>% map(filter, rk == 1) %>% 
  # bind_rows() %>% 
  bind_rows(tmp %>% mutate(status = "updated")) %>% 
  mutate(relative = relative - 1) %>% 
  ggplot(., aes(x = scenarios, y = relative, color = status)) +
  geom_boxplot(outlier.shape = NA)+
  # geom_line(aes(group = interaction(status, population)), alpha = 0.1) +
  # geom_jitter() +
  # geom_point(size = 2, alpha = 0.5) +
  ggsci::scale_color_lancet(labels = c("Baseline VE",
                                       "Dosing-interval Sensitive VE")) +
  geom_vline(xintercept = 6, linetype = 2) +
  scale_y_continuous(breaks = seq(0,0.45,0.15),
                     limits = c(-0.05, 0.5),
                     labels = c("0%","15%","30%","45%")) +
  geom_hline(yintercept = 0, linetype = 2 ) +
  theme_cowplot() +
  theme(axis.title = element_text(size = 16),
        legend.position = "top",
        legend.text = element_text(size = 16),
        # title = element_text(size = 18),
        plot.subtitle = element_text(size = 16)) +
  labs(x = "Dosing Interval Strategy",
       subtitle = "w/ VOC\n1st Dose wane at 360 days",
       color = "",
       #  title = "Relative Differences in Cumulative Outcome Compared to B1",
       y = "%Differences Compared to B1")  -> p1

p <- plot_grid(p1, p2)
p <- plot_grid(title, p1,p2,ncol = 1, rel_heights = c(0.1, 1,1), labels = c("","(A)","(B)"),
               label_fontface = "bold", label_size = 20)
title <- ggdraw() + draw_label("Relative Differences in Cumulative COVID-19 Mortality Compared to B1\nwv = 360 days",
                               fontface = "bold", x = 0, hjust = 0, size = 18)
ggsave("figs/supplemental/res_3_dynamic_VE.png", p, width = 10, height = 14)

#### SA around VE ####
# SA_VE_360 %>% bind_rows() %>% filter(grepl("death_o", compartment)) %>% 
#   group_by(ve_set, scenarios, population) %>% 
#   summarise(value = sum(V1)) %>% mutate(status = "w/o VOC") %>%  
#   bind_rows(
SA_VE_360_VOC  %>% 
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
  left_join(ve_new, by = "ve_set") %>%   
  # dplyr::select(-B1, -value) %>% 
  # pivot_wider(names_from = status, values_from = relative) %>% 
  # filter(`w/o VOC` != `w/ VOC`) %>% View()
  ggplot(., aes(x = scenarios, y = relative-1, group = ve_set)) +
  geom_line(alpha = 0.01) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_hline(yintercept = -0.01, linetype = 2, color = "orange") +
  geom_hline(yintercept = -0.05, linetype = 2, color = "red") +
  geom_hline(yintercept = 0.01, linetype = 2, color = "orange") +
  geom_hline(yintercept = 0.05, linetype = 2, color = "red") +
  facet_wrap(~population, scales = "free", ncol = 3) +
  theme_cowplot() +
  theme(axis.title = element_text(size = 16),
        legend.position = "none",
        strip.background = element_rect(fill = NA, colour = "black"),
        # title = element_text(size = 18),
        plot.subtitle = element_text(size = 16)) +
  labs(x = "Dosing Interval Strategy",
       subtitle = "w/ VOC; wv = 360",
       color = "",
       title = "Sensitivity analysis around first and second doses\ninfection-blocking and disease-reducing VEs",
       y = "Proportional Differences Compared to B1\nRCM") -> p

ggsave("figs/supplemental/SA_VE_VOC_360.png", p, width = 10, height = 15)
  
  
