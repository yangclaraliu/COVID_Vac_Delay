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

#### extract non_S ####
lapply(1:13, function(x){
  params_3_vp[[x]]$res[[1]] %>% 
    cm_simulate%>% 
    .[["dynamics"]] %>% 
    # bind_rows(.id = "scenarios") %>% 
    filter(compartment == "S") %>% 
    group_by(t, population) %>% 
    summarise(value = sum(value), .groups = "drop") %>% 
    left_join(model_selected_ur[,c("country_name", "t")],
              by = c("population" = "country_name")) %>% 
    mutate(date = ymd("2019-12-01") + t.x + t.y) %>% 
    filter(t.x == 0 | date == "2021-03-01") %>% 
    bind_cols(tag = c("start","end")) %>% 
    dplyr::select(population, tag, value) %>% 
    pivot_wider(names_from = tag, values_from = value) %>% 
    mutate(non_S = 1 - end/start)
  
}) %>% bind_rows() -> non_S

non_S %>% arrange(non_S)


all_euro <- c("Albania, Andorra, Armenia, Austria, Azerbaijan, Belarus, Belgium, Bosnia and Herzegovina, Bulgaria, Croatia, Cyprus, Czech Republic, Denmark, Estonia, Finland, France, Georgia, Germany, Greece, Hungary, Iceland, Ireland, Israel, Italy, Kazakhstan, Kyrgyzstan, Latvia, Lithuania, Luxembourg, Malta, Monaco, Montenegro, Netherlands, Norway, Poland, Portugal, Republic of Moldova, Romania, Russian Federation, San Marino, Serbia, Slovakia, Slovenia, Spain, Sweden, Switzerland, Tajikistan, The former Yugoslav Republic of Macedonia, Turkey, Turkmenistan, Ukraine, United Kingdom, Uzbekistan") %>% 
  strsplit(., ", ") %>% unlist %>% 
  data.frame(cn = .) %>% 
  mutate(iso3c = countrycode(cn, "country.name", "iso3c"))

cm_populations %>% 
  mutate(iso3c = countrycode(name, "country.name", "iso3c"),
         continent = countrycode(name, "country.name", "continent")) -> pop_struct


pop_struct %>% 
  filter(continent == "Europe",
         location_type == 4) %>% 
  # filter(iso3c %in% euro_all) %>% 
  separate(age, sep = "-", into = c("LL", "UL")) %>% 
  mutate(LL = parse_number(LL),
         ag = case_when(LL < 20 ~ "c",
                        LL >= 20 & LL < 60 ~ "a",
                        LL >= 60 ~ "oa"),
         MIC = if_else(iso3c %in% euro_lmic,
                       T, F)) %>% 
  group_by(name, iso3c) %>% 
  mutate(tot = sum(f) + sum(m)) %>% 
  group_by(name, iso3c, ag, tot, MIC) %>% 
  summarise(tot_group = sum(f) + sum(m)) %>% 
  mutate(p = tot_group/tot) %>% 
  group_by(ag, MIC) %>% summarise(p_mu = mean(p),
                                  LL = min(p),
                                  UL = max(p)) %>% 
  pivot_longer(cols = c(p_mu, LL, UL)) %>% 
  pivot_wider(names_from = MIC, values_from = value) %>% View()
  filter(name == "p_mu")

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

load("data/intermediate/params_3_vp_additional_add_supply_delay.rdata")
params_3_vp_add_delays <- params_3_vp
params_3_VOC_vp_add_delays <- params_3_VOC_vp

# params_3_vp_add_delays[[2]][[2]]$res %>% 
#   map(cm_simulate) %>% 
#   map(~.[["dynamics"]]) %>% 
#   map(filter, compartment == "death_o") %>% 
#   map(group_by, t) %>% 
#   map(summarise, value = sum(value)) %>% 
#   bind_rows(.id = "scenario") %>% 
#   ggplot(., aes(x = t, y = value, group = scenario, color = scenario)) +
#   geom_line()

# get_res_SDelay <- function(m, VOC = F){
# 
#   if(m == 1 & VOC == F) tmp <- params_3_vp_add_delays[[1]]
#   if(m == 2 & VOC == F) tmp <- params_3_vp_add_delays[[2]]
#   
#   if(m == 1 & VOC == T) tmp <- params_3_VOC_vp_add_delays[[1]]
#   if(m == 2 & VOC == T) tmp <- params_3_VOC_vp_add_delays[[2]]
# 
#   res <- list()
# 
#   for(i in 1:length(tmp)){
#     tmp[[i]]$res %>%
#       map(cm_simulate) %>%
#       map(~.[["dynamics"]]) %>%
#       bind_rows(.id = "scenario") %>%
#       filter(grepl("death|hosp|case", compartment)) %>% 
#       group_by(scenario, t, compartment) %>%
#       summarise(value = sum(value), .groups = "drop") %>%
#       mutate(iso3c = euro_inuse[i]) %>%
#       left_join(model_selected_ur %>%
#                   dplyr::select(-date, -r, -ur, -to_replace),
#                 by = "iso3c") %>%
#       mutate(date = ymd("2019-12-01") + t.x + t.y) %>%
#       dplyr::select(-t.x, -t.y)  -> res[[i]]
#   }
#   
#   if(VOC == F){
#     res %<>% 
#       map(filter, (!compartment %in% c("hosp_p", "cases_reported"))) %>% 
#       map(group_by, scenario, compartment, iso3c, country_name) %>% 
#       map(summarise, value = sum(value)) %>% 
#       bind_rows()
#   }
# 
#   if(VOC == T){
#     res %<>% 
#       map(mutate, VOC = if_else(date > "2021-04-15", T, F)) %>% 
#       map(filter, !(compartment %in% c("hosp_p", "hosp_i", "death_o") & VOC == T)) %>% 
#       map(filter, !(compartment %in% c("hosp_voc_p", "hosp_voc_i", "death_voc_o") & VOC == F)) %>% 
#       map(filter, compartment != "cases_reported") %>% 
#       map(filter, !grepl("hosp_p|hosp_voc_p", compartment)) %>% 
#       map(mutate, tag = substr(compartment, 1,4)) %>% 
#       map(group_by, scenario, tag, iso3c, country_name) %>% 
#       map(summarise, value = sum(value), .groups = "drop") %>% 
#       bind_rows()
#   }
# 
# 
#   return(res)
# 
#   }
# 
# 
# SDelay_outcome <- list()
# SDelay_outcome[["12_wk"]] <- get_res_SDelay(m =1, VOC = F)
# SDelay_outcome[["12_wk_VOC"]] <- get_res_SDelay(m =1, VOC = T)
# SDelay_outcome[["52_wk"]] <- get_res_SDelay(m = 2, VOC = F)
# SDelay_outcome[["52_wk_VOC"]] <- get_res_SDelay(m = 2, VOC = T)
# write_rds(SDelay_outcome, "data/R1_SDelay_outcome.rds")



# 
i = 1
res_baseline_v4[[i]] %>% filter(compartment %in% c("death_o", "death")) %>%
  group_by(scenario, population) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  mutate(scenario = factor(scenario, levels = c(1:5, 7, 6), labels = strategy_labels)) %>%
  group_by(population) %>% mutate(rk = rank(value)) %>%
  left_join(res_baseline_v4[[i]] %>% filter(compartment %in% c("death_o", "death")) %>%
              group_by(scenario, population) %>%
              summarise(value = sum(value), .groups = "drop") %>%
              filter(scenario == 7) %>%
              rename(B1 = value) %>% dplyr::select(-scenario),
            by = "population") %>%
  mutate(relative = value/B1 - 1) %>%
  # filter(scenario %in% c("B1", "B2")) %>% select(-B1, -value, -rk) %>% pivot_wider(names_from = scenario, values_from = relative) %>% arrange(B2)
  #filter(rk == 1)
  filter(scenario == "A1") %>% pull(relative) %>% mean



# 
# 
# 
SDelay_outcome <- read_rds("data/R1_SDelay_outcome.rds")

SDelay_outcome %>%
  bind_rows(.id = "set") %>%
  mutate(tag = if_else(is.na(tag), as.character(compartment), tag),
         tag = substr(tag, 1,4)) %>%
  separate(set, into = c("SDelay", "dep","VOC")) %>%
  dplyr::select(-dep) %>%
  mutate(VOC = if_else(is.na(VOC), F, T)) -> p_tab

# 
# 
p_tab %>%
  filter(VOC == T, tag == "deat", SDelay == 12) %>%
  group_by(SDelay, VOC, iso3c, country_name) %>% mutate(rk = rank(value)) %>%
  filter(rk == 1) %>% arrange(iso3c) %>%
  left_join(p_tab
            %>%
              filter(VOC == T, tag == "deat", scenario == 7, SDelay == 52) %>%
              dplyr::select(-scenario) %>%
              rename(B1 = value) %>%
              ungroup %>%
              select(iso3c, B1),
            by = "iso3c") %>%
  mutate(relative = value/B1 - 1)

# 
# 
# p_tab %>% 
#   filter(tag == "deat") %>% 
#   left_join(p_tab %>% 
#               filter(tag == "deat", scenario == 7) %>% 
#               group_by(SDelay, VOC, iso3c) %>% 
#               rename(B1 = scenario,
#                      B1_value = value) %>% 
#               dplyr::select(SDelay, VOC, tag, iso3c, B1_value, B1),
#             by = c("SDelay", "VOC", "tag", "iso3c")) %>% 
#   group_by(SDelay, VOC, scenario, iso3c) %>% 
#   mutate(scenario = factor(scenario,
#                            levels = c(1:5,7,6),
#                            labels = strategy_labels),
#          relative = value/B1_value) -> p_tab
# 
# p_tab %>% 
#   mutate(scenario_t = case_when(scenario == "A1" ~ 4,
#                                 scenario == "A2" ~ 8,
#                                 scenario == "A3" ~ 12,
#                                 scenario == "A4" ~ 16,
#                                 scenario == "A5" ~ 20,
#                                 scenario == "B1" & SDelay == 12 ~ 32,
#                                 scenario == "B2" & SDelay == 12 ~ 49,
#                                 scenario == "B1" & SDelay == 52 ~ 31,
#                                 scenario == "B2" & SDelay == 52 ~ 66),
#          SDelay = factor(SDelay, levels = c(12, 52), 
#                          labels = c("Supply Delay = 12 wks",
#                                     "Supply Delay = 52 wks")),
#          VOC = factor(VOC, levels = c(T,F), labels = c("w/ VOC", "w/o VOC")),
#          relative = relative - 1) %>% 
#   group_by(SDelay, VOC, iso3c, country_name) %>% 
#   mutate(rk = rank(relative)) %>% 
#   filter(rk == 1) %>% 
#   group_by(SDelay, VOC, scenario) %>% tally
# 
# 
# 
# p_tab %>% 
#   mutate(scenario_t = case_when(scenario == "A1" ~ 4,
#                                 scenario == "A2" ~ 8,
#                                 scenario == "A3" ~ 12,
#                                 scenario == "A4" ~ 16,
#                                 scenario == "A5" ~ 20,
#                                 scenario == "B1" & SDelay == 12 ~ 32,
#                                 scenario == "B2" & SDelay == 12 ~ 49,
#                                 scenario == "B1" & SDelay == 52 ~ 31,
#                                 scenario == "B2" & SDelay == 52 ~ 66),
#          SDelay = factor(SDelay, levels = c(12, 52), 
#                          labels = c("Supply Delay = 12 wks",
#                                     "Supply Delay = 52 wks")),
#          VOC = factor(VOC, levels = c(T,F), labels = c("w/ VOC", "w/o VOC")),
#          relative = relative - 1) %>% 
#   ggplot(., aes(x = scenario_t, y = relative, group = iso3c, color = iso3c)) +
#   geom_line() +
#   geom_point(size = 2, alpha = 0.5) +
#   ggsci::scale_color_lancet() +
#   facet_grid(SDelay~VOC)  +
#   scale_y_continuous(breaks = c(0, 0.05, 0.15, 0.3), labels = c("0%", "5%","15%", "30%")) +
#   
#   # geom_dl(aes(label = wb), 
#   #         method = list(dl.combine("first.points", "last.points")))
#   geom_vline(xintercept = 25, linetype = 2) +
#   geom_hline(yintercept = c(0), linetype = 2 ) +
#   theme_cowplot() +
#   theme(axis.title = element_text(size = 16),
#         legend.position = "right",
#         # title = element_text(size = 18),
#         plot.subtitle = element_text(size = 16)) +
#   guides(colour = guide_legend(nrow = 1)) +
#   labs(x = "Effective Dosing Intervals",
#        color = "",
#        #  title = "Relative Differences in Cumulative Outcome Compared to B1",
#        y = "%Differences Compared to B1") # -> p
# 
#  ggsave("figs/R1_SA_SDelay.png", p,
#         width = 6, height = 6)

# require(ggridges)
# unroll <- function(tmp){
#   tmp2 <- list()
#   for(p in 1:7){
#     tmp1 <- tmp %>% 
#       filter(scenario == p)
#     roll_sheet <- list()
#     for(n in 1:nrow(tmp1)){
#       roll_sheet[[n]] <- rep(tmp1$elapse, tmp1$value/1000)
#     }
#     tmp2[[p]] <- roll_sheet %>% unlist() 
#   }
#   
#   tmp2 %<>% 
#     map(data.table::data.table) %>% 
#     bind_rows(.id = "scenario")
#   
#   return(tmp2)
# }
# 
# # lapply(1:length(params_3_vp_add_delays[[1]]), function(x){
# #   lapply(params_3_vp_add_delays[[1]][[x]]$scenarios, "[[", "pending") %>%
# lapply(1:length(params_3_vp), function(x){
#   lapply(params_3_vp[[x]]$scenarios, "[[", "pending") %>%
#     bind_rows(.id = "scenario") %>%
#     filter(elapse >= 28) %>%
#     filter(!is.na(t_dose2)) %>%
#     pivot_longer(cols = starts_with("Y")) %>% 
#     group_by(status, scenario, elapse) %>% 
#     summarise(value = sum(value), .groups = "drop")
# }) %>% 
#   setNames(euro_inuse) -> tab
# 
# x <- tab %>% map(unroll) %>% bind_rows(.id = "iso3c")
# 
# x %>% group_by(scenario) %>% group_split() %>% map(mutate, V1 = V1/7) %>% map(pull, V1) %>% map(~summary(.)) %>% bind_rows() %>% .[,c(1,3,4,6)]
# 
# x %>% 
#   mutate(scenario = factor(scenario, levels = c(2,1),
#                            labels = c("B1","B2")),
#          country_name = countrycode::countrycode(iso3c, "iso3c", "country.name")) %>%
#   ggplot(., aes(x = V1, y = country_name, fill = scenario)) +
#   geom_density_ridges(alpha = 0.75) +
#   labs(y = "", color = "", x = "Dosing Interval (Days)") +
#   theme_cowplot() +
#   theme(legend.position = "top",
#         legend.text = element_text(size = 14),
#         axis.text = element_text(size = 14),
#         axis.title = element_text(size = 14)) +
#   ggsci::scale_fill_lancet() +
#   labs(fill = "") +
#   scale_y_discrete(limits = rev) +
#   geom_vline(xintercept = c(4*7, 20*7)) +
#   scale_x_continuous(breaks = c(7*4, 7*20, 7*30, 40*7,50*7,60*7),
#                      labels = c("A1", "A5", "30wk", "40wk", "50wk","60wk")) -> p
# 
# ggsave("figs/EI_R1_52wkDelay.png", p, width = 12, height = 8)



# 

# 
# tab %>% 
#   filter(scenario %in% 1:5) %>%
#   ungroup %>% 
#   dplyr::select(-iso3c, -cn) %>% 
#   group_by(scenario, elapse) %>% 
#   tally()

# params_3_vp <-  params_3_VOC_vp <- list()
# for(m in 1:length(additional_delay)){
#   params_3_vp[[m]] <- params_3_VOC_vp[[m]] <- list()
#   # delay = 12, i = 12
#   # delay = 52, i = 5, 11, 13
#   for(i in 12:length(euro_inuse)){
#     which(model_selected_ur$iso3c ==  euro_inuse[i]) -> j
#     params_3_vp[[m]][[i]] <- add_vp(params = params_3[[j]], 
#                                     delay = additional_delay[m])
#     params_3_VOC_vp[[m]][[i]] <- add_vp(params_3_VOC[[j]], 
#                                         delay = additional_delay[m])
#     print(round(i*100/length(euro_inuse),0))
#     save(params_3_vp, params_3_VOC_vp,
#          file = "data/intermediate/params_3_vp_additional_add_supply_delay.rdata")
#   }
# }

# w1 <- 1:13#; w1 <- w1[!w1 %in% c(12)]
# w2 <- 1:13# ; w2 <- w2[!w2 %in% c(13)]
# 
# for(i in w1){
#   lapply(params_3_VOC_vp[[1]][[i]]$scenarios, "[[", "daily_vac_scenarios") %>% 
#     map(arrange, date) %>% 
#     map(mutate_at, vars(starts_with("Y", ignore.case = T)), cumsum) %>% 
#     bind_rows(.id = "scenario") %>%
#     select(-starts_with("supply")) %>% 
#     pivot_longer(starts_with("Y")) %>% 
#     separate(name, into = c("ag", "dose")) %>% 
#     group_by(scenario, date, dose) %>% 
#     summarise(value = sum(value), .groups = "drop") -> x
#  
#   
#   lapply(params_3_VOC_vp[[1]][[i]]$scenarios, "[[", "daily_vac_scenarios") %>% 
#     map(arrange, date) %>% 
#     map(mutate_at, vars(starts_with("supply", ignore.case = T)), cumsum) %>% 
#     bind_rows(.id = "scenario") %>% 
#     dplyr::select(date, starts_with("supply"), -supply_cum) -> y
#   
#   y %>% filter(supply > 10000) %>% .[1,] %>% pull(date) -> t1
#   y %>% filter(supply_2 > 10000) %>% .[1,] %>% pull(date) -> t2
#   
#   print(t2-t1)
#   
#   x %>% 
#     left_join(y, by = "date") %>% 
#     filter(date >= "2021-01-01") %>%
#     mutate(pop =  vac_denom %>% filter(iso3c == euro_inuse[i]) %>% pull(tot),
#            value = value/pop,
#            supply = supply/pop,
#            supply_2 = supply_2/pop) %>% 
#     ggplot(., aes(x = date, y = value, group = dose, color = dose)) +
#     geom_line() +
#     geom_line(aes(y = supply), color = "black", linetype = 2) +
#     geom_line(aes(y = supply_2), color = "black", linetype = 3) +
#     facet_wrap(~scenario, nrow = 1) +
#     labs(title = euro_inuse[i]) -> p
#   
#   ggsave(paste0("figs/supplemental/Delay_12wk_Sanity/",i,
#                 "_",euro_inuse[i],".png"),
#          p, width = 20, height = 6)
# }
# 
# for(i in w2){
#   lapply(params_3_VOC_vp[[2]][[i]]$scenarios, "[[", "daily_vac_scenarios") %>% 
#     map(arrange, date) %>% 
#     map(mutate_at, vars(starts_with("Y", ignore.case = T)), cumsum) %>% 
#     bind_rows(.id = "scenario") %>%
#     select(-starts_with("supply")) %>% 
#     pivot_longer(starts_with("Y")) %>% 
#     separate(name, into = c("ag", "dose")) %>% 
#     group_by(scenario, date, dose) %>% 
#     summarise(value = sum(value), .groups = "drop") -> x
#   
#   
#   lapply(params_3_VOC_vp[[2]][[i]]$scenarios, "[[", "daily_vac_scenarios") %>% 
#     map(arrange, date) %>% 
#     map(mutate_at, vars(starts_with("supply", ignore.case = T)), cumsum) %>% 
#     bind_rows(.id = "scenario") %>% 
#     dplyr::select(date, starts_with("supply"), -supply_cum) -> y
#   
#   y %>% filter(supply > 10000) %>% .[1,] %>% pull(date) -> t1
#   y %>% filter(supply_2 > 10000) %>% .[1,] %>% pull(date) -> t2
#   
#   print(t2-t1)
#   
#   x %>% 
#     left_join(y, by = "date") %>% 
#     filter(date >= "2021-01-01") %>% 
#     mutate(pop =  vac_denom %>% filter(iso3c == euro_inuse[i]) %>% pull(tot),
#            value = value/pop,
#            supply = supply/pop,
#            supply_2 = supply_2/pop) %>% 
#     ggplot(., aes(x = date, y = value, group = dose, color = dose)) +
#     geom_line() +
#     geom_line(aes(y = supply), color = "black", linetype = 2) +
#     geom_line(aes(y = supply_2), color = "black", linetype = 3) +
#     facet_wrap(~scenario, nrow = 1) +
#     labs(title = euro_inuse[i]) -> p
#   
#   ggsave(paste0("figs/supplemental/Delay_52wk_Sanity/",i,
#                 "_",euro_inuse[i],".png"),
#          p, width = 20, height = 6)
# }
 
  
#### vac cov comparison ####
owid_vac %>% 
    group_by(wb) %>% 
    mutate(date_max = max(date)) %>% 
    filter(date == date_max,
           wb %in% all_euro$iso3c) %>% 
    mutate(con = countrycode(iso_code , "iso3c", "continent")) %>% 
    filter(con == "Europe") %>% 
    mutate(MIC = if_else(iso_code %in% euro_lmic, T, F)) %>% 
    dplyr::select(location, iso_code, date, people_vaccinated_per_hundred, MIC) %>% 
    rename(cov = people_vaccinated_per_hundred) %>% 
    group_by(MIC) %>% 
    summarise(t = mean(cov, na.rm = T))
