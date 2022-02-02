path_dropbox <- "C:/Users/eideyliu/Dropbox/Github_Data/COVID-Vac_Delay/"

#### load packages ####
pacman::p_load(
  tidyverse, sf, countrycode, rnaturalearth, magrittr, data.table,
  ggsflabel, mgcv, pspline, viridis, ggsci, mgcv, imputeTS, ggpattern,
  ggpubr, gridExtra, grid,   tidyverse, sf, countrycode, rnaturalearth, 
  magrittr, data.table, ggsflabel, mgcv, pspline, viridis, ggsci, mgcv, 
  imputeTS, cowplot, qs, testthat
)

# imported fitted results
# old fitting table, needs to be updated
# model_selected <- read_rds(paste0(path_dropbox, "DEoptim2_selected_debug.rds")) %>% 
#   mutate(iso = countrycode(wb, "iso3c","wb"))

model_selected_ur <- read_rds(paste0("data/intermediate/DEoptim3_selected.rds"))

# read_rds(paste0(path_dropbox, "DEoptim3_selected_debug.rds")) %>% 
#   mutate(iso = countrycode(wb, "iso3c","wb")) %>% 
#   dplyr::select(country_name, iso) %>% 
#   write_rds(., "data/country_dictionary.rds")

country_dictionary <- read_rds("data/country_dictionary.rds")

members_all <- model_selected_ur$iso3c
members_remove <- read_rds(paste0(path_dropbox, 
                                 "/intermediate/members_remove.rds"))

## World bank list as of June 2020
# https://databank.worldbank.org/data/download/site-content/CLASS.xls
euro_lmic <- c("ALB","ARM","AZE","BLR","BIH","BGR","GEO","KAZ","XKX","KGZ",
               "MDA","MNE","MKD","RUS","SRB","TJK","TUR","TKM","UKR","UZB")
euro_inuse <- setdiff(euro_lmic, members_remove)

# updated contact matrices
load(paste0(path_dropbox, "contact_all.rdata"))
load(paste0(path_dropbox, "contact_work.rdata"))
load(paste0(path_dropbox, "contact_home.rdata"))
load(paste0(path_dropbox, "contact_school.rdata"))
load(paste0(path_dropbox, "contact_others.rdata"))

model_selected_ur %<>% 
  mutate(to_replace = iso3c %in% names(contact_all)) %>% 
  left_join(country_dictionary, c("iso3c" = "iso",
                                  "country_name"))

tmp <- cm_parameters_SEI3R("Thailand")
ag_labels <- tmp$pop[[1]]$group_names; rm(tmp)

for(i in 1:nrow(model_selected_ur)){
  if(model_selected_ur$to_replace[i]){
    cm_matrices[[model_selected_ur$country_name[i]]]$home <-
      as.matrix(contact_home[[model_selected_ur$iso3c[i]]]) %>% 
      set_colnames(ag_labels) %>% 
      set_rownames(ag_labels)
    
    cm_matrices[[model_selected_ur$country_name[i]]]$work <-
      as.matrix(contact_work[[model_selected_ur$iso3c[i]]]) %>% 
      set_colnames(ag_labels) %>% 
      set_rownames(ag_labels)
    
    cm_matrices[[model_selected_ur$country_name[i]]]$school <-
      as.matrix(contact_school[[model_selected_ur$iso3c[i]]]) %>% 
      set_colnames(ag_labels) %>% 
      set_rownames(ag_labels)
    
    cm_matrices[[model_selected_ur$country_name[i]]]$other <-
      as.matrix(contact_others[[model_selected_ur$iso3c[i]]]) %>% 
      set_colnames(ag_labels) %>% 
      set_rownames(ag_labels)
  }
}

# fix discrepancies among country names
names(cm_matrices)[names(cm_matrices) == "TFYR of Macedonia"] <- 
  "North Macedonia"

# proxy for contact matrices (currently unavailable ?)
cm_matrices[["Republic of Moldova"]] <- cm_matrices$Romania
cm_matrices[["Turkmenistan"]] <- cm_matrices$Uzbekistan


#### load epidemic parameters ####
#####  Clinical Fraction #####
# (based on Davies et al, Nature paper) 
cf <- c(
  0.2904047, 0.2904047, 0.2070468, 0.2070468, 0.2676134,
  0.2676134, 0.3284704, 0.3284704, 0.3979398, 0.3979398,
  0.4863355, 0.4863355, 0.6306967, 0.6306967, 0.6906705, 0.6906705
)

##### susceptibility #####
# (based on Davies et al, Nature paper)
sus <- c(
  0.3956736, 0.3956736, 0.3815349, 0.3815349, 0.7859512,
  0.7859512, 0.8585759, 0.8585759, 0.7981468, 0.7981468,
  0.8166960, 0.8166960, 0.8784811, 0.8784811, 0.7383189, 0.7383189
)

#### Google mobility ####
gm_type <- c("retail", "grocery", "parks", "transit", "work", "residential")
# gm <- fread("https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv") %>%
# .[sub_region_1 == "" & sub_region_2 == "" & metro_area == ""] %>%
# .[, wb := countrycode::countrycode(
#   country_region,
#   "country.name",
#   "wb"
# )] %>%
#   .[, !c(
#     "sub_region_1", "sub_region_2", "metro_area", "iso_3166_2_code",
#     "census_fips_code", "country_region_code", "place_id"
#   )] %>%
#   setnames(., c(
#     "country_name", "date",
#     gm_type, "wb"
#   ))
# 
# write_rds(gm, "data/gm.rds")
gm <- read_rds(paste0(path_dropbox, "gm.rds"))
print(paste0("The Google Mobility Report in used was downloaded on ", 
             file.info(paste0(path_dropbox, "gm.rds"))$mtime, "."))

# mobility scalers 
curves <- data.table(
  work_scaler = c(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0.008, 0.021, 0.033, 0.046, 0.058, 0.071, 0.083, 0.096, 0.108, 0.121, 0.133,
    0.146, 0.158, 0.171, 0.183, 0.196, 0.208, 0.221, 0.233, 0.246, 0.258, 0.271,
    0.283, 0.296, 0.308, 0.321, 0.334, 0.346, 0.359, 0.371, 0.384, 0.397, 0.41,
    0.422, 0.435, 0.448, 0.461, 0.474, 0.487, 0.5, 0.513, 0.526, 0.539, 0.552,
    0.566, 0.579, 0.592, 0.606, 0.619, 0.633, 0.646, 0.66, 0.674, 0.687, 0.701,
    0.715, 0.729, 0.743, 0.757, 0.771, 0.785, 0.799, 0.813, 0.828, 0.842, 0.856,
    0.87, 0.885, 0.899, 0.914, 0.928, 0.942, 0.957, 0.971, 0.986, 1, 1.014, 1.029,
    1.043, 1.058, 1.072, 1.087, 1.101, 1.115, 1.13, 1.144, 1.159, 1.173, 1.188,
    1.202, 1.216, 1.231, 1.245, 1.26, 1.274, 1.289, 1.303, 1.317, 1.332, 1.346, 1.361
  ),
  other_scaler = c(
    0.064, 0.066, 0.067, 0.068, 0.069, 0.071, 0.072, 0.073, 0.075, 0.076, 0.077, 0.078,
    0.08, 0.081, 0.082, 0.084, 0.085, 0.086, 0.087, 0.089, 0.09, 0.091, 0.092, 0.094,
    0.095, 0.096, 0.098, 0.099, 0.1, 0.101, 0.103, 0.104, 0.105, 0.106, 0.108, 0.109,
    0.11, 0.112, 0.113, 0.114, 0.116, 0.118, 0.119, 0.121, 0.123, 0.125, 0.128, 0.13,
    0.132, 0.135, 0.137, 0.14, 0.143, 0.146, 0.15, 0.154, 0.159, 0.164, 0.169, 0.175,
    0.182, 0.19, 0.198, 0.207, 0.217, 0.228, 0.24, 0.252, 0.266, 0.28, 0.295, 0.31,
    0.327, 0.344, 0.361, 0.379, 0.398, 0.418, 0.438, 0.459, 0.48, 0.502, 0.525, 0.549,
    0.572, 0.597, 0.621, 0.647, 0.672, 0.698, 0.725, 0.751, 0.778, 0.805, 0.833, 0.86,
    0.888, 0.916, 0.944, 0.972, 1, 1.028, 1.056, 1.084, 1.112, 1.14, 1.168, 1.196, 1.224,
    1.252, 1.28, 1.308, 1.337, 1.365, 1.393, 1.421, 1.449, 1.477, 1.505, 1.533, 1.561,
    1.589, 1.617, 1.645, 1.673, 1.701
  ),
  perc = round(seq(0, 1.25, 0.01), 2)
)

#### OXCGRT ####
##### school related parameters #####
# oxcgrt_raw <- fread("https://github.com/OxCGRT/covid-policy-tracker/raw/master/data/OxCGRT_latest.csv")
# wb_dic <- oxcgrt_raw[,"CountryName"] %>% 
#   distinct() %>% 
#   mutate(wb = countrycode(CountryName, "country.name", "wb"))
# 
# oxcgrt_raw %>%
#   filter(Jurisdiction == "NAT_TOTAL") %>% 
#   mutate(date = lubridate::ymd(Date)) %>% 
#   left_join(wb_dic,
#             by = "CountryName"
#   ) %>% 
#   dplyr::filter(C1_Flag == 1 | is.na(C1_Flag)) %>%
#   dplyr::select(wb, date, `C1_School closing`) %>% 
#   rename(C1 = `C1_School closing`) %>%
#   pivot_wider(
#     names_from = date,
#     values_from = C1
#   ) %>%
#   ungroup()  %>%
#   pivot_longer(
#     cols = starts_with("202"),
#     names_to = "date",
#     values_to = "C1"
#   ) %>%
#   pivot_wider(names_from = wb, values_from = C1) %>% 
#   group_by(date) %>% 
#   pivot_longer(cols = wb_dic$wb,
#                names_to = "wb",
#                values_to = "C1") %>% 
#   mutate(C1 = if_else(is.na(C1), 0, C1)) -> oxcgrt
# 
# ##### stringency index #####
# si <- oxcgrt_raw %>%
#   dplyr::select(CountryName, Date, StringencyIndex) %>%
#   mutate(date = lubridate::ymd(Date)) %>%
#   left_join(wb_dic, by = "CountryName") %>% 
#   dplyr::select(
#     -Date,
#     -CountryName
#   ) %>%
#   distinct()
# 
# si %>%
#   filter(!is.na(StringencyIndex)) %>%
#   group_by(date) %>%
#   tally() %>%
#   mutate(
#     n_max = max(n),
#     missing = (n_max - n) / n_max
#   ) %>%
#   filter(missing > 0.1,
#          date > "2021-01-01") %>%
#   pull(date) %>%
#   min() -> si_stopdate
# 
# si %<>%
#   filter(date < si_stopdate) %>%
#   group_by(wb) %>%
#   arrange(date) %>%
#   group_split() %>%
#   map(mutate, StringencyIndex = zoo::na.locf(StringencyIndex)) %>%
#   bind_rows() %>%
#   as.data.table()

# qsave(si, "data/si.qs")
# qsave(oxcgrt, "data/oxcgrt.qs")

si <- qread(paste0(path_dropbox, "si.qs"))
oxcgrt <- qread(paste0(path_dropbox,"oxcgrt.qs"))
print(paste0("The Stringency Index and School Closure Data in used was downloaded on ", 
             file.info(paste0(path_dropbox, "si.qs"))$mtime, "."))

si %>%
  filter(!is.na(StringencyIndex)) %>%
  group_by(date) %>%
  tally() %>%
  mutate(
    n_max = max(n),
    missing = (n_max - n) / n_max
  ) %>%
  filter(missing > 0.1,
         date > "2021-01-01") %>%
  pull(date) %>%
  min() -> si_stopdate

#### our world in data ####
# epi data
# download.file(,
#               paste0("data/owid.csv"))
# owid_epi <- read_csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv")
# 
# owid_epi[,"location"] %>% distinct() %>%
#   mutate(iso3c = countrycode(location, "country.name", "iso3c")) %>%
#   filter(!is.na(iso3c)) %>%
#   left_join(owid_epi, by = "location") %>%
#   .[, c("location", "iso3c","date", "new_deaths_smoothed", "new_cases_smoothed")] %>%
#   setNames(c("loc", "iso3c","date", "deaths", "cases")) %>%
#   mutate_at(vars(c("deaths", "cases")), ~if_else(is.na(.), 0, .))%>%
#   mutate_at(vars(c("deaths", "cases")), ~if_else(.<0, 0, .)) %>%
#   data.table %>%
#   .[,date := lubridate::ymd(date)] -> epi
# 
# qsave(epi, "data/epi.qs")
# owid_epi <- qread(paste0(path_dropbox, "epi.qs")) 
owid_epi <- qread("data/epi.qs") 

# vaccination data
# owid_vac <- read_csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/vaccinations.csv") 
# owid_vac[,"location"] %>% distinct() %>% 
#   mutate(wb = countrycode(location, "country.name", "wb")) %>% 
#   filter(!is.na(wb)) %>%
#   left_join(owid_vac, by = "location") -> owid_vac
# qsave(owid_vac, "data/owid_vac.qs")

owid_vac <- qread(paste0(path_dropbox, "owid_vac.qs"))
print(paste0("Data from Our World in Data was downloaded on ", 
             file.info(paste0(path_dropbox, "owid_vac.qs"))$mtime,
             "."))

cm_populations %>% 
  mutate(tot = f + m) %>% 
  filter(location_type == 4) %>% 
  group_by(name) %>% 
  summarise(tot = sum(tot) * 1000) %>% 
  mutate(iso3c = countrycode(name, "country.name", "iso3c")) -> vac_denom

# vaccine efficacy tested
data.table(ve_i_o = c(0.67, 0.68),
           ve_d_o = c(0.67, 0.78)) %>% 
  mutate(
    ve_d = exp_ve(ve_d_o, ve_i_o),
    ve_h = c(0.845, 0.9),
    ve_mort = c(0.845, 0.95)) -> ve

##### healthcare parameters ####
critical2 <- 0
picu_cocin_func <- function(age) {
  x <- c(-0.1309118, 0, 17.2398874, 65.7016492, 100)
  y <- c(-2.1825091, -2.1407043, -1.3993552, -1.2344361, -8.8191062)
  p <- splinefun(x, y)(age)
  exp(p) / (1 + exp(p))
}
picu_cocin <- picu_cocin_func(0:85)

# Infection fatality rate (derived from Levin et al., preprint)
ifr_levin <- 100 * exp(-7.56 + 0.121 * 0:85) / (100 + exp(-7.56 + 0.121 * 0:85)) / 100
# Infection hospitalisation rate (derived from Salje et al., Science)
ihr_salje <- exp(-7.37 + 0.068 * 0:85) / (1 + exp(-7.37 + 0.068 * 0:85))
# Amalgamate probabilities
probabilities <- data.table(age = 0:85, ihr = ihr_salje, ifr = ifr_levin, picu = picu_cocin)
probabilities[, age_group := pmin(15, age %/% 5)]
probabilities <- probabilities[, lapply(.SD, mean), by = age_group, .SDcols = 2:4]

# Create model burden processes
P.critical <- probabilities[, ihr * picu]
P.severe <- probabilities[, ihr * (1 - picu)]
P.death <- probabilities[, ifr]
P.hosp <- P.critical + P.severe

delay_2death <- cm_delay_gamma(26, 5, 60, 0.25)$p
delay_2severe <- cm_delay_gamma(8.5, 5, 60, 0.25)$p
delay_2hosp <- cm_delay_gamma(14.6, 5, 60, 0.25)$p

burden_processes <- list(
  cm_multinom_process("E",       data.frame(death = P.death),                   delays = data.frame(death = delay_2death), report = "o"),
  cm_multinom_process("Ev",      data.frame(death = P.death*(1-ve$ve_mort[1])), delays = data.frame(death = delay_2death), report = "o"),
  cm_multinom_process("Ev2",     data.frame(death = P.death*(1-ve$ve_mort[2])), delays = data.frame(death = delay_2death), report = "o"),
  

  cm_multinom_process("E",       data.frame(to_hosp = P.hosp),                  delays = data.frame(to_hosp = delay_2severe)),
  cm_multinom_process("Ev",      data.frame(to_hosp = P.hosp*(1-ve$ve_h[1])),   delays = data.frame(to_hosp = delay_2severe)),
  cm_multinom_process("Ev2",     data.frame(to_hosp = P.hosp*(1-ve$ve_h[2])),   delays = data.frame(to_hosp = delay_2severe)),
  
  cm_multinom_process("to_hosp", data.frame(hosp = rep(1,16)),                  delays = data.frame(hosp = delay_2hosp),   report = "ip")
)

#### impute contact ####
# simulate end dates, determine the end of mobility imputation
sim_start <- "2020-02-15"
sim_end <- "2022-12-31"
si_post <- 10 # generally should be lower than si_post after examining the data
recovery_period <- 365 # recovery_period days to recover to si_post level of stringency
asc = function(x, y0, y1, s0, s1)
{
  xx = s0 + x * (s1 - s0)
  h0 = exp(s0) / (1 + exp(s0))
  h1 = exp(s1) / (1 + exp(s1))
  h = (exp(xx) / (1 + exp(xx)) - h0) / (h1 - h0)
  y0 + (y1 - y0) * h
}

#### process SI  #### 
si_stopvalue <- si[date == as.Date(si_stopdate)-1]

# draw gradual changes
# si_changes <- list()
# for(i in 1:nrow(si_stopvalue)){
#   asc(seq(from = 0, 
#           to = 1, 
#           length.out = recovery_period),
#       si_stopvalue$StringencyIndex[i],
#       si_post,
#       -5,
#       5) -> si_changes[[i]]
# }
# 
# si_changes %>% 
#   setNames(si_stopvalue$wb) %>%
#   map(data.table) %>%
#   map(~.[, date := seq(as.Date(si$date %>% range %>% .[2]),
#                        as.Date(si$date %>% range %>% .[2]) + recovery_period - 1,
#                        1)]) %>% 
#   bind_rows(.id = "wb") %>% 
#   setNames(c("wb", "StringencyIndex", "date")) %>%
#   bind_rows(si[,list(wb, StringencyIndex, date)]) %>% 
#   right_join(CJ(date = seq(as.Date(sim_start), 
#                            as.Date(sim_end), 
#                            1),
#                 wb = unique(si$wb)), 
#              by = c("wb","date")) %>% 
#   mutate(StringencyIndex = if_else(is.na(StringencyIndex),
#                                    si_post,
#                                    StringencyIndex)) %>% 
#   mutate(status = case_when(date <= si$date %>% range %>% .[2] ~ "empirical",
#                             date %in% seq(as.Date(si$date %>% range %>% .[2]),
#                                           as.Date(si$date %>% range %>% .[2]) + 
#                                             recovery_period - 1,
#                                           1) ~ "transition",
#                             TRUE ~ "post-pandemic")) %>% 
#   mutate(d = as.numeric(date)) -> si_imputed

# qsave(si_imputed, "data/si_imputed.qs")
# si_imputed <- qread(paste0(path_dropbox, "si_imputed.qs"))
# 
# # ggplot(si_imputed, aes(x = date,
# #                        color = status,
# #                        y = StringencyIndex)) +
# #   # geom_jitter() +
# #   geom_line(aes(group = wb)) +
# #   geom_point() +
# #   theme_bw() +
# #   scale_color_lancet() +
# #   facet_wrap(~wb) +
# #   labs(x = "Date", color = "",
# #        title = "OXCGRT - COVID19 Stringecy Index") +
# #   theme(axis.text.x = element_text(angle = 90),
# #         axis.text = element_text(size = 16),
# #         strip.text = element_text(size = 16),
# #         legend.position = "bottom",
# #         legend.text = element_text(size = 16),
# #         title = element_text(size = 16))  -> p
# # 
# # ggsave("figs/figureSX_si.png", p, width = 20, height = 10); rm(p)
# 
# 
# #### process mobility ####
# gm %>%
#   melt(.,
#        id.vars = c("country_name", "wb", "date"),
#        measures.var = gm_type) %>%
#   .[, c("m",
#         "dow",
#         "doy",
#         "d",
#         "variable",
#         "value") := list(factor(month(date), levels = 1:12),
#                          factor(lubridate::wday(date)),
#                          lubridate::yday(date),
#                          as.numeric(date),
#                          factor(variable),
#                          (value + 100)/100)] %>%
#   left_join(., si %>%
#               mutate(d = as.numeric(date)) %>%
#               .[,list(d,wb,StringencyIndex)],
#             by = c("wb", "d")) %>%
#   filter(date < si_stopdate) %>%
#   filter(variable %in% c("retail",
#                          "transit",
#                          "grocery",
#                          "work")) %>%
#   mutate(wb = factor(wb)) -> tab
#  
# fit <- gam(formula = value ~  dow + s(wb, bs = "re") + variable +
#              variable*dow + StringencyIndex + variable*m,
#            data = tab,
#            na.action = na.omit)
# # 
# CJ(date = seq(range(tab$date)[1],
#               as.Date(sim_end),
#               by = 1),
#    wb = unique(tab$wb),
#    variable = c("retail",
#                 "transit",
#                 "grocery",
#                 "work")) %>%
#   .[, c("dow",
#         "doy",
#         "m",
#         "d") := list(lubridate::wday(date) %>% factor,
#                      lubridate::yday(date),
#                      month(date),
#                      as.numeric(date))]  %>%
#   .[, m := if_else(m == 12, 11, m)] %>%
#   .[, m := if_else(m == 1, 2, m)] %>%
#   .[, m := factor(m)] %>%
#   left_join(si_imputed[,c("wb","d", "StringencyIndex","status")],
#             by = c("wb", "d")) %>%
#   split(by = "variable") %>%
#   map(arrange, date) %>%
#   bind_rows() %>%
#   mutate(d = if_else(d >= max(tab$d), max(tab$d), d),
#          country_missing = if_else(is.na(status), T, F)) -> pre_tab
# # 
# val <- predict(fit, pre_tab)
# pre_tab %<>% mutate(predicted = val)
# pre_tab %>%
#   left_join(tab[,c("wb", "date", "variable","value")],
#             by = c("wb","date", "variable")) %>%
#   mutate(imputed = if_else(is.na(value),
#                            predicted,
#                            value)) -> gm_forecast
# # 
# gm_forecast %<>%
#   data.table() %>%
#   dplyr::select(date, wb, variable, status, imputed) %>%
#   distinct() %>%
#   pivot_wider(names_from = variable,
#               values_fn = function(x) mean(x, na.rm = T),
#               values_from = imputed) %>%
#   arrange(date, wb)
# # 
# gm_forecast %>%
#   mutate(work = if_else(work > 1.25, 1.25, work),
#          othx = 0.345*retail + 0.445*transit + 0.21*grocery,
#          othx = if_else(othx > 1.25, 1.25, othx),
#          work = round(work, 2),
#          othx = round(othx, 2)) %>%
#   left_join(curves[,c("perc","work_scaler")], by = c("work" = "perc")) %>%
#   left_join(curves[,c("perc", "other_scaler")], by = c("othx" = "perc")) %>%
#   dplyr::select(-c(grocery, retail, transit, work, othx)) %>%
#   rename(work = work_scaler,
#          other = other_scaler) -> gm_scaled
# # 
# country_data_length <-
#   gm_scaled %>% group_by(wb) %>% group_split() %>% map(nrow) %>% unlist()
# # 
# schedule_raw <- gm_scaled %>%
#   mutate(home = 1,
#          date = as.character(date)) %>%
#   left_join(oxcgrt %>%
#               dplyr::select(date, C1, wb) %>%
#               setNames(c("date", "school", "wb")), # %>%
#             # mutate(school = case_when(school == 0 ~ 1,
#             #                           is.na(school) ~ 1,
#             #                           school == 3 ~ 0,
#             #                           TRUE ~ 0.5),
#             by = c("date", "wb"))  %>%
#   # filter(is.na(school)) %>% pull(date) %>% table
#   # filter(is.na(school))
#   # dplyr::filter(!is.na(school)) %>%
#   mutate(school = case_when(school == 0 ~ 1,
#                             school == 3 ~ 0,
#                             is.na(school) ~ 1,
#                             TRUE ~ 0.5))
# 
# CJ(date = seq(as.Date("2019-12-01"), as.Date("2020-02-14"),1),
#    wb = unique(si$wb)) %>%
#   .[,status := "assumed"] %>%
#   .[,c("work",
#        "other",
#        "home",
#        "school",
#        "date") :=
#       list(1,1,1,1,
#            as.character(date))] -> schedule_pre
# 
# # school holidays
# schedule_raw %<>%
#   bind_rows(schedule_pre) %>%
#   arrange(date) %>%
#   mutate(date = lubridate::ymd(date),
#          month = lubridate::month(date),
#          day = lubridate::day(date),
#          year = lubridate::year(date)) %>%
#   mutate(holiday = if_else(
#     #winter holiday,
#     (year > 2020 & month == 12 & day >=  15) |
#       (year > 2020 & month == 1 & day < 5) |
#       # summer holiday
#       month %in% c(7,8),
#     T,
#     F),
#     school = if_else(holiday, 0, school)) %>%
#   dplyr::select(-holiday) %>%
#   mutate(status = if_else(is.na(status), "averaged", status))
# 
# # qsave(schedule_raw, "data/schedule_raw.qs")
# schedule_raw <- qread(paste0(path_dropbox, "schedule_raw.qs")) %>%   
#   dplyr::select(date, wb, status, home, work, school, other, month, day, year)
# 
# # imputation
# shp <- sf::read_sf(paste0(path_dropbox, "CNTR_RG_60M_2020_4326.shp"))
# nb <- spdep::poly2nb(shp)
# wb_missing <- schedule_raw %>% group_by(wb) %>% tally() %>% 
#   mutate(n_max = max(n)) %>% filter(n != n_max) %>% 
#   filter(wb %in% model_selected$wb) %>% 
#   mutate(iso = countrycode(wb, "wb","country.name"))
# 
# nb_id <- sapply(1:nrow(wb_missing), function(x) {nb[[which(shp$ISO3_CODE == wb_missing$wb[x])]]}) 
# nb_id <- lapply(1:length(nb_id), function(y) nb_id[[y]] %>% map(~shp[.,]) %>% bind_rows() %>% pull(ISO3_CODE))
# tmp <- lapply(1:length(nb_id), function(x) {
#       schedule_raw %>% filter(wb %in% nb_id[[x]]) %>% 
#     group_by(date, month, day, year) %>% 
#     summarise(home = mean(home),
#               work = mean(work),
#               school = mean(school),
#               other = mean(other),
#               .groups = "drop") %>% 
#     mutate(status = "imputed",
#            wb = wb_missing$wb[x]) %>% 
#     dplyr::select(colnames(schedule_raw))
#     }
#   )
# 
# tmp_2replace <- schedule_raw %>% filter(wb %in% wb_missing$wb) %>% 
#   group_by(wb) %>% group_split()
# 
# # fix element 1
# existing_dates <- (bind_rows(tmp_2replace[[1]], tmp[[1]]) %>% group_by(date) %>% tally() %>% 
#                     filter(n > 1) %>% pull(date))
# tmp[[1]] %<>% 
#   filter(!date %in% existing_dates) %>% 
#   bind_rows(tmp_2replace[[1]]) %>% distinct()
# 
# # fix element 2
# existing_dates <- (bind_rows(tmp_2replace[[2]], tmp[[2]]) %>% group_by(date) %>% tally() %>% 
#                      filter(n > 1) %>% pull(date))
# tmp[[2]] %<>% 
#   filter(!date %in% existing_dates) %>% 
#   bind_rows(tmp_2replace[[2]]) %>% distinct()
#   
# # fix element 3
# existing_dates <- (bind_rows(tmp_2replace[[3]], tmp[[3]]) %>% group_by(date) %>% tally() %>% 
#                      filter(n > 1) %>% pull(date))
# 
# tmp[[3]] %<>% 
#   filter(!date %in% existing_dates) %>% 
#   bind_rows(tmp_2replace[[3]]) %>% distinct()
# 
# # merge all imputations back
# schedule_raw %<>% 
#   filter(!wb %in% wb_missing$wb) %>% 
#   bind_rows(tmp)

schedule_raw <- read_rds(paste0(path_dropbox,"schedule_raw.rds"))



