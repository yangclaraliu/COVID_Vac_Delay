pacman::p_load(
  tidyverse, covidm
)

source("code/util_data.R")
source("code/util_functions.R")


##### load covidm #####
cm_path <- "code/covidm_for_fitting/"
cm_force_rebuild <- F
cm_build_verbose <- T
cm_version <- 2
source(paste0(cm_path, "/R/covidm.R"))


# ve = RCT results, disease blocking
# ei_v = infection-blocking VE
# ed_vi = internal parameter that involves conditionality
i = 12
params <- gen_country_basics("Thailand") %>% 
  update_vac_char(., 
                  ve_i = 0.6, ve_i2 = 0.9,
                  ve_d = 0, ve_d2 = 0)

n_age_groups = length(params$pop[[1]]$v)

params$pop[[1]]$v = rep(10000, n_age_groups);
params$pop[[1]]$v12 = rep(8000, n_age_groups);
run2 = cm_simulate(params)

# params$schedule[[1]] = list(
#   list(
#     parameter = "v",
#     pops = 0,
#     mode = "assign",
#     values = list(rep(c(500000, 0), times = c(10, n_age_groups - 10)), rep(0, n_age_groups)),
#     times = c(0, 20)
#   )
# )

# params$schedule[["v2"]] = list(
#   list(
#     parameter = "v2",
#     pops = 0,
#     mode = "assign",
#     values = list(rep(c(500000, 0), times = c(10, n_age_groups - 10)), rep(0, n_age_groups)),
#     times = c(0, 40)
#   )
# )

# res <- cm_simulate(params)

run2$dynamics %>%
  group_by(t, compartment) %>%
  summarise(tot = sum(value)) %>%
  ggplot(., aes(x = t, y = tot)) +
  geom_line() +
  facet_wrap(~compartment, scales = "free")
