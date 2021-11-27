require(DEoptim)

country_dictionary %>% 
  filter(iso %in% euro_lmic) %>% 
  pull(country_name) -> fit_yet

fit_func_3 <- function(input){
  gen_country_basics(country = cn,
                     date_start = as.character(ymd("2019-12-01") + input[1]),
                     waning_nat = 52*7*3,
                     date_end = "2021-03-01",
                     R0_assumed = input[2]) -> params
  
  # lapply(params$processes, "[[", "names")
  
  params$processes[[1]]$prob[1,] <- params$processes[[1]]$prob[1,]*input[3] 
  
  cm_simulate(params) %>%
    .[["dynamics"]] %>% 
    .[compartment == "death_o", sum(value), by = t] %>%
    rename(t_internal = t,
           deaths_predicted = V1) %>%
    mutate(date = ymd("2019-12-01") + floor(input[1] + as.numeric(t_internal))) %>% 
    full_join(., data_tmp, by = "date") %>% 
    mutate(deaths = if_else(is.na(deaths), 
                            0, 
                            deaths) %>% round,
           deaths_predicted = if_else(is.na(deaths_predicted), 
                                      0, 
                                      deaths_predicted),
           ll = dpois(deaths, deaths_predicted, log = T),
           ll = if_else(is.infinite(ll), -12, ll)) %>% 
    pull(ll) %>% sum  -> a
  
  return(-a)
}

out_3 <-  list()
cal_start <- "2019-12-01"

for(i in 18:length(fit_yet)){
  
  cn <- fit_yet[i]
  iso3c_tmp <- countrycode(cn, "country.name","iso3c")
  data_tmp <- owid_epi[iso3c == iso3c_tmp] %>% .[order(date)] %>% .[date <= "2021-03-01"]
  
  owid_epi[iso3c == iso3c_tmp & deaths > 1] %>% 
    .[order(date)] %>% 
    pull(date) %>% .[1] - 90 -> intro_LL
  
  intro_LL + 90 -> intro_UL
  
  t_LL <- as.numeric(ymd(intro_LL) - ymd(cal_start))
  t_UL <- as.numeric(ymd(intro_UL) - ymd(cal_start))
  
  controlDE <- list(reltol=1e-8, steptol=20, itermax = 100, trace = 10,
                    parallelType = 2)
  
  # three parameter fitting
  DEoptim(fn = fit_func_3,
          lower = c(t_LL, 1.5, 0.01),
          upper = c(t_UL, 5, 1),
          control = controlDE) -> out_3[[i]]
  res_3 <- out_3[[i]]$member$bestmemit %>% tail(1)
  res_3[1] <- round(res_3[1])
  
  # plot results
  gen_country_basics(country = cn,
                     date_start = as.character(ymd("2019-12-01") + res_3[1]),
                     date_end = "2021-03-01",
                     waning_nat = 52*7*3,
                     R0_assumed = res_3[2]) %>% 
    cm_simulate() %>% 
    .[["dynamics"]] %>% 
    .[compartment == "death_o", sum(value), by = list(t)] %>% 
    mutate(n_param = 3,
           t_external = res_3[1])-> dyna_3
  
  dyna_3 %>% 
    mutate(status = "scaled",
           V1 = V1*res_3[3]) %>% 
    bind_rows(dyna_3) %>% 
    mutate(status = if_else(is.na(status), "raw", status)) %>% 
    rename(t_internal = t,
           deaths_predicted = V1) %>%
    mutate(date = ymd("2019-12-01") + t_external + as.numeric(t_internal)) %>% 
    full_join(., data_tmp, by = "date") %>% 
    mutate(deaths = if_else(is.na(deaths),
                            0,
                            deaths) %>% round,
           n_param = paste0("# of varying parameters = ", n_param)) %>%
    ggplot(., aes(x = date, color = status)) +
    geom_point((aes(y = deaths)), color = "black") +
    geom_point(aes(y = deaths_predicted)) +
    scale_color_manual(values = c("red", "orange")) +
    facet_wrap(~n_param, ncol = 1) +
    labs(title = cn, x = "") +
    # geom_point(aes(y = deaths_predicted*res_3[3]), color = "blue") +
    theme_bw() ->  p_tmp_3
  
  ggsave(paste0("figs/intermediate/fitted/", cn, ".png"),
         plot  = p_tmp_3,
         width = 10, height = 10)
  
  print(paste0(cn, " completed!"))
  save( out_3, file = "data/intermediate/out_combined.rdata")
  rm(cn)
}

lapply(lapply(out_3, "[[", "member"), "[[", "bestmemit") %>%
  map(tail, 1) %>%
  map(data.frame) %>%
  bind_rows(., .id = "country_index") %>%
  left_join(data.frame(iso3c = euro_lmic) %>% 
              filter(iso3c != "XKX") %>% 
              rownames_to_column(var = "country_index"), 
            by = "country_index") %>%
  rename(t = par1, r = par2, ur = par3) %>%
  mutate(t = round(t),
         date = ymd("2020-12-01") + t) %>%
  write_rds(., "data/intermediate/DEoptim3_selected.rds")

