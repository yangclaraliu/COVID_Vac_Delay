# WHO AEFI project
# FGS, 13.09.2021


# load librarires
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# load simulated data from Yang
FGS_VacAllocated = read_rds("data/FGS_VacAllocated.rds")
res_baseline <- read_rds("data/intermediate/res_baseline_v3.rds")
# load("data/intermediate/res_baseline.rdata")

# load resources of (QA)LE estimates
load("data/LE_estimates.rda")

# set path
econ_path = "C:/Users/Frank.Sandmann/Desktop/WHO-AEFI/" 


# add back strategies for each dataset ; Yang: "in res_baseline c(1:5,7,6) correspond to c(A1-5,B1-2)"
res_baseline2 = lapply(1:length(res_baseline), function(x) {
  
  res_baseline[[ names(res_baseline)[x] ]] %>% 
          mutate(strategy = ifelse(scenario == 1, "A1",
                            ifelse(scenario == 2, "A2",
                            ifelse(scenario == 3, "A3",
                            ifelse(scenario == 4, "A4",
                            ifelse(scenario == 5, "A5",
                            ifelse(scenario == 6, "B2",
                            ifelse(scenario == 7, "B1",
                                   NA))))))),
                 strategy = factor(strategy, levels = levels(FGS_VacAllocated$strategy) )) %>% 
    mutate(iso3c = countrycode(population, "country.name", "iso3c")) %>% 
    filter(iso3c %in% euro_inuse)
  
  }) 

# change list names back
names(res_baseline2) = names(res_baseline)
res_baseline = data.table::rbindlist(res_baseline2, fill=TRUE, use.names=TRUE, idcol=TRUE)
rm(res_baseline2)


#res_2 only fitted for t0 and r0
#res_3 fitted for t0, r0 AND underreporting
#there are countries where we weren't able to get a result for you - for example for serbia, 
#res_2 and res_3 is identical cuz underreporting landed on 1

# thus just going with res3 and res3_VOC


# add QALY losses due to mortality
res_mort = res_baseline %>% 
  #filter deaths and res_3
  dplyr::filter(compartment %in% c("death_o", "death"), 
                .id %in% c("res_3", "res_3_VOC")) %>% 
  mutate(compartment = "death") %>%  # rename to align.. only works as long as only looking at mortality alone
  #clean
  dplyr::select(-scenario, -t_start) %>%
  # add LE and QALYs
  left_join(
    LE_estimates %>% 
      rename("population" = "country_name",
             "group" = "AgeGroup") %>%
      dplyr::select(population, group, LEdisc, adjQALEdisc) 
    ) %>%
  mutate(LEdisc      = value*LEdisc,
         adjQALEdisc = value*adjQALEdisc )




# figure 9 in the report (without VOC); old figure, use with VOC instead
res_baseline %>% 
  dplyr::group_by(.id, population, strategy, compartment) %>% 
  summarise_if(is.numeric, sum, na.rm=TRUE) %>% 
  dplyr::group_by(.id, population, compartment) %>% 
  mutate(value_scaled = value/max(value)) %>% #View
  dplyr::filter(.id == "res_3") %>%
  ggplot(aes(x=strategy, y=value_scaled, colour = compartment)) + 
  geom_line(aes(group = compartment, linetype = compartment), size = 1.5) +
  facet_wrap(~population, scales = "free", ncol = 3) +
  labs(x = "Dosing Interval Strategy",
       y = "Proportion scaled to maximum value",
       subtitle = "w/o VOC") +
  theme(text = element_text(size = 20),
        #axis.title.y = element_text(size = 16),
        strip.background =element_rect(fill="white", colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank(),
        axis.line = element_line()) + 
  scale_colour_manual(name = "",
                      values = c("death_o" = "grey20", "hosp_i" = "grey50", "cases" = "grey80"),
                      labels = c("cases", "hospitalisations", "deaths"), #c("deaths", "hospitalisations", "cases"),
                      breaks=c("cases", "hosp_i", "death_o")) + #c("death_o", "hosp_i", "cases")) + 
  scale_linetype_manual(name = "",
                        values = c("death_o" = "solid", "hosp_i" = "longdash", "cases" = "dashed"),
                        labels = c("cases", "hospitalisations", "deaths"), #c("deaths", "hospitalisations", "cases"),
                        breaks=c("cases", "hosp_i", "death_o")) #c("death_o", "hosp_i", "cases"))


# save
ggsave(paste0("Fig9_", format(Sys.time(), "%d%b%Y"), "_v1.png"), 
       path = "figs/supplemental/", width = 40, height = 40, units = "cm")


# figure 9 in the report (with VOC)
res_baseline %>% 
  dplyr::group_by(.id, population, strategy, compartment) %>% 
  summarise_if(is.numeric, sum, na.rm=TRUE) %>% 
  dplyr::group_by(.id, population, compartment) %>% 
  mutate(value_scaled = value/max(value)) %>% #View
  dplyr::filter(.id == "res_3_VOC") %>%
  ggplot(aes(x=strategy, y=value_scaled, colour = compartment)) + 
  geom_line(aes(group = compartment, linetype = compartment), size = 1.5) +
  facet_wrap(~population, scales = "free", ncol = 3) +
  labs(x = "Dosing Interval Strategy",
       y = "Proportion scaled to maximum value",
      subtitle = "w/ VOC") +
  theme(text = element_text(size = 20),
        #axis.title.y = element_text(size = 16),
        strip.background =element_rect(fill="white", colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank(),
        axis.line = element_line()) + 
  scale_colour_manual(name = "",
                      values = c("death" = "grey20", "hosp" = "grey50", "cases" = "grey80"),
                      labels = c("cases", "hospitalisations", "deaths"), #c("deaths", "hospitalisations", "cases"),
                      breaks=c("cases", "hosp", "death")) + #c("death_o", "hosp_i", "cases")) + 
  scale_linetype_manual(name = "",
                        values = c("death" = "solid", "hosp" = "longdash", "cases" = "dashed"),
                        labels = c("cases", "hospitalisations", "deaths"), #c("deaths", "hospitalisations", "cases"),
                        breaks=c("cases", "hosp", "death")) #c("death_o", "hosp_i", "cases")) 
  


# save
ggsave(paste0("Fig9_", format(Sys.time(), "%d%b%Y"), 
              "_v1b.png"), path = "figs/supplemental/", 
       width = 40, height = 40, units = "cm")

# total deaths, and different visuatlisation (not used)

plot_res_mort = res_mort %>%
  dplyr::group_by(.id, population, strategy, compartment) %>% 
  summarise_if(is.numeric, sum, na.rm=TRUE)
  
# deaths per country
plot_res_mort %>%
  dplyr::filter(.id == "res_3") %>%
  #dplyr::filter(.id == "res_3_VOC") %>%
  ggplot(aes(x=forcats::fct_rev(strategy), y=value, fill = strategy)) + 
  geom_bar(stat = "identity") +
  geom_line(aes(group = population), size = 1.5, colour = "darkblue") +
  coord_flip() +
  scale_fill_grey() +
  facet_wrap(~population, scales = "free", ncol = 3) +
  labs(x = "Dosing interval strategy",
       y = "Deaths per country")


# calculate incremental change of deaths
# there's not a no-vac scenario but B1 is used as the baseline scenario - which is 
# (older adults 1st dose) --> (older adults 2nd dose) --> (younger adults first dose) --> (younger adult second dose)

plot_res_mort_incr = plot_res_mort %>%
  dplyr::filter(strategy == "A1") %>%
  dplyr::rename(value_base = value, 
                LEdisc_base = LEdisc, 
                adjQALEdisc_base = adjQALEdisc) %>%
  ungroup() %>%
  dplyr::select(-strategy) %>%
  full_join(plot_res_mort) %>%
  # incremental
  mutate(value = value_base - value, 
         LEdisc = LEdisc_base - LEdisc, 
         adjQALEdisc = adjQALEdisc_base - adjQALEdisc) %>%
  dplyr::select(-value_base, -LEdisc_base, -adjQALEdisc_base)


# plot incremental change of deaths
p1 = plot_res_mort_incr %>%
  dplyr::filter(.id == "res_3") %>%
  ggplot(aes(x=forcats::fct_rev(strategy), y=value, fill = strategy)) + 
  geom_bar(stat = "identity") +
  geom_line(aes(group = population), size = 1.5, colour = "darkblue") +
  coord_flip() +
  scale_fill_grey() +
  facet_wrap(~population, scales = "free", ncol = 3) +
  labs(x = "Dosing interval strategy",
       y = "Change (benefit) in prevented mortality (compared to strategy A1")

p1


# could estimate QALYs lost based on: 
# cases = 0.0307 QALYs # (not used currently!)



## harm from AEFIs

# Based on UK data for AZ; updated 09. Sep 2021
# https://www.gov.uk/government/publications/coronavirus-covid-19-vaccine-adverse-reactions/coronavirus-vaccine-summary-of-yellow-card-reporting

#Up to 1 September 2021, the MHRA had received Yellow Card reports of 416 cases of major thromboembolic events (blood clots) 
#with concurrent thrombocytopenia (low platelet counts) in the UK following vaccination with COVID-19 Vaccine AstraZeneca. 
#Forty five of the 416 reports have been reported after a second dose. Of the 416 reports, 210 occurred in women, and 202 
#occurred in men aged from 18 to 93 years. The overall case fatality rate was 17% with 72 deaths, six of which occurred after 
#the second dose.

#Taking into account the different numbers of patients vaccinated with COVID-19 Vaccine AstraZeneca in different age groups, 
#the data shows that there is a higher reported incidence rate in the younger adult age groups following the first dose 
#compared to the older groups (20.5 per million doses in those aged 18-49 years compared to 10.9 per million doses in those aged 
#50 years and over). The number of first doses given to those in the 18-49 years age group is estimated to be 8.5 million while 
#an estimated 16.3 million first doses have been given to patients aged 50+ years.

#Taking into account the different numbers of patients vaccinated with COVID-19 Vaccine AstraZeneca in different age groups, 
#the data shows that there is a lower reported incidence rate in younger adult age groups following the second dose 
#compared to the older groups (0.9 per million doses in those aged 18-49 years compared to 1.9 per million doses in those aged 
#50 years and over). The number of second doses given to those in the 18-49 years age group is estimated to be 8.1 million while 
#an estimated 15.9 million second doses have been given to patients aged 50+ years


rates_AEFI = data.frame("group" = levels(res_baseline$group),
                        "AEFI_d1" = c(rep(20.5/1000000, 10), rep(10.9/1000000, 6)),
                        "AEFI_d2" = c(rep(0.9/1000000, 10), rep(1.9/1000000, 6)),
                        "AEFImort_d1" = c(rep(20.5/1000000, 10), rep(10.9/1000000, 6)) * 0.178,
                        "AEFImort_d2" = c(rep(0.9/1000000, 10), rep(1.9/1000000, 6)) * 0.133)


# vaccine doses per country and strategy
# ignore mild AEFIs; estimate severe AEFIs only
res_AEFI = FGS_VacAllocated %>% 
  dplyr::group_by(country_name, strategy) %>% 
  summarise_if(is.numeric, sum, na.rm=TRUE) %>% 
  tidyr::pivot_longer(-c(country_name, strategy), names_to = "AgeGroup", values_to = "doses_ttl") %>%
  mutate(VaccDose = sub(".*[A-Z][0-9]+_d", "", AgeGroup),
         AgeGroup = sub("_.*", "", AgeGroup)) %>%
  # get numeric ranges of age groups again
  full_join(data.frame("group" = levels(res_baseline$group),
                       "AgeGroup" = paste0("Y", 1:16))) %>%
  # assign AEFIs
  full_join(rates_AEFI) %>%
  mutate(AEFI_nr     = ifelse(VaccDose == "1", doses_ttl*AEFI_d1,     doses_ttl*AEFI_d2),
         AEFImort_nr = ifelse(VaccDose == "1", doses_ttl*AEFImort_d1, doses_ttl*AEFImort_d2) ) %>%
  dplyr::select(-AEFI_d1, -AEFI_d2, -AEFImort_d1, -AEFImort_d2)


# add QALY losses due to mortality
res_AEFI = res_AEFI %>% 
  left_join(
    LE_estimates %>% 
      dplyr::select(country_name, AgeGroup, LEdisc, adjQALEdisc), 
    by = c("country_name", "group" = "AgeGroup")  ) %>%
  #mutate_at(vars(matches("LE", "^QALE", "^QALEdisc", "^adjQALEdisc" )), .funs = funs(. * AEFImort_nr))
  mutate(LEdisc      = AEFImort_nr*LEdisc,
         adjQALEdisc = AEFImort_nr*adjQALEdisc )

# print total loss per country and strategy
plot_res_AEFI = res_AEFI %>% 
  dplyr::group_by(country_name, strategy) %>% 
  summarise_if(is.numeric, sum, na.rm=TRUE) %>%
  mutate(AEFI_rate        = AEFI_nr/doses_ttl*1000000,
         AEFImort_rate    = AEFImort_nr/doses_ttl*1000000,
         LEdisc_rate      = LEdisc/doses_ttl*1000000,
         adjQALEdisc_rate = adjQALEdisc/doses_ttl*1000000 )

# number of doses per strategy and country
plot_res_AEFI %>%
  ggplot(aes(x=forcats::fct_rev(strategy), y=doses_ttl, fill = strategy)) + 
  geom_bar(stat = "identity") +
  geom_line(aes(group = country_name), size = 1.5, colour = "darkblue") +
  coord_flip() +
  scale_fill_grey() +
  facet_wrap(~country_name, scales = "free", ncol = 3) +
  labs(x = "Dosing interval strategy",
       y = "Rate of life years lost (per 1 million doses)")



# calculate incremental change of deaths
plot_res_AEFI_incr = plot_res_AEFI %>%
  dplyr::filter(strategy == "A1") %>%
  dplyr::select(country_name, strategy,
                AEFImort_nr_base =  AEFImort_nr, 
                LEdisc_base = LEdisc, 
                adjQALEdisc_base = adjQALEdisc) %>%
  ungroup() %>%
  dplyr::select(-strategy) %>%
  full_join(plot_res_AEFI %>%
              dplyr::select(country_name, strategy,
                            AEFImort_nr, LEdisc, adjQALEdisc) ) %>%
  # incremental
  mutate(AEFImort_nr = AEFImort_nr - AEFImort_nr_base, 
         LEdisc = LEdisc - LEdisc_base, 
         adjQALEdisc = adjQALEdisc - adjQALEdisc_base ) %>%
  dplyr::select(-AEFImort_nr_base, -LEdisc_base, -adjQALEdisc_base)


# plot incremental change of deaths
p2 = plot_res_AEFI_incr %>%
  ggplot(aes(x=forcats::fct_rev(strategy), y=AEFImort_nr, fill = strategy)) + 
  geom_bar(stat = "identity") +
  geom_line(aes(group = country_name), size = 1.5, colour = "darkblue") +
  coord_flip() +
  scale_fill_grey() +
  facet_wrap(~country_name, scales = "free", ncol = 3) +
  labs(x = "Dosing interval strategy",
       y = "Change (harm) in additional mortality (compared to strategy A1")

p2



# benefit-risk ratios
p3 = plot_res_AEFI_incr %>% 
  rename("LEdisc_AEFI"      = "LEdisc",
         "adjQALEdisc_AEFI" = "adjQALEdisc") %>%
  full_join( plot_res_mort_incr %>%   
               dplyr::filter(.id == "res_3") %>%
               rename("country_name" = "population" ) ) %>%
  dplyr::group_by(country_name, strategy, .id) %>% 
  summarise_if(is.numeric, sum, na.rm=TRUE) %>%
  mutate(AEFImort_nrRatio = ifelse(AEFImort_nr == 0, 0, value/AEFImort_nr),  #/doses_ttl*1000,
         LEdiscRatio      = ifelse(LEdisc_AEFI == 0, 0, LEdisc/LEdisc_AEFI),  #/doses_ttl*1000,
         adjQALEdiscRatio = ifelse(adjQALEdisc_AEFI == 0, 0, adjQALEdisc/adjQALEdisc_AEFI)) %>%  #/doses_ttl*1000) 
  # only look at non-VOC; results only magnified
  dplyr::filter(.id == "res_3") %>%
  
  ggplot(aes(x=forcats::fct_rev(strategy), y=AEFImort_nrRatio, fill = strategy)) + 
  geom_bar(stat = "identity") +
  geom_line(aes(group = country_name), size = 1.5, colour = "darkblue") +
  coord_flip() +
  scale_fill_grey() +
  facet_wrap(~country_name, scales = "free", ncol = 3) +
  labs(x = "Dosing interval strategy",
       y = "Ratio of change in prevented deaths by the vaccine vs. change in extra deaths due to AEFIs")


p3

p1 / p2 / p3

p2 / p3 + 
  plot_layout(heights = c(13/5, 1)) +
  plot_annotation(tag_levels = 'A', tag_prefix = "(", tag_suffix = ")") & 
  theme(plot.tag = element_text(size = 35))



# figure 8 in the report

# plot ratios with VOC
# p1 = 
plot_res_AEFI_incr %>% 
  rename("LEdisc_AEFI"      = "LEdisc",
         "adjQALEdisc_AEFI" = "adjQALEdisc") %>%
  full_join( plot_res_mort_incr %>% 
               rename("country_name" = "population" ) ) %>%
  dplyr::group_by(country_name, strategy, .id) %>% 
  summarise_if(is.numeric, sum, na.rm=TRUE) %>%
  mutate(AEFImort_nrRatio = ifelse(AEFImort_nr == 0, 0, value/AEFImort_nr),  #/doses_ttl*1000,
         LEdiscRatio      = ifelse(LEdisc_AEFI == 0, 0, LEdisc/LEdisc_AEFI),  #/doses_ttl*1000,
         adjQALEdiscRatio = ifelse(adjQALEdisc_AEFI == 0, 0, adjQALEdisc/adjQALEdisc_AEFI)) %>%  #/doses_ttl*1000) 
  # only look at non-VOC; results only magnified
  dplyr::filter(.id == "res_3_VOC") %>%
  filter(AEFImort_nrRatio < 0 | AEFImort_nrRatio < 0 | LEdiscRatio < 0 | adjQALEdiscRatio < 0)
  
  
  ggplot(aes(x=strategy, y=AEFImort_nrRatio)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group = country_name), alpha =0.5, colour = "grey") +
  geom_point(size = 3, alpha =0.5, colour = "black") +
  geom_vline(xintercept=1, linetype = "dashed") +
  geom_hline(yintercept=0, linetype = "dashed") +
  theme(text = element_text(size = 20),
        #axis.title.y = element_text(size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line()) +
  labs(x = "Dosing Interval Strategy",
       y = "Benefit-risk ratios of change in deaths prevented by\nthe vaccines vs. caused by AEFIs, compared to A1") + 
  geom_text(data= data.frame(xpos = -Inf,ypos =  Inf, annotateText = c("w/ VOC"),
                             hjustvar = c(-0.1), vjustvar = c(1)), size=6.5, 
            aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText))

# plot ratios without VOC
# p2 = 
  plot_res_AEFI_incr %>% 
  rename("LEdisc_AEFI"      = "LEdisc",
         "adjQALEdisc_AEFI" = "adjQALEdisc") %>%
  full_join( plot_res_mort_incr %>% 
               rename("country_name" = "population" ) ) %>%
  dplyr::group_by(country_name, strategy, .id) %>% 
  summarise_if(is.numeric, sum, na.rm=TRUE) %>%
  mutate(AEFImort_nrRatio = ifelse(AEFImort_nr == 0, 0, value/AEFImort_nr),  #/doses_ttl*1000,
         LEdiscRatio      = ifelse(LEdisc_AEFI == 0, 0, LEdisc/LEdisc_AEFI),  #/doses_ttl*1000,
         adjQALEdiscRatio = ifelse(adjQALEdisc_AEFI == 0, 0, adjQALEdisc/adjQALEdisc_AEFI)) %>%  #/doses_ttl*1000) 
  # only look at non-VOC; results only magnified
  dplyr::filter(.id == "res_3") %>%
  filter(AEFImort_nrRatio < 0 | AEFImort_nrRatio < 0 | LEdiscRatio < 0 | adjQALEdiscRatio < 0)
  
  ggplot(aes(x=strategy, y=AEFImort_nrRatio)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group = country_name), alpha =0.5, colour = "grey") +
  geom_point(size = 3, alpha =0.5, colour = "black") +
  geom_vline(xintercept=1, linetype = "dashed") +
  geom_hline(yintercept=0, linetype = "dashed") +
  theme(text = element_text(size = 20),
        #axis.title.y = element_text(size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line()) +
  labs(x = "Dosing Interval Strategy",
       y = "Benefit-risk ratios of change in deaths prevented by\nthe vaccines vs. caused by AEFIs, compared to A1") + 
  geom_text(data= data.frame(xpos = -Inf,ypos =  Inf, annotateText = c("w/o VOC"),
                             hjustvar = c(-0.1), vjustvar = c(1)), size=6.5, 
            aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText))


p1 / p2 + 
  plot_annotation(tag_levels = 'a', tag_prefix = "(", tag_suffix = ")") & 
  theme(plot.tag = element_text(size = 30))


# save
ggsave(paste0("fig5", format(Sys.time(), "%d%b%Y"), "_v1.png"), 
       path = "figs/", width = 40, height = 40, units = "cm")




# countries with negative ratios
plot_res_mort_incr %>% dplyr::filter(value<0)

# countries with ratios <1.0
plot_res_AEFI_incr %>% 
  rename("LEdisc_AEFI"      = "LEdisc",
         "adjQALEdisc_AEFI" = "adjQALEdisc") %>%
  full_join( plot_res_mort_incr %>% 
               rename("country_name" = "population" ) ) %>%
  dplyr::group_by(country_name, strategy, .id) %>% 
  summarise_if(is.numeric, sum, na.rm=TRUE) %>%
  mutate(AEFImort_nrRatio = ifelse(AEFImort_nr == 0, 0, value/AEFImort_nr),  #/doses_ttl*1000,
         LEdiscRatio      = ifelse(LEdisc_AEFI == 0, 0, LEdisc/LEdisc_AEFI),  #/doses_ttl*1000,
         adjQALEdiscRatio = ifelse(adjQALEdisc_AEFI == 0, 0, adjQALEdisc/adjQALEdisc_AEFI)) %>% 
  dplyr::filter(AEFImort_nrRatio>0 & AEFImort_nrRatio<1.0)

# nr positive
(1 - ( nrow(plot_res_mort_incr %>% dplyr::filter(value<0)) / nrow(plot_res_mort_incr) ))*100

# nr total
length(unique(plot_res_mort_incr$population))


(1 - ( nrow(plot_res_AEFI_incr %>% 
       rename("LEdisc_AEFI"      = "LEdisc",
              "adjQALEdisc_AEFI" = "adjQALEdisc") %>%
       full_join( plot_res_mort_incr %>% 
                    rename("country_name" = "population" ) ) %>%
       dplyr::group_by(country_name, strategy, .id) %>% 
       summarise_if(is.numeric, sum, na.rm=TRUE) %>%
       mutate(AEFImort_nrRatio = ifelse(AEFImort_nr == 0, 0, value/AEFImort_nr),  #/doses_ttl*1000,
              LEdiscRatio      = ifelse(LEdisc_AEFI == 0, 0, LEdisc/LEdisc_AEFI),  #/doses_ttl*1000,
              adjQALEdiscRatio = ifelse(adjQALEdisc_AEFI == 0, 0, adjQALEdisc/adjQALEdisc_AEFI)) %>% 
       dplyr::filter(AEFImort_nrRatio>0 & AEFImort_nrRatio<1.0))  / 
         nrow(plot_res_mort_incr %>% dplyr::filter(!value<0) ) ))*100



# few more possible visualisations (not used)


# plot AEFIs
plot_res_AEFI %>%
  ggplot(aes(x=forcats::fct_rev(country_name), y=AEFI_rate, fill = strategy)) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~strategy, ncol = 2)


# rates of AEFI, fatal AEFI, discounted LEs lost, and discounted QALEs lost (adjusted for comorbidity)
plot_res_AEFI %>%
  gather(parm, vals, AEFI_rate:adjQALEdisc_rate) %>%
  dplyr::group_by(strategy, parm) %>% 
  summarise(mean_rate = mean(vals),
            lowerUI = quantile(vals, probs = 0.025),
            upperUI = quantile(vals, probs = 0.975) ) %>% 
  mutate(parm = factor(parm, levels = c("AEFI_rate", "AEFImort_rate", 
                                        "LEdisc_rate", "adjQALEdisc_rate"))) %>%
  ggplot(aes(x=strategy, y=mean_rate)) + #, fill = strategy, colour = strategy
  geom_errorbar(aes(ymin=lowerUI, ymax=upperUI), width=.1, position=position_dodge(0.1)) +
  geom_point(position=position_dodge(0.1)) +
  facet_wrap(~parm, scales = "free", ncol = 2) 



# rates of fatal AEFI per country
plot_res_AEFI %>%
  ggplot(aes(x=forcats::fct_rev(strategy), y=AEFImort_rate, fill = strategy)) + 
  geom_bar(stat = "identity") +
  geom_line(aes(group = country_name), size = 1.5, colour = "darkblue") +
  coord_flip() +
  scale_fill_grey() +
  facet_wrap(~country_name, ncol = 3) +
  labs(x = "Dosing interval strategy",
       y = "Rate of fatal adverse events (per 1 million doses)")


# rates of LEs lost per country
plot_res_AEFI %>%
  ggplot(aes(x=forcats::fct_rev(strategy), y=LEdisc_rate, fill = strategy)) + 
  geom_bar(stat = "identity") +
  geom_line(aes(group = country_name), size = 1.5, colour = "darkblue") +
  coord_flip() +
  scale_fill_grey() +
  facet_wrap(~country_name, ncol = 3) +
  labs(x = "Dosing interval strategy",
       y = "Rate of life years lost (per 1 million doses)")


res_AEFI %>%
  ungroup() %>%
  #mutate(AEFImort_nr_rate = AEFImort_nr/doses_ttl*1000000) %>%
  dplyr::filter(!group %in% c("0-4", "10-14", "15-19")) %>%
  ggplot(aes(x=forcats::fct_rev(strategy), y=AEFImort_nr, fill = group)) + 
  geom_bar(position="stack", stat = "identity") +
  coord_flip() +
  scale_fill_grey() +
  facet_wrap(~country_name, scales = "free", ncol = 3) +
  labs(x = "Dosing interval strategy",
       y = "Rate of life years lost (per 1 million doses)")
