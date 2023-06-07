library( tidyverse )
library( tidybayes )
library( patchwork )
library( here )
library( furrr )
library( runner )

if(!exists("functions_sourced")){
  source( "0. functions.R" )
}
if(!exists("settings_sourced")){
  source( "0. settings.R" )
}

load_if_needed( "df_variants", here( runname, "output", "model_data", variant_filename) )	
load_if_needed( "df_vaccination", here( runname, "output", "model_data", vaccin_filename) )	
load_if_needed( "df_municipality", here( runname, "output", "model_data", viralload_hospital_filename) )	

load_if_needed( "df_nederland", here(runname,"output","Netherlands","df_netherlands_hosp.RData"))
load_if_needed( "fit_hospitalization", here( runname, "output", "model_data", "fit_hosp.RData"))

##### Time line for the variants (used to colour the data) ####

time_variants <- bind_rows(bind_rows(
  bind_rows(df_variants %>%filter(date <= as.Date("2021-05-01") & variant == "Alpha"),
            df_variants %>%filter(date > as.Date("2021-05-01") & date <= as.Date("2021-11-01") & 
                                    variant == "Delta") %>% mutate(percentage_variant = percentage_variant + 1)),
  df_variants %>%filter(date > as.Date("2021-11-01") & variant == "Omicron") %>% mutate(percentage_variant = percentage_variant + 2)),
  tibble(date = as.Date("2022-12-12"), percentage_variant = 3)) # For a nice legend

##### We plot the fitted hospitalizations for each municipality #####

list.files(here(runname,"output","municipality"), full.names = T) %>%
  as.list() %>% future_walk(function(name){
    x <- read.csv(name) %>% mutate(date = as.Date(date), 
                                   age_group = if_else(age_group == "Totaal",municipality,age_group))
    
    ggplot(left_join(x,time_variants, by = "date"),
           mapping = aes(x = date, y = expected_hospitalizations,
                         ymin = simulated_hospitalizations.lower, 
                         ymax = simulated_hospitalizations.upper)) +
      geom_point(aes(y = 1), alpha = 0) + # Ensure a nice y-axis
      geom_ribbon(alpha = 0.25) +
      geom_text(aes(label = age_group), size = 5, hjust = "right", data = x %>%
                  group_by(age_group) %>% 
                  mutate(expected_hospitalizations = 0.9*max(1,hospitalizations,
                                                         simulated_hospitalizations.upper,
                                                         expected_hospitalizations.upper)) %>%
                  filter(date == lastday -35) %>% 
                  ungroup()) +
      # We plot the expected number of hospitalizations as well as the 95%-intervals
      geom_ribbon(aes(ymin = expected_hospitalizations.lower,
                      ymax = expected_hospitalizations.upper), color = cbPalette[6], alpha = .5) +
      geom_point(aes(y = hospitalizations, color = percentage_variant) ,size = 2.5, alpha = .5) +
      geom_line(aes(y = expected_hospitalizations.upper), color = cbPalette[6], linetype = "dotted") +
      geom_line(aes(y = expected_hospitalizations.lower), color = cbPalette[6], linetype = "dashed") + 
      geom_line(color = cbPalette[6], linetype = "solid") + 
      scale_x_date( "Date", breaks = seq.Date(startday,lastday, by = "3 months"), date_labels = "%m/%y") +      
      scale_y_continuous("Hospitalisations", breaks = ~round(unique(pretty(.)))) + 
      scale_colour_gradientn(colours = cbPalette[c(5,7,6,3)], guide = "none") +
      coord_cartesian(xlim = c(startday + 14,lastday - 21)) + 
      theme_bw(base_size = 20) +
      facet_wrap(vars(age_group),ncol = 1, scales = "free_y") +
      theme(
        plot.title = element_text(color = cbPalette[6]),
        panel.grid.major = element_line(size = 0.7),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text.x = element_blank())
    
    ggsave( here(runname,"figures", "municipality", str_c("hosp_", pull(x,municipality)[1], ".png")),
            width = 0.75*297, height = 300, units = "mm")
    })

#### We plot the largest municipalities as well as 2 smaller ones ####

x <- list.files(here(runname,"output","municipality"), full.names = T, 
                pattern =  df_municipality %>% 
                  filter(date == first(date)) %>%
                  group_by(municipality) %>% 
                  summarize(population = sum(population)) %>%
                  arrange(desc(population)) %>% 
                  pull(municipality) %>% 
                  .[1:4] %>% 
                  paste(collapse ="|") %>% 
                  paste0("|Smallingerland|Vaals") %>%
                  paste0("(",.,")\\.csv")) %>%
  as.list() %>%
  lapply(function(name){x <- read.csv(name) %>% mutate(date = as.Date(date)) %>%
    filter(age_group == "Totaal")}) %>%
  bind_rows() %>%
  mutate(municipality= if_else(municipality == "s-Gravenhage","The Hague",municipality),
    municipality = factor(municipality, 
                               levels = df_municipality %>%
                                 mutate(municipality = if_else(municipality == "s-Gravenhage",
                                                              "The Hague",
                                                              municipality %>% as.character())) %>%
                                 filter(date == first(date)) %>%
                                 group_by(municipality) %>% 
                                 summarize(population = sum(population)) %>%
                                 arrange(desc(population)) %>% 
                                 pull(municipality) %>% 
                                 .[1:4] %>% 
                                 as.character() %>%
                                 c("Smallingerland","Vaals")))

ggplot(left_join(x,time_variants, by = "date"),
       mapping = aes(x = date, y = expected_hospitalizations,
                     ymin = simulated_hospitalizations.lower, 
                     ymax = simulated_hospitalizations.upper)) +
  geom_point(aes(y = 1), alpha = 0) + # Ensures a nice y-axis
  geom_ribbon(alpha = 0.25) +
  geom_text(aes(label = municipality), size = 5, hjust = "right", data = x %>%
              group_by(municipality) %>% 
              mutate(expected_hospitalizations = 0.9*max(1,hospitalizations,
                                                         simulated_hospitalizations.upper,
                                                         expected_hospitalizations.upper)) %>%
              filter(date == lastday -35) %>% 
              ungroup()) +
  # We plot the expected number of hospitalizations as well as the 95%-intervals
  geom_ribbon(aes(ymin = expected_hospitalizations.lower,
                  ymax = expected_hospitalizations.upper), color = cbPalette[6], alpha = .5) +
  geom_point(aes(y = hospitalizations, color = percentage_variant) ,size = 2.5, alpha = .5) +
  geom_line(aes(y = expected_hospitalizations.upper), color = cbPalette[6], linetype = "dotted") +
  geom_line(aes(y = expected_hospitalizations.lower), color = cbPalette[6], linetype = "dashed") + 
  geom_line(color = cbPalette[6], linetype = "solid") + 
  scale_x_date( "Date", breaks = seq.Date(startday,lastday, by = "3 months"), date_labels = "%m/%y") +      
  scale_y_continuous("Hospitalisations", breaks = ~round(unique(pretty(.)))) + 
  scale_colour_gradientn(colours = cbPalette[c(5,7,6,3)], guide = "none") +
  coord_cartesian(xlim = c(startday + 14,lastday - 21)) + 
  theme_bw(base_size = 20) +
  facet_wrap(vars(municipality),ncol = 1, scales = "free_y") +
  theme(
    plot.title = element_text(color = cbPalette[6]),
    panel.grid.major = element_line(size = 0.7),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text.x = element_blank())

ggsave( here(runname,"figures", "manuscript", str_c("hosp_grote_gemeentes.png")),
        width = .75*297, height = 300, units = "mm")


#### A plot of the Netherlands #####

x <- df_nederland %>%
  mutate(simulated_hospitalizations = rpois(n(),expected_hospitalizations)) %>%
  group_by(date,hospitalizations,load) %>%
  median_qi(expected_hospitalizations,simulated_hospitalizations)

ggplot(x, mapping = aes(x = date, y = expected_hospitalizations,
                        ymin = simulated_hospitalizations.lower, 
                        ymax = simulated_hospitalizations.upper)) +
  geom_ribbon(alpha = 0.25) +
  # We plot the expected number of hospitalizations as well as the 95%-intervals
  geom_ribbon(aes(ymin = expected_hospitalizations.lower,
                  ymax = expected_hospitalizations.upper), color = cbPalette[6], alpha = .5) +
  geom_point(aes(y = hospitalizations), color = cbPalette[6] ,size = 2.5, alpha = .5) +
  geom_line(aes(y = expected_hospitalizations.upper), color = cbPalette[6], linetype = "dotted") +
  geom_line(aes(y = expected_hospitalizations.lower), color = cbPalette[6], linetype = "dashed") + 
  geom_line(color = cbPalette[6], linetype = "solid") + 
  geom_line(aes(y = 10^(load-11.7)), color = cbPalette[2]) +
  scale_x_date( "Date", breaks = seq.Date(startday,lastday, by = "3 months"), date_labels = "%m/%y") +      
  scale_y_continuous("Hospitalisations", breaks = ~round(unique(pretty(.))),
                     sec.axis = sec_axis(trans = function(x){x*10^(11.7)}, 
                                         name = expression(paste("Load (10"^"-5"," persons"^"-1",")")))) + 
  coord_cartesian(xlim = c(startday + 14,lastday - 21)) + 
  theme_bw(base_size = 20) +
  theme(
    plot.title = element_text(color = cbPalette[6]),
    panel.grid.major = element_line(size = 0.7),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text.x = element_blank())

ggsave( here(runname,"figures", "manuscript", str_c("Nederland.png")),
        width = .75*297, height = 100, units = "mm")

#### Plot hyper distribution of the hosprates #####

# Hyper distributions
hyper_distribution <- fit_hospitalization %>% 
  recover_types(df_municipality,df_variants) %>%
  spread_draws(mean_hosp_rate[age_group],sigma_hosp_rate[age_group],hosp_rate_variant[variant]) %>% 
  group_by(age_group,variant) %>%
  mutate(shape = mean_hosp_rate * sigma_hosp_rate,
         rate = sigma_hosp_rate,
         xmax = mean(hosp_rate_variant) * mean(shape)/mean(rate)*3) %>%
  group_by(.draw,shape,rate,hosp_rate_variant,age_group,variant) %>%
  summarize(x = seq(0,xmax,length.out = 50),
            .groups = "drop") %>%
  mutate(chance = dgamma(x/hosp_rate_variant,shape,rate),
         chance = if_else(is.infinite(chance),NA_real_,chance)) %>%
  group_by(age_group,variant,x) %>%
  mean_qi(chance, na.rm = T) %>%
  # Scale everything to approximately one
  group_by(age_group,variant) %>% 
  mutate(across(c(".lower","chance",".upper"), function(x){x/max(.upper)})) %>%
  ungroup() %>%
  mutate(age_group = paste0("Age group:\n",age_group),
         variant = paste0("Variant:\n",variant))

# Distribution of the 4 largest municipalities
city_distribution <- fit_hospitalization %>%
  recover_types(df_municipality,df_variants) %>%
  spread_draws(hosp_rate[age_group,municipality], hosp_rate_variant[variant]) %>%
  mutate(hosp_rate = hosp_rate * hosp_rate_variant) %>% 
  select( - hosp_rate_variant) %>%
  filter(municipality %in% c("Amsterdam","Rotterdam","s-Gravenhage","Utrecht")) %>%
  mutate(municipality = if_else(municipality == "s-Gravenhage",
                                "The Hague",
                                municipality %>% as.character())) %>%
  left_join(hyper_distribution %>%
              group_by(age_group,variant) %>%
              summarize(x = max(x), .groups = "drop")) %>%
  group_by(age_group,variant) %>%
  mutate(resolutie = .9*(max(hosp_rate) - min(hosp_rate))/(15-1),
         hosp_rate = cut(hosp_rate,
                         breaks = seq(min(hosp_rate),
                                      max(hosp_rate),
                                      length.out= 15),
                         labels = seq(min(hosp_rate),
                                      max(hosp_rate),
                                      length.out = 15)%>%
                           mean_run(k = 2) %>% .[-1],
                         include.lowest = T) %>%
           as.character() %>% as.numeric()) %>%
  rename("Municipality" = "municipality")  %>% 
  # Scale everything to max 1
  group_by(Municipality,age_group,variant,hosp_rate) %>% 
  summarize(count = n(),
            resolutie = first(resolutie),
            .groups = "drop_last") %>% 
  mutate(count = count / max(count)) %>%
  mutate(age_group = paste0("Age group:\n",age_group),
         variant = paste0("Variant:\n",variant))

# Plot the hyper distribution for the different age groups during wildtype
ggplot() + 
  geom_line(aes(x = x, y = chance),
            color = cbPalette[6],
            data = hyper_distribution %>% filter(str_detect(variant,"Wildtype",negate = F))) +
  geom_ribbon(aes(x= x, ymin = .lower, ymax = .upper),
              alpha = .5,
              fill = cbPalette[6],
              data = hyper_distribution %>% filter(str_detect(variant,"Wildtype",negate = F))) +
  geom_col(aes(x = hosp_rate, y = count, fill = Municipality, width = resolutie),
           alpha = .6,
           position = "identity",
           data = city_distribution %>% filter(str_detect(variant,"Wildtype",negate = F))) +
  scale_y_continuous("Density (scaled to [0,1])") +
  scale_x_continuous(expression(paste("Hospitalisations (10"^"-13"," virus load"^"-1"," 10"^"-6"," persons"^"-1",")"))) +
  facet_grid(cols = vars(age_group),
             scales = "free")

ggsave(here(runname,"figures","manuscript","hosp_rates.png"), width = 9, height = 3, units = "in")

# Plot the hyper distribution for the different variants
hyper_distribution %>%
  mutate(variant = factor(variant, levels = paste0("Variant:\n",variant_names))) %>%
  ggplot() + geom_line(aes(x = x, y = chance, color = variant)) +
  geom_ribbon(aes(x= x, ymin = .lower, ymax = .upper, fill = variant),
              alpha = .5) +
  scale_y_continuous("Density (scaled to [0,1])") +
  scale_x_continuous(expression(paste("Hospitalisations (10"^"-13"," virus load"^"-1"," 10"^"-6"," persons"^"-1",")"))) +
  facet_grid(cols = vars(age_group), rows = vars(variant), scales = "free") +
  guides(colour = "none", fill = "none") 

ggsave(here(runname,"figures","manuscript","hosp_rates_variants.png"), width = 9, height = 8, units = "in")


#### Plot the different effects of vaccination ####

fit_hospitalization %>% 
  recover_types(df_municipality,df_variants,df_vaccination) %>%
  spread_draws(prevention_vax[age_group,vaccination,variant]) %>%
  filter(variant != "Wildtype", 
         vaccination != "Unvaccinated") %>%
  rename("Age group" = "age_group",
         "Vaccination" = "vaccination",
         "Variant"  = "variant") %>%
  mutate(Vaccination = factor(Vaccination, levels = c("Partially","Fully","Boosted"),
                              labels = paste0("Vaccination:\n",c("Partially","Fully","Boosted"))),
         Variant = paste0("Variant:\n",Variant)) %>%
  ggplot() + 
  geom_density(aes(x = x,after_stat(scaled)), data = tibble(x = rbeta(10^5,2,8)),
               alpha = .4,
               linetype = "dotted") +
  geom_density(aes(x = prevention_vax, color = `Age group`, after_stat(scaled))) +
  xlim(0,1) + 
  facet_grid(rows = vars(Vaccination), 
             cols = vars(Variant),
             labeller = labeller(.cols = )) +
  xlab("Effect of vaccination on propotionality constant") +
  ylab("Density (scaled to [0,1])") +
  scale_color_discrete(name = "Age group", aesthetics = c("colour","fill")) 

ggsave(here(runname,"figures","manuscript","effect_vaccins.png"), width = 8, height = 4, units = "in")


#### Pearson residue's of the largest municipalities as well as the Netherlands ####

x <- list.files(here(runname,"output","municipality"), full.names = T, 
                pattern = c("(Amsterdam|Rotterdam|s-Gravenhage)\\.csv")) %>%
  as.list() %>%
  lapply(function(name){x <- read.csv(name) %>% mutate(date = as.Date(date)) %>%
    filter(age_group == "Totaal")}) %>%
  bind_rows() %>%
  mutate(municipality= if_else(municipality == "s-Gravenhage","The Hague",municipality),
         municipality = as.factor(municipality),
         pearson_res = (hospitalizations - expected_hospitalizations) / 
           sqrt(expected_hospitalizations))

ggplot(left_join(x,time_variants, by = "date"),
       mapping = aes(x = expected_hospitalizations, y = pearson_res)) +
  geom_text(aes(label = municipality), size = 5, hjust = "right", data = x %>%
              group_by(municipality) %>% 
              summarize(pearson_res = 0.9*max(1,pearson_res),
                        expected_hospitalizations = 1 * max(expected_hospitalizations)) %>% 
              ungroup()) +
  # We plot the expected number of hospitalizations as well as the 95%-intervals
  geom_point(aes( color = percentage_variant) ,size = 2.5, alpha = .5) +
  scale_y_continuous("Pearson's residuals", breaks = ~round(unique(pretty(.)))) + 
  xlab("Median of modelled hospitalisations") +
  scale_colour_gradientn(colours = cbPalette[c(5,7,6,3)], guide = "legend",
                         name = "Dominant\nvariant", breaks = 0:3,
                         labels = c("Wildtype","Alpha","Delta","Omicron"),
                         limits = c(0,3)) +
  theme_bw(base_size = 20) +
  facet_wrap(vars(municipality),ncol = 1, scales = "free") +
  theme(
    plot.title = element_text(color = cbPalette[6]),
    panel.grid.major = element_line(size = 0.7),
    panel.grid.minor = element_blank(),
    #   legend.position = "none",
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text.x = element_blank())

ggsave( here(runname,"figures", "manuscript", str_c("pearson_residuals.png")),
        width = .75*297, height = 150, units = "mm")


ggplot(left_join(df_nederland,time_variants, by = "date") %>%
         mutate(pearson_res = (hospitalizations - expected_hospitalizations) / 
                  sqrt(expected_hospitalizations)) %>% 
         group_by(date) %>%
         summarize(pearson_res = median(pearson_res),
                   expected_hospitalizations = median(expected_hospitalizations),
                   percentage_variant = first(percentage_variant)),
       mapping = aes(x = expected_hospitalizations, y = pearson_res)) +
  # We plot the expected number of hospitalizations as well as the 95%-intervals
  geom_point(aes( color = percentage_variant) ,size = 2.5, alpha = .5) +
  scale_y_continuous("Pearson's residuals", breaks = ~round(unique(pretty(.)))) + 
  xlab("Median of modelled hospitalisations") +
  scale_colour_gradientn(colours = cbPalette[c(5,7,6,3)], guide = "legend",
                         name = "Dominant\nvariant", breaks = 0:3,
                         labels = c("Wildtype","Alpha","Delta","Omicron"),
                         limits = c(0,3)) +
  theme_bw(base_size = 20) +
  theme(
    plot.title = element_text(color = cbPalette[6]),
    panel.grid.major = element_line(size = 0.7),
    panel.grid.minor = element_blank(),
    #   legend.position = "none",
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text.x = element_blank())

ggsave( here(runname,"figures", "manuscript", str_c("pearson_residuals_national.png")),
        width = .75*297, height = 100, units = "mm")