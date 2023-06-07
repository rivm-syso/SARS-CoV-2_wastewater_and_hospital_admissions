### Poisson regression for hospital data

library( tidyverse )
library( here )
library( rstan )
library( tidybayes )
library( loo )
library( furrr )
library( runner )

setwd( here() )

if(!exists("functions_sourced")){
  source( "0. functions.R" )
}
if(!exists("settings_sourced")){
  source( "0. settings.R" )
}

#options(buildtools.check = NULL)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#
# Load files
# 

load_if_needed( "df_variants", variant_filename)
load_if_needed( "df_vaccination", vaccin_filename)
load_if_needed( "df_municipality", viralload_hospital_filename)

# Create factors from the necessary data and order df_municipality

df_variants <- mutate(df_variants,
                      date = as.Date(date),
                      variant = factor(variant, levels = variant_names))

df_vaccination <- mutate(df_vaccination,
                         date = as.Date(date),
                         vaccination = factor(vaccination, levels = vaccination_status))

df_municipality <- mutate(df_municipality,
                          date = as.Date(date),
                          municipality = as.factor(as.character(municipality)),
                          age_group = as.factor(as.character(age_group))) %>%
                          arrange(date,municipality,age_group)

# Save the dataframes
save(df_variants, file = here(runname, "output", "model_data",variant_filename))
save(df_vaccination, file = here(runname, "output", "model_data",vaccin_filename))
save(df_municipality, file = here(runname, "output", "model_data",viralload_hospital_filename))

# run Stan model
fit_hospitalization = stan(
  "hospitalizations.stan",
  model_name = "wastewater_model",
  data = compose_data(
    ref_load = ref_load,
    create_stan_data(df_municipality, df_variants, df_vaccination)),
  init = initials_hosp,
  chains = 10,
  warmup = 200,
  iter = 400,
  refresh = 10,
  thin = 2,
  control = list(adapt_delta = .95, max_treedepth = 12)
)

save(fit_hospitalization, file = here( runname, "output", "model_data", "fit_hosp.RData"))

print(traceplot(fit_hospitalization,"mean_hosp_rate"))
print(traceplot(fit_hospitalization,"sigma_hosp_rate"))
print(traceplot(fit_hospitalization,"prevention_vax"))

# Determine the outcome for each municipality and write as a .csv
future_walk(df_municipality %>% 
              pull(municipality) %>% 
              as.character() %>% 
              unique() %>% 
              as.list(),function(name){
  # First sum over the different vaccinations
  df_vaccination %>% 
    filter(municipality == name) %>%
    left_join(fit_hospitalization %>%
                recover_types(df_municipality,df_variants,df_vaccination) %>%
                spread_draws(prevention_vax[age_group,vaccination,variant]),
              by = c("age_group", "vaccination")) %>%
    group_by(.draw,date,municipality,age_group,variant) %>%
    summarize(hosp_rate_vaccins = sum(prevention_vax * percentage_vaccination), .groups = "drop") %>%
    # Then sum over the different variants
    left_join(fit_hospitalization %>%
                recover_types(df_variants) %>%
                spread_draws(hosp_rate_variant[variant]),
              by = c(".draw","variant")) %>%
    left_join(df_variants, by = c("date","variant")) %>%
    group_by(.draw,date,age_group,municipality) %>%
    summarize(hosp_rate_variant = sum(hosp_rate_variant * percentage_variant * hosp_rate_vaccins), 
              .groups = "drop") %>%
    # Finally compute the expected hospitalizations per municipality
    left_join(fit_hospitalization %>%
                recover_types(df_municipality) %>%
                spread_draws(hosp_rate[age_group,municipality]),
              by = c(".draw","age_group","municipality")) %>%
    left_join(df_municipality, by = c("municipality","age_group","date")) %>%
    mutate(expected_hospitalizations = hosp_rate * hosp_rate_variant * 10^(load - ref_load) * population) %>%
    select(.draw,date,municipality,age_group,load,hospitalizations,expected_hospitalizations,population) %>%
    # Compute the total for the municipality
    bind_rows(group_by(.,.draw,date,municipality) %>%
                summarize(across(contains("hosp"),sum),
                          population = sum(population),
                          load = first(load),
                          age_group = "Totaal",
                          .groups = "drop")) %>%
    # Perform a randow draw for 95%-intervals
    mutate(simulated_hospitalizations = rpois(n(),expected_hospitalizations)) %>%
    group_by(age_group,date,municipality,load,hospitalizations) %>%
    median_qi(expected_hospitalizations,simulated_hospitalizations) %>%
    # Write as a csv
    write.csv(here(runname,"output","municipality",paste0("hosp_",name,".csv")),
              row.names = F)
  
  gc()
}, .options = furrr_options(seed = T))

# Determine the outcome for the Netherlands

df_posteriors_hosp <- fit_hospitalization %>%
  recover_types(df_municipality,df_variants,df_vaccination) %>%
  spread_draws(hosp_rate[age_group,municipality],
               hosp_rate_variant[variant],
               prevention_vax[age_group,vaccination,variant]) %>%
  mutate(hosp_rate = hosp_rate * hosp_rate_variant) %>%
  select( - hosp_rate_variant)

df_nederland <- list(NULL)

# Compute the number of hospitalisation per day for each chain seperately
for(i in 1:10){
  
  df_nederland[[i]] <- df_posteriors_hosp %>% 
    filter(.chain == i) %>%
    group_by(.draw) %>% 
    group_split() %>%
    future_map(function(x){
      left_join(x, df_municipality, by = c("municipality", "age_group")) %>%
        left_join(df_variants, by = c("date","variant")) %>%
        left_join(df_vaccination, by = c("date", "municipality","age_group","vaccination")) %>%
        # Calculate the effect of vaccination
        group_by(date,municipality,age_group,.draw,variant) %>%
        summarize(hospitalizations = first(hospitalizations),
                  population= first(population),
                  load = first(load),
                  hosp_rate = first(hosp_rate) * sum(prevention_vax * percentage_vaccination),
                  percentage_variant = first(percentage_variant),
                  .groups = "drop_last") %>%
        summarize(hospitalizations = first(hospitalizations),
                  population= first(population),
                  load = first(load),
                  hosp_rate = sum( hosp_rate * percentage_variant),
                  .groups = "drop") %>%
        mutate(expected_hospitalizations = hosp_rate * 10^(load - ref_load) * population) %>%
        group_by(date) %>%
        summarize(hospitalizations = sum(hospitalizations),
                  expected_hospitalizations = sum(expected_hospitalizations),
                  load = log10( sum(10^load*population)/sum(population)))
    }) %>%
    bind_rows()
  
}

df_nederland <- bind_rows(df_nederland)

save(df_nederland, file = here(runname,"output","Netherlands","df_netherlands_hosp.RData"))