load_if_needed <- function( object, filename ){
# Object becomes a list in case we need to load multiple file from the same file
    if( object %>% 
      lapply(function(object){
        length(ls( pattern=object, envir = .GlobalEnv)) == 0}) %>%
      c(recursive = T) %>% any())
    load( filename, envir = .GlobalEnv )
}

# initialize variables
initials_hosp = function() {
  
  # We have 5 age groups, 3 different vaccination statuses besides unvaccinated,
  # and 4 different variants including wildtype for the effect of vaccination
  prevention_vax <- array(.3, dim = c(5,3,4))

  hosp_rate_small <- array(c(7.5/4, 7.5/2, 7.5, 7.5*2, 7.5*4))
  
  sigma_hosp <- array(10, 5)
  
  # The hosp rate is the same for each age group across the municipalities 
  hosp_rate <- array(hosp_rate_small, dim = c(5, df_municipality %>% 
                                                pull(municipality) %>% 
                                                nlevels()))
  
  return(list(
    sigma_hosp_rate = sigma_hosp,
    mean_hosp_rate = hosp_rate_small,
    hosp_rate = hosp_rate,
    prevention_vax = prevention_vax
  ))
}

# colorblind-friendly scheme
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 


# Additional check to ensure that the variant, vaccination, hospitalisation, and
# viral load data is all ordered in the same way without creating unnecessary large
# data frames.
create_stan_data <- function(df_muni, df_variants, df_vaccins){
  
  # We arrange all data frames in the same way, namely:
    # Chronological
    # Municipality on the same day
    # Age groups within municipality
    # vaccination status
    
  # Order vaccination data
  df_vaccins <- df_vaccins %>% 
    # Furthermore, we only consider the dates which are present in df_muni
    filter(date %in% (df_muni %>% pull(date) %>% unique())) %>%
    mutate(municipality = as.factor(municipality),
           age_group = as.factor(age_group)) %>%
    arrange(date,municipality,age_group,vaccination)
  
  for(counter in 1:nlevels(pull(df_vaccins,vaccination))){
    
    # As df_vaccins has an extra level (vaccination) we only select the entries
    # which should be the same vaccination level.
    # We then check if these are ordered the same way as df_muni and if they are
    # indeed all the same vaccination level.
    size <- seq(from = counter,
                to = nlevels(pull(df_vaccins,vaccination)) * nrow(df_muni),
                by = nlevels(pull(df_vaccins,vaccination)))
    # We check if the order of df_vaccins and df_muni match
    if(any( (df_vaccins %>% select(date,municipality,age_group) %>% .[size,]) !=
            (df_muni %>% select(date,municipality,age_group)))){
      stop("The rows of df_muni and df_vaccins are not ordered correctly")
    }
    
    if(any( (df_vaccins %>% pull(vaccination) %>% .[size]) != 
            (df_vaccins %>% pull(vaccination) %>% levels() %>% .[counter]))){
      stop("The vaccinations are not ordered correctly in df_vaccins")
    }
    
  }
  
  # We do the same thing for the variant data
  # As we do not have municipality nor age groups in the variant data, we manually
  # repeat the variant data in order to "arrange" by municipality and age.
  df_variants <- df_variants %>% 
    filter(between(date,df_muni %>% pull(date) %>% min(),
                   df_muni %>% pull(date) %>% max())) %>%
    group_by(variant,date) %>%
    summarize(percentage_variant = rep(percentage_variant, 
                                    df_muni %>% pull(age_group) %>% factor() %>% nlevels() *
                                      df_muni %>% pull(municipality) %>% factor() %>% nlevels()),
              .groups = "drop") %>%
    arrange(variant,date) 
  
  for(counter in 1:nlevels(pull(df_variants,variant))){
    
    size <- (counter-1)*nrow(df_muni) + (1:nrow(df_muni))
    
    if(any( (df_variants %>% select(date) %>% .[size,]) !=
            (df_muni %>% select(date)))){
      stop("The rows of df_muni and df_variants are not ordered correctly")
    }
    
    if(any( (df_variants %>% pull(variant) %>% .[size]) != 
            (df_variants %>% pull(variant) %>% levels() %>% .[counter]))){
      stop("The variants are not ordered correctly in df_variants")
    }
    
  }
  
  # Once the different data frames are all ordered in the same way, we combine
  # them to create a single object for stan
  stan_data <- c(compose_data(df_muni),
                 compose_data(df_vaccins %>% select(vaccination, percentage_vaccination)) %>%
                   within(rm("n")),
                 compose_data(df_variants %>% select(variant, percentage_variant)) %>%
                   within(rm("n")))
  
  return(stan_data)
  
}


functions_sourced <- T
