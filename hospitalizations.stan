/* Stan program to determine the relation between viral loads of          */ 
/* SARS-CoV-2 in wastewater and hospitalisations on the municipality      */
/* level. Both the influence of variants and vaccination on the relation  */
/* is taken into account.                                                 */
/* (c) Wouter Hetebrij, licensed with BSD3 clause (use but mention)       */
/* Created: 2021                                                          */

data {
  int<lower=1> n;
  int<lower = 1> n_municipality;  // number of regions (province, safety region, municipality)
  int<lower = 1> n_age_group;     // number of different age-groups  
  
  int<lower = 1> n_variant;       // number of different variants
  int<lower = 1> n_vaccination;   // number of different vaccinations

  /* Classifications of the data */
  int<lower=1, upper=n_municipality> municipality[n];
  int<lower=1, upper=n_age_group> age_group[n];

  /* For vector-based calculations, we work exclusively with vectors */
  vector<lower=1>[n] population;
  vector<lower=0>[n] load;
  
  row_vector<lower=0,upper=1>[n * n_vaccination] percentage_vaccination;
  vector<lower=0,upper=1>[n * n_variant] percentage_variant;
  
  int hospitalizations[n];        // now included at the level of sewage plants

  int ref_load;                   // reference sewage load for hospitalization rates
}

transformed data {
  
  /* The data with the different variants and vaccinations will be made a matrix */
  matrix<lower=0,upper=1>[n,n_variant] percentage_variant_matrix;
  matrix<lower=0,upper=1>[n,n_vaccination * n_age_group] percentage_vaccination_matrix = 
                                        rep_matrix(0, n, n_vaccination * n_age_group);

  /* Replace the vector with the ratios of VOC with a matrix where each column 
  is a single VOC */
  for(i in 1:n_variant){
    percentage_variant_matrix[,i] = percentage_variant[((i-1)*n + 1):((i-1)*n + n)];
  }
  
  /* Replace the vector with the different vaccinations with a blockmatrix. In
  essence we wrap the vaccination row vector as a matrix. However, depending
  on the correspondent each group, we stagger the row. */
  for(i in 1:n){
    
    percentage_vaccination_matrix[i,((age_group[i]-1)*n_vaccination + 1):
                              ((age_group[i]-1)*n_vaccination + n_vaccination)] =
                          percentage_vaccination[((i-1)*n_vaccination + 1):
                                          ((i-1)*n_vaccination + n_vaccination)];
    
  }
  
}

parameters {
  
  /* hospitalization rates and hospitalization delay */
  vector<lower = 0>[n_age_group] mean_hosp_rate;             // hospitalization rates with random effect 
  vector<lower = 0>[n_age_group] sigma_hosp_rate;
  
  matrix<lower = 0, upper = 1>[n_vaccination-1,n_variant] prevention_vax_small[n_age_group];
  
  matrix<lower = 0>[n_age_group,n_municipality] hosp_rate;
  
  vector<lower = 0>[n_variant-1] hosp_rate_variant_small;
  
}

transformed parameters {
  
  matrix<lower = 0, upper = 1>[n_vaccination,n_variant] prevention_vax[n_age_group];
  vector<lower = 0>[n_variant] hosp_rate_variant;

  /*  The baseline is unvaccinated and the first variant */
  prevention_vax[,1,] = rep_array(rep_row_vector(1,n_variant),n_age_group);
  hosp_rate_variant[1] = 1;
  
  /* Add the other vaccination statuses */
  prevention_vax[,2:n_vaccination,] = prevention_vax_small;
  hosp_rate_variant[2:n_variant] = hosp_rate_variant_small;
  
}

model {
  vector[n] sum_load = exp( log(10) *(load - ref_load));
  vector[n] hosp_parameter;
  
  /* For a slight increase in speed, we want to use vector-valued calculations */
  vector[n] hosp_rate_matrix;
  matrix[n_vaccination * n_age_group, n_variant] prevention_vax_matrix;
  
  /* Priors on the effectiveness of vaccination. To make sure the problem is 
  identifiable, we ensure that in a fully unvaccinated population the hosp rate
  is the coefficient between load and hospitalizations */
  for(i in 1:n_vaccination){
    for(j in 1:n_variant){
      if(i > 1){
        prevention_vax[,i,j] ~ beta(2,8);
      }
    }
  }
  
  /* To avoid that sigma_hosp_rate trails to infinity, we put a slightly restrictive
  prior on sigma_hosp_rate. As large values of sigma_hosp_rate mean that all hosp_rates
  will be more or less the same, this shouldn't effect the outcome. */

  sigma_hosp_rate ~ gamma(10 , 1);
  
  
  /* The hospitalizations ratios for each municipality is governed by the overall
  hosp_rate for the concerning municipality */
  for(i in 1:n_age_group){
      hosp_rate[i,] ~ gamma(sigma_hosp_rate[i] * mean_hosp_rate[i], sigma_hosp_rate[i]);
  }
  
  /* We fill out the matrices with the hosp-rates for the different variants*/
  for(i in 1:n){
    hosp_rate_matrix[i] = hosp_rate[age_group[i],municipality[i]];
  }


  /* For each variant, we want to sum the effectiveness of the different vaccinations,
   in other words, we want to multiply the prevention with the percentage of the vaccinations.
   Due to the different age_groups, we do this seperately for each age group */
  for(i in 1:n_age_group){
        prevention_vax_matrix[((i-1)*n_vaccination + 1):((i-1)*n_vaccination + n_vaccination)] =
                                                                                prevention_vax[i,,];
  }
  

  /* The chance of hospitalization is load (variant -independent )
                                      hosp_rate (variant-dependent)
                                      percentage_vaccination (variant-independent)
                                      prevention of vaccination (variant-dependent)
                                      population (variant-independent)
                                      percentage varianten ( variant-dependent) */
  hosp_parameter = (percentage_variant_matrix .*
                      (percentage_vaccination_matrix * prevention_vax_matrix)) * hosp_rate_variant;
                      
  hosp_parameter =  sum_load .* hosp_rate_matrix .* hosp_parameter .* population;

  hospitalizations ~ poisson(hosp_parameter);

}
