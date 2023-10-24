# Inferring Hospital Admissions from SARS-CoV-2 Virus Loads in Wastewater in the Netherlands, August 2020 – February 2022
[Wouter A. Hetebrij](mailto://wouter.hetebrij@rivm.nl), Ana Maria de Roda Husman, Erwin Nagelkerke, Rudolf F.H.J. van der Beek, Senna C.J.L. van Iersel, Titus G.V. Breuning, Willemijn J. Lodder, Michiel van Boven
Under review by Science of the Total Enviroment

# How to use with own data

To use the scripts with your own data sets, the following three data sets are needed:

  - _df_variants_ with columns:
      -  date: The date
      -  variant: The variant name, which should include for all variants in _variant_names_ in 0. settings, and only include those variants.
      -  percentage_variant: The relative presence on day _x_ of variant _y_ and should lie between 0 and 1.
  
  - _df_vaccination_ with columns:
      - date: The date
      - municipality: The municipalities which are included in the model. (Should be the same municipalities as in _df_municipality_).
      - age_group: The different age groups which are considered. (Should be the age groups as in _df_municipality_). If not 5, then _initials_hosp_ in 0. function should be updated to reflect the number of age groups.
      - vaccination: The different vaccination statuses included in the model, which should include the statuses in _vaccination_status_ in 0. settings.
      - percentage_vaccination: The vaccination coverage on day _x-t_ of vaccination status _y_, with _t_ the time it takes for a vaccine to become fully effective after administration. Vaccination coverage should lie between 0 and 1 and for a given day, municipality, and age group it should add to 1.
    
  - _df_municipality_ with columns:
      - date: The date
      - municipality: The municipalities which are included in the model. (Should be the same municipalities as in _df_vaccination_).
      - age_group: The different age groups which are considered. (Should be the age groups as in _df_vaccination_). If not 5, then _initials_hosp_ in 0. function should be updated to reflect the number of age groups.
      - load: The 10-log of the virus load at day _x_ in municipality _y_
      - hospitalizations: The number of hospitalisations at day _x_ in municipality _y_ of age group _z_. Should be an integer greater than or equal to 0.
      - population: The population at day _x_ of municipality _y_ of age group _z_, which should be at least 1.

  # Data availability
Contact [Wouter Hetebrij](mailto://wouter.hetebrij@rivm.nl) for questions regarding the data set used in the manuscript.
