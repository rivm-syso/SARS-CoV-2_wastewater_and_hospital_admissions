# Needed for finding the most recent datafiles from EPI
library(tidyverse)
library(here)
library(furrr)

# Under this subfolder, all relevant files will be written
runname <- "results/The_Netherlands_2020_2022"
dir.create(here(runname), showWarnings = FALSE, recursive = TRUE ) # allows all runs to be in a single subdirectory if runname consists of stacked folders

startday <- as.Date("2020-08-01")
lastday  <- as.Date("2022-02-08")

# File names for viral load, hospital, vaccination, and variant data
vaccin_filename <- "Vaccination_coverage.RData"
viralload_hospital_filename <- "Hospitalisations_viralload.RData"
variant_filename <- "Variant.RData"

# Create folders to save outcomes
dir.create.if.needed <- function(x){
  if( !dir.exists(x)) dir.create(x)
}

# directories
dir.create.if.needed(here(runname, "model_code_backup"))
file.copy( from = list.files( pattern = "^[0-9].*R|^.*stan"), 
           to   = here(runname, "model_code_backup"))

walk( c("figures", "output"),
      function(x){
        outdir <- here(runname, x)
        dir.create.if.needed(outdir)
        dir.create.if.needed(here(outdir, "model_data" ))
        dir.create.if.needed(here(outdir, "Netherlands" ))
        dir.create.if.needed(here(outdir, "municipality" ))
        dir.create.if.needed(here(outdir, "manuscript" ))})

# Settings for hospitalizations-model
ref_load      <- 19
variant_names <- c("Wildtype","Alpha","Delta","Omicron")
vaccination_status <- c("Unvaccinated","Partially","Fully","Boosted")

# Number of cores for parallelization
plan(multisession, workers = 10)

settings_sourced <- T