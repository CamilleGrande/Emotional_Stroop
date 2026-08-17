## Driver_Descriptives.R
## Author: Camille Grandé
## Description:
##      This script runs the descriptives stats for the emotional stroop task analysis


## Libraries needed to run the script
    library(tidyverse)
    library(parameters)

## Utility script
    source("./1-Descriptives/Utilities_Descriptives.R")
    
## Paths
    EST.Tidy = "./0-Data_Prep/Outputs/Tidy_EST.rds"
    Demog.Tidy = "./0-Data_Prep/Outputs/Tidy_demog.rds"


EST.Descriptives(EST.Tidy)
Demographics.Descriptives(Demog.Tidy)