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
    Data.Tidy = "./0-Data_Prep/Outputs/Tidy_all.rds"

EST.Descriptives(Data.Tidy)
Demographics.Descriptives(Data.Tidy)