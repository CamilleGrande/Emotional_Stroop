## Driver_DataPrep.R
## Author: Camille Grandé
## Description:
##      This script preps the data for the emotional stroop task analysis


## Libraries needed to run the script
    library(tidyverse)

## Utility script
    source("./0-Data_Prep/Utilities_DataPrep.R")
    
## Paths
    EST.Raw = "./Working_Data/data/cognitive/emotional_stroop_raw_score.csv"

## Clean data
    Tidy.Data(EST.Raw)