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
    Demographics.Raw = "./Working_Data/datacut_demographics_core_data.csv"
    NIH.Raw = "./Working_Data/data/cognitive/nih_toolbox_assessment_scores.csv"
    RAVLT.Raw = "./Working_Data/data/cognitive/ravlt_item_level_score.csv"
    Matrices.Raw = "./Working_Data/data/cognitive/ravlt_subtest_data.csv"
    LMT.Raw = "./Working_Data/data/cognitive/little_man_task_summary_score.csv"
    CSSRS.Raw = "./Working_Data/data/cssrs_baseline_data.csv"
    ASR.Raw = "./Working_Data/data/aseba/aseba_adult_selfreport_asr_data.csv"
    ASR.GI.Raw = "./Working_Data/data/aseba/aseba_adult_selfreport_asr_gender_inclusive_tscore_data.csv"
    YSR.Raw = "./Working_Data/data/aseba/aseba_youth_selfreport_ysr_data.csv"
    YSR.GI.Raw = "./Working_Data/data/aseba/aseba_youth_selfreport_ysr_gender_inclusive_tscore_data.csv"


## Clean data EST and demographics
    Tidy.Data(EST.Raw, Demographics.Raw)

## Clean data for Cognitive Factor analysis
    Tidy.Cognition(NIH.Raw, RAVLT.Raw, Matrices.Raw, LMT.Raw)

## Clean clinical data
    Tidy.Clinical(CSSRS.Raw, ASR.Raw, ASR.GI.Raw, YSR.Raw, YSR.GI.Raw)