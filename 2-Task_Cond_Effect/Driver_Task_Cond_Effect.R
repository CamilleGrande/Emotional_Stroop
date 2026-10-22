## Driver_Task_Cond_Effect.R
## Author: Camille Grandé
## Description:
##      This script runs the stats to test the task conditions effect for the emotional stroop task analysis


## Libraries needed to run the script
    library(tidyverse)
    library(lme4)
    library(DHARMa)
    library(performance)
    library(emmeans)
    library(ggeffects)

## Utility script
    source("./2-Task_Cond_Effect/Utilities_Task_Cond_Effect.R")

## Paths
    data.all = "./0-Data_Prep/Outputs/Tidy_all.rds"
    ## These files will get created when running this script
    model.rt = "./2-Task_Cond_Effect/Outputs/RT_Model_results.rds"
    model.acc = "./2-Task_Cond_Effect/Outputs/Acc_Model_results.rds"

## Main analysis
    ANOVA.RT(data.all)
    Percent.Change.RT(model.rt)
    Plots.RT(model.rt, data.all)
    LRT.Accuracy(data.all)
    Odds.Ratio.Acc(model.acc)
    Plots.Accuracy(model.acc, data.all)

## Supplements
    Model.RT.RandomSlope(data.all)
    Model.Acc.RandomSlope(data.all)