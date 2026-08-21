## Driver_Task_Cond_Effect.R
## Author: Camille Grandé
## Description:
##      This script runs the stats to test the task conditions effect for the emotional stroop task analysis


## To do:
##    Influential points for loop

library(tidyverse)
library(lmerTest)
library(effectsize)
library(DHARMa)
library(performance)


source("./2-Task_Cond_Effect/Utilities_Task_Cond_Effect.R")

data.all = "./0-Data_Prep/Outputs/Tidy_all.rds"

ANOVA.RT(data.all)
ChiSq.Accuracy(data.all)