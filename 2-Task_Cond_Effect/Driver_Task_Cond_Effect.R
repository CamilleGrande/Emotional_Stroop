## Driver_Task_Cond_Effect.R
## Author: Camille Grandé
## Description:
##      This script runs the stats to test the task conditions effect for the emotional stroop task analysis


## To do:
##    Influential points for loop
##    test percentage accuracy
##      still ORs or not needed? OR for RT too?
##      coefficients RT?

library(tidyverse)
library(lmerTest)
library(effectsize)
library(DHARMa)
library(performance)
library(broom.mixed)
library(afex)
library(paletteer)
## might need lme4


source("./2-Task_Cond_Effect/Utilities_Task_Cond_Effect.R")

data.all = "./0-Data_Prep/Outputs/Tidy_all.rds"
## These files will get created when running this script
model.rt = "./2-Task_Cond_Effect/Outputs/RT_Final_Model_results.rds"

ANOVA.RT(data.all)
Percent.Change.RT(model.rt)
LRT.Accuracy(data.all)