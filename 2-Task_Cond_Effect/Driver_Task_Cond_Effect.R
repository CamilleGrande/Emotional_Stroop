## Driver_Task_Cond_Effect.R
## Author: Camille Grandé
## Description:
##      This script runs the stats to test the task conditions effect for the emotional stroop task analysis

## Wipe and recreate the Outputs folder so old files don't interfere with this run
    unlink("./2-Task_Cond_Effect/Outputs", recursive = TRUE)
    dir.create("./2-Task_Cond_Effect/Outputs", recursive = TRUE)

## Libraries needed to run the script
    library(tidyverse)
    library(lme4)
    library(DHARMa)
    library(performance)
    library(emmeans)
    library(afex)
    library(effectsize)

## Utility script
    source("./2-Task_Cond_Effect/Utilities_Task_Cond_Effect.R")

## Paths
    data.all = "./0-Data_Prep/Outputs/Tidy_all.rds"
    ## These files will get created when running this script
    model.rt = "./2-Task_Cond_Effect/Outputs/RT_Model_results.rds"
    model.acc = "./2-Task_Cond_Effect/Outputs/Acc_Model_results.rds"
    model.acc.age = "./2-Task_Cond_Effect/Outputs/Age_Acc_Model_results.rds"

## Main analysis

## --- RT ---
    ## Models RT and computes ANOVA
    ANOVA.RT(data.all)
    ## Computes % change for each effect/interaction
    Percent.Change.RT(model.rt)
    ## Plots the three way interaction and the age interaction plots
    Plots.RT(model.rt, data.all)

## --- Accuracy ---
    ## Models Accuracy and computes Likelihood Ratio Test
    LRT.Accuracy(data.all)
    ## Computes OR for each effect/interaction
    Odds.Ratio.Acc(model.acc)
    ## Plots % predicted accuracy for congruence and age
    Plots.Accuracy(model.acc, data.all)

## --- Supplements ---
    ## To test if results with random intercept only are similar to random slope + intercept
    Model.RT.RandomSlope(data.all)
    Model.Acc.RandomSlope(data.all)

    ## Model RT interaction between age and task conditions
    Model.RT.Interaction.Age(data.all)
    ## Model Acc interaction between age and task conditions
    Model.Acc.Interaction.Age(data.all)
    ## Plot modulation frequency effect by age
    Plot.Acc.Interaction.Age(data.all, model.acc.age)

    ## Model RT interaction between sex and task conditions
    Model.RT.Interaction.Sex(data.all)
    ## Model Acc interaction between sex and task conditions
    Model.Acc.Interaction.Sex(data.all)