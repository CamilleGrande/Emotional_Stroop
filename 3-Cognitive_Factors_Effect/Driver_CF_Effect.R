## Driver_CF_Effect.R
## Author: Camille Grandé
## Description:
##      This script runs the cognition CFA analysis in the TAY sample and the modeling 
##      part for the EST task with the cognitive factors as covariates

## Wipe and recreate the Outputs folder so old files don't interfere with this run
    unlink("./3-Cognitive_Factors_Effect/Outputs", recursive = TRUE)
    dir.create("./3-Cognitive_Factors_Effect/Outputs", recursive = TRUE)

## Libraries needed to run the script
    library(tidyverse)
    library(lavaan)
    library(lme4)
    library(DHARMa)
    library(performance)
    library(emmeans)
    library(afex)
    library(effectsize)

## Utility scripts
    source("./3-Cognitive_Factors_Effect/Utilities_CFA.R")
    source("./3-Cognitive_Factors_Effect/Utilities_CF_Effect.R")

## Paths
    EST.data = "./0-Data_Prep/Outputs/Tidy_all.rds"
    cognition.data = "./0-Data_Prep/Outputs/Tidy_cognition.rds"
    FactorAnalysis.output = "./3-Cognitive_Factors_Effect/Outputs/Factor_Scores.csv"
    model.rt = "./3-Cognitive_Factors_Effect/Outputs/RT_Model_results.rds"
    model.rt.f2 = "./3-Cognitive_Factors_Effect/Outputs/f2_RT_Model_results.rds"
    model.acc = "./3-Cognitive_Factors_Effect/Outputs/Acc_Model_results.rds"
    model.acc.f3 = "./3-Cognitive_Factors_Effect/Outputs/f3_Acc_Model_results.rds"


## Confirmatory Factor Analysis on whole sample
## Output is a csv file in Outputs folder
    CFA.Analysis(cognition.data)

## Modelling RT and Accuracy with factors from CFA

## --- RT ---
    ## Models RT and computes ANOVA
    ANOVA.RT(EST.data, FactorAnalysis.output)
    ## Computes % change for each effect/interaction, controlling for cognitive abilities
    Percent.Change.RT(model.rt)
    ## Computes % change for each effect/interaction, this time 
    ## adding interaction term for f2
    Percent.Change.RT(model.rt.f2, f2.interaction = TRUE)
    ## Plot interaction of f2 with frequency on RT
    Plot.RT.f2(model.rt.f2)

## --- Accuracy ---
    ## Models Accuracy and computes Likelihood Ratio Test
    LRT.Accuracy(EST.data, FactorAnalysis.output)
    ## Computes OR for each effect/interaction, controlling for cognitive abilities
    Odds.Ratio.Acc(model.acc)
    ## Computes OR for each effect/interaction, this time
    ## adding interaction term for f3
    Odds.Ratio.Acc(model.acc.f3, f3.interaction = TRUE)
    ## Plot interaction of f3 with congruency and frequency on accuracy
    Plot.Acc.f3(model.acc.f3)

## --- Supplements ---
    ## To test if results with random intercept only are similar to random slope + intercept
    ## Here, no test on accuracy expended model as does not converge, cannot trust results
    Model.RT.RandomSlope(EST.data, FactorAnalysis.output)