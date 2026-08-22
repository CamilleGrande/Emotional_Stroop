## Utilities_CFA.R
## Author: Camille Grandé
## Description:
##      This script runs the cognition CFA analysis in the TAY sample and the modeling 
##      part for the EST task with the cognitive factors as covariates


## Libraries needed to run the script
    library(tidyverse)
    library(lavaan)
    library(lmerTest)
    library(effectsize)
    library(DHARMa)
    library(performance)

## Utility script
    source("./3-Cognitive_Factors_Effect/Utilities_CFA.R")
    source("./3-Cognitive_Factors_Effect/Utilities_CF_Effect.R")

## Paths
    EST.data = "./0-Data_Prep/Outputs/Tidy_all.rds"
    cognition.data = "./0-Data_Prep/Outputs/Tidy_cognition.rds"
    FactorAnalysis.output = "./3-Cognitive_Factors_Effect/Outputs/Factor_Scores.csv"


## Confirmatory Factor Analysis on whole sample
## Output is a csv file in Outputs folder
    CFA.Analysis(cognition.data)

## Modelling RT and Accuracy with factors from CFA
    ANOVA.RT(EST.data, FactorAnalysis.output)
    ChiSq.Accuracy(EST.data, FactorAnalysis.output)