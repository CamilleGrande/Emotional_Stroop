## Driver_Clinical_Prediction.R
## Author: Camille Grandé
## Description:
##      This script runs the prediction of task performance using 
##      clinical metrics (syndrome scales and suicide attempts)

## Libraries needed to run the script
    library(tidyverse)
    library(lme4)
    library(lmerTest)
    library(broom.mixed)
    library(emmeans)

## Utility script
    source("./4-Clinical_Sympt_Prediction/Utilities_Clinical_Prediction.R")

## Paths
    EST.data = "./0-Data_Prep/Outputs/Tidy_all.rds"
    clinical.data = "./0-Data_Prep/Outputs/Tidy_clinical.rds"
    suicide.data = "./0-Data_Prep/Outputs/Tidy_suicide_att.rds"

## Model Congruent trials only
## Answers the question: does X clinical subscale predicts accuracy/RT
##                          on congruent trials?
Congruent.RT(EST.data, clinical.data, suicide.data)
Congruent.Accuracy(EST.data, clinical.data, suicide.data)

## Models Incongruency effect, face emotion effect, and their interaction,
## i.e. incongruency effect depending on face emotion (e.g. higher incongruency
## effect for angry faces)
## Answers the question: does X clinical subscale predicts accuracy/RT
##                          difference on cong vs incong trials, happy vs angry faces,
##                          and incongruency effect depending on face emotion
Interference.FaceEmotion.RT(EST.data, clinical.data, suicide.data)
Interference.FaceEmotion.Accuracy(EST.data, clinical.data, suicide.data)