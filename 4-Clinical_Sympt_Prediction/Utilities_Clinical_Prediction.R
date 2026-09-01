## Utilities_Clinical_Prediction.R
## Author: Camille Grandé
## Description:
##      This script holds the utility function to run the task performance prediction
##      using clinical metrics


## Helper function to create one big summary of our models
Create.Summary.Models <- function(list.models) {

    results <- imap_dfr(
                list.models, ~ tidy(.x, effects = "fixed") |>
                    mutate(model = .y) |>
                        relocate(model) |>
                            ## Adjusting with FDR because clinical constructs are correlated
                            mutate(p_fdr = p.adjust(p.value, method = "fdr"))
                )
    
    return(results)
}


## Runs the analysis for the models on interference effect, face emotion effect
## and their interaction for RT
Interference.FaceEmotion.RT <- function(EST.data, clinical.data, suicide.data) {

    EST <- read_rds(EST.data) |>
                ## Keeping only RT for correct trials
                subset(correct == "correct") 
    
    EST$log_rt <- log(EST$reaction_time)
    
    clinical <- read_rds(clinical.data)

    suicide <- read_rds(suicide.data)

    ## ----- Suicide analysis -----
    ## don't bind all 3 (EST, clinical, suicide) together bc would loose a lot of suicide data
    ## EST + clinical: n = 690
    ## EST + suicide: n = 640
    ## ptcp w EST + clinical + suicide data: n = 616
    dat_suicide <- inner_join(EST, suicide, by = "SID")

    model_suicide <- lmer(log_rt ~ congruence * face_emotion * cssrs_attempt_lt + log_age + sex + ethnicity 
                    + (1 | SID),
                    data = dat_suicide, REML = T)

    ## ----- Clinical sympt analysis -----
    dat_clinical <- inner_join(EST, clinical, by = "SID")

    model_anxdepr <- lmer(log_rt ~ congruence * face_emotion * anxdepr_t + log_age + sex + ethnicity 
                    + (1 | SID),
                    data = dat_clinical, REML = T)

    model_withdrdepr <- lmer(log_rt ~ congruence * face_emotion * withdrdepr_t + log_age + sex + ethnicity 
                    + (1 | SID),
                    data = dat_clinical, REML = T)

    model_somatic <- lmer(log_rt ~ congruence * face_emotion * somatic_t + log_age + sex + ethnicity 
                    + (1 | SID),
                    data = dat_clinical, REML = T)

    model_thought <- lmer(log_rt ~ congruence * face_emotion * thought_t + log_age + sex + ethnicity 
                    + (1 | SID),
                    data = dat_clinical, REML = T)
    
    model_attention <- lmer(log_rt ~ congruence * face_emotion * attention_t + log_age + sex + ethnicity 
                    + (1 | SID),
                    data = dat_clinical, REML = T)

    model_aggressive <- lmer(log_rt ~ congruence * face_emotion * aggressive_t + log_age + sex + ethnicity 
                    + (1 | SID),
                    data = dat_clinical, REML = T)

    model_rulebreaking <- lmer(log_rt ~ congruence * face_emotion * rulebreaking_t + log_age + sex + ethnicity 
                    + (1 | SID),
                    data = dat_clinical, REML = T)

    all_models <- list(
                    suicide = model_suicide,
                    anxdepr = model_anxdepr,
                    withdrdepr = model_withdrdepr,
                    somatic = model_somatic,
                    thought = model_thought,
                    attention = model_attention,
                    aggressive = model_aggressive,
                    rulebreaking = model_rulebreaking
                )

    results <- Create.Summary.Models(all_models)
    write_csv(results, "./4-Clinical_Sympt_Prediction/Outputs/RT_Interference-FaceEmotion_Models.csv")
}


## Runs the analysis for the models on interference effect, face emotion effect
## and their interaction for accuracy
Interference.FaceEmotion.Accuracy <- function(EST.data, clinical.data, suicide.data) {

    EST <- read_rds(EST.data) |> 
                mutate(correct = case_when(
                    correct == "correct"   ~ 1,
                    correct == "incorrect" ~ 0
                ), .keep = "unused")

    clinical <- read_rds(clinical.data)

    suicide <- read_rds(suicide.data)

    ## ----- Suicide analysis -----
    ## don't bind all 3 (EST, clinical, suicide) together bc would loose a lot of suicide data
    ## EST + clinical: n = 690
    ## EST + suicide: n = 640
    ## ptcp w EST + clinical + suicide data: n = 616
    dat_suicide <- inner_join(EST, suicide, by = "SID")
    
    model_suicide <- glmer(correct ~ congruence * face_emotion * cssrs_attempt_lt + log_age + sex + ethnicity
                    + (1 | SID),
                    data = dat_suicide, family = binomial(link = "logit"))
    
    ## ----- Clinical sympt analysis -----
    dat_clinical <- inner_join(EST, clinical, by = "SID")

    model_anxdepr <- glmer(correct ~ congruence * face_emotion * anxdepr_t + log_age + sex + ethnicity
                    + (1 | SID),
                    data = dat_clinical, family = binomial(link = "logit"))

    model_withdrdepr <- glmer(correct ~ congruence * face_emotion * withdrdepr_t + log_age + sex + ethnicity
                    + (1 | SID),
                    data = dat_clinical, family = binomial(link = "logit"))

    model_somatic <- glmer(correct ~ congruence * face_emotion * somatic_t + log_age + sex + ethnicity
                    + (1 | SID),
                    data = dat_clinical, family = binomial(link = "logit"))

    model_thought <- glmer(correct ~ congruence * face_emotion * thought_t + log_age + sex + ethnicity
                    + (1 | SID),
                    data = dat_clinical, family = binomial(link = "logit"))

    model_attention <- glmer(correct ~ congruence * face_emotion * attention_t + log_age + sex + ethnicity
                    + (1 | SID),
                    data = dat_clinical, family = binomial(link = "logit"))

    model_aggressive <- glmer(correct ~ congruence * face_emotion * aggressive_t + log_age + sex + ethnicity
                    + (1 | SID),
                    data = dat_clinical, family = binomial(link = "logit"))
    
    model_rulebreaking <- glmer(correct ~ congruence * face_emotion * rulebreaking_t + log_age + sex + ethnicity
                    + (1 | SID),
                    data = dat_clinical, family = binomial(link = "logit"))

    all_models <- list(
                    suicide = model_suicide,
                    anxdepr = model_anxdepr,
                    withdrdepr = model_withdrdepr,
                    somatic = model_somatic,
                    thought = model_thought,
                    attention = model_attention,
                    aggressive = model_aggressive,
                    rulebreaking = model_rulebreaking
                )

    results <- Create.Summary.Models(all_models)
    write_csv(results, "./4-Clinical_Sympt_Prediction/Outputs/Accuracy_Interference-FaceEmotion_Models.csv")
}


## Runs the analysis for the models on congruent trials only
## for RT
Congruent.RT <- function(EST.data, clinical.data, suicide.data) {

    EST <- read_rds(EST.data) |>
                ## Keeping only congruent trials
                subset(congruence == "congruent") |>
                    ## Keeping only RT for correct trials
                    subset(correct == "correct") 
    
    EST$log_rt <- log(EST$reaction_time)
    
    clinical <- read_rds(clinical.data)

    suicide <- read_rds(suicide.data)

    ## ----- Suicide analysis -----
    ## don't bind all 3 (EST, clinical, suicide) together bc would loose a lot of suicide data
    ## EST + clinical: n = 690
    ## EST + suicide: n = 640
    ## ptcp w EST + clinical + suicide data: n = 616
    dat_suicide <- inner_join(EST, suicide, by = "SID")

    model_suicide <- lmer(log_rt ~ cssrs_attempt_lt + log_age + sex + ethnicity 
                    + (1 | SID),
                    data = dat_suicide, REML = T)

    ## ----- Clinical sympt analysis -----
    dat_clinical <- inner_join(EST, clinical, by = "SID")

    model_anxdepr <- lmer(log_rt ~ anxdepr_t + log_age + sex + ethnicity 
                    + (1 | SID),
                    data = dat_clinical, REML = T)

    model_withdrdepr <- lmer(log_rt ~ withdrdepr_t + log_age + sex + ethnicity 
                    + (1 | SID),
                    data = dat_clinical, REML = T)

    model_somatic <- lmer(log_rt ~ somatic_t + log_age + sex + ethnicity 
                    + (1 | SID),
                    data = dat_clinical, REML = T)

    model_thought <- lmer(log_rt ~ thought_t + log_age + sex + ethnicity 
                    + (1 | SID),
                    data = dat_clinical, REML = T)
    
    model_attention <- lmer(log_rt ~ attention_t + log_age + sex + ethnicity 
                    + (1 | SID),
                    data = dat_clinical, REML = T)

    model_aggressive <- lmer(log_rt ~ aggressive_t + log_age + sex + ethnicity 
                    + (1 | SID),
                    data = dat_clinical, REML = T)

    model_rulebreaking <- lmer(log_rt ~ rulebreaking_t + log_age + sex + ethnicity 
                    + (1 | SID),
                    data = dat_clinical, REML = T)

    all_models <- list(
                    suicide = model_suicide,
                    anxdepr = model_anxdepr,
                    withdrdepr = model_withdrdepr,
                    somatic = model_somatic,
                    thought = model_thought,
                    attention = model_attention,
                    aggressive = model_aggressive,
                    rulebreaking = model_rulebreaking
                )

    results <- Create.Summary.Models(all_models)
    write_csv(results, "./4-Clinical_Sympt_Prediction/Outputs/RT_Congruent_Models.csv")
}


## Runs the analysis for the models on congruent trials only
## for Accuracy
Congruent.Accuracy <- function(EST.data, clinical.data, suicide.data) {

    EST <- read_rds(EST.data) |> 
                ## Keeping only congruent trials
                subset(congruence == "congruent") |>
                    mutate(correct = case_when(
                        correct == "correct"   ~ 1,
                        correct == "incorrect" ~ 0
                    ), .keep = "unused")

    clinical <- read_rds(clinical.data)

    suicide <- read_rds(suicide.data)

    ## ----- Suicide analysis -----
    ## don't bind all 3 (EST, clinical, suicide) together bc would loose a lot of suicide data
    ## EST + clinical: n = 690
    ## EST + suicide: n = 640
    ## ptcp w EST + clinical + suicide data: n = 616
    dat_suicide <- inner_join(EST, suicide, by = "SID")
    
    model_suicide <- glmer(correct ~ cssrs_attempt_lt + log_age + sex + ethnicity
                    + (1 | SID),
                    data = dat_suicide, family = binomial(link = "logit"))
    
    ## ----- Clinical sympt analysis -----
    dat_clinical <- inner_join(EST, clinical, by = "SID")

    model_anxdepr <- glmer(correct ~ anxdepr_t + log_age + sex + ethnicity
                    + (1 | SID),
                    data = dat_clinical, family = binomial(link = "logit"))

    model_withdrdepr <- glmer(correct ~ withdrdepr_t + log_age + sex + ethnicity
                    + (1 | SID),
                    data = dat_clinical, family = binomial(link = "logit"))

    model_somatic <- glmer(correct ~ somatic_t + log_age + sex + ethnicity
                    + (1 | SID),
                    data = dat_clinical, family = binomial(link = "logit"))

    model_thought <- glmer(correct ~ thought_t + log_age + sex + ethnicity
                    + (1 | SID),
                    data = dat_clinical, family = binomial(link = "logit"))

    model_attention <- glmer(correct ~ attention_t + log_age + sex + ethnicity
                    + (1 | SID),
                    data = dat_clinical, family = binomial(link = "logit"))

    model_aggressive <- glmer(correct ~ aggressive_t + log_age + sex + ethnicity
                    + (1 | SID),
                    data = dat_clinical, family = binomial(link = "logit"))
    
    model_rulebreaking <- glmer(correct ~ rulebreaking_t + log_age + sex + ethnicity
                    + (1 | SID),
                    data = dat_clinical, family = binomial(link = "logit"))

    all_models <- list(
                    suicide = model_suicide,
                    anxdepr = model_anxdepr,
                    withdrdepr = model_withdrdepr,
                    somatic = model_somatic,
                    thought = model_thought,
                    attention = model_attention,
                    aggressive = model_aggressive,
                    rulebreaking = model_rulebreaking
                )

    results <- Create.Summary.Models(all_models)
    write_csv(results, "./4-Clinical_Sympt_Prediction/Outputs/Accuracy_Congruent_Models.csv")
}