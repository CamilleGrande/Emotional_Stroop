## Utilities_Clinical_Prediction.R
## Author: Camille Grandé
## Description:
##      This script holds the utility function to run the task performance prediction
##      using clinical metrics


## Helper function to create one big summary of our models
## and FDR correct p values across scales being tested for each outcome
## We don't correct across outcomes, we correct across scales for each outcome
## i.e., we see if for one outcome, testing multiple scales gave rise to chance errors
Create.Summary.Models <- function(list.models) {

    predictors <- "cssrs_attempt_lt\\d*|anxdepr_t|withdrdepr_t|somatic_t|thought_t|attention_t|aggressive_t|rulebreaking_t"

    results <- imap_dfr(
                list.models, ~ tidy(.x, effects = "fixed") |>
                    mutate(model = .y) |>
                        relocate(model)
                ) |>
        mutate(is_predictor = str_detect(term, predictors),
               ## label the effect without the predictor name, so the same effect matches across models
               effect = str_remove(term, predictors)) |>
        group_by(effect) |>
        ## Adjusting with FDR across the 8 models (clinical constructs are correlated)
        mutate(p_fdr = ifelse(is_predictor, p.adjust(p.value, "fdr"), NA)) |>
        ungroup()

    return(results)
}


## Helper function to center log age for each population
## It centers based on unique SID, not for each row
## so participants with more/less rows (e.g., more/less correct trials)
## are not down or up-weighted in mean log age calculation
Center.Age <- function(d) {
    d |> mutate(log_age = log(age) - mean(log(age[!duplicated(SID)])))
}


## Centers clinical scales at the participant level (one value per SID)
## like with age
## Negative values = lower symptoms than mean; positive values = higher symptoms than mean
Center.Clinical <- function(d) {

    clinical_vars <- c("anxdepr_t", "withdrdepr_t", "somatic_t", "thought_t",
                   "attention_t", "aggressive_t", "rulebreaking_t")

    means <- d |>
        distinct(SID, across(all_of(clinical_vars))) |>
        summarise(across(all_of(clinical_vars), \(x) mean(x, na.rm = TRUE)))

    d |>
        mutate(across(all_of(clinical_vars), \(x) x - means[[cur_column()]]))
}


## This function computes the confidence intervals for our ratios
## to get CI for percent faster/slower RT
Compute.CI <- function(my.emm, label) {

    ## get the CI on the contrasts
    ci <- as.data.frame(confint(my.emm$contrasts, adjust = "none"))
    lcl <- grep("LCL|lower.CL", names(ci), value = TRUE)
    ucl <- grep("UCL|upper.CL", names(ci), value = TRUE)

    ci$analysis <- label
    ## negative = faster
    ci$pct_change <- round((100 * (ci$ratio - 1)), 2)
    ## order preserved, no flip as ratio - 1 is an increasing function
    ## before had 1 - ratio, which is a decreasing function so had to reverse lcl as uci and ucl as lci
    ci$pct_lo <- round((100 * (ci[[lcl]] - 1)), 2)
    ci$pct_hi <- round((100 * (ci[[ucl]] - 1)), 2)
    
    clean_ci <- ci[, setdiff(names(ci), c("SE", "df", "null", "z.ratio", "p.value", lcl, ucl, "ratio"))]
    return(clean_ci)
}


## OR computation for suicide attempt (binomial)
Suicide.OR <- function(model, contr, spec) {

    emm <- emmeans(model, as.formula(paste(contr, "~", spec)), type = "response")
    OR <- Compute.OR(emm, label = spec)

    print(emm)
    return(OR)
}


## OR computation for clinical scales (continuous)
Clinical.OR <- function(model, spec, baseline = FALSE) {

    ## compute sd scale
    mf <- model.frame(model)
    sd_scale <- sd(mf[[spec]][!duplicated(mf$SID)])

    at_list <- setNames(list(c(-sd_scale, 0, sd_scale)), spec)

    if (baseline == FALSE) {

        ## How big is the congruency/face emotion/ congruency per face emotion effect at -1 SD, mean, and +1 SD of scale?
        emm_cong <- emmeans(model, as.formula(paste("revpairwise ~ congruence | ", spec)),
                            at = at_list, type = "response")
        print(emm_cong)
        
        emm_emo <- emmeans(model, as.formula(paste("revpairwise ~ face_emotion | ", spec)),
                            at = at_list, type = "response")
        print(emm_emo)

        emm_cong_emo <- emmeans(model, as.formula(paste("revpairwise ~ congruence | face_emotion * ", spec)),
                            at = at_list, type = "response")
        print(emm_cong_emo)

        ## 1st line is -1SD, second one mean, third one +1SD
        res_cong <- Compute.OR(emm_cong, label = paste0("congruence | ", spec, " (-1 SD, mean, +1 SD)"))
        res_emo <- Compute.OR(emm_emo, label = paste0("face emotion | ", spec, " (-1 SD, mean, +1 SD)"))
        res_cong_emo <- Compute.OR(emm_cong_emo, label = paste0("congruence | face emotion * ", spec, " (-1 SD, mean, +1 SD)"))

        all_res <- bind_rows(res_cong, res_emo, res_cong_emo)
    }

    if (baseline == TRUE) {

        emm_baseline <- emmeans(model, as.formula(paste("revpairwise ~ ", spec)),
                            at = at_list, type = "response")
        print(emm_baseline)

        all_res <- Compute.OR(emm_baseline, label = paste0(spec, ", (-1 SD, mean, +1 SD)"))
    }
    
    return(all_res)
}


## Computes the % change for +1 SD for each scale
Percent.Change.per.SD <- function(model, spec, by = NULL) {
    
    emm_options(lmer.df = "asymptotic")

    mf <- model.frame(model)
    sd_spec <- sd(mf[[spec]][!duplicated(mf$SID)])

    ## two points -> a single contrast per group: +1 SD vs mean
    ## (single contrast = no multiplicity adjustment, CI matches the Wald CI)
    at_list <- setNames(list(c(0, sd_spec)), spec)

    rhs <- if (is.null(by)) spec else paste(spec, "|", paste(by, collapse = " * "))

    emm <- emmeans(model, as.formula(paste("revpairwise ~", rhs)),
                   at = at_list, tran = "log", type = "response")

    label <- paste0(spec, if (!is.null(by)) paste0(" | ", paste(by, collapse = " * ")),
                    " (+1 SD vs mean)")

    return(Compute.CI(emm, label = label))
}


## Computes odds ratios + CIs from emmeans contrasts (binomial GLMM)
Compute.OR <- function(my.emm, label) {

    ci  <- as.data.frame(confint(my.emm$contrasts))
    lcl <- grep("LCL|lower.CL", names(ci), value = TRUE)
    ucl <- grep("UCL|upper.CL", names(ci), value = TRUE)

    ci$analysis <- label
    ci$OR <- round(ci$odds.ratio, 2)
    ci$OR_lo <- round(ci[[lcl]], 2)
    ci$OR_hi <- round(ci[[ucl]], 2)

    return(ci[, setdiff(names(ci), c("SE", "df", "null", "z.ratio", "p.value", "odds.ratio", lcl, ucl))])
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

    ## Set contrasts to sum-to-zero coding
    ## Age is centered in data prep
    options(contrasts = c("contr.sum", "contr.poly"))

    emm_options(lmer.df = "asymptotic")

    ## ----- Suicide analysis -----
    ## don't bind all 3 (EST, clinical, suicide) together bc would loose a lot of suicide data
    ## EST + clinical: n = 699
    ## EST + suicide: n = 646
    ## ptcp w EST + clinical + suicide data: n = 616
    dat_suicide <- inner_join(EST, suicide, by = "SID") |>
                        ## Centering log age for easier interpretation
                        ## The intercept will now refer to an average-aged participant
                        Center.Age()

    model_suicide <- lmer(log_rt ~ congruence * face_emotion * cssrs_attempt_lt + log_age + sex + ethnicity 
                    + (1 | SID),
                    data = dat_suicide, REML = T)

    ## Get RT difference between no history of pior SA and history of prior SA
    ## to get ratio (diff RT btw the two): exp(response0 - response1)
    emm_sui_cong <- emmeans(model_suicide, "pairwise" ~ "cssrs_attempt_lt | congruence", tran = "log", type = "response")
    emm_sui_emo <- emmeans(model_suicide, "pairwise" ~ "cssrs_attempt_lt | face_emotion", tran = "log", type = "response")
    emm_sui_cong_emo <- emmeans(model_suicide, "pairwise" ~ "cssrs_attempt_lt | congruence * face_emotion", tran = "log", type = "response")
    ## Get % change in RT between the two and CI
    res_suicide_cong <- Compute.CI(emm_sui_cong, label = "cssrs_attempt_lt")
    res_suicide_emo <- Compute.CI(emm_sui_emo, label = "cssrs_attempt_lt")
    res_suicide_cong_emo <- Compute.CI(emm_sui_cong_emo, label = "cssrs_attempt_lt")

    ## ----- Clinical sympt analysis -----
    dat_clinical <- inner_join(EST, clinical, by = "SID") |>
                        ## Centering log age for easier interpretation
                        ## The intercept will now refer to an average-aged participant
                        Center.Age() |>
                            ## Centering clinical scales for easier interpretation based on the sample mean
                            ## Negative values = low t-score; positive values = high t-score
                            Center.Clinical()

    ## Fitting all models first
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


    ## % change for each model
    res_anx_cong <- Percent.Change.per.SD(model_anxdepr, "anxdepr_t", by = "congruence")
    res_anx_emo  <- Percent.Change.per.SD(model_anxdepr, "anxdepr_t", by = "face_emotion")
    res_anx_3way <- Percent.Change.per.SD(model_anxdepr, "anxdepr_t", by = c("congruence", "face_emotion"))

    res_with_cong <- Percent.Change.per.SD(model_withdrdepr, "withdrdepr_t", by = "congruence")
    res_with_emo  <- Percent.Change.per.SD(model_withdrdepr, "withdrdepr_t", by = "face_emotion")
    res_with_3way <- Percent.Change.per.SD(model_withdrdepr, "withdrdepr_t", by = c("congruence", "face_emotion"))

    res_som_cong <- Percent.Change.per.SD(model_somatic, "somatic_t", by = "congruence")
    res_som_emo  <- Percent.Change.per.SD(model_somatic, "somatic_t", by = "face_emotion")
    res_som_3way <- Percent.Change.per.SD(model_somatic, "somatic_t", by = c("congruence", "face_emotion"))

    res_thou_cong <- Percent.Change.per.SD(model_thought, "thought_t", by = "congruence")
    res_thou_emo  <- Percent.Change.per.SD(model_thought, "thought_t", by = "face_emotion")
    res_thou_3way <- Percent.Change.per.SD(model_thought, "thought_t", by = c("congruence", "face_emotion"))

    res_att_cong <- Percent.Change.per.SD(model_attention, "attention_t", by = "congruence")
    res_att_emo  <- Percent.Change.per.SD(model_attention, "attention_t", by = "face_emotion")
    res_att_3way <- Percent.Change.per.SD(model_attention, "attention_t", by = c("congruence", "face_emotion"))

    res_agg_cong <- Percent.Change.per.SD(model_aggressive, "aggressive_t", by = "congruence")
    res_agg_emo  <- Percent.Change.per.SD(model_aggressive, "aggressive_t", by = "face_emotion")
    res_agg_3way <- Percent.Change.per.SD(model_aggressive, "aggressive_t", by = c("congruence", "face_emotion"))

    res_rul_cong <- Percent.Change.per.SD(model_rulebreaking, "rulebreaking_t", by = "congruence")
    res_rul_emo  <- Percent.Change.per.SD(model_rulebreaking, "rulebreaking_t", by = "face_emotion")
    res_rul_3way <- Percent.Change.per.SD(model_rulebreaking, "rulebreaking_t", by = c("congruence", "face_emotion"))

    ## List with all models primary results (look at p FDR to know if should look at % change)
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

    ## All % change results
    all_pct_change <- bind_rows(
                        res_suicide_cong, res_suicide_emo, res_suicide_cong_emo,
                        res_anx_cong, res_anx_emo, res_anx_3way,
                        res_with_cong, res_with_emo, res_with_3way,
                        res_som_cong, res_som_emo, res_som_3way,
                        res_thou_cong, res_thou_emo, res_thou_3way,
                        res_att_cong, res_att_emo, res_att_3way,
                        res_agg_cong, res_agg_emo, res_agg_3way,
                        res_rul_cong, res_rul_emo, res_rul_3way)

    write_csv(all_pct_change, "./4-Clinical_Sympt_Prediction/Outputs/RT_Interference-FaceEmotion_Pct_Change.csv")
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

    ## Set contrasts to sum-to-zero coding
    ## Age is centered in data prep
    options(contrasts = c("contr.sum", "contr.poly"))

    ## ----- Suicide analysis -----
    ## don't bind all 3 (EST, clinical, suicide) together bc would loose a lot of suicide data
    ## EST + clinical: n = 699
    ## EST + suicide: n = 646
    ## ptcp w EST + clinical + suicide data: n = 616
    dat_suicide <- inner_join(EST, suicide, by = "SID") |>
                        ## Centering log age for easier interpretation
                        ## The intercept will now refer to an average-aged participant
                        Center.Age()
    
    model_suicide <- glmer(correct ~ congruence * face_emotion * cssrs_attempt_lt + log_age + sex + ethnicity
                    + (1 | SID),
                    data = dat_suicide, family = binomial(link = "logit"))

    ## OR for each outcome
    res_sui_cong <- Suicide.OR(model_suicide, "revpairwise", "congruence | cssrs_attempt_lt")
    res_sui_emo <- Suicide.OR(model_suicide, "revpairwise", "face_emotion | cssrs_attempt_lt")
    res_sui_cong_emo <- Suicide.OR(model_suicide, "revpairwise", "congruence | face_emotion * cssrs_attempt_lt")
    
    ## ----- Clinical sympt analysis -----
    dat_clinical <- inner_join(EST, clinical, by = "SID") |>
                        ## Centering log age for easier interpretation
                        ## The intercept will now refer to an average-aged participant
                        Center.Age() |>
                            ## Centering clinical scales for easier interpretation based on the sample mean
                            ## Negative values = low t-score; positive values = high t-score
                            Center.Clinical()

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

    ## OR for each model
    res_anxdep <- Clinical.OR(model_anxdepr, "anxdepr_t")
    res_withdep <- Clinical.OR(model_withdrdepr, "withdrdepr_t")
    res_som <- Clinical.OR(model_somatic, "somatic_t")
    res_thou <- Clinical.OR(model_thought, "thought_t")
    res_att <- Clinical.OR(model_attention, "attention_t")
    res_agg <- Clinical.OR(model_aggressive, "aggressive_t")
    res_ruleb <- Clinical.OR(model_rulebreaking, "rulebreaking_t")

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
    write_csv(results, "./4-Clinical_Sympt_Prediction/Outputs/Acc_Interference-FaceEmotion_Models.csv")


    ## All OR results
    all_OR <- bind_rows(res_sui_cong, res_sui_emo, res_sui_cong_emo,
                        res_anxdep, res_withdep, res_som, res_thou, res_att, res_agg, res_ruleb)

    write_csv(all_OR, "./4-Clinical_Sympt_Prediction/Outputs/Acc_Interference-FaceEmotion_OR.csv")

    ## Plots for significant findings
    

    y_lim <- c(0.95, 1)

    no_y <- theme(axis.title.y = element_blank(),
              axis.text.y  = element_blank(),
              axis.ticks.y = element_blank(),
              axis.line.y  = element_blank())

    p1 <- Plot.Accuracy.by.Scale(model_somatic, "somatic_t", "Somatic complaints")
    p2 <- Plot.Accuracy.by.Scale(model_thought, "thought_t", "Thought problems") + no_y
    p3 <- Plot.Accuracy.by.Scale(model_attention, "attention_t", "Attention problems") + no_y
    
    final_p <- (p1 | p2 | p3) + 
                plot_annotation(tag_levels = "A") + 
                plot_layout(guides = "collect") &
                coord_cartesian(ylim = y_lim) &
                theme(legend.position = "bottom")

    ggsave("./4-Clinical_Sympt_Prediction/Outputs/Acc_all_by_cong_emotion.pdf", final_p, width = 15)
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

    ## Set contrasts to sum-to-zero coding
    ## Age is centered in data prep
    options(contrasts = c("contr.sum", "contr.poly"))

    emm_options(lmer.df = "asymptotic")

    ## ----- Suicide analysis -----
    ## don't bind all 3 (EST, clinical, suicide) together bc would loose a lot of suicide data
    ## EST + clinical: n = 699
    ## EST + suicide: n = 646
    ## ptcp w EST + clinical + suicide data: n = 616
    dat_suicide <- inner_join(EST, suicide, by = "SID") |>
                        ## Centering log age for easier interpretation
                        ## The intercept will now refer to an average-aged participant
                        Center.Age()

    model_suicide <- lmer(log_rt ~ cssrs_attempt_lt + log_age + sex + ethnicity 
                    + (1 | SID),
                    data = dat_suicide, REML = T)
    
    ## Get RT difference between no history of pior SA and history of prior SA
    ## to get ratio (diff RT btw the two): exp(response0 - response1)
    emm_suicide <- emmeans(model_suicide, "pairwise" ~ "cssrs_attempt_lt", tran = "log", type = "response")
    ## Get % change in RT between the two and CI
    res_suicide <- Compute.CI(emm_suicide, label = "cssrs_attempt_lt")

    ## ----- Clinical sympt analysis -----
    dat_clinical <- inner_join(EST, clinical, by = "SID") |>
                        ## Centering log age for easier interpretation
                        ## The intercept will now refer to an average-aged participant
                        Center.Age() |>
                            ## Centering clinical scales for easier interpretation based on the sample mean
                            ## Negative values = low t-score; positive values = high t-score
                            Center.Clinical()

    ## Fit all models first
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
    
    ## % change for each model
    res_anxdepr <- Percent.Change.per.SD(model_anxdepr, "anxdepr_t")
    res_withdrdepr <- Percent.Change.per.SD(model_withdrdepr, "withdrdepr_t")
    res_somatic <- Percent.Change.per.SD(model_somatic, "somatic_t")
    res_thought <- Percent.Change.per.SD(model_thought, "thought_t")
    res_attention <- Percent.Change.per.SD(model_attention, "attention_t")
    res_aggressive <- Percent.Change.per.SD(model_aggressive, "aggressive_t")
    res_rulebreaking <- Percent.Change.per.SD(model_rulebreaking, "rulebreaking_t")

    ## List with all models primary results (look at p FDR to know if should look at % change)
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

    ## All % change results
    all_pct_change <- bind_rows(res_suicide, res_anxdepr, res_withdrdepr, res_somatic, res_thought, res_attention, res_aggressive, res_rulebreaking)
    write_csv(all_pct_change, "./4-Clinical_Sympt_Prediction/Outputs/RT_Congruent_Pct_Change.csv")
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

    ## Set contrasts to sum-to-zero coding
    ## Age is centered in data prep
    options(contrasts = c("contr.sum", "contr.poly"))

    ## ----- Suicide analysis -----
    ## don't bind all 3 (EST, clinical, suicide) together bc would loose a lot of suicide data
    ## EST + clinical: n = 699
    ## EST + suicide: n =646
    ## ptcp w EST + clinical + suicide data: n = 616
    dat_suicide <- inner_join(EST, suicide, by = "SID") |>
                        ## Centering log age for easier interpretation
                        ## The intercept will now refer to an average-aged participant
                        Center.Age()
    
    model_suicide <- glmer(correct ~ cssrs_attempt_lt + log_age + sex + ethnicity
                    + (1 | SID),
                    data = dat_suicide, family = binomial(link = "logit"))
    
    ## OR for each outcome
    res_sui_base <- Suicide.OR(model_suicide, "revpairwise", "cssrs_attempt_lt")

    ## ----- Clinical sympt analysis -----
    dat_clinical <- inner_join(EST, clinical, by = "SID") |>
                        ## Centering log age for easier interpretation
                        ## The intercept will now refer to an average-aged participant
                        Center.Age() |>
                            ## Centering clinical scales for easier interpretation based on the sample mean
                            ## Negative values = low t-score; positive values = high t-score
                            Center.Clinical()

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

    ## OR for each model
    res_anxdep <- Clinical.OR(model_anxdepr, "anxdepr_t", baseline = TRUE)
    res_withdep <- Clinical.OR(model_withdrdepr, "withdrdepr_t", baseline = TRUE)
    res_som <- Clinical.OR(model_somatic, "somatic_t", baseline = TRUE)
    res_thou <- Clinical.OR(model_thought, "thought_t", baseline = TRUE)
    res_att <- Clinical.OR(model_attention, "attention_t", baseline = TRUE)
    res_agg <- Clinical.OR(model_aggressive, "aggressive_t", baseline = TRUE)
    res_ruleb <- Clinical.OR(model_rulebreaking, "rulebreaking_t", baseline = TRUE)


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
    write_csv(results, "./4-Clinical_Sympt_Prediction/Outputs/Acc_Congruent_Models.csv")

    ## All OR results
    all_OR <- bind_rows(res_sui_base, res_anxdep, res_withdep, res_som, res_thou, res_att, res_agg, res_ruleb)

    write_csv(all_OR, "./4-Clinical_Sympt_Prediction/Outputs/Acc_Congruent_OR.csv")
}


## Predicted probability of a correct response as a function of a centered
## clinical scale. Lines span the observed range of the scale (participant level)
## Works for any glmer accuracy model of the form correct ~ congruence * face_emotion * <scale>
Plot.Accuracy.by.Scale <- function(model, scale, scale_label = scale) {

    ## participant-level SD and observed range of the centered scale
    mf   <- model.frame(model)
    pt   <- mf[!duplicated(mf$SID), ]
    sd_s <- sd(pt[[scale]])
    at_list <- setNames(list(seq(min(pt[[scale]]), max(pt[[scale]]), length.out = 60)), scale)

    ## model-based probabilities (response scale) with 95% CI
    ## averaged over sex and ethnicity; age at its centered mean (0)
    pred <- emmeans(model, as.formula(paste("~", scale, "* congruence * face_emotion")),
                    at = at_list, type = "response") |>
        as.data.frame()

    ## CI column names differ by df method (asymp.* vs lower.CL/upper.CL)
    lcl <- grep("LCL|lower.CL", names(pred), value = TRUE)
    ucl <- grep("UCL|upper.CL", names(pred), value = TRUE)

    p <- ggplot(pred, aes(x = .data[[scale]], y = prob, colour = congruence, fill = congruence)) +
            geom_ribbon(aes(ymin = .data[[lcl]], ymax = .data[[ucl]]), alpha = 0.2, colour = NA) +
            geom_line(linewidth = 1) +
            facet_wrap(~ face_emotion, labeller = labeller(face_emotion = str_to_title)) +
            scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
            scale_color_manual(values = c("#C2DDB2FF", "#B74954FF"), name = "Congruence", labels = c("Congruent", "Incongruent")) +
            scale_fill_manual(values = c("#C2DDB2FF", "#B74954FF"), name = "Congruence", labels = c("Congruent", "Incongruent")) +
            labs(x = paste0(scale_label, " (centered)"),
                y = "Predicted accuracy (%)") +
            theme_classic() +
            theme(legend.position = "bottom",
                    axis.text = element_text(size = 12),
                    axis.title = element_text(size = 14, face = "bold"),
                    strip.text.x = element_text(size = 12),
                    legend.text = element_text(size = 12),
                    panel.grid.minor = element_blank(),
                    strip.background = element_rect(fill = "grey95"))
}