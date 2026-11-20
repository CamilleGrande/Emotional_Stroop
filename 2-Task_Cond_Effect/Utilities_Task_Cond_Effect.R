## Utilities_Task_Cond_Effect.R
## Author: Camille Grandé
## Description:
##      This script holds the utility functions to run the stats to test the task conditions effect for the emotional stroop task analysis

## Helper function to center log age for each population
## It centers based on unique SID, not for each row
## so participants with more/less rows (e.g., more/less correct trials)
## are not down or up-weighted in mean log age calculation
Center.Age <- function(d) {
    d |> mutate(log_age = log(age) - mean(log(age[!duplicated(SID)])))
}

## Helper function to prep the RT data
Prep.Data.RT <- function(data) {

    dat <- read_rds(data) |>
                ## Keeping only RT for correct trials
                subset(correct == "correct") |>
                    ## Centering log age for easier interpretation
                    ## The intercept will now refer to an average-aged participant
                    Center.Age()
    
    return(dat)
}

## Helper function to prep the accuracy data
Prep.Data.Acc <- function(data) {

    dat <- read_rds(data) |>
                mutate(correct = case_when(
                    correct == "correct"   ~ 1,
                    correct == "incorrect" ~ 0
                ), .keep = "unused") |>
                    ## Centering log age for easier interpretation
                    ## The intercept will now refer to an average-aged participant
                    Center.Age()
    
    return(dat)
}

## Function for initial RT model and model assumptions
Modelling.RT <- function(data) {

    dat <- Prep.Data.RT(data)
    
    ## Set contrasts to sum-to-zero coding so that it matches mixed() output later
    options(contrasts = c("contr.sum", "contr.poly"))

    ## Initial model
    model1 <- lmer(reaction_time ~ congruence * frequency * face_emotion + log_age + sex + ethnicity 
                    + (1 | SID),
                    data = dat, REML = T)
    
    cat("\nResiduals are distributed in fan-shape; log-scoring the RT\n")
    ## Compute RT log-scores
    dat$log_rt <- log(dat$reaction_time)

    ## Excluded trials + RT Log-scored
    model2 <- lmer(log_rt ~ congruence * frequency * face_emotion + log_age + sex + ethnicity 
                    + (1 | SID),
                    data = dat, REML = T)

    write_rds(model2, "./2-Task_Cond_Effect/Outputs/RT_Model_results.rds")

    ## Residuals are linearly distributed (similar deviations)
    ## Homoskedasticity OK (no funnel plot)
    ## Normality of residuals OK 
    ## (see the two plots below)

    ## Check collinearity among fixed effects
    collinearity <- check_collinearity(model2)
    cat("\nLook at csv to assess collinearity\n")
    write_csv(collinearity, "./2-Task_Cond_Effect/Outputs/RT_Collinearity.csv")

    cat("\nAssumptions for LME all met; look at plots!\n\n")
    pdf("./2-Task_Cond_Effect/Outputs/RT_Assumptions.pdf")
        # model 1
        plot(fitted(model1), residuals(model1),
                ylab = "Residuals", xlab = "Reaction Time (ms)",
                main = "Residuals against fitted RT (Initial model)")
        # model 2
        plot(fitted(model2), residuals(model2),
                ylab = "Residuals", xlab = "log-scored Reaction Time",
                main = "Residuals against fitted log-scored RT (log model)")
        # normality residuals
        hist(residuals(model2), 
                ylab = "Frequency", xlab = "Residuals",
                main = "Histogram of residuals")
        qqnorm(residuals(model2),
                main = "QQ-plot of residuals")
    dev.off()

    return(model2)
}


## Function for ANOVA part and effect size
ANOVA.RT <- function(my.data) {

    ## Calls function for initial model and assumptions
    my.model <- Modelling.RT(my.data)

    ## Reformatting data for mixed() models like in Modelling.RT function
    dat <- Prep.Data.RT(data)

    ## Compute RT log-scores
    dat$log_rt <- log(dat$reaction_time)

    RT_Anova <- afex::mixed(log_rt ~ congruence * frequency * face_emotion + log_age + sex + ethnicity
                     + (1 | SID), data = dat, method = "S")

    my.ef <- effectsize::eta_squared(RT_Anova$full_model, partial = TRUE) |>
                dplyr::rename(term = Parameter)

    results_RT <- afex::nice(RT_Anova) |>
                    dplyr::rename(term = Effect) |>
                        dplyr::left_join(my.ef, by = "term")

    write.csv(results_RT, "./2-Task_Cond_Effect/Outputs/RT_Anova_results.csv")
}


## This function plots RTs for the three way interaction and age
Plots.RT <- function(model, EST.data) {

    emm_options(lmer.df = "asymptotic")

    ## ---- Plot three way interaction ----
    ## Fetch our model (saved it so no need to re run analysis for plotting)
    RT_perf <- read_rds(model)

    ## Plot three way interaction
    emm <- emmeans(RT_perf, ~ congruence * frequency * face_emotion)
    ## Save as df for plotting, otherwise won't work
    emm_df <- as.data.frame(emm)

    ## Exponentiate coefficients as are on the log scale (log RT)
    emm_df$RT_ms <- exp(emm_df$emmean)
    ## Lower and Upper CI
    emm_df$lower_ms <- exp(emm_df$asymp.LCL)
    emm_df$upper_ms <- exp(emm_df$asymp.UCL)

    three_way_plot <- ggplot(emm_df, aes(x = congruence, y = RT_ms, color = frequency, group = frequency)) +
                        geom_point(size = 4) +
                        scale_y_continuous(breaks = seq(800, 1000, by = 25)) +
                        geom_line(linewidth = 1) +
                        geom_errorbar(aes(ymin = lower_ms, ymax = upper_ms), linewidth = 1, width = .1) +
                        facet_wrap(~ face_emotion) +
                        labs(x = "Congruence", y = "Predicted RT (ms)", color = "Frequency") +
                        theme_classic() +
                        theme(axis.text = element_text(size = 18),
                                axis.title = element_text(size = 20, face = "bold"),
                                strip.text.x = element_text(size = 18),
                                panel.grid.minor = element_blank(),
                                strip.background = element_rect(fill = "grey95"),
                                legend.text = element_text(size = 18),
                                legend.title = element_text(size = 20, face = "bold"),
                                legend.position = "bottom") +
                        scale_color_manual(values = c("#E3C454FF", "#A27CBAFF"), labels = c("Equal", "More congruent"))
   
    ggsave("./2-Task_Cond_Effect/Outputs/RT_Three_way_interaction.pdf", three_way_plot, width = 7, height = 7)


    ## ---- Plot predicted RT for age ----
    EST <- read_rds(EST.data) |> 
            subset(correct == "correct")

    ## Same centering constant as Center.Age() used when fitting the model
    age_center <- mean(log(EST$age[!duplicated(EST$SID)]))

    ## Get predicted values for log age from our model
    pred_age <- ggpredict(RT_perf, terms = "log_age [all]")

    ## Back-transform centered log age to age in years
    pred_age$age <- exp(pred_age$x + age_center)

    ## Exponentiate coefficients as are on the log scale (log RT)
    pred_age$RT_ms <- exp(pred_age$predicted)
    pred_age$lower_ms <- exp(pred_age$conf.low)
    pred_age$upper_ms <- exp(pred_age$conf.high)

    age_plot <- ggplot(pred_age, aes(age, RT_ms)) +
            geom_line(colour = "#339999", linewidth = 1) +
            geom_ribbon(aes(ymin = lower_ms, ymax = upper_ms), fill = "#339999", alpha = .2) +
            labs(x = "Age (years)",y = "Predicted RT (ms)") +
            theme_classic() +
            theme(axis.text = element_text(size = 18),
                  axis.title = element_text(size = 20, face = "bold"))
    
    ggsave("./2-Task_Cond_Effect/Outputs/RT_Age.pdf", age_plot, width = 8, height = 5)
}


## This function computes the confidence intervals for our ratios
## to get CI for percent faster/slower RT
Compute.CI <- function(my.emm, label) {

    ## get the CI on the contrasts
    ci <- as.data.frame(confint(my.emm$contrasts))
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


Percent.Change.RT <- function(model) {

    ## Model already fitted
    RT_perf <- read_rds(model)

    emm_options(lmer.df = "asymptotic")

    ## all the comparisons we want to test
    specs <- c(
        "congruence",
        "frequency",
        "face_emotion",
        "congruence | frequency",
        "congruence | face_emotion",
        "face_emotion | frequency",
        "congruence | frequency * face_emotion",
        "sex",
        "ethnicity"
    )

    ## one list with all comparisons and CIs
    res <- lapply(specs, function(s) {
        ## revpairwise gives incongruent / congruent (positive = slower for incongruent)
        contr <- if (grepl("^congruence", s)) "revpairwise" else "pairwise"
        emm <- emmeans(RT_perf, as.formula(paste(contr, "~", s)),
                    tran = "log", type = "response")
        Compute.CI(emm, label = s)
    })

    ## For age, % change in RT per 1 SD of log-age
    ## beta
    b  <- unname(fixef(RT_perf)["log_age"])
    ## confidence interval
    ci <- confint(RT_perf, parm = "log_age", method = "Wald")

    ## SD at the participant level (age is between-subject, so we don't use trial-level rows,
    ## which would weight people by how many trials they contributed)
    mf <- model.frame(RT_perf)
    sd_logage <- sd(mf$log_age[!duplicated(mf$SID)])

    res_age <- data.frame(
        analysis   = "age",
        contrast   = "+1 SD of log(age)",
        pct_change = round((100 * (exp(b * sd_logage) - 1)), 2),
        pct_lo     = round((100 * (exp(ci[1, 1] * sd_logage) - 1)), 2),
        pct_hi     = round((100 * (exp(ci[1, 2] * sd_logage) - 1)), 2)
    )

    all_res <- dplyr::bind_rows(c(res, list(res_age)))
    write_csv(all_res, "./2-Task_Cond_Effect/Outputs/RT_Percent_Change.csv")
    return(all_res)
}


## Function for initial accuracy model and model assumptions
Modelling.Accuracy <- function(data) {

    dat <- Prep.Data.Acc(data)

    ## Set contrasts to sum-to-zero coding so that it matches mixed() output later
    ## Age is centered in data prep
    options(contrasts = c("contr.sum", "contr.poly"))

    ## random intercept only as random slope increases dispersion
    accuracy <- glmer(correct ~ congruence * frequency * face_emotion + log_age + sex + ethnicity
                    + (1 | SID),
                    data = dat, family = binomial(link = "logit"))

    write_rds(accuracy, "./2-Task_Cond_Effect/Outputs/Acc_Model_results.rds")

    cat("\nChecking assumptions for our Binomial Generalized Linear Mixed Model\n")

    ## Check collinearity among fixed effects
    ## Computes VIF, 0-5 = low, 5-10 = moderate, > 10 = high multicollinearity
    collinearity <- check_collinearity(accuracy)
    cat("\nLook at csv to assess collinearity\n")
    write_csv(collinearity, "./2-Task_Cond_Effect/Outputs/Acc_Collinearity.csv")


    ## Simulation-based residual diagnostics (DHARMa)
    ## Calculates randomized quantile residuals 
    sim_res <- simulateResiduals(fittedModel = accuracy, n = 1000, plot = F)

    pdf("./2-Task_Cond_Effect/Outputs/Acc_Assumptions.pdf")
        ## QQ plot + residual plot residuals fitted vs. simulated
        plot(sim_res)
        ## Dispersion check
        testDispersion(sim_res)
        ## plot outliers
        testOutliers(sim_res)
        testUniformity(sim_res)
        ## Residuals against predictors
        plotResiduals(sim_res, form = dat$congruence)
        plotResiduals(sim_res, form = dat$frequency)
        plotResiduals(sim_res, form = dat$face_emotion)
        plotResiduals(sim_res, form = dat$log_age)
        plotResiduals(sim_res, form = dat$sex)
        plotResiduals(sim_res, form = dat$ethnicity)
    dev.off()

    return(accuracy)
}


## Function for Likelihood Ratio test part
LRT.Accuracy <- function(data) {

    ## Calls function for initial model and assumptions
    Modelling.Accuracy(data)

    ## Reformatting data for mixed() models like in Modelling.Accuracy function
    dat <- Prep.Data.Acc(data)

    ## re-fit model with afex to get LRT
    ## check_contrast = TRUE will do a sum-to-zero coding for each factor,
    ## i.e., each factor's coefficients are deviations from the grand mean (e.g., congruent = +1, incongruent = -1)
    ## rather than tratment coding (R default; e.g., congruent = 0, incongruent = 1)
    ## However, it only affects factors, so we still need to center age (continuous variable)
    accuracy_LRT <- afex::mixed(correct ~ congruence * frequency * face_emotion + log_age + sex + ethnicity
                     + (1 | SID), family = binomial(), data = dat, check_contrasts = TRUE, method = "LRT")

    ## Likelihood ratio test chi square, one row per model term
    ## Each test compares the full model against one with that term's column removed (Type 3)
    ## A factor with 3+ levels gets a single omnibus test here, rather than several coefficient rows
    print(accuracy_LRT)
    write_rds(accuracy_LRT, "./2-Task_Cond_Effect/Outputs/Acc_LRT_results.rds")

    result_accuracy <- afex::nice(accuracy_LRT)
    write_csv(result_accuracy, "./2-Task_Cond_Effect/Outputs/Acc_LRT_results.csv")
}


## Plots % predicted accuracy for congruence and age
Plots.Accuracy <- function(model, EST.data) {

    acc_perf <- read_rds(model)

    ## ---- Congruency ----
    ## Model-based accuracy per congruence level, averaged over other factors
    ## type = "response" back-transforms from log-odds to probability
    emm <- emmeans(acc_perf, ~ congruence, type = "response")
    emm_df <- as.data.frame(emm)   # columns: congruence, prob, SE, df, asymp.LCL, asymp.UCL

    cong_acc <- ggplot(emm_df, aes(x = congruence, y = prob * 100)) +
            geom_point(size = 4) +
            geom_errorbar(aes(ymin = asymp.LCL * 100, ymax = asymp.UCL * 100),
                          linewidth = 1, width = .1) +
            labs(x = "Congruence", y = "Predicted accuracy (%)") +
            theme_classic() +
            theme(axis.text = element_text(size = 18),
                  axis.title = element_text(size = 20, face = "bold"))

    ggsave("./2-Task_Cond_Effect/Outputs/Acc_Congruence.pdf", cong_acc, width = 5, height = 5)


    ## ---- Age ----

    EST <- read_rds(EST.data)
    ## Same centering constant as Center.Age() used when fitting the model
    age_center <- mean(log(EST$age[!duplicated(EST$SID)]))

    ## Grid of ages across the observed range (centered log age is what the model uses)
    ages <- seq(min(EST$age), max(EST$age), length.out = 100)
    log_age <- log(ages) - age_center

    ## Predicted probability of a correct response by age, averaged over other factors
    ## type = "response" back-transforms from log-odds to probability
    emm_age <- emmeans(acc_perf, ~ log_age, at = list(log_age = log_age), type = "response")
    ## columns: log_age, prob, SE, df, asymp.LCL, asymp.UCL
    pred_age <- as.data.frame(emm_age)

    ## Back to age in years
    pred_age$age <- exp(pred_age$log_age + age_center)

    age_plot <- ggplot(pred_age, aes(age, prob * 100)) +
            geom_line(colour = "#339999", linewidth = 1) +
            geom_ribbon(aes(ymin = asymp.LCL * 100, ymax = asymp.UCL * 100), fill = "#339999", alpha = .2) +
            labs(x = "Age (years)", y = "Predicted accuracy (%)") +
            theme_classic() +
            theme(axis.text = element_text(size = 18),
                  axis.title = element_text(size = 20, face = "bold"))

    ggsave("./2-Task_Cond_Effect/Outputs/Acc_Age.pdf", age_plot, width = 6, height = 4)
}


Compute.OR <- function(my.emm, label) {

    ci  <- as.data.frame(confint(my.emm$contrasts))
    lcl <- grep("LCL|lower.CL", names(ci), value = TRUE)
    ucl <- grep("UCL|upper.CL", names(ci), value = TRUE)

    ci$analysis <- label
    ci$OR <- round(ci$odds.ratio, 2)
    ci$OR_lo <- round(ci[[lcl]], 2)
    ci$OR_hi <- round(ci[[ucl]], 2)

    ci[, setdiff(names(ci), c("SE", "df", "odds.ratio", lcl, ucl))]
}

Odds.Ratio.Acc <- function(model) {

    acc_perf <- read_rds(model)

    specs <- c(
        "congruence",
        "frequency",
        "face_emotion",
        "congruence | frequency",
        "congruence | face_emotion",
        "face_emotion | frequency",
        "congruence | frequency * face_emotion",
        "sex",
        "ethnicity"
    )

    res <- lapply(specs, function(s) {
        ## revpairwise gives incongruent / congruent odds of a correct response
        contr <- if (grepl("^congruence", s)) "revpairwise" else "pairwise"
        emm <- emmeans(acc_perf, as.formula(paste(contr, "~", s)), type = "response")
        Compute.OR(emm, label = s)
    })

    ## For age, OR per 1 SD of log age
    ## If want per 1 log age, just exponentiate the coefficient
    b  <- unname(fixef(acc_perf)["log_age"])
    ci <- confint(acc_perf, parm = "log_age", method = "Wald")

    ## Participant-level SD (age is between-subject, so avoid trial-level rows)
    mf <- model.frame(acc_perf)
    sd_logage <- sd(mf$log_age[!duplicated(mf$SID)])

    res_age <- data.frame(
        analysis = "age",
        contrast = "+1 SD of log(age)",
        OR       = round(exp(b * sd_logage), 2),
        OR_lo    = round(exp(ci[1, 1] * sd_logage), 2),
        OR_hi    = round(exp(ci[1, 2] * sd_logage), 2)
    )

    all_res <- dplyr::bind_rows(c(res, list(res_age)))
    write_csv(all_res, "./2-Task_Cond_Effect/Outputs/Acc_Odds_Ratios.csv")
    return(all_res)
}


## To test if results with random intercept only are similar to random slope + intercept
Model.RT.RandomSlope <- function(data) {

    dat <- Prep.Data.RT(data)
    
    ## Set contrasts to sum-to-zero coding so that it matches mixed() output later
    options(contrasts = c("contr.sum", "contr.poly"))

    ## Compute RT log-scores
    dat$log_rt <- log(dat$reaction_time)

    ## Excluded trials + RT Log-scored
    RT_random_slope <- lmer(log_rt ~ congruence * frequency * face_emotion + log_age + sex + ethnicity 
                    + (1 + congruence + frequency + face_emotion | SID),
                    data = dat, REML = T)

    write_rds(RT_random_slope, "./2-Task_Cond_Effect/Outputs/Suppl_RT_RS_Model_results.rds")

    RT_Random_slope_Anova <- afex::mixed(log_rt ~ congruence * frequency * face_emotion + log_age + sex + ethnicity
                     + (1 + congruence + frequency + face_emotion | SID), data = dat, method = "S")

    my.ef <- effectsize::eta_squared(RT_Random_slope_Anova$full_model, partial = TRUE) |>
                dplyr::rename(term = Parameter)

    results_RT <- afex::nice(RT_Random_slope_Anova) |>
                    dplyr::rename(term = Effect) |>
                        dplyr::left_join(my.ef, by = "term")

    write.csv(results_RT, "./2-Task_Cond_Effect/Outputs/Suppl_RT_RS_Anova_results.csv")

    ## Three way interaction plot does not change between random intercept only model
    ## and random intercept + random slope model; got rid of the plotting part here
}


## To test if results with random intercept only are similar to random slope + intercept
Model.Acc.RandomSlope <- function(data) {

    dat <- Prep.Data.Acc(data)

    ## Set contrasts to sum-to-zero coding so that it matches mixed() output later
    ## Age is centered in data prep
    options(contrasts = c("contr.sum", "contr.poly"))

    ## Accuracy glmm model with random slope and random intercept
    accuracy_random_slope <- glmer(correct ~ congruence * frequency * face_emotion + log_age + sex + ethnicity
                    + (1 + congruence + face_emotion + frequency | SID),
                    data = dat, family = binomial(link = "logit"))
    
    ## save model as rds for model comp
    write_rds(accuracy_random_slope, "./2-Task_Cond_Effect/Outputs/Suppl_Acc_RS_Model_results.rds")

    ## Likelihood ratio part
    accuracy_LRT <- afex::mixed(correct ~ congruence * frequency * face_emotion + log_age + sex + ethnicity
                    + (1 + congruence + face_emotion + frequency | SID), 
                    family = binomial(), data = dat, check_contrasts = TRUE, method = "LRT")

    ## save results as csv for reporting
    result_accuracy <- afex::nice(accuracy_LRT)
    write_csv(result_accuracy, "./2-Task_Cond_Effect/Outputs/Suppl_Acc_RS_LRT_results.csv")


    ## ---- Plot three way ----
    ## Three way interaction between accuracy and congruence/frequency/face emotion
    ## type = "response" back-transforms from log-odds to probability
    emm <- emmeans(accuracy_random_slope, ~ congruence * frequency * face_emotion, type = "response")
    emm_df <- as.data.frame(emm)

    cong_acc <- ggplot(emm_df, aes(x = congruence, y = prob * 100, color = frequency, group = frequency)) +
                    geom_point(size = 4) +
                    geom_line(linewidth = 1) +
                    geom_errorbar(aes(ymin = asymp.LCL * 100, ymax = asymp.UCL * 100),
                                linewidth = 1, width = .1) +
                    facet_wrap(~ face_emotion) +
                        ylab("Predicted RT (ms)") +
                        labs(x = "Congruence", y = "Predicted accuracy (%)", color = "Frequency") +
                    theme_classic() +
                    theme(legend.position = "bottom",
                            axis.text = element_text(size = 18),
                            axis.title = element_text(size = 20, face = "bold"),
                            strip.text.x = element_text(size = 20),
                            legend.text = element_text(size = 18),
                            legend.title = element_text(size = 20, face = "bold"),
                            panel.grid.minor = element_blank(),
                            strip.background = element_rect(fill = "grey95")) +
                    scale_color_manual(values = c("#E3C454FF", "#A27CBAFF"), labels = c("Equal", "More congruent"))

    ggsave("./2-Task_Cond_Effect/Outputs/Suppl_Acc_Congruence.pdf", cong_acc, width = 7, height = 7)
}

OR.Acc.Suppl <- function(model_RS) {

    accuracy_random_slope <- read_rds(model_RS)

    emm <- emmeans(accuracy_random_slope, ~ congruence | frequency * face_emotion)

    cong_or <- contrast(emm, "revpairwise") 
    summary(cong_or, type = "response", infer = TRUE)

    pairs(cong_or, by = "congruence")
}


## This functions fits and tests the model for the interaction 
## between age and task effects
Model.RT.Interaction.Age <- function(data) {

    ## Model fitting part
    dat <- Prep.Data.RT(data)
    
    ## Set contrasts to sum-to-zero coding so that it matches mixed() output later
    options(contrasts = c("contr.sum", "contr.poly"))

    dat$log_rt <- log(dat$reaction_time)

    rt_age <- lmer(log_rt ~ congruence * frequency * face_emotion * log_age + sex + ethnicity 
                    + (1 | SID),
                    data = dat, REML = T)

    write_rds(rt_age, "./2-Task_Cond_Effect/Outputs/Age_RT_Model_results.rds")

    ## ANOVA part
    RT_Anova <- afex::mixed(log_rt ~ congruence * frequency * face_emotion * log_age + sex + ethnicity
                     + (1 | SID), data = dat, method = "S")

    my.ef <- effectsize::eta_squared(RT_Anova$full_model, partial = TRUE) |>
                dplyr::rename(term = Parameter)

    results_RT <- afex::nice(RT_Anova) |>
                    dplyr::rename(term = Effect) |>
                        dplyr::left_join(my.ef, by = "term")

    write.csv(results_RT, "./2-Task_Cond_Effect/Outputs/Age_RT_Anova_results.csv")

    ## no follow up or plot because no results significant
}


Model.Acc.Interaction.Age <- function(data) {

    dat <- Prep.Data.Acc(data)

    ## Set contrasts to sum-to-zero coding so that it matches mixed() output later
    ## Age is centered in data prep
    options(contrasts = c("contr.sum", "contr.poly"))

    accuracy_age <- glmer(correct ~ congruence * frequency * face_emotion * log_age + sex + ethnicity
                    + (1 | SID),
                    data = dat, family = binomial(link = "logit"))

    write_rds(accuracy_age, "./2-Task_Cond_Effect/Outputs/Age_Acc_Model_results.rds")

    ## Calls function for initial model and assumptions
    Modelling.Accuracy(data)

    ## re-fit model with afex to get LRT
    ## check_contrast = TRUE will do a sum-to-zero coding for each factor,
    ## i.e., each factor's coefficients are deviations from the grand mean (e.g., congruent = +1, incongruent = -1)
    ## rather than tratment coding (R default; e.g., congruent = 0, incongruent = 1)
    ## However, it only affects factors, so we still need to center age (continuous variable)
    accuracy_LRT <- afex::mixed(correct ~ congruence * frequency * face_emotion * log_age + sex + ethnicity
                     + (1 | SID), family = binomial(), data = dat, check_contrasts = TRUE, method = "LRT")

    ## Likelihood ratio test chi square, one row per model term
    ## Each test compares the full model against one with that term's column removed (Type 3)
    ## A factor with 3+ levels gets a single omnibus test here, rather than several coefficient rows
    print(accuracy_LRT)
    write_rds(accuracy_LRT, "./2-Task_Cond_Effect/Outputs/Age_Acc_LRT_results.rds")

    result_accuracy <- afex::nice(accuracy_LRT)
    write_csv(result_accuracy, "./2-Task_Cond_Effect/Outputs/Age_Acc_LRT_results.csv")

    ## Odds ratio part
    ## compute sd age
    mf <- model.frame(accuracy_age)
    sd_age <- sd(mf$log_age[!duplicated(mf$SID)])

    ## How big is the frequency effect at -1 SD, mean, and +1 SD of age?
    emm_freq_age <- emmeans(accuracy_age, pairwise ~ frequency | log_age,
                            at = list(log_age = c(-sd_age, 0, sd_age)), type = "response")

    ## 1st line is -1SD, second one mean, third one +1SD
    res_freq_age <- Compute.OR(emm_freq_age, label = "frequency | log age (-1 SD, mean, +1 SD)")
    write_csv(res_freq_age, "./2-Task_Cond_Effect/Outputs/Age_Acc_OR_results.csv")
}


Plot.Accuracy.Age <- function(EST.data, model.acc.age) {

    EST <- read_rds(EST.data)
    ## Same centering constant as Center.Age() used when fitting the model
    age_center <- mean(log(EST$age[!duplicated(EST$SID)]))

    ## Grid of ages across the observed range (centered log age is what the model uses)
    ages <- seq(min(EST$age), max(EST$age), length.out = 100)
    log_age <- log(ages) - age_center

    acc_perf <- read_rds(model.acc.age)

    ## Predicted probability of a correct response by age, averaged over other factors
    ## type = "response" back-transforms from log-odds to probability
    emm_age <- emmeans(acc_perf, ~ frequency | log_age, at = list(log_age = log_age), type = "response")
    ## columns: log_age, prob, SE, df, asymp.LCL, asymp.UCL
    pred_age <- as.data.frame(emm_age)

    ## Back to age in years
    pred_age$age <- exp(pred_age$log_age + age_center)

    freq_age <- ggplot(pred_age, aes(x = age, y = prob * 100, colour = frequency, fill = frequency)) + 
                geom_line(linewidth = 1) +
                geom_ribbon(aes(ymin = asymp.LCL * 100, ymax = asymp.UCL * 100), alpha = .2, colour = NA) +
                labs(x = "Age (years)", y = "Predicted accuracy (%)") +
                theme_classic() +
                theme(axis.text = element_text(size = 18),
                        axis.title = element_text(size = 20, face = "bold"),
                        legend.text = element_text(size = 18),
                        legend.title = element_text(size = 20, face = "bold"),
                        legend.position = "bottom") +
                scale_color_manual(values = c("#E3C454FF", "#A27CBAFF"), name = "Frequency", labels = c("Equal", "More congruent")) +
                scale_fill_manual(values = c("#E3C454FF", "#A27CBAFF"), name = "Frequency", labels = c("Equal", "More congruent"))
    
    ggsave("./2-Task_Cond_Effect/Outputs/Age_Acc_Freq.pdf", freq_age, width = 7, height = 6)
}


Model.RT.Interaction.Sex <- function(data) {

    ## Model fitting part
    dat <- Prep.Data.RT(data)
    
    ## Set contrasts to sum-to-zero coding so that it matches mixed() output later
    options(contrasts = c("contr.sum", "contr.poly"))

    dat$log_rt <- log(dat$reaction_time)

    rt_sex <- lmer(log_rt ~ congruence * frequency * face_emotion * sex + log_age + ethnicity 
                    + (1 | SID),
                    data = dat, REML = T)

    write_rds(rt_sex, "./2-Task_Cond_Effect/Outputs/Sex_RT_Model_results.rds")

    ## ANOVA part
    RT_Anova <- afex::mixed(log_rt ~ congruence * frequency * face_emotion * sex + log_age + ethnicity
                     + (1 | SID), data = dat, method = "S")

    my.ef <- effectsize::eta_squared(RT_Anova$full_model, partial = TRUE) |>
                dplyr::rename(term = Parameter)

    results_RT <- afex::nice(RT_Anova) |>
                    dplyr::rename(term = Effect) |>
                        dplyr::left_join(my.ef, by = "term")

    write.csv(results_RT, "./2-Task_Cond_Effect/Outputs/Sex_RT_Anova_results.csv")

    ## % change part
    emm_options(lmer.df = "asymptotic")
    emm <- emmeans(rt_sex, pairwise ~ face_emotion | sex, tran = "log", type = "response")

    cat("\nRT EMM result face emotion | sex:\n")
    print(emm)

    ci <- Compute.CI(emm, label = s)
    write.csv(ci, "./2-Task_Cond_Effect/Outputs/Sex_RT_Percent_Change.csv")

    ## Plot sex * face emotion interaction
    ## Save as df for plotting, otherwise won't work
    emm_df <- as.data.frame(emm$emmeans)

    p_sex_emo <- ggplot(emm_df, aes(x = face_emotion, y = response, color = sex, group = sex)) +
                        geom_point(size = 4) +
                        scale_y_continuous(breaks = seq(700, 1100, by = 25)) +
                        geom_line(linewidth = 1) +
                        geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), linewidth = 1, width = .1) +
                        labs(x = "Face emotion", y = "Predicted RT (ms)", color = "Sex") +
                        theme_classic() +
                        theme(axis.text = element_text(size = 18),
                                axis.title = element_text(size = 20, face = "bold"),
                                legend.text = element_text(size = 18),
                                legend.title = element_text(size = 20, face = "bold"),
                                legend.position = "bottom") +
                        scale_color_manual(values = c("#A27CBAFF", "#E3C454FF"), labels = c("Female", "Male"))
   
    ggsave("./2-Task_Cond_Effect/Outputs/Sex_RT_Face_emotion.pdf", p_sex_emo, width = 7, height = 7)
}


Model.Acc.Interaction.Sex <- function(data) { 

    ## Accuracy * sex interaction (no need for followup?)
    ## dat already has log age??
    dat <- Prep.Data.Acc("./0-Data_Prep/Outputs/Tidy_all.rds")
    options(contrasts = c("contr.sum", "contr.poly"))

    acc_sex <- glmer(correct ~ congruence * frequency * face_emotion * sex + log_age + ethnicity
                     + (1 | SID),
                     data = dat, family = binomial(link = "logit"))
    write_rds(acc_sex, "./2-Task_Cond_Effect/Outputs/Sex_Acc_Model_results.rds")

    LRT_sex <- afex::mixed(correct ~ congruence * frequency * face_emotion * sex + log_age + ethnicity
                      + (1 | SID), family = binomial(), data = dat, check_contrasts = TRUE, method = "LRT")

    print(LRT_sex)
    result_accuracy <- afex::nice(LRT_sex)
    write_csv(result_accuracy, "./2-Task_Cond_Effect/Outputs/Sex_Acc_LRT_results.csv")

    emm <- emmeans(acc_sex, pairwise ~ face_emotion | frequency * sex, type = "response")

    cat("\nAccuracy EMM result sex | frequency * face emotion:\n")
    print(emm)

    or <- Compute.OR(emm, label = "face_emotion | frequency * sex")
    write.csv(or, "./2-Task_Cond_Effect/Outputs/Sex_Acc_OR.csv")

    emm_df <- as.data.frame(emm$emmeans)

    p_sex_acc <- ggplot(emm_df, aes(x = face_emotion, y = prob * 100, color = frequency, group = frequency)) +
                    geom_point(size = 4) +
                    geom_line(linewidth = 1) +
                    geom_errorbar(aes(ymin = asymp.LCL * 100, ymax = asymp.UCL * 100),
                                linewidth = 1, width = .1) +
                    facet_wrap(~ sex) +
                        labs(x = "Face emotion", y = "Predicted accuracy (%)", color = "Frequency") +
                    theme_classic() +
                    theme(legend.position = "bottom",
                            axis.text = element_text(size = 18),
                            axis.title = element_text(size = 20, face = "bold"),
                            strip.text.x = element_text(size = 20),
                            legend.text = element_text(size = 18),
                            legend.title = element_text(size = 20, face = "bold"),
                            panel.grid.minor = element_blank(),
                            strip.background = element_rect(fill = "grey95")) +
                    scale_color_manual(values = c("#E3C454FF", "#A27CBAFF"), labels = c("Equal", "More congruent"))

    ggsave("./2-Task_Cond_Effect/Outputs/Sex_Acc_Freq_Face_emotion.pdf", p_sex_acc, width = 7, height = 8)
}