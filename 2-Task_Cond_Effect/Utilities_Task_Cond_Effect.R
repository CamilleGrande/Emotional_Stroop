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

## Function for initial RT model and model assumptions
Modelling.RT <- function(data) {

    dat <- read_rds(data) |>
                ## Keeping only RT for correct trials
                subset(correct == "correct") |>
                    ## Centering log age for easier interpretation
                    ## The intercept will now refer to an average-aged participant
                    Center.Age()
    
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

    write_rds(model2, "./2-Task_Cond_Effect/Outputs/RT_Final_Model_results.rds")


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
ANOVA.RT <- function(data) {

    ## Calls function for initial model and assumptions
    my.model <- Modelling.RT(data)

    ## Reformatting data for mixed() models like in Modelling.RT function
    dat <- read_rds(data) |>
                ## Keeping only RT for correct trials
                subset(correct == "correct") |>
                    ## Centering log age for easier interpretation
                    ## The intercept will now refer to an average-aged participant
                    Center.Age()

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

    ## Plotting
    data_plot <- read_rds(data)
    Plots.RT(my.model, data_plot)
}


Plots.RT <- function(model, EST.data) {

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

    pdf("./2-Task_Cond_Effect/Outputs/RT_Three_way_interaction.pdf")
        ggplot(emm_df, aes(
                    x = congruence,
                    y = RT_ms,
                    color = frequency,
                    group = frequency
                )) +
                geom_point(size = 4) +
                scale_y_continuous(breaks = seq(800, 1000, by = 25)) +
                geom_line(linewidth = 1) +
                geom_errorbar(aes(ymin = lower_ms, ymax = upper_ms), linewidth = 1, width = .1) +
                facet_wrap(~ face_emotion) +
                ylab("Predicted RT (ms)") +
                labs(color = "Frequency") +
                theme_classic() +
                theme(axis.text = element_text(size = 12),
                        axis.title = element_text(size = 14, face = "bold"),
                        strip.text.x = element_text(size = 12),
                        legend.text = element_text(size = 12)) +
                scale_color_manual(values = c("#E3C454FF", "#A27CBAFF"), labels = c("Equal", "More frequent"))
    dev.off()


    ## ---- Plot predicted RT for age ----
    EST <- read_rds(EST.data)

    ## Get predicted values for log age from our model
    pred_age <- ggpredict(RT_perf, terms = "log_age [all]")

    ## Predicted age from log centered age in our model
    pred_age$age <- exp(pred_age$x + mean(log(EST$age)))

    ## Exponentiate coefficients as are on the log scale (log RT)
    pred_age$RT_ms <- exp(pred_age$predicted)
    pred_age$lower_ms <- exp(pred_age$conf.low)
    pred_age$upper_ms <- exp(pred_age$conf.high)

    pdf("./2-Task_Cond_Effect/Outputs/RT_Age.pdf")
        ggplot(pred_age, aes(age, RT_ms)) +
            geom_line() +
            geom_ribbon(
                aes(ymin = lower_ms,
                    ymax = upper_ms),
                alpha = .2
            ) +
            labs(
                x = "Age (years)",
                y = "Predicted RT (ms)"
            ) +
            theme_classic()
    dev.off()
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
    ci$pct_change <- 100 * (ci$ratio - 1)
    ## order preserved, no flip as ratio - 1 is an increasing function
    ## before had 1 - ratio, which is a decreasing function so had to reverse lcl as uci and ucl as lci
    ci$pct_lo <- 100 * (ci[[lcl]] - 1)
    ci$pct_hi <- 100 * (ci[[ucl]] - 1)
    
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
        emm <- emmeans(RT_perf, as.formula(paste("pairwise ~", s)),
                      tran = "log", type = "response")
        Compute.CI(emm, label = s)
    })

    ## For age, % change in RT per doubling of age (assumes log_age = natural log)
    b  <- unname(fixef(RT_perf)["log_age"])
    ci <- confint(RT_perf, parm = "log_age", method = "Wald")
    res_age <- data.frame(
        analysis   = "age",
        contrast   = "doubling of age",
        pct_change = 100 * (2^b - 1),
        pct_lo     = 100 * (2^ci[1, 1] - 1),
        pct_hi     = 100 * (2^ci[1, 2] - 1)
    )

    all_res <- dplyr::bind_rows(c(res, list(res_age)))
    write_csv(all_res, "./2-Task_Cond_Effect/Outputs/RT_Percent_Change.csv")
    return(all_res)
}


## Function for initial accuracy model and model assumptions
Modelling.Accuracy <- function(data) {

    dat <- read_rds(data) |>
                mutate(correct = case_when(
                    correct == "correct"   ~ 1,
                    correct == "incorrect" ~ 0
                ), .keep = "unused") |>
                    ## Centering log age for easier interpretation
                    ## The intercept will now refer to an average-aged participant
                    Center.Age()

    ## Set contrasts to sum-to-zero coding so that it matches mixed() output later
    ## Age is centered in data prep
    options(contrasts = c("contr.sum", "contr.poly"))

    ## random intercept only as random slope increases dispersion
    accuracy <- glmer(correct ~ congruence * frequency * face_emotion + log_age + sex + ethnicity
                    + (1 | SID),
                    data = dat, family = binomial(link = "logit"))

    write_rds(accuracy, "./2-Task_Cond_Effect/Outputs/Accuracy_Model_results.rds")


    cat("\nChecking assumptions for our Binomial Generalized Linear Mixed Model\n")

    ## Check collinearity among fixed effects
    ## Computes VIF, 0-5 = low, 5-10 = moderate, > 10 = high multicollinearity
    collinearity <- check_collinearity(accuracy)
    cat("\nLook at csv to assess collinearity\n")
    write_csv(collinearity, "./2-Task_Cond_Effect/Outputs/Accuracy_Collinearity.csv")


    ## Simulation-based residual diagnostics (DHARMa)
    ## Calculates randomized quantile residuals 
    sim_res <- simulateResiduals(fittedModel = accuracy, n = 1000, plot = F)

    pdf("./2-Task_Cond_Effect/Outputs/Accuracy_Assumptions.pdf")
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
    dat <- read_rds(data) |>
                mutate(correct = case_when(
                    correct == "correct"   ~ 1,
                    correct == "incorrect" ~ 0
                ), .keep = "unused") |>
                    ## Centering log age for easier interpretation
                    ## The intercept will now refer to an average-aged participant
                    Center.Age()

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
    write_rds(accuracy_LRT, "./2-Task_Cond_Effect/Outputs/Accuracy_LRT_results.rds")

    result_accuracy <- afex::nice(accuracy_LRT)
    write_csv(result_accuracy, "./2-Task_Cond_Effect/Outputs/Accuracy_LRT_results.csv")
}