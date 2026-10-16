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
    Modelling.RT(data)

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