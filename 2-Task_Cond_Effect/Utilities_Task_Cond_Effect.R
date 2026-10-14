## Utilities_Task_Cond_Effect.R
## Author: Camille Grandé
## Description:
##      This script holds the utility functions to run the stats to test the task conditions effect for the emotional stroop task analysis


Modelling.RT <- function(data) {

    dat <- read_rds(data) |>
                ## Keeping only RT for correct trials
                subset(correct == "correct") 

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
                main = "Residuals against observed RT (Initial model)")
        # model 2
        plot(fitted(model2), residuals(model2),
                ylab = "Residuals", xlab = "log-scored Reaction Time",
                main = "Residuals against log-scored RT (log model)")
        # normality residuals
        hist(residuals(model2), 
                ylab = "Frequency", xlab = "Residuals",
                main = "Histogram of residuals")
        qqnorm(residuals(model2),
                main = "QQ-plot of residuals")
    dev.off()

    return(model2)
}


ANOVA.RT <- function(data) {

    model <- Modelling.RT(data)

    my.anova <- anova(model)
    my.ef <- eta_squared(my.anova, partial = TRUE) |> rename(term = Parameter)
    
    summary <- my.anova |> rownames_to_column("term") |> left_join(my.ef, by = "term")

    write.csv(summary, "./2-Task_Cond_Effect/Outputs/RT_Anova_results.csv")
}



Modelling.Accuracy <- function(data) {

    dat <- read_rds(data)

    dat <- dat |> 
                mutate(correct = case_when(
                    correct == "correct"   ~ 1,
                    correct == "incorrect" ~ 0
                ), .keep = "unused")
    
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


ChiSq.Accuracy <- function(data) {

    model <- Modelling.Accuracy(data)

    ## computes odds ratios and their CI
    ## reports coefficient comparing a specific level to the reference level
    OR_table <- tidy(model, effects = "fixed") |>
                    mutate(
                        OR = exp(estimate),
                        OR_low = exp(estimate - 1.96 * std.error),
                        OR_high = exp(estimate + 1.96 * std.error))

    write.csv(OR_table, "./2-Task_Cond_Effect/Outputs/Accuracy_OR_results.csv")

    ## Wald tests (χ² distribution)
    ## tests that the tested coefficient is zero, given everything else in the model?


    ## tests overall effect of the factor
    ## here we are doing a likelihood ratio test where we perform model comparison
    ## on the full model vs alternative model which excludes each predictor
    ## is is LR or chisquare or same thing?
    ## Wald test that the relevant coefficients are jointly zero, given everything else in the model.
    chisq_table <- car::Anova(model, type = 3) |>
        as.data.frame() |>
            tibble::rownames_to_column("term")    
    
    write.csv(chisq_table, "./2-Task_Cond_Effect/Outputs/Accuracy_ChiSq_results.csv")
}


## dat$var <- scale(dat$var) to center variables? only age?
library(afex)
accuracy_afex <- afex::mixed(correct ~ congruence * frequency * face_emotion + log_age + sex + ethnicity
                     + (1 | SID), family = binomial(), data = dat, method = "LRT")

## Wald z test
## one row per coeff
## uses dummy coding for contrasts, each coeff is a difference from the reference level
summary(accuracy)
## report result of accuracy_afex or summary(accuracy_afex) ?

## LRT chi square
## one row per model term
## Each test compares the full model against one with that term’s columns removed (Type 3)
## A factor with 3+ levels gets a single omnibus test here, rather than several coefficient rows.
## this is the table we use for reporting
accuracy_afex

## Wald z test, not LRT!!
## summary of accuracy_afex$full_model
## but refits the model with sum to zero contrasts
## Each coefficient is a deviation from the grand mean (e.g., congruence1, not congruenceIncongruent)
summary(accuracy_afex)



## to do:
##      still ORs or not needed? OR for RT too?
##      get rid of car::Anova and implement afex::mixed for accuracy
##      afex:mixed for RT too?
##      can save output as csv?
##      contrasts?
##      scaling (centering) of variables? Age only?
##      coefficients RT?