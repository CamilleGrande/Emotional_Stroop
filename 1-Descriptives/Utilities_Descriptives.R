## Utilities_Descriptives.R
## Author: Camille Grandé
## Description:
##      This script holds the utility functions to run the descriptives stats for the emotional stroop task data


Subject.Level.Descr <- function(my.df, ...) {

    summary <- my.df |> 
                    subset(reaction_time > 200 & reaction_time < 2000) |>
                    group_by(...) |>
                        summarise(Mean_Accuracy = sum(correct == "correct") / sum(!is.na(correct)),
                                    Mean_RT = sum(reaction_time) / sum(!is.na(reaction_time)),
                                    .groups = "drop")
    return(summary)
}


Cohort.Level.Descr <- function(subj.level) {

    summary <- as.data.frame(subj.level) |>
                    summarise(Accuracy = round(mean(Mean_Accuracy) * 100, 2),
                                SD_Accuracy = round(sd(Mean_Accuracy) *100, 2),
                                RT = round(mean(Mean_RT), 2),
                                SD_RT = round(sd(Mean_RT), 2),
                                .groups = "drop")
    
    return(summary)
}



EST.Descriptives <- function(data) {

    my.data <- read_rds(data)

    ## ---- Global mean accuracy and RT ----
    all_subj <- Subject.Level.Descr(my.data, SID)
    all_cohort <- Cohort.Level.Descr(all_subj)


    ## ---- Mean accuracy and RT depending on congruency ----              
    congruency <- Subject.Level.Descr(my.data, SID, congruence)

    ## Separating congruent VS incongruent subj summary levels
    congruent_subj <- subset(congruency, congruence == "congruent")
    incongruent_subj <- subset(congruency, congruence == "incongruent")

    congruent_cohort <- Cohort.Level.Descr(congruent_subj)
    incongruent_cohort <- Cohort.Level.Descr(incongruent_subj)


    ## ---- Mean accuracy and RT depending on frequency ----               
    frequency <- Subject.Level.Descr(my.data, SID, frequency)

    ## Separating equal VS more congruent subj summary levels
    equal_subj <- subset(frequency, frequency == "equal")
    more_cong_subj <- subset(frequency, frequency == "more_cong")

    equal_cohort <- Cohort.Level.Descr(equal_subj)
    more_cong_cohort <- Cohort.Level.Descr(more_cong_subj)


    ## ---- Mean accuracy and RT depending on distractor valence ----               
    face_emotion <- Subject.Level.Descr(my.data, SID, face_emotion)

    ## Separating angry VS happy distractors (faces) subj summary levels
    angry_subj <- subset(face_emotion, face_emotion == "angry")
    happy_subj <- subset(face_emotion, face_emotion == "happy")

    angry_cohort <- Cohort.Level.Descr(angry_subj)
    happy_cohort <- Cohort.Level.Descr(happy_subj)


    ## ---- Mean accuracy and RT depending on all 8 combinations ----               
    across <- Subject.Level.Descr(my.data, SID, congruence, frequency, face_emotion)

    ## Separating all 8 combinations
    Cong_Eq_Angry <- subset(across, congruence == "congruent" & frequency == "equal" & face_emotion == "angry")
    Incong_Eq_Angry <- subset(across, congruence == "incongruent" & frequency == "equal" & face_emotion == "angry")
    Cong_Eq_Happy <- subset(across, congruence == "congruent" & frequency == "equal" & face_emotion == "happy")
    Incong_Eq_Happy <- subset(across, congruence == "incongruent" & frequency == "equal" & face_emotion == "happy")
    Cong_MC_Angry <- subset(across, congruence == "congruent" & frequency == "more_cong" & face_emotion == "angry")
    Incong_MC_Angry <- subset(across, congruence == "incongruent" & frequency == "more_cong" & face_emotion == "angry")
    Cong_MC_Happy <- subset(across, congruence == "congruent" & frequency == "more_cong" & face_emotion == "happy")
    Incong_MC_Happy <- subset(across, congruence == "incongruent" & frequency == "more_cong" & face_emotion == "happy")

    Cong_Eq_Angry_cohort <- Cohort.Level.Descr(Cong_Eq_Angry)
    Incong_Eq_Angry_cohort <- Cohort.Level.Descr(Incong_Eq_Angry)
    Cong_Eq_Happy_cohort <- Cohort.Level.Descr(Cong_Eq_Happy)
    Incong_Eq_Happy_cohort <- Cohort.Level.Descr(Incong_Eq_Happy)
    Cong_MC_Angry_cohort <- Cohort.Level.Descr(Cong_MC_Angry)
    Incong_MC_Angry_cohort <- Cohort.Level.Descr(Incong_MC_Angry)
    Cong_MC_Happy_cohort <- Cohort.Level.Descr(Cong_MC_Happy)
    Incong_MC_Happy_cohort <- Cohort.Level.Descr(Incong_MC_Happy)
 

    ## ---- Final summaries ----
    my_summaries <- c("All", "Congruent", "Incongruent", "Equal", "More_Congruent", "Angry", "Happy", 
                        "Cong_Eq_Angry", "Incong_Eq_Angry", "Cong_Eq_Happy", "Incong_Eq_Happy",
                                "Cong_MC_Angry", "Incong_MC_Angry", "Cong_MC_Happy", "Incong_MC_Happy")

    summary_cohort <- bind_rows(all_cohort, congruent_cohort, incongruent_cohort, equal_cohort, more_cong_cohort, angry_cohort, happy_cohort,
                                Cong_Eq_Angry_cohort, Incong_Eq_Angry_cohort, Cong_Eq_Happy_cohort, Incong_Eq_Happy_cohort,
                                Cong_MC_Angry_cohort, Incong_MC_Angry_cohort, Cong_MC_Happy_cohort, Incong_MC_Happy_cohort)
    
    summary_cohort <- bind_cols(my_summaries, summary_cohort)

    write.csv(summary_cohort, "./1-Descriptives/Outputs/Summary_Cohort_Mean_Acc_RT.csv", row.names = F)
}


## This function creates a numeric summary for all our variables (keep numeric only first)
Numeric.Summary <- function(vars) {

  numeric_summary <- vars |>
                        summarise(across(everything(),
                                    list(
                                            n = ~ sum(!is.na(.x), na.rm = T),
                                            mean = ~ round(mean(.x, na.rm = T), 2),
                                            sd = ~ round(sd(.x, na.rm = T), 2),
                                            median = ~ round(median(.x, na.rm = T), 2),
                                            iqr = ~ round(IQR(.x, na.rm = T), 2),
                                            min = ~ round(min(.x, na.rm = T), 2),
                                            max = ~ round(max(.x, na.rm = T), 2),
                                            skewness = ~ round(as.numeric(skewness(.x, na.rm = T)), 2), 
                                            kurtosis = ~ round(as.numeric(kurtosis(.x, na.rm = T)), 2)
                                        ), .names = "{.col}.{.fn}"
                                    )) |>
                                        pivot_longer(cols = everything(), names_to = c("variable", ".value"), names_sep = "\\.")

  return(numeric_summary)
}


## This function creates a frequency summary for all our variables (keep factor only first)
Frequency.Summary <- function(vars) {

  frequency_summary <- vars |>
                        pivot_longer(cols = -SID, names_to = "variable", values_to = "level") |>
                            group_by(variable) |>
                                mutate(total_n = n(),
                                        valid_n = sum(!is.na(level))) |>
                                    filter(!is.na(level)) |>
                                        count(variable, level, total_n, valid_n, name = "n") |>
                                            group_by(variable) |>
                                                mutate(pct_valid = round((n / valid_n * 100), 2),
                                                        pct_total = round((n / total_n * 100), 2),
                                                        pct_valid_cum = round(cumsum(pct_valid), 2),
                                                        pct_total_cum = round(cumsum(pct_total), 2)) |>
                                                    ungroup()
  
  return(frequency_summary)
}


## This function computes numeric/frequency summaries for the descriptive stats
Descriptives <- function(data.demog, data.cog, data.clinical, data.suicide) {

    ## ---- Demographics ----
    demog <- read_rds(data.demog) |>
                distinct(SID, .keep_all = TRUE)

    numeric_demog <- select(demog, all_of(c("age")))
    numeric_demog_summary <- Numeric.Summary(numeric_demog)
    
    frequency_demog <- select(demog, all_of(c("SID", "sex", "gender_diversity", "ethnicity", "family_income")))
    frequency_demog_summary <- Frequency.Summary(frequency_demog)


    ## ---- Cognition ----
    cog <- read_rds(data.cog) |>
                rename(SID = "PIN")
    
    ## Keep only subjects used in analysis
    cog <- filter(cog, cog$SID %in% demog$SID)

    numeric_cog <- select(cog, -all_of(c("SID")))
    numeric_cog_summary <- Numeric.Summary(numeric_cog)


    ## ---- Clinical symptoms ----
    clinical <- read_rds(data.clinical)

    ## Keep only subjects used in analysis
    clinical <- filter(clinical, clinical$SID %in% demog$SID)

    numeric_clinical <- select(clinical, -all_of(c("SID")))
    numeric_clinical_summary <- Numeric.Summary(numeric_clinical)


    ## ---- Suicide attempt ----
    suicide <- read_rds(data.suicide) |>
                mutate(cssrs_attempt_lt = as.factor(cssrs_attempt_lt), .keep = "unused")

    ## Keep only subjects used in analysis
    suicide <- filter(suicide, suicide$SID %in% demog$SID)

    frequency_suicide_summary <- Frequency.Summary(suicide)


    ## Bind all summaries
    all_numeric <- bind_rows(numeric_demog_summary, numeric_cog_summary, numeric_clinical_summary)
    all_frequency <- bind_rows(frequency_demog_summary, frequency_suicide_summary)

    write.csv(all_numeric, "./1-Descriptives/Outputs/Numeric_Descriptives.csv", row.names = F)
    write.csv(all_frequency, "./1-Descriptives/Outputs/Frequency_Descriptives.csv", row.names = F)
}