## Utilities_DataPrep.R
## Author: Camille Grandé
## Description:
##      This script holds the utility functions to prep the emotional stroop task data


## This function cleans the EST data from the TAY dataset to keep only the variables we are interested in
Tidy.EST <- function(data) {

    my.data <- read_csv(data) |>
                    select(all_of(c("kcni_participant_id", "congruence", "blockcode", "faceemotion", "stimulusitem3", "correctButton", "response", "latency", 
                                    "word_y", "wordposition"))) |>
                                        filter(blockcode == "testMC" | blockcode == "testE") |>
                                            drop_na() |>
                                                mutate(SID = kcni_participant_id,
                                                    ## 1 = Congruent, 2 = Incongruent
                                                    congruence = fct_recode(as.factor(congruence), "congruent" = "1", "incongruent" = "2"),
                                                    frequency = fct_recode(blockcode, "equal" = "testE", "more_cong" = "testMC"),
                                                    face_emotion = as.factor(faceemotion),
                                                    word = as.factor(stimulusitem3),
                                                    word_emotion = as.factor(case_when(stimulusitem3 == "joyful" | stimulusitem3 == "happy" | stimulusitem3 == "glad" ~ "positive",
                                                                                stimulusitem3 == "mad" | stimulusitem3 == "angry" | stimulusitem3 == "upset"  ~ "negative")),
                                                    correct_button = as.factor(correctButton),
                                                    response = as.factor(response),
                                                    correct = as.factor(case_when(as.character(correct_button) == as.character(response) ~ "correct",
                                                                            TRUE ~ "incorrect")),
                                                    reaction_time = as.numeric(latency),
                                                    word_y = as.factor(word_y),
                                                    word_position = as.factor(wordposition),
                                                    .keep = "unused") |>
                                                        ## One subject has 0% accuracy and another 16% so excluding them
                                                        subset(SID != "TAY01_CMH_00000248" & SID != "TAY01_CMH_00000543")
    
    ## Excluded anticipated and timed-out trials
    cat("\n\nExcluding anticipated trials (RT < 200ms) and time-out sessions (RT ceiling 2000ms)\n")

    to_exclude <- data.frame(excluded_timeout = sum(my.data$reaction_time >= 2000), excluded_anticipatory = sum(my.data$reaction_time < 200))
    write.csv(to_exclude, "./0-Data_Prep/Outputs/Sum_Excluded_Trials.csv", row.names = F)

    ## 822 trials should be excluded total (ceiling + floor)
    ## ceiling were also all incorrect
    dat_clean <- subset(my.data, my.data$reaction_time >= 200 & my.data$reaction_time < 2000)

        cat(paste0("\nCheck all trials where RT < 200ms or >= 2000ms (time-out) were correctly excluded:\n\nSum trials < 200ms: ",
        sum(dat_clean$reaction_time < 200),
        "\nSum trials >= 2000ms: ", 
        sum(dat_clean$reaction_time >= 2000), 
        "\n\nIf both not 0, there is a problem\n"))
                                                    
    ## Ordering variables
    dat_clean <- dat_clean[, c("SID", "frequency", "congruence", "face_emotion", "word_emotion", "word", "correct_button", "response", "correct", "reaction_time",
                            "word_y", "word_position")]

    return(dat_clean)
}


Tidy.Data <- function(data.EST, data.demog) {

    EST <- Tidy.EST(data.EST)

    ids <- as.data.frame(unique(EST$SID))
    colnames(ids) <- "SID"
    
    demog <- read_csv(data.demog, na = c("ASKU", "NA", "UNK", "PNA", "NI", "MI", "NASK", "Prefer not to answer", "Do not know")) |>
                    mutate(
                        SID = subject_id,
                        age = as.numeric(age_enrolment),
                        sex = as.factor(assigned_sex_at_birth),
                        gender_diversity = as.factor(gender_diversity),
                        ethnicity = as.factor(ethnicity),
                        family_income = as.factor(family_income),
                        .keep = "none"
                    ) |> 
                        inner_join(ids, by = "SID")

    demog_after_exclusions <- demog |>
                                drop_na(sex, ethnicity)

    sum_exclusions <- setdiff(demog, demog_after_exclusions) |>
                        mutate(
                            missing_sex = sum(is.na(sex)),
                            missing_ethnicity = sum(is.na(ethnicity)),
                            .keep = "none")

    write_csv(sum_exclusions, "./0-Data_Prep/Outputs/Sum_Excluded_Demographics.csv")


    all <- inner_join(EST, demog_after_exclusions, by = "SID")
    write_rds(all, "./0-Data_Prep/Outputs/Tidy_all.rds")
}