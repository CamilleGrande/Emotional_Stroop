## Utilities_DataPrep.R
## Author: Camille Grandé
## Description:
##      This script holds the utility functions to prep the emotional stroop task data


## This function cleans the EST data from the TAY dataset to keep only the variables we are interested in
Tidy.Data <- function(data) {

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
                                                    .keep = "unused")

    ## Ordering variables
    my.data <- my.data[, c("SID", "frequency", "congruence", "face_emotion", "word_emotion", "word", "correct_button", "response", "correct", "reaction_time",
                            "word_y", "word_position")]

    write_rds(my.data, "./0-Data_Prep/Outputs/Tidy_EST.rds")
}