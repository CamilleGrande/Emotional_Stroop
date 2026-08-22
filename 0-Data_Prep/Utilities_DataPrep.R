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
                        log_age = log(age),
                        sex = as.factor(assigned_sex_at_birth),
                        gender_diversity = as.factor(gender_diversity),
                        ethnicity = as.factor(ethnicity),
                        family_income = as.factor(family_income),
                        .keep = "none"
                    ) |> 
                        inner_join(ids, by = "SID")

    demog_after_exclusions <- demog |>
                                drop_na(sex, ethnicity) |>
                                    subset(sex != "Intersex") |>
                                        ## need to exclude "intersex" factor
                                        mutate(sex = as.character(sex), .keep = "unused") |>
                                            mutate(sex = as.factor(sex), .keep = "unused")

    sum_exclusions <- setdiff(demog, demog_after_exclusions) |>
                        summarise(
                            missing_sex = sum(is.na(sex)),
                            missing_ethnicity = sum(is.na(ethnicity)),
                            sex_intersex = sum(sex == "Intersex", na.rm = TRUE))

    write_csv(sum_exclusions, "./0-Data_Prep/Outputs/Sum_Excluded_Demographics.csv")


    all <- inner_join(EST, demog_after_exclusions, by = "SID")
    write_rds(all, "./0-Data_Prep/Outputs/Tidy_all.rds")
}


## Functions used in NIH data cleaning
Load.Data <- function(data) {
     dat <- as.data.frame(read.csv(data, na = c("NA", "PNA", "UNK", "Prefer not to answer", "MI", "NASK", "ASKU", " ", "", "NCO", "NI", "OTH")))
     return(dat)
}

Unite.Data <- function(data, new.name, col1, col2) {
    data |> 
        unite({{ new.name }}, {{ col1 }}, {{ col2 }}, na.rm = T, remove = T) |>
            mutate({{ new.name }} := as.integer({{ new.name }}))
}


## NIH toolbox data cleaning
Tidy.NIH <- function(data) {

    NIH <- Load.Data(data) |> 
            filter(str_detect(Inst, "NIH Toolbox") & !is.na(Uncorrected.Standard.Score) & InstrumentBreakoff == 2) |>
                select(all_of(c("PIN", "Inst", "Uncorrected.Standard.Score", "Age.Corrected.Standard.Score", "National.Percentile.age.adjusted"))) |>
                    arrange(PIN)

    ## Used to rename columns
    lookup <- c(NIH_Picture_Vocab = "NIH Toolbox Picture Vocabulary Test Age 3+ v2.1", 
                NIH_Reading_Recog = "NIH Toolbox Oral Reading Recognition Test Age 3+ v2.1",
                NIH_Pattern_Comp = "NIH Toolbox Pattern Comparison Processing Speed Test Age 7+ v2.1",
                NIH_List_Sort = "NIH Toolbox List Sorting Working Memory Test Age 7+ v2.1",
                NIH_Picture_Seq = "NIH Toolbox Picture Sequence Memory Test Age 8+ Form A v2.1")


    NIH <- NIH |>
            pivot_wider(id_cols = "PIN", names_from = "Inst", values_from ="Uncorrected.Standard.Score") |>
                rename(any_of(lookup)) |>
                    Unite.Data(NIH_Flanker, "NIH Toolbox Flanker Inhibitory Control and Attention Test Age 12+ v2.1", 
                            "NIH Toolbox Flanker Inhibitory Control and Attention Test Ages 8-11 v2.1") |>
                                Unite.Data(NIH_Card_Sort, "NIH Toolbox Dimensional Change Card Sort Test Age 12+ v2.1",
                                                "NIH Toolbox Dimensional Change Card Sort Test Ages 8-11 v2.1")

    write.csv(NIH, "./0-Data_Prep/Outputs/NIH_Raw_7.csv", row.names = F)

    return(NIH)
}


## RAVLT data cleaning
Tidy.RAVLT <- function(filename) {
    
    RAVLT <- Load.Data(filename)|>
                arrange(kcni_participant_id) |>
                    rename(PIN = "kcni_participant_id") |>
                        pivot_wider(id_cols = "PIN", names_from = "Item", values_from = "Raw.score") |>
                            mutate(across(starts_with("RAVLT"), as.integer)) |>
                                select(all_of(c("PIN", "RAVLT Short Delay Trial I", "RAVLT Short Delay Trial II", "RAVLT Short Delay Trial III", "RAVLT Short Delay Trial IV", 
                                        "RAVLT Short Delay Trial V", "RAVLT Short Delay List B Trial", "RAVLT Short Delay Trial VI", "RAVLT Long Delay Trial VII"))) |>
                                    rename(
                                        RAVLT_Short_1 = "RAVLT Short Delay Trial I", 
                                        RAVLT_Short_2 = "RAVLT Short Delay Trial II", 
                                        RAVLT_Short_3 = "RAVLT Short Delay Trial III",
                                        RAVLT_Short_4 = "RAVLT Short Delay Trial IV", 
                                        RAVLT_Short_5 = "RAVLT Short Delay Trial V", 
                                        RAVLT_List_B = "RAVLT Short Delay List B Trial",
                                        RAVLT_Short_PostB = "RAVLT Short Delay Trial VI", 
                                        RAVLT_Long = "RAVLT Long Delay Trial VII"
                                    )
    
    RAVLT_Total <- RAVLT |> 
                    select(all_of(c("RAVLT_Short_1", "RAVLT_Short_2", "RAVLT_Short_3", "RAVLT_Short_4", "RAVLT_Short_5", "RAVLT_Short_PostB", "RAVLT_Long"))) |>
                        apply(1, sum)

    RAVLT$RAVLT_Total <- RAVLT_Total

    ## Participants with no missing data:
    RAVLT |>
        drop_na() |>
            write.csv("./0-Data_Prep/Outputs/RAVLT.csv", row.names = F)
    
    return(RAVLT)
}


## Matrices test data cleaning
Tidy.Matrices <- function(filename, score.type, path) {

    lookup <- paste0("MATRIX_", score.type)

    Matrices <- Load.Data(filename) |>
                    arrange(kcni_participant_id) |>
                        rename(PIN = "kcni_participant_id") |>
                            select(all_of(c("PIN", "Subtest", {{ score.type }} ))) |>
                                filter(str_detect(Subtest, "Reasoning")) |>
                                    distinct(.keep_all = F) |>
                                        pivot_wider(id_cols = "PIN", names_from = "Subtest", values_from = {{ score.type }}) |>
                                            Unite.Data(MATRICES, "WISC-V-CDN Matrix Reasoning", "WAIS-IV-CDN Matrix Reasoning")

    ## Participants with no missing data:
    Matrices |>
        drop_na() |>
            write.csv(paste0("./0-Data_Prep/Outputs/Matrix_", {{ path }}, ".csv"), row.names = F)

    return(Matrices)
}


## Little man task data cleaning
Tidy.LMT <- function(filename) {

    LMT <- Load.Data(filename) |>
                select(all_of(c("kcni_participant_id", "propCorrect", "meanRT"))) |>
                    arrange(kcni_participant_id) |>
                        rename(
                            PIN = "kcni_participant_id", 
                            LMT_Prop_Accurate = "propCorrect", 
                            LMT_Mean_RT = "meanRT"
                        ) |>
                                distinct(PIN, .keep_all = T) |>
                                    ## computes efficiency score and creates new column to hold the values
                                    mutate(LMT_Efficiency = (LMT_Prop_Accurate / LMT_Mean_RT) * 1000)

    write.csv(LMT, "./0-Data_Prep/Outputs/LittleManTask.csv", row.names = F)

    return(LMT)
}


## Summarizes n participants with missing data for each test
Missing.Data.Cog <- function(NIH.Tidy, RAVLT.Tidy, Matrices.Tidy, LMT.Tidy) {

    summarise_missing <- function(data, test_name) {
        data |>
            summarise(
                test = test_name,
                n_total = n(),
                n_missing = sum(if_any(everything(), is.na)),
                pct_missing = round(100 * n_missing / n_total, 1)
            )
    }

    summary_table <- bind_rows(
        summarise_missing(NIH.Tidy, "NIH"),
        summarise_missing(RAVLT.Tidy, "RAVLT"),
        summarise_missing(Matrices.Tidy, "Matrices"),
        summarise_missing(LMT.Tidy, "LMT")
    )

    write.csv(summary_table, "./0-Data_Prep/Outputs/Sum_Missing_Cognition.csv", row.names = F)
    
}


## Tidies cognitive data
Tidy.Cognition <- function(NIH.Raw, RAVLT.Raw, Matrices.Raw, LMT.Raw) {

    NIH <- Tidy.NIH(NIH.Raw)
    RAVLT <- Tidy.RAVLT(RAVLT.Raw)
    Matrices <- Tidy.Matrices(Matrices.Raw, "Raw.score", "Raw")
    LMT <- Tidy.LMT(LMT.Raw)

    Missing.Data.Cog(NIH, RAVLT, Matrices, LMT)

    cog <- full_join(NIH, Matrices, by = "PIN") |>
            full_join(RAVLT, by = "PIN") |>
                full_join(LMT, by = "PIN")
    
    write_rds(cog, "./0-Data_Prep/Outputs/Tidy_cognition.rds")
}


## This function cleans the item "previous suicide attempt" of the C-SSRS
Tidy.CSSRS <- function(filename) {
   
    suicide <- Load.Data(filename) |>
                    rename(SID = "subject_id") |>
                            distinct(SID, .keep_all = T) |>
                                arrange(SID) |>
                                    select(all_of(c("SID", "cssrs_attempt_lt"))) |>
                                        mutate(as.integer(cssrs_attempt_lt), .keep = "unused")

    ## Participants with no missing data:
    suicide |>
        drop_na() |>
            write.csv("./0-Data_Prep/Outputs/CSSRS.csv", row.names = F)

    return(suicide)
}


## This function selects the appropriate score for each A/YSR file 
## because of coding difference in column names
AYSR.Select.Score <- function(data, score.type) {

    ASR_t <- c("asr_syndr_anxdepr_t", "asr_syndr_withdrdep_t", "asr_syndr_somatic_t", "asr_syndr_thought_t",
                    "asr_syndr_attention_t", "asr_syndr_aggressive_t", "asr_syndr_rulebreaking_t")
    
    ASR_t_GI <- c("asr_syndr_anxdepr_t", "asr_syndr_withdrdepr_t", "asr_syndr_somatic_t", "asr_syndr_thought_t",
                    "asr_syndr_attention_t", "asr_syndr_aggressive_t", "asr_syndr_rulebreaking_t")

    YSR_t <- c("ysr_syndr_anxdepr_t", "ysr_syndr_withdrdepr_t", "ysr_syndr_somatic_t", "ysr_syndr_thought_t", 
                    "ysr_syndr_attention_t", "ysr_syndr_aggressive_t", "ysr_syndr_rulebreaking_t")
    
    lookup <- c("withdrdepr_t" = "withdrdep_t", "withdrdepr_p" = "withdrdep_p")

    Load.Data(data) |>
        rename(SID = "subject_id") |>
            arrange(SID) |>
                ## Select our columns
                select(all_of(c("SID", get(score.type)))) |>
                    ## 9999 also = NA, so put all as NAs
                    mutate(across(-SID, ~ na_if(.x, 9999))) |>
                        ## cleans our names for next step: merging!
                        rename_with(~str_remove(.x, "asr_syndr_")) |>
                        rename_with(~str_remove(.x, "ysr_syndr_")) |>
                        rename(any_of(lookup))
    
}


## This function tidies the ASR and YSR in one clean dataset
Tidy.AYSR <- function(data.ASR, data.ASR.GI, data.YSR, data.YSR.GI) {

    ASR <- AYSR.Select.Score(data.ASR, "ASR_t")
    ASR_GI <- AYSR.Select.Score(data.ASR.GI, "ASR_t_GI")
    YSR <- AYSR.Select.Score(data.YSR, "YSR_t")
    YSR_GI <- AYSR.Select.Score(data.YSR.GI, "YSR_t")

    ## put together YSR and ASR data
    YASR <- rows_patch(ASR, YSR, by = "SID")

    ## add gender inclusive data
    GI <- rows_patch(YSR_GI, ASR_GI, by = "SID")
        YASR <- rows_upsert(YASR, GI, by = "SID") |>
                arrange(SID)

    YASR |>
        drop_na() |>
            write.csv("./0-Data_Prep/Outputs/YASR_t.csv", row.names = F)

    return(YASR)
}


## Summarizes n participants with missing data for each test
Missing.Data.Clin <- function(suicide.Tidy, aysr.Tidy) {

    summarise_missing <- function(data, test_name) {
        data |>
            summarise(
                test = test_name,
                n_total = n(),
                n_missing = sum(if_any(everything(), is.na)),
                pct_missing = round(100 * n_missing / n_total, 1)
            )
    }

    summary_table <- bind_rows(
        summarise_missing(suicide.Tidy, "Suicide_Attempt"),
        summarise_missing(aysr.Tidy, "AYSR")
    )

    write.csv(summary_table, "./0-Data_Prep/Outputs/Sum_Missing_Clinical.csv", row.names = F)
}


Tidy.Clinical <- function(CSSRS.Raw, ASR.raw, ASR.GI.raw, YSR.raw, YSR.GI.raw) {

    suicide <- Tidy.CSSRS(CSSRS.Raw)
    AYSR <- Tidy.AYSR(ASR.raw, ASR.GI.raw, YSR.raw, YSR.GI.raw)
   
    Missing.Data.Clin(suicide, AYSR)
    
    write_rds(AYSR, "./0-Data_Prep/Outputs/Tidy_clinical.rds")
    write_rds(suicide, "./0-Data_Prep/Outputs/Tidy_suicide_att.rds")
}