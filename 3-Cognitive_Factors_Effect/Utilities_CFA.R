## Utilities_CFA.R
## Author: Camille Grandé
## Description:
##      This script holds the utility functions to run the cognition CFA analysis in the TAY sample


CFA.Analysis <- function(cog.data) {

    cog <- read_rds(cog.data)

    efa_vars <- c("NIH_Picture_Vocab", "NIH_Reading_Recog", "NIH_Flanker", "NIH_Card_Sort", "NIH_Pattern_Comp",
                    "LMT_Efficiency", "NIH_Picture_Seq", "RAVLT_Total", "NIH_List_Sort")

    #efa_vars <- c("NIH_Picture_Vocab", "NIH_Reading_Recog", "NIH_Flanker", "NIH_Card_Sort", "NIH_Pattern_Comp",
    #                "LMT_Efficiency", "NIH_Picture_Seq", "RAVLT_Total", "NIH_List_Sort", "MATRICES")
    
    ## Get all subjects with no NA
    tay_pca_full <- cog[complete.cases(cog[, efa_vars]), ]

    ## Some missing data is coded as 0
    tay_pca_full$RAVLT_Total[tay_pca_full$RAVLT_Total == 0] <- NA
    tay_pca_full$LMT_Efficiency[tay_pca_full$LMT_Efficiency == 0] <- NA

    ## should be 879 for july 2026 datacut
    nrow(tay_pca_full)
    ## 2 missing for RAVLT
    sum(is.na(tay_pca_full$RAVLT_Total))
    varTable(tay_pca_full[, efa_vars])

    #Specifiy and run CFA Model for three factor structure consisting of Verbal Ability, Executive Functioning and Learning/Memory# - Note that Matrix cross loads onto F1 and F2# 
    model_3factor <- '
        f1 =~ NIH_Picture_Vocab + NIH_Reading_Recog
        f2 =~ NIH_Flanker + NIH_Card_Sort + NIH_Pattern_Comp + LMT_Efficiency
        f3 =~ NIH_Picture_Seq + RAVLT_Total + NIH_List_Sort'
    
    #model_3factor <- '
    #    f1 =~ NIH_Picture_Vocab + NIH_Reading_Recog + MATRICES
    #    f2 =~ NIH_Flanker + NIH_Card_Sort + NIH_Pattern_Comp + LMT_Efficiency
    #    f3 =~ NIH_Picture_Seq + RAVLT_Total + NIH_List_Sort + MATRICES'

    fit_3factor <- cfa(model_3factor, data = tay_pca_full, missing = "fiml")
    fs <- as.data.frame(lavPredict(fit_3factor))

    fitMeasures(fit_3factor, c("chisq", "df", "cfi", "tli", "rmsea", "srmr", "aic", "bic"))

    #check N - should be 878 
    lavInspect(fit_3factor, "nobs")


    fs <- as.data.frame(lavPredict(fit_3factor))
    fs$PIN <- tay_pca_full$PIN[lavInspect(fit_3factor, "case.idx")]
    nrow(fs) 

    write.csv(fs, "./3-Cognitive_Factors_Effect/Outputs/Factor_Scores.csv", row.names = FALSE)
}
