mcnum <- 2L # Number of iterations
source("./global.R")
IMPACTncd <- Simulation$new("./inputs/sim_design.yaml")

IMPACTncd$
    del_logs()$
    del_outputs()$
    run(1:mcnum, multicore = TRUE, "sc0")

# F&V ----
IMPACTncd$update_primary_prevention_scn(
    function(synthpop) {
        setyear <- 2001L
        tbl <-
            read_fst("./inputs/exposure_distributions/Table_Fruit_vege.fst",
                as.data.table = TRUE
            )[Year == setyear & between(Age, self$design$sim_prm$ageL - self$design$sim_prm$maxlag, self$design$sim_prm$ageH)]
        tbl[, Year := NULL]
        setnames(tbl, tolower(names(tbl)))
        tbl[, sex := factor(sex, 0:1, c("men", "women"))]

        col_nam <-
            setdiff(names(tbl), intersect(names(synthpop$pop), names(tbl)))
        absorb_dt(synthpop$pop, tbl)
        if (anyNA(synthpop$pop[, ..col_nam])) stop("NA in the exposure distribution")
        synthpop$pop[year > setyear, Fruit_vege_curr_xps := qZINBI(rank_Fruit_vege, mu, sigma, nu)]
        synthpop$pop[, c(col_nam) := NULL]
        NULL
    }
)

IMPACTncd$
    run(1:mcnum, multicore = TRUE, "FV")

# Smoking ----
IMPACTncd$update_primary_prevention_scn(
    function(synthpop) {
        setyear <- 2001L
        rs <- .Random.seed
        on.exit(set.seed(rs))
        set.seed(2121870L + synthpop$mc) # Not mc_aggr # same as in Synthpop class

        synthpop$pop[, tax_tabaco := 0L] # true for 2001
        # synthpop$pop[, tax_tabaco := fcase(
        #     year < 2006L,                 0L,
        #     year >= 2006L & year < 2010L, 1L,
        #     year >= 2010L & year < 2018L, 2L,
        #     year >= 2018L,                3L
        # )]
        synthpop$pop[, tax_tabaco := factor(tax_tabaco, 0:3, 0:3)]
        synthpop$pop[, Smoking_curr_xps := as.integer(Smoking_curr_xps) - 1L]

        tbl <-
            read_fst("./inputs/exposure_distributions/Table_Smoking_NevEx_vs_current.fst",
                as.data.table = TRUE
            )[Year == setyear & between(Age, self$design$sim_prm$ageL - self$design$sim_prm$maxlag, self$design$sim_prm$ageH)]
        tbl[, Year := NULL]
        setnames(tbl, tolower(names(tbl)))
        tbl[, sex := factor(sex, 0:1, c("men", "women"))]

        col_nam <-
            setdiff(names(tbl), intersect(names(synthpop$pop), names(tbl)))
        absorb_dt(synthpop$pop, tbl)
        if (anyNA(synthpop$pop[, ..col_nam])) stop("NA in the exposure distribution")
        synthpop$pop[year > setyear, Smoking_curr_xps := as.integer(rankstat_Smoking_act < mu) * 2L] # 0 = never smoker or ex, 2 = current
        synthpop$pop[, c(col_nam) := NULL]


        tbl <-
            read_fst("./inputs/exposure_distributions/Table_Smoking_never_vs_ex.fst",
                as.data.table = TRUE
            )[Year == setyear & between(Age, self$design$sim_prm$ageL - self$design$sim_prm$maxlag, self$design$sim_prm$ageH)]
        tbl[, Year := NULL]
        setnames(tbl, tolower(names(tbl)))
        tbl[, sex := factor(sex, 0:1, c("men", "women")), ]
        col_nam <-
            setdiff(names(tbl), intersect(names(synthpop$pop), names(tbl)))
        absorb_dt(synthpop$pop, tbl)
        if (anyNA(synthpop$pop[, ..col_nam])) stop("NA in the exposure distribution")
        synthpop$pop[year > setyear & Smoking_curr_xps == 0L, Smoking_curr_xps := as.integer(rankstat_Smoking_ex < mu)] # 0 = never smoker, 1=ex, 2=current

        synthpop$pop[, Smoking_curr_xps := factor(Smoking_curr_xps + 1L)]
        if (anyNA(synthpop$pop$Smoking_curr_xps)) stop("NA in Smoking_curr_xps")

        synthpop$pop[, c(col_nam, "tax_tabaco") := NULL]


        tbl <-
            read_fst("./inputs/exposure_distributions/Table_Smoking_number.fst",
                as.data.table = TRUE
            )[Year == setyear & between(Age, self$design$sim_prm$ageL - self$design$sim_prm$maxlag, self$design$sim_prm$ageH)]
        tbl[, Year := NULL]
        setnames(tbl, tolower(names(tbl)))
        tbl[, sex := factor(sex, 0:1, c("men", "women")), ]

        col_nam <-
            setdiff(names(tbl), intersect(names(synthpop$pop), names(tbl)))
        absorb_dt(synthpop$pop, tbl)
        if (anyNA(synthpop$pop[, ..col_nam])) stop("NA in the exposure distribution")
        synthpop$pop[year > setyear, Smoking_number_curr_xps := NA_integer_] # Clean slate
        synthpop$pop[
            year > setyear & Smoking_curr_xps == "3",
            Smoking_number_grp := (rankstat_Smoking_number > pa0) +
                (rankstat_Smoking_number > pa1) +
                (rankstat_Smoking_number > pa2) +
                (rankstat_Smoking_number > pa3) +
                (rankstat_Smoking_number > pa4) +
                (rankstat_Smoking_number > pa5) +
                (rankstat_Smoking_number > pa6) +
                (rankstat_Smoking_number > pa7)
        ]

        synthpop$pop[Smoking_number_grp == 0L, Smoking_number_curr_xps := 5L]
        synthpop$pop[Smoking_number_grp == 1L, Smoking_number_curr_xps := 10L]
        synthpop$pop[Smoking_number_grp == 2L, Smoking_number_curr_xps := 15L]
        synthpop$pop[Smoking_number_grp == 3L, Smoking_number_curr_xps := 20L]
        synthpop$pop[Smoking_number_grp == 4L, Smoking_number_curr_xps := 25L]
        synthpop$pop[Smoking_number_grp == 5L, Smoking_number_curr_xps := 30L]
        synthpop$pop[Smoking_number_grp == 6L, Smoking_number_curr_xps := 35L]
        synthpop$pop[Smoking_number_grp == 7L, Smoking_number_curr_xps := 40L]
        # I do not explicitly set.seed because I do so at the beginning of the scenario
        synthpop$pop[Smoking_number_grp == 8L, Smoking_number_curr_xps := sample(c(50L, 60L, 80L), .N, TRUE, prob = c(0.4, 0.45, 0.15))]

        synthpop$pop[, c(col_nam, "Smoking_number_grp") := NULL]
        NULL
    }
)

IMPACTncd$
    run(1:mcnum, multicore = TRUE, "smoking")


# Physical activity ----
IMPACTncd$update_primary_prevention_scn(
    function(synthpop) {
        setyear <- 2001L
        tbl <-
            read_fst("./inputs/exposure_distributions/Table_PA_days.fst",
                as.data.table = TRUE
            )[between(Age, self$design$sim_prm$ageL - self$design$sim_prm$maxlag, self$design$sim_prm$ageH)]
        setnames(tbl, tolower(names(tbl)))
        tbl[, sex := factor(sex, 0:1, c("men", "women"))]

        col_nam <-
            setdiff(names(tbl), intersect(names(synthpop$pop), names(tbl)))

        synthpop$pop[, trueyear := year]
        synthpop$pop[trueyear > setyear, year := setyear]
        synthpop$pop[age >= 70 & year > 2025L, year := 2025L] # for setyear = 2001L, none
        synthpop$pop[age >= 50 & sex == "men" & year < 2010L, year := 2010L]
        synthpop$pop[age >= 70 & sex == "women" & year < 2015L, year := 2015L]

        absorb_dt(synthpop$pop, tbl)
        if (anyNA(synthpop$pop[, ..col_nam])) stop("NA in the exposure distribution")
        synthpop$pop[, `:=`(year = trueyear, trueyear = NULL)]

        synthpop$pop[
            year > setyear,
            PA_days_curr_xps := factor(
                (rank_PA_days > pa0) +
                    (rank_PA_days > pa1) +
                    (rank_PA_days > pa2) +
                    (rank_PA_days > pa3) +
                    (rank_PA_days > pa4) +
                    (rank_PA_days > pa5) +
                    (rank_PA_days > pa6),
                levels = 0:7, labels = 0:7, ordered = TRUE
            )
        ]
        synthpop$pop[, c(col_nam) := NULL]
        NULL
    }
)

IMPACTncd$
    run(1:mcnum, multicore = TRUE, "PA")


# BMI ----
IMPACTncd$update_primary_prevention_scn(
    function(synthpop) {
        setyear <- 2001L
        tbl <-
            read_fst("./inputs/exposure_distributions/Table_BMI.fst",
                as.data.table = TRUE
            )[Year == setyear & between(Age, self$design$sim_prm$ageL - self$design$sim_prm$maxlag, self$design$sim_prm$ageH)]
        tbl[, Year := NULL]
        setnames(tbl, tolower(names(tbl)))
        tbl[, sex := factor(sex, 0:1, c("men", "women"))]

        ### Make PA days category
        synthpop$pop[, PA_3cat := fifelse(
            PA_days_curr_xps %in% as.character(0:1), 1L,
            fifelse(
                PA_days_curr_xps %in% as.character(2:4), 2L,
                fifelse(PA_days_curr_xps %in% as.character(5:7), 3L, NA_integer_)
            )
        )]
        synthpop$pop[, PA_3cat := factor(PA_3cat)]


        col_nam <-
            setdiff(names(tbl), intersect(names(synthpop$pop), names(tbl)))
        absorb_dt(synthpop$pop, tbl)
        if (anyNA(synthpop$pop[, ..col_nam])) stop("NA in the exposure distribution")
        synthpop$pop[year > setyear, BMI_curr_xps := qBCTo(rank_BMI, mu, sigma, nu, tau)]
        synthpop$pop[BMI_curr_xps < 10, BMI_curr_xps := 10] # Truncate BMI predictions to avoid unrealistic values.
        synthpop$pop[BMI_curr_xps > 70, BMI_curr_xps := 70] # Truncate BMI predictions to avoid unrealistic values.

        synthpop$pop[, c(col_nam, "PA_3cat") := NULL]
        NULL
    }
)

IMPACTncd$
    run(1:mcnum, multicore = TRUE, "BMI")

# HbA1c + Med_DM ----
IMPACTncd$update_primary_prevention_scn(
    function(synthpop) {
        setyear <- 2001L

        # first do the med_DM
        tbl <-
            read_fst("./inputs/exposure_distributions/Table_Med_DM.fst",
                as.data.table = TRUE
            )[Year == setyear & between(Age, self$design$sim_prm$ageL - self$design$sim_prm$maxlag, self$design$sim_prm$ageH)]
        tbl[, Year := NULL]
        setnames(tbl, tolower(names(tbl)))
        tbl[, sex := factor(sex, 0:1, c("men", "women")), ]
        col_nam <-
            setdiff(names(tbl), intersect(names(synthpop$pop), names(tbl)))
        absorb_dt(synthpop$pop, tbl)
        if (anyNA(synthpop$pop[, ..col_nam])) stop("NA in the exposure distribution")
        synthpop$pop[year > setyear, Med_DM_curr_xps := qbinom(rankstat_Med_DM, 1L, mu)] # , n_cpu = design_$sim_prm$n_cpu)]
        synthpop$pop[, c(col_nam) := NULL]

        # then do the HbA1c
        tbl <-
            read_fst("./inputs/exposure_distributions/Table_HbA1c.fst",
                as.data.table = TRUE
            )[Year == setyear & between(Age, self$design$sim_prm$ageL - self$design$sim_prm$maxlag, self$design$sim_prm$ageH)]
        tbl[, Year := NULL]
        setnames(tbl, c("Age", "Sex", "BMI", "Med_DM"), c("age", "sex", "BMI_round", "Med_DM_curr_xps"))
        tbl[, sex := factor(sex, 0:1, c("men", "women"))]
        tbl[, BMI_round := as.integer(10 * BMI_round)]
        synthpop$pop[, BMI_round := as.integer(round(10 * BMI_curr_xps, 0))]
        col_nam <-
            setdiff(names(tbl), intersect(names(synthpop$pop), names(tbl)))
        absorb_dt(synthpop$pop, tbl)
        if (anyNA(synthpop$pop[, ..col_nam])) stop("NA in the exposure distribution")
        synthpop$pop[year > setyear, HbA1c_curr_xps := qBCT(rank_HbA1c, mu, sigma, nu, tau)]
        synthpop$pop[, c(col_nam, "BMI_round") := NULL]
        NULL
    }
)

IMPACTncd$
    run(1:mcnum, multicore = TRUE, "HbA1c")

# LDLc + Med_HL ----
IMPACTncd$update_primary_prevention_scn(
    function(synthpop) {
        setyear <- 2001L

        # first do the med_HL
        tbl <-
            read_fst("./inputs/exposure_distributions/Table_Med_HL.fst",
                as.data.table = TRUE
            )[Year == setyear & between(Age, self$design$sim_prm$ageL - self$design$sim_prm$maxlag, self$design$sim_prm$ageH)]
        tbl[, Year := NULL]
        setnames(tbl, tolower(names(tbl)))
        tbl[, sex := factor(sex, 0:1, c("men", "women")), ]

        col_nam <-
            setdiff(names(tbl), intersect(names(synthpop$pop), names(tbl)))
        absorb_dt(synthpop$pop, tbl)
        if (anyNA(synthpop$pop[, ..col_nam])) stop("NA in the exposure distribution")
        synthpop$pop[year > setyear, Med_HL_curr_xps := qbinom(rankstat_Med_HL, 1L, mu)] # , n_cpu = design_$sim_prm$n_cpu)]
        synthpop$pop[, c(col_nam) := NULL]

        # then do the LDLc
        tbl <-
            read_fst("./inputs/exposure_distributions/Table_LDLc.fst",
                as.data.table = TRUE
            )[Year == setyear & between(Age, self$design$sim_prm$ageL - self$design$sim_prm$maxlag, self$design$sim_prm$ageH)]
        tbl[, Year := NULL]
        setnames(tbl, c("Age", "Sex", "BMI", "Med_HL"), c("age", "sex", "BMI_round", "Med_HL_curr_xps"))
        tbl[, sex := factor(sex, 0:1, c("men", "women"))]
        tbl[, BMI_round := as.integer(10 * BMI_round)]
        synthpop$pop[, BMI_round := as.integer(round(10 * BMI_curr_xps, 0))]
        col_nam <-
            setdiff(names(tbl), intersect(names(synthpop$pop), names(tbl)))
        absorb_dt(synthpop$pop, tbl)
        if (anyNA(synthpop$pop[, ..col_nam])) stop("NA in the exposure distribution")
        synthpop$pop[year > setyear, LDLc_curr_xps := qBCT(rank_LDLc, mu, sigma, nu, tau)]
        synthpop$pop[, c(col_nam, "BMI_round") := NULL]
        NULL
    }
)

IMPACTncd$
    run(1:mcnum, multicore = TRUE, "LDLc")

# SBP + Med_HT ----
IMPACTncd$update_primary_prevention_scn(
    function(synthpop) {
        setyear <- 2001L

        # first do the med_HT
        tbl <-
            read_fst("./inputs/exposure_distributions/Table_Med_HT.fst",
                as.data.table = TRUE
            )[Year == setyear & between(Age, self$design$sim_prm$ageL - self$design$sim_prm$maxlag, self$design$sim_prm$ageH)]
        tbl[, Year := NULL]
        setnames(tbl, tolower(names(tbl)))
        tbl[, sex := factor(sex, 0:1, c("men", "women")), ]

        col_nam <-
            setdiff(names(tbl), intersect(names(synthpop$pop), names(tbl)))
        absorb_dt(synthpop$pop, tbl)
        if (anyNA(synthpop$pop[, ..col_nam])) stop("NA in the exposure distribution")
        synthpop$pop[year > setyear, Med_HT_curr_xps := qbinom(rankstat_Med_HT, 1L, mu)] # , n_cpu = design_$sim_prm$n_cpu)]
        synthpop$pop[, c(col_nam) := NULL]

        # then do the SBPs
        tbl <-
            read_fst("./inputs/exposure_distributions/Table_SBP.fst",
                as.data.table = TRUE
            )[Year == setyear & between(Age, self$design$sim_prm$ageL - self$design$sim_prm$maxlag, self$design$sim_prm$ageH)]
        tbl[, Year := NULL]
        setnames(tbl, c("Age", "Sex", "BMI", "Med_HT", "Smoking"), c("age", "sex", "BMI_round", "Med_HT_curr_xps", "smoking_tmp"))
        tbl[, `:=`(
            sex = factor(sex, 0:1, c("men", "women")),
            smoking_tmp = as.integer(smoking_tmp),
            BMI_round = as.integer(BMI_round)
        )]
        synthpop$pop[, `:=`(
            BMI_round = as.integer(round(BMI_curr_xps)), # TODO consider Rfast::Round to speedup
            smoking_tmp = as.integer(Smoking_curr_xps == "3")
        )]
        col_nam <-
            setdiff(names(tbl), intersect(names(synthpop$pop), names(tbl)))
        absorb_dt(synthpop$pop, tbl)
        if (anyNA(synthpop$pop[, ..col_nam])) stop("NA in the exposure distribution")
        synthpop$pop[year > setyear, SBP_curr_xps := qBCPE(rank_SBP, mu, sigma, nu, tau)]
        synthpop$pop[, c(col_nam, "BMI_round", "smoking_tmp") := NULL]
        NULL
    }
)

IMPACTncd$
    run(1:mcnum, multicore = TRUE, "SBP")



# All xps ----
IMPACTncd$update_primary_prevention_scn(
    function(synthpop) {
        setyear <- 2001L

        # F&V
        tbl <-
            read_fst("./inputs/exposure_distributions/Table_Fruit_vege.fst",
                as.data.table = TRUE
            )[Year == setyear & between(Age, self$design$sim_prm$ageL - self$design$sim_prm$maxlag, self$design$sim_prm$ageH)]
        tbl[, Year := NULL]
        setnames(tbl, tolower(names(tbl)))
        tbl[, sex := factor(sex, 0:1, c("men", "women"))]

        col_nam <-
            setdiff(names(tbl), intersect(names(synthpop$pop), names(tbl)))
        absorb_dt(synthpop$pop, tbl)
        if (anyNA(synthpop$pop[, ..col_nam])) stop("NA in the exposure distribution")
        synthpop$pop[year > setyear, Fruit_vege_curr_xps := qZINBI(rank_Fruit_vege, mu, sigma, nu)]
        synthpop$pop[, c(col_nam) := NULL]

        # Smoking
        rs <- .Random.seed
        set.seed(2121870L + synthpop$mc) # Not mc_aggr # same as in Synthpop class
        synthpop$pop[, tax_tabaco := 0L] # true for 2001
        # synthpop$pop[, tax_tabaco := fcase(
        #     year < 2006L,                 0L,
        #     year >= 2006L & year < 2010L, 1L,
        #     year >= 2010L & year < 2018L, 2L,
        #     year >= 2018L,                3L
        # )]
        synthpop$pop[, tax_tabaco := factor(tax_tabaco, 0:3, 0:3)]
        synthpop$pop[, Smoking_curr_xps := as.integer(Smoking_curr_xps) - 1L]

        tbl <-
            read_fst("./inputs/exposure_distributions/Table_Smoking_NevEx_vs_current.fst",
                as.data.table = TRUE
            )[Year == setyear & between(Age, self$design$sim_prm$ageL - self$design$sim_prm$maxlag, self$design$sim_prm$ageH)]
        tbl[, Year := NULL]
        setnames(tbl, tolower(names(tbl)))
        tbl[, sex := factor(sex, 0:1, c("men", "women"))]

        col_nam <-
            setdiff(names(tbl), intersect(names(synthpop$pop), names(tbl)))
        absorb_dt(synthpop$pop, tbl)
        if (anyNA(synthpop$pop[, ..col_nam])) stop("NA in the exposure distribution")
        synthpop$pop[year > setyear, Smoking_curr_xps := as.integer(rankstat_Smoking_act < mu) * 2L] # 0 = never smoker or ex, 2 = current
        synthpop$pop[, c(col_nam) := NULL]


        tbl <-
            read_fst("./inputs/exposure_distributions/Table_Smoking_never_vs_ex.fst",
                as.data.table = TRUE
            )[Year == setyear & between(Age, self$design$sim_prm$ageL - self$design$sim_prm$maxlag, self$design$sim_prm$ageH)]
        tbl[, Year := NULL]
        setnames(tbl, tolower(names(tbl)))
        tbl[, sex := factor(sex, 0:1, c("men", "women")), ]
        col_nam <-
            setdiff(names(tbl), intersect(names(synthpop$pop), names(tbl)))
        absorb_dt(synthpop$pop, tbl)
        if (anyNA(synthpop$pop[, ..col_nam])) stop("NA in the exposure distribution")
        synthpop$pop[year > setyear & Smoking_curr_xps == 0L, Smoking_curr_xps := as.integer(rankstat_Smoking_ex < mu)] # 0 = never smoker, 1=ex, 2=current

        synthpop$pop[, Smoking_curr_xps := factor(Smoking_curr_xps + 1L)]
        if (anyNA(synthpop$pop$Smoking_curr_xps)) stop("NA in Smoking_curr_xps")
        synthpop$pop[, c(col_nam, "tax_tabaco") := NULL]


        tbl <-
            read_fst("./inputs/exposure_distributions/Table_Smoking_number.fst",
                as.data.table = TRUE
            )[Year == setyear & between(Age, self$design$sim_prm$ageL - self$design$sim_prm$maxlag, self$design$sim_prm$ageH)]
        tbl[, Year := NULL]
        setnames(tbl, tolower(names(tbl)))
        tbl[, sex := factor(sex, 0:1, c("men", "women")), ]

        col_nam <-
            setdiff(names(tbl), intersect(names(synthpop$pop), names(tbl)))
        absorb_dt(synthpop$pop, tbl)
        if (anyNA(synthpop$pop[, ..col_nam])) stop("NA in the exposure distribution")
        synthpop$pop[year > setyear, Smoking_number_curr_xps := NA_integer_] # Clean slate
        synthpop$pop[
            year > setyear & Smoking_curr_xps == "3",
            Smoking_number_grp := (rankstat_Smoking_number > pa0) +
                (rankstat_Smoking_number > pa1) +
                (rankstat_Smoking_number > pa2) +
                (rankstat_Smoking_number > pa3) +
                (rankstat_Smoking_number > pa4) +
                (rankstat_Smoking_number > pa5) +
                (rankstat_Smoking_number > pa6) +
                (rankstat_Smoking_number > pa7)
        ]

        synthpop$pop[Smoking_number_grp == 0L, Smoking_number_curr_xps := 5L]
        synthpop$pop[Smoking_number_grp == 1L, Smoking_number_curr_xps := 10L]
        synthpop$pop[Smoking_number_grp == 2L, Smoking_number_curr_xps := 15L]
        synthpop$pop[Smoking_number_grp == 3L, Smoking_number_curr_xps := 20L]
        synthpop$pop[Smoking_number_grp == 4L, Smoking_number_curr_xps := 25L]
        synthpop$pop[Smoking_number_grp == 5L, Smoking_number_curr_xps := 30L]
        synthpop$pop[Smoking_number_grp == 6L, Smoking_number_curr_xps := 35L]
        synthpop$pop[Smoking_number_grp == 7L, Smoking_number_curr_xps := 40L]
        # I do not explicitly set.seed because I do so at the beginning of the scenario
        synthpop$pop[Smoking_number_grp == 8L, Smoking_number_curr_xps := sample(c(50L, 60L, 80L), .N, TRUE, prob = c(0.4, 0.45, 0.15))]

        synthpop$pop[, c(col_nam, "Smoking_number_grp") := NULL]
        set.seed(rs)

        # PA
        tbl <-
            read_fst("./inputs/exposure_distributions/Table_PA_days.fst",
                as.data.table = TRUE
            )[between(Age, self$design$sim_prm$ageL - self$design$sim_prm$maxlag, self$design$sim_prm$ageH)]
        setnames(tbl, tolower(names(tbl)))
        tbl[, sex := factor(sex, 0:1, c("men", "women"))]

        col_nam <-
            setdiff(names(tbl), intersect(names(synthpop$pop), names(tbl)))

        synthpop$pop[, trueyear := year]
        synthpop$pop[trueyear > setyear, year := setyear]
        synthpop$pop[age >= 70 & year > 2025L, year := 2025L] # for setyear = 2001L, none
        synthpop$pop[age >= 50 & sex == "men" & year < 2010L, year := 2010L]
        synthpop$pop[age >= 70 & sex == "women" & year < 2015L, year := 2015L]

        absorb_dt(synthpop$pop, tbl)
        if (anyNA(synthpop$pop[, ..col_nam])) stop("NA in the exposure distribution")
        synthpop$pop[, `:=`(year = trueyear, trueyear = NULL)]

        synthpop$pop[
            year > setyear,
            PA_days_curr_xps := factor(
                (rank_PA_days > pa0) +
                    (rank_PA_days > pa1) +
                    (rank_PA_days > pa2) +
                    (rank_PA_days > pa3) +
                    (rank_PA_days > pa4) +
                    (rank_PA_days > pa5) +
                    (rank_PA_days > pa6),
                levels = 0:7, labels = 0:7, ordered = TRUE
            )
        ]
        synthpop$pop[, c(col_nam) := NULL]

        # BMI
        tbl <-
            read_fst("./inputs/exposure_distributions/Table_BMI.fst",
                as.data.table = TRUE
            )[Year == setyear & between(Age, self$design$sim_prm$ageL - self$design$sim_prm$maxlag, self$design$sim_prm$ageH)]
        tbl[, Year := NULL]
        setnames(tbl, tolower(names(tbl)))
        tbl[, sex := factor(sex, 0:1, c("men", "women"))]

        ### Make PA days category
        synthpop$pop[, PA_3cat := fifelse(
            PA_days_curr_xps %in% as.character(0:1), 1L,
            fifelse(
                PA_days_curr_xps %in% as.character(2:4), 2L,
                fifelse(PA_days_curr_xps %in% as.character(5:7), 3L, NA_integer_)
            )
        )]
        synthpop$pop[, PA_3cat := factor(PA_3cat)]


        col_nam <-
            setdiff(names(tbl), intersect(names(synthpop$pop), names(tbl)))
        absorb_dt(synthpop$pop, tbl)
        if (anyNA(synthpop$pop[, ..col_nam])) stop("NA in the exposure distribution")
        synthpop$pop[year > setyear, BMI_curr_xps := qBCTo(rank_BMI, mu, sigma, nu, tau)]
        synthpop$pop[BMI_curr_xps < 10, BMI_curr_xps := 10] # Truncate BMI predictions to avoid unrealistic values.
        synthpop$pop[BMI_curr_xps > 70, BMI_curr_xps := 70] # Truncate BMI predictions to avoid unrealistic values.
        synthpop$pop[, c(col_nam, "PA_3cat") := NULL]


        # med_DM
        tbl <-
            read_fst("./inputs/exposure_distributions/Table_Med_DM.fst",
                as.data.table = TRUE
            )[Year == setyear & between(Age, self$design$sim_prm$ageL - self$design$sim_prm$maxlag, self$design$sim_prm$ageH)]
        tbl[, Year := NULL]
        setnames(tbl, tolower(names(tbl)))
        tbl[, sex := factor(sex, 0:1, c("men", "women")), ]
        col_nam <-
            setdiff(names(tbl), intersect(names(synthpop$pop), names(tbl)))
        absorb_dt(synthpop$pop, tbl)
        if (anyNA(synthpop$pop[, ..col_nam])) stop("NA in the exposure distribution")
        synthpop$pop[year > setyear, Med_DM_curr_xps := qbinom(rankstat_Med_DM, 1L, mu)] # , n_cpu = design_$sim_prm$n_cpu)]
        synthpop$pop[, c(col_nam) := NULL]

        synthpop$pop[, BMI_round := as.integer(round(10 * BMI_curr_xps, 0))]

        # then do the HbA1c
        tbl <-
            read_fst("./inputs/exposure_distributions/Table_HbA1c.fst",
                as.data.table = TRUE
            )[Year == setyear & between(Age, self$design$sim_prm$ageL - self$design$sim_prm$maxlag, self$design$sim_prm$ageH)]
        tbl[, Year := NULL]
        setnames(tbl, c("Age", "Sex", "BMI", "Med_DM"), c("age", "sex", "BMI_round", "Med_DM_curr_xps"))
        tbl[, sex := factor(sex, 0:1, c("men", "women"))]
        tbl[, BMI_round := as.integer(10 * BMI_round)]
        col_nam <-
            setdiff(names(tbl), intersect(names(synthpop$pop), names(tbl)))
        absorb_dt(synthpop$pop, tbl)
        if (anyNA(synthpop$pop[, ..col_nam])) stop("NA in the exposure distribution")
        synthpop$pop[year > setyear, HbA1c_curr_xps := qBCT(rank_HbA1c, mu, sigma, nu, tau)]
        synthpop$pop[, c(col_nam) := NULL]

        # first do the med_HL
        tbl <-
            read_fst("./inputs/exposure_distributions/Table_Med_HL.fst",
                as.data.table = TRUE
            )[Year == setyear & between(Age, self$design$sim_prm$ageL - self$design$sim_prm$maxlag, self$design$sim_prm$ageH)]
        tbl[, Year := NULL]
        setnames(tbl, tolower(names(tbl)))
        tbl[, sex := factor(sex, 0:1, c("men", "women")), ]

        col_nam <-
            setdiff(names(tbl), intersect(names(synthpop$pop), names(tbl)))
        absorb_dt(synthpop$pop, tbl)
        if (anyNA(synthpop$pop[, ..col_nam])) stop("NA in the exposure distribution")
        synthpop$pop[year > setyear, Med_HL_curr_xps := qbinom(rankstat_Med_HL, 1L, mu)] # , n_cpu = design_$sim_prm$n_cpu)]
        synthpop$pop[, c(col_nam) := NULL]

        # then do the LDLc
        tbl <-
            read_fst("./inputs/exposure_distributions/Table_LDLc.fst",
                as.data.table = TRUE
            )[Year == setyear & between(Age, self$design$sim_prm$ageL - self$design$sim_prm$maxlag, self$design$sim_prm$ageH)]
        tbl[, Year := NULL]
        setnames(tbl, c("Age", "Sex", "BMI", "Med_HL"), c("age", "sex", "BMI_round", "Med_HL_curr_xps"))
        tbl[, sex := factor(sex, 0:1, c("men", "women"))]
        tbl[, BMI_round := as.integer(10 * BMI_round)]
        col_nam <-
            setdiff(names(tbl), intersect(names(synthpop$pop), names(tbl)))
        absorb_dt(synthpop$pop, tbl)
        if (anyNA(synthpop$pop[, ..col_nam])) stop("NA in the exposure distribution")
        synthpop$pop[year > setyear, LDLc_curr_xps := qBCT(rank_LDLc, mu, sigma, nu, tau)]
        synthpop$pop[, c(col_nam) := NULL]


        # first do the med_HT
        tbl <-
            read_fst("./inputs/exposure_distributions/Table_Med_HT.fst",
                as.data.table = TRUE
            )[Year == setyear & between(Age, self$design$sim_prm$ageL - self$design$sim_prm$maxlag, self$design$sim_prm$ageH)]
        tbl[, Year := NULL]
        setnames(tbl, tolower(names(tbl)))
        tbl[, sex := factor(sex, 0:1, c("men", "women")), ]

        col_nam <-
            setdiff(names(tbl), intersect(names(synthpop$pop), names(tbl)))
        absorb_dt(synthpop$pop, tbl)
        if (anyNA(synthpop$pop[, ..col_nam])) stop("NA in the exposure distribution")
        synthpop$pop[year > setyear, Med_HT_curr_xps := qbinom(rankstat_Med_HT, 1L, mu)] # , n_cpu = design_$sim_prm$n_cpu)]
        synthpop$pop[, c(col_nam) := NULL]

        # then do the SBPs
        tbl <-
            read_fst("./inputs/exposure_distributions/Table_SBP.fst",
                as.data.table = TRUE
            )[Year == setyear & between(Age, self$design$sim_prm$ageL - self$design$sim_prm$maxlag, self$design$sim_prm$ageH)]
        tbl[, Year := NULL]
        setnames(tbl, c("Age", "Sex", "BMI", "Med_HT", "Smoking"), c("age", "sex", "BMI_round", "Med_HT_curr_xps", "smoking_tmp"))
        tbl[, `:=`(
            sex = factor(sex, 0:1, c("men", "women")),
            smoking_tmp = as.integer(smoking_tmp),
            BMI_round = as.integer(BMI_round)
        )]
        synthpop$pop[, `:=`(
            BMI_round = as.integer(round(BMI_curr_xps)), # TODO consider Rfast::Round to speedup
            smoking_tmp = as.integer(Smoking_curr_xps == "3")
        )]
        col_nam <-
            setdiff(names(tbl), intersect(names(synthpop$pop), names(tbl)))
        absorb_dt(synthpop$pop, tbl)
        if (anyNA(synthpop$pop[, ..col_nam])) stop("NA in the exposure distribution")
        synthpop$pop[year > setyear, SBP_curr_xps := qBCPE(rank_SBP, mu, sigma, nu, tau)]
        synthpop$pop[, c(col_nam, "BMI_round", "smoking_tmp") := NULL]

        NULL
    }
)

IMPACTncd$
    run(1:mcnum, multicore = TRUE, "All")

IMPACTncd$export_summaries(
    multicore = TRUE,
    type = c(
        "le", "hle", "dis_char", "prvl",
        "incd", "dis_mrtl", "mrtl",
        "allcause_mrtl_by_dis", "cms"
    )
)

source("./auxil/process_out_for_HF.R")


# sp$pop colnames
#  [1] "pid"                     "year"                    "age"                     "sex"                     "type"
#  [6] "rank_Fruit_vege"         "rankstat_Smoking_act"    "rankstat_Smoking_ex"     "rankstat_Med_HT"         "rankstat_Med_HL"
# [11] "rankstat_Med_DM"         "rank_PA_days"            "rank_BMI"                "rank_HbA1c"              "rank_LDLc"
# [16] "rank_SBP"                "rankstat_Smoking_number" "Fruit_vege_curr_xps"     "Smoking_curr_xps"        "Smoking_number_curr_xps"
# [21] "Med_HT_curr_xps"         "Med_HL_curr_xps"         "Med_DM_curr_xps"         "PA_days_curr_xps"        "BMI_curr_xps"
# [26] "HbA1c_curr_xps"          "LDLc_curr_xps"           "SBP_curr_xps"            "pid_mrk"                 "wt_immrtl"
# [31] "all_cause_mrtl"          "cms_score"               "cms_count"
