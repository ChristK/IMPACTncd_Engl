source("global.R")

IMPACTncd <- Simulation$new("./scenarios/Ineq/sim_design_ineq.yaml")

n_runs <- 100L

IMPACTncd$
del_logs()$
del_outputs()$
# del_synthpops()$
run(1:n_runs, multicore = TRUE, "sc0")

# Everyone goes to 5 ----
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    sc_year <- 26L
    sp$pop[
      year >= sc_year,
      `:=`(dimd = "10 least deprived", qimd = "5 least deprived")
    ]

    # first delete existing columns to avoid confusion
    xps_nam <- c(
      "active_days",
      "met",
      "fruit",
      "veg",
      "smok_status",
      "smok_quit_yrs",
      "smok_dur",
      "smok_cig",
      "ets",
      "alcohol",
      "bmi",
      "sbp",
      "bpmed",
      "tchol",
      "statin_px"
    )
    sp$pop[, (c(paste0(xps_nam, "_curr_xps"), "income", "education")) := NULL] # to be regenerated

    if (max(sp$pop$age) > 90L) {
      sp$pop[, age100 := age]
      sp$pop[age > 90L, age := 90L]
    }

    # Generate education
    tbl <- self$design$exposures$education$get_table()
    nam <- intersect(names(sp$pop), names(tbl))
    tt <- tbl[age == min(age)]
    tt <- clone_dt(tt, self$design$sim_prm$sim_horizon_max) # TODO adding design_$sim_prm$sim_horizon_max + design_$sim_prm$maxlag for longer projections
    tt[, ':='(age = age - .id, year = year - .id)] # as the sim progress these will become 30 yo
    tt[, .id := NULL]
    tbl <- rbind(tt, tbl)

    idx <- sp$pop[, .I[1], by = pid]$V1 # index of first row for each individual (where education is generated)
    tt <- sp$pop[
      idx,
      .(pid, year, age, ethnicity, qimd, sex, sha, rankstat_education)
    ]
    tt[
      tbl,
      education := (rankstat_education > ed1) +
        (rankstat_education > ed2) +
        (rankstat_education > ed3) +
        (rankstat_education > ed4) +
        (rankstat_education > ed5) +
        (rankstat_education > ed6) +
        1L,
      on = nam
    ]
    tt[,
      education := factor(
        education,
        levels = 1:7,
        labels = c(
          "NVQ4/NVQ5/Degree or equiv",
          "Higher ed below degree",
          "NVQ3/GCE A Level equiv",
          "NVQ2/GCE O Level equiv",
          "NVQ1/CSE other grade equiv",
          "Foreign/other",
          "No qualification"
        )
      )
    ]
    sp$pop[tt, on = "pid", education := i.education]

    sp$pop[, rankstat_education := NULL]
    rm(tt, idx, tbl, nam)

    # Generate income ----
    self$design$exposures$income$generate(sp$pop, self$design)

    # Generate active days ----
    self$design$exposures$active_days$generate(sp$pop, self$design)

    sp$pop[,
      met := as.integer(floor(
        active_days *
          (3L + qbinom(rankstat_pa_met, 8, 3 / 11)) *
          (30 + qexp(rankstat_pa_dur, 1 / 7)) /
          100
      ))
    ]

    sp$pop[, c("rankstat_pa_met", "rankstat_pa_dur") := NULL]

    # Generate fruit consumption (ZISICHEL) ----
    self$design$exposures$fruit$generate(sp$pop, self$design)

    # Generate veg consumption (DEL) ----
    self$design$exposures$veg$generate(sp$pop, self$design)

    # Smoking simulation ----

    # Assign smok_status when pid_mrk == true (the first year an individual enters the simulation (with lags))
    self$design$exposures$smok_status$generate(sp$pop, self$design)
    sp$pop[, smok_status_ref := smok_status]

    # Initialize smoking duration variables
    set(sp$pop, NULL, "smok_quit_yrs", 0L)
    set(sp$pop, NULL, "smok_dur", 0L)

    # Assign smok_quit_yrs and smok_dur for ex-smokers (pid_mrk & smok_status %in% 2:3)
    # Only generate values for rows that need them for efficiency
    idx <- sp$pop[, which((pid_mrk) & smok_status %in% 2:3)]
    self$design$exposures$smok_quit_yrs$generate(sp$pop, self$design, idx)
    self$design$exposures$smok_dur_ex$generate(sp$pop, self$design, idx)

    # Assign smok_dur for current smokers (pid_mrk & smok_status == 4)
    idx <- sp$pop[, which((pid_mrk) & smok_status == 4L)]
    self$design$exposures$smok_dur_curr$generate(sp$pop, self$design, idx)

    # Clean up rank variables if not keeping
    if (!self$design$sim_prm$keep_simulants_rn) {
      for (rn in c(
        "rankstat_smok_quit_yrs",
        "rankstat_smok_dur_ex",
        "rankstat_smok_dur_curr"
      )) {
        if (rn %in% names(sp$pop)) sp$pop[, (rn) := NULL]
      }
    }

    # Ensure smoking histories start from age 12
    sp$pop[age - smok_quit_yrs < 12L, smok_quit_yrs := age - 12L]
    sp$pop[age - smok_dur < 12L, smok_dur := age - 12L]
    sp$pop[
      age - smok_dur - smok_quit_yrs < 12L,
      `:=`(
        smok_dur = as.integer(
          smok_dur / ((smok_dur + smok_quit_yrs) / (age - 12L))
        ),
        smok_quit_yrs = as.integer(
          smok_quit_yrs / ((smok_dur + smok_quit_yrs) / (age - 12L))
        )
      )
    ]

    # Assign smok_incid probabilities
    self$design$exposures$smok_incid$generate(sp$pop, self$design)

    # Assign smok_cessation probabilities
    self$design$exposures$smok_cess$generate(sp$pop, self$design)

    # Handle smok_relapse probabilities
    tbl <-
      read_parquet_dt("./inputs/exposure_distributions/smok_relapse")
    tbl <-
      dcast(tbl, sex + qimd ~ smok_quit_yrs, value.var = "pr")
    nam <- tbl[, paste0(sex, " ", qimd)]
    tbl <-
      as.matrix(tbl[, mget(paste0(1:15))], rownames = nam)

    simsmok(sp$pop, tbl, self$design$sim_prm$smoking_relapse_limit)
    # sp$pop[!(pid_mrk), table(smok_status)]
    # sp$pop[pid == 1, plot(year, smok_status, ylim = c(0, 4))]
    # sp$pop[pid == 10, .(age, smok_status, smok_quit_yrs, smok_dur)]
    # sp$pop[, sum(smok_status == 4)/.N, keyby = year]

    if (self$design$sim_prm$simsmok_calibration) {
      # calculate dif between ref (multinom) and simsmok
      # I will further calibrate to better match HSE
      # Deterministic cell-level random value for calibration adjustments
      cell_rn <- sp$pop[,
        .(rn_clb = rank_smok_clb[1L]),
        keyby = .(year, age, sex, qimd)
      ]
      obs <-
        sp$pop[smok_status == 1L, .(nsa = .N), keyby = .(year, age, sex, qimd)] # ? add sha/ethn
      ref <-
        sp$pop[
          smok_status_ref == 1L,
          .(nsr = .N),
          keyby = .(year, age, sex, qimd)
        ] # ? add sha/ethn
      absorb_dt(ref, obs)
      setnafill(ref, "c", 0L, cols = "nsa")
      ref[, `:=`(dif = nsr - nsa, nsr = NULL, nsa = NULL)]
      ref[cell_rn, on = .(year, age, sex, qimd), rn_clb := i.rn_clb]
      # Further calibrate to match better with HSE (deterministic via rank_smok_clb)
      ref[
        sex == "men" &
          rn_clb < 0.9,
        dif := dif - 2L
      ]
      ref[
        sex == "men" &
          rn_clb < (1 * clamp((year - min(year)) / 10, 0, 1)),
        dif := dif + 2L
      ]
      ref[
        sex == "women" &
          rn_clb < 0.2,
        dif := dif - 1L
      ]
      ref[
        age < 49 &
          rn_clb < 0.5,
        dif := dif + 1L
      ]
      ref[
        age < 49 &
          sex == "men" &
          qimd == "3" &
          rn_clb < 0.5,
        dif := dif - 1L
      ]
      ref[
        age < 49 &
          sex == "women" &
          qimd %in% c("4", "5 least deprived") &
          rn_clb < 0.4,
        dif := dif - 1L
      ]
      ref[
        between(age, 50, 69) &
          rn_clb < 0.2,
        dif := dif + 1L
      ]
      ref[
        between(age, 50, 69) &
          sex == "women" &
          qimd == "5 least deprived" &
          rn_clb < 0.4,
        dif := dif + 1L
      ]
      ref[
        between(age, 70, 89) &
          rn_clb < 0.4,
        dif := dif + 1L
      ]
      ref[
        between(age, 70, 89) &
          sex == "men" &
          qimd %in% c("1 most deprived", "2") &
          rn_clb < 0.4,
        dif := dif - 1L
      ]
      ref[
        between(age, 70, 89) &
          sex == "women" &
          qimd %in% c("4d", "2") &
          rn_clb < 0.4,
        dif := dif + 1L
      ]
      ref[, rn_clb := NULL]
      absorb_dt(sp$pop, ref)

      # when not enough never smokers convert those ex smokers with the longer quit years
      tt <-
        sp$pop[
          smok_status %in% 2:3,
          .(year, age, sex, qimd, pid, smok_quit_yrs, dif)
        ]
      setnafill(tt, "c", 0, cols = "dif")
      tt[dif < 0, dif := 0L]
      setkey(tt, year, age, sex, qimd, smok_quit_yrs)
      pid_to_conv <-
        tt[
          dif > 0,
          .(pid = tail(pid, max(dif))),
          keyby = .(year, age, sex, qimd)
        ]
      sp$pop[
        pid_to_conv,
        on = .(year, pid),
        `:=`(
          smok_status = 1L,
          smok_quit_yrs = 0L,
          smok_dur = 0L,
          smok_cig = 0L
        )
      ]

      # when too many never smokers convert to smok status 2 (occasional)
      tt <-
        sp$pop[
          smok_status %in% 1,
          .(year, age, sex, qimd, pid, rank_smok_clb, dif)
        ]
      setnafill(tt, "c", 0, cols = "dif")
      tt[dif > 0, dif := 0L]
      tt[, dif := -dif]
      setkey(tt, year, age, sex, qimd, rank_smok_clb)
      # Ensure there are enough people to sample from
      ttt <-
        tt[,
          .(lpid = length(pid), mdif = max(dif)),
          by = .(year, age, sex, qimd)
        ][mdif > lpid, ]
      tt[ttt, on = .NATURAL, dif := i.lpid]
      pid_to_conv <-
        tt[
          dif > 0,
          .(pid = head(pid, max(dif))),
          keyby = .(year, age, sex, qimd)
        ]
      sp$pop[
        pid_to_conv,
        on = .(year, pid),
        `:=`(
          smok_status = 2L,
          smok_quit_yrs = sp$pop[smok_status == 2L, median(smok_quit_yrs)],
          smok_dur = sp$pop[smok_status == 2L, median(smok_dur)],
          smok_cig = 1L
        )
      ]
      sp$pop[, dif := NULL]

      # Same logic for active smokers
      # I will further calibrate to better match HSE
      # calculate dif between ref (multinom) and simsmok
      obs <-
        sp$pop[smok_status == 4L, .(nsa = .N), keyby = .(year, age, sex, qimd)] # ? add sha/ethn
      ref <-
        sp$pop[
          smok_status_ref == 4L,
          .(nsr = .N),
          keyby = .(year, age, sex, qimd)
        ] # ? add sha/ethn
      absorb_dt(ref, obs)
      setnafill(ref, "c", 0, cols = "nsa")
      ref[, `:=`(dif = nsr - nsa, nsr = NULL, nsa = NULL)]
      ref[cell_rn, on = .(year, age, sex, qimd), rn_clb := i.rn_clb]
      # Further calibrate to match better with HSE (deterministic via rank_smok_clb)
      ref[
        sex == "men" &
          rn_clb < 0.3,
        dif := dif - 1L
      ]
      ref[
        sex == "women" &
          qimd != "1 most deprived" &
          qimd != "5 least deprived" &
          rn_clb < 0.8,
        dif := dif - 1L
      ]
      ref[
        qimd == "1 most deprived" &
          rn_clb < 0.5,
        dif := dif + 1L
      ]
      ref[
        qimd == "5 least deprived" &
          rn_clb < 0.6,
        dif := dif - 1L
      ] # - reduces
      ref[, rn_clb := NULL]
      absorb_dt(sp$pop, ref)

      # when not enough active smokers convert those ex smokers with the shortest quit years
      tt <-
        sp$pop[
          smok_status == 3L,
          .(year, age, sex, qimd, pid, smok_quit_yrs, dif)
        ]
      setnafill(tt, "c", 0, cols = "dif")
      tt[dif < 0, dif := 0L]
      setkey(tt, year, age, sex, qimd, smok_quit_yrs)
      pid_to_conv <-
        tt[
          dif > 0,
          .(pid = head(pid, max(dif))),
          keyby = .(year, age, sex, qimd)
        ]
      sp$pop[
        pid_to_conv,
        on = .(year, pid),
        `:=`(
          smok_status = 4L,
          smok_quit_yrs = 0L,
          smok_dur = smok_dur + 1L
        )
      ] # TODO fix smoking duration

      # when too many never smokers convert to smok status 3
      tt <-
        sp$pop[
          smok_status == 4L,
          .(year, age, sex, qimd, pid, rank_smok_clb, dif)
        ]
      setnafill(tt, "c", 0, cols = "dif")
      tt[dif > 0, dif := 0L]
      tt[, dif := -dif]
      setkey(tt, year, age, sex, qimd, rank_smok_clb)
      # Ensure there are enough people to sample from
      ttt <-
        tt[,
          .(lpid = length(pid), mdif = max(dif)),
          by = .(year, age, sex, qimd)
        ][mdif > lpid, ]
      tt[ttt, on = .NATURAL, dif := i.lpid]
      pid_to_conv <-
        tt[
          dif > 0,
          .(pid = head(pid, max(dif))),
          keyby = .(year, age, sex, qimd)
        ]
      sp$pop[
        pid_to_conv,
        on = .(year, pid),
        `:=`(smok_status = 3L, smok_quit_yrs = 1L)
      ]
      sp$pop[, dif := NULL]
      sp$pop[, rank_smok_clb := NULL]
      rm(tt, ttt, obs, ref, pid_to_conv, cell_rn)
    }

    # Assign smok_cig
    set(sp$pop, NULL, "smok_cig", 0L)

    # smok_cig_curr for current smokers (smok_status == 4)
    idx <- sp$pop[, which(smok_status == 4L)]

    self$design$exposures$smok_cig_curr$generate(sp$pop, self$design, idx)

    # smok_cig_ex for ex-smokers at initialization (pid_mrk & smok_status == 3)

    idx <- sp$pop[, which((pid_mrk) & smok_status == 3L)]

    self$design$exposures$smok_cig_ex$generate(sp$pop, self$design, idx)

    simsmok_cig(sp$pop) # carry forward smok_cig if smok_status == 3

    sp$pop[smok_cig == 0L & smok_status != 1L, smok_cig := 1L]

    if (self$design$sim_prm$simsmok_calibration) {
      simsmok_postcalibration(sp$pop)
    } # need to be post cig simulation

    sp$pop[, smok_status := factor(smok_status)]

    colnam <- c("prb_smok_incid", "prb_smok_cess", "smok_status_ref")
    if (!self$design$sim_prm$keep_simulants_rn) {
      colnam <- c(
        colnam,
        "rankstat_simsmok"
      )
    }
    sp$pop[, (colnam) := NULL]


    # Generate ETS (BI) ----

    # Note at the moment this is independent of smoking prevalence TODO
    # calculate how many each smoker pollutes by year, SHA (not qimd) to
    # be used in scenarios. Ideally correct for mortality
    self$design$exposures$ets$generate(sp$pop, self$design)
    # NOTE for line above invert = true in sim_design yaml so 1 - mu to be
    # equivalent to qbinom(rank_ets, 1, mu). Otherwise rank_ets < mu is
    # equivalent to qbinom(rank_ets, 1, mu, lower.tail = FALSE)). The two
    # have the same prevalence, but correlations are captured correctly only
    # with the 1-mu variant.

    # View(sp$pop[, prop_if(ets == 1)/prop_if(smok_status == "4"), keyby = .(year, sha)])

    # Generate alcohol (ZINBI) ----
    self$design$exposures$alcohol$generate(sp$pop, self$design)

    # Generate BMI (BCPEo) ----
    self$design$exposures$bmi$generate(sp$pop, self$design)

    # Generate SBP (BCPEo) ----
    self$design$exposures$sbp$generate(sp$pop, self$design)

    # Generate BP medication (BI) -----
    # Temporarily round sbp for join, then restore
    sp$pop[, `:=`(
      sbp_acc = sbp,
      sbp = as.integer(round(clamp(sbp, 110, 200), -1))
    )]
    self$design$exposures$bp_med$generate(sp$pop, self$design)
    sp$pop[, `:=`(sbp = sbp_acc, sbp_acc = NULL)]

    # TODO calculate probability of dgn HTN

    # Generate tchol (BCT) ----
    self$design$exposures$tchol$generate(sp$pop, self$design)

    # Generate HDL (to tchol ratio) (GB1) ----
    self$design$exposures$hdl_to_tchol$generate(sp$pop, self$design)


    # NOTE this very highly correlated with hdl level (~0.76) and
    #  highly to tchol (~-0.47). The latter is captured by the correlated RNs

    # Generate statins medication (BI) -----
    # Temporarily round tchol for join, then restore
    sp$pop[, `:=`(tchol_acc = tchol, tchol = round(clamp(tchol, 2, 12), 0))]
    self$design$exposures$statin_px$generate(sp$pop, self$design)
    sp$pop[, `:=`(tchol = tchol_acc, tchol_acc = NULL)]

    sp$pop[,
      statin_adherence := qBE(
        rankstat_statin_adherence,
        self$design$sim_prm$statin_adherence,
        0.2
      )
    ]
    sp$pop[,
      bpmed_adherence := qBE(
        rankstat_bpmed_adherence,
        self$design$sim_prm$bpmed_adherence,
        0.2
      )
    ]
    if (!self$design$sim_prm$keep_simulants_rn) {
      sp$pop[,
        c("rankstat_bpmed_adherence", "rankstat_statin_adherence") := NULL
      ]
    }

    if ("age100" %in% names(sp$pop)) {
      sp$pop[, age := NULL]
      setnames(sp$pop, "age100", "age")
    }

    setnames(sp$pop, xps_nam, c(paste0(xps_nam, "_curr_xps")))
  }
)

IMPACTncd$run(1:n_runs, multicore = TRUE, "future_mediated_scn")


IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    sc_year <- 26L
    sp$pop[
      year >= sc_year,
      `:=`(dimd = "10 least deprived", qimd = "5 least deprived")
    ]

        lapply(self$design$diseases, function(x) {
          x$gen_parf(sp, self$design, self$design$diseases)
        })

    # first delete existing columns to avoid confusion
    xps_nam <- c(
      "active_days",
      "met",
      "fruit",
      "veg",
      "smok_status",
      "smok_quit_yrs",
      "smok_dur",
      "smok_cig",
      "ets",
      "alcohol",
      "bmi",
      "sbp",
      "bpmed",
      "tchol",
      "statin_px"
    )
    sp$pop[, (c(paste0(xps_nam, "_curr_xps"), "income", "education")) := NULL] # to be regenerated

    if (max(sp$pop$age) > 90L) {
      sp$pop[, age100 := age]
      sp$pop[age > 90L, age := 90L]
    }

    # Generate education
    tbl <- self$design$exposures$education$get_table()
    nam <- intersect(names(sp$pop), names(tbl))
    tt <- tbl[age == min(age)]
    tt <- clone_dt(tt, self$design$sim_prm$sim_horizon_max) # TODO adding design_$sim_prm$sim_horizon_max + design_$sim_prm$maxlag for longer projections
    tt[, ':='(age = age - .id, year = year - .id)] # as the sim progress these will become 30 yo
    tt[, .id := NULL]
    tbl <- rbind(tt, tbl)

    idx <- sp$pop[, .I[1], by = pid]$V1 # index of first row for each individual (where education is generated)
    tt <- sp$pop[
      idx,
      .(pid, year, age, ethnicity, qimd, sex, sha, rankstat_education)
    ]
    tt[
      tbl,
      education := (rankstat_education > ed1) +
        (rankstat_education > ed2) +
        (rankstat_education > ed3) +
        (rankstat_education > ed4) +
        (rankstat_education > ed5) +
        (rankstat_education > ed6) +
        1L,
      on = nam
    ]
    tt[,
      education := factor(
        education,
        levels = 1:7,
        labels = c(
          "NVQ4/NVQ5/Degree or equiv",
          "Higher ed below degree",
          "NVQ3/GCE A Level equiv",
          "NVQ2/GCE O Level equiv",
          "NVQ1/CSE other grade equiv",
          "Foreign/other",
          "No qualification"
        )
      )
    ]
    sp$pop[tt, on = "pid", education := i.education]

    sp$pop[, rankstat_education := NULL]
    rm(tt, idx, tbl, nam)

    # Generate income ----
    self$design$exposures$income$generate(sp$pop, self$design)

    # Generate active days ----
    self$design$exposures$active_days$generate(sp$pop, self$design)

    sp$pop[,
      met := as.integer(floor(
        active_days *
          (3L + qbinom(rankstat_pa_met, 8, 3 / 11)) *
          (30 + qexp(rankstat_pa_dur, 1 / 7)) /
          100
      ))
    ]

    sp$pop[, c("rankstat_pa_met", "rankstat_pa_dur") := NULL]

    # Generate fruit consumption (ZISICHEL) ----
    self$design$exposures$fruit$generate(sp$pop, self$design)

    # Generate veg consumption (DEL) ----
    self$design$exposures$veg$generate(sp$pop, self$design)

    # Smoking simulation ----

    # Assign smok_status when pid_mrk == true (the first year an individual enters the simulation (with lags))
    self$design$exposures$smok_status$generate(sp$pop, self$design)
    sp$pop[, smok_status_ref := smok_status]

    # Initialize smoking duration variables
    set(sp$pop, NULL, "smok_quit_yrs", 0L)
    set(sp$pop, NULL, "smok_dur", 0L)

    # Assign smok_quit_yrs and smok_dur for ex-smokers (pid_mrk & smok_status %in% 2:3)
    # Only generate values for rows that need them for efficiency
    idx <- sp$pop[, which((pid_mrk) & smok_status %in% 2:3)]
    self$design$exposures$smok_quit_yrs$generate(sp$pop, self$design, idx)
    self$design$exposures$smok_dur_ex$generate(sp$pop, self$design, idx)

    # Assign smok_dur for current smokers (pid_mrk & smok_status == 4)
    idx <- sp$pop[, which((pid_mrk) & smok_status == 4L)]
    self$design$exposures$smok_dur_curr$generate(sp$pop, self$design, idx)

    # Clean up rank variables if not keeping
    if (!self$design$sim_prm$keep_simulants_rn) {
      for (rn in c(
        "rankstat_smok_quit_yrs",
        "rankstat_smok_dur_ex",
        "rankstat_smok_dur_curr"
      )) {
        if (rn %in% names(sp$pop)) sp$pop[, (rn) := NULL]
      }
    }

    # Ensure smoking histories start from age 12
    sp$pop[age - smok_quit_yrs < 12L, smok_quit_yrs := age - 12L]
    sp$pop[age - smok_dur < 12L, smok_dur := age - 12L]
    sp$pop[
      age - smok_dur - smok_quit_yrs < 12L,
      `:=`(
        smok_dur = as.integer(
          smok_dur / ((smok_dur + smok_quit_yrs) / (age - 12L))
        ),
        smok_quit_yrs = as.integer(
          smok_quit_yrs / ((smok_dur + smok_quit_yrs) / (age - 12L))
        )
      )
    ]

    # Assign smok_incid probabilities
    self$design$exposures$smok_incid$generate(sp$pop, self$design)

    # Assign smok_cessation probabilities
    self$design$exposures$smok_cess$generate(sp$pop, self$design)

    # Handle smok_relapse probabilities
    tbl <-
      read_parquet_dt("./inputs/exposure_distributions/smok_relapse")
    tbl <-
      dcast(tbl, sex + qimd ~ smok_quit_yrs, value.var = "pr")
    nam <- tbl[, paste0(sex, " ", qimd)]
    tbl <-
      as.matrix(tbl[, mget(paste0(1:15))], rownames = nam)

    simsmok(sp$pop, tbl, self$design$sim_prm$smoking_relapse_limit)
    # sp$pop[!(pid_mrk), table(smok_status)]
    # sp$pop[pid == 1, plot(year, smok_status, ylim = c(0, 4))]
    # sp$pop[pid == 10, .(age, smok_status, smok_quit_yrs, smok_dur)]
    # sp$pop[, sum(smok_status == 4)/.N, keyby = year]

    if (self$design$sim_prm$simsmok_calibration) {
      # calculate dif between ref (multinom) and simsmok
      # I will further calibrate to better match HSE
      # Deterministic cell-level random value for calibration adjustments
      cell_rn <- sp$pop[,
        .(rn_clb = rank_smok_clb[1L]),
        keyby = .(year, age, sex, qimd)
      ]
      obs <-
        sp$pop[smok_status == 1L, .(nsa = .N), keyby = .(year, age, sex, qimd)] # ? add sha/ethn
      ref <-
        sp$pop[
          smok_status_ref == 1L,
          .(nsr = .N),
          keyby = .(year, age, sex, qimd)
        ] # ? add sha/ethn
      absorb_dt(ref, obs)
      setnafill(ref, "c", 0L, cols = "nsa")
      ref[, `:=`(dif = nsr - nsa, nsr = NULL, nsa = NULL)]
      ref[cell_rn, on = .(year, age, sex, qimd), rn_clb := i.rn_clb]
      # Further calibrate to match better with HSE (deterministic via rank_smok_clb)
      ref[
        sex == "men" &
          rn_clb < 0.9,
        dif := dif - 2L
      ]
      ref[
        sex == "men" &
          rn_clb < (1 * clamp((year - min(year)) / 10, 0, 1)),
        dif := dif + 2L
      ]
      ref[
        sex == "women" &
          rn_clb < 0.2,
        dif := dif - 1L
      ]
      ref[
        age < 49 &
          rn_clb < 0.5,
        dif := dif + 1L
      ]
      ref[
        age < 49 &
          sex == "men" &
          qimd == "3" &
          rn_clb < 0.5,
        dif := dif - 1L
      ]
      ref[
        age < 49 &
          sex == "women" &
          qimd %in% c("4", "5 least deprived") &
          rn_clb < 0.4,
        dif := dif - 1L
      ]
      ref[
        between(age, 50, 69) &
          rn_clb < 0.2,
        dif := dif + 1L
      ]
      ref[
        between(age, 50, 69) &
          sex == "women" &
          qimd == "5 least deprived" &
          rn_clb < 0.4,
        dif := dif + 1L
      ]
      ref[
        between(age, 70, 89) &
          rn_clb < 0.4,
        dif := dif + 1L
      ]
      ref[
        between(age, 70, 89) &
          sex == "men" &
          qimd %in% c("1 most deprived", "2") &
          rn_clb < 0.4,
        dif := dif - 1L
      ]
      ref[
        between(age, 70, 89) &
          sex == "women" &
          qimd %in% c("4d", "2") &
          rn_clb < 0.4,
        dif := dif + 1L
      ]
      ref[, rn_clb := NULL]
      absorb_dt(sp$pop, ref)

      # when not enough never smokers convert those ex smokers with the longer quit years
      tt <-
        sp$pop[
          smok_status %in% 2:3,
          .(year, age, sex, qimd, pid, smok_quit_yrs, dif)
        ]
      setnafill(tt, "c", 0, cols = "dif")
      tt[dif < 0, dif := 0L]
      setkey(tt, year, age, sex, qimd, smok_quit_yrs)
      pid_to_conv <-
        tt[
          dif > 0,
          .(pid = tail(pid, max(dif))),
          keyby = .(year, age, sex, qimd)
        ]
      sp$pop[
        pid_to_conv,
        on = .(year, pid),
        `:=`(
          smok_status = 1L,
          smok_quit_yrs = 0L,
          smok_dur = 0L,
          smok_cig = 0L
        )
      ]

      # when too many never smokers convert to smok status 2 (occasional)
      tt <-
        sp$pop[
          smok_status %in% 1,
          .(year, age, sex, qimd, pid, rank_smok_clb, dif)
        ]
      setnafill(tt, "c", 0, cols = "dif")
      tt[dif > 0, dif := 0L]
      tt[, dif := -dif]
      setkey(tt, year, age, sex, qimd, rank_smok_clb)
      # Ensure there are enough people to sample from
      ttt <-
        tt[,
          .(lpid = length(pid), mdif = max(dif)),
          by = .(year, age, sex, qimd)
        ][mdif > lpid, ]
      tt[ttt, on = .NATURAL, dif := i.lpid]
      pid_to_conv <-
        tt[
          dif > 0,
          .(pid = head(pid, max(dif))),
          keyby = .(year, age, sex, qimd)
        ]
      sp$pop[
        pid_to_conv,
        on = .(year, pid),
        `:=`(
          smok_status = 2L,
          smok_quit_yrs = sp$pop[smok_status == 2L, median(smok_quit_yrs)],
          smok_dur = sp$pop[smok_status == 2L, median(smok_dur)],
          smok_cig = 1L
        )
      ]
      sp$pop[, dif := NULL]

      # Same logic for active smokers
      # I will further calibrate to better match HSE
      # calculate dif between ref (multinom) and simsmok
      obs <-
        sp$pop[smok_status == 4L, .(nsa = .N), keyby = .(year, age, sex, qimd)] # ? add sha/ethn
      ref <-
        sp$pop[
          smok_status_ref == 4L,
          .(nsr = .N),
          keyby = .(year, age, sex, qimd)
        ] # ? add sha/ethn
      absorb_dt(ref, obs)
      setnafill(ref, "c", 0, cols = "nsa")
      ref[, `:=`(dif = nsr - nsa, nsr = NULL, nsa = NULL)]
      ref[cell_rn, on = .(year, age, sex, qimd), rn_clb := i.rn_clb]
      # Further calibrate to match better with HSE (deterministic via rank_smok_clb)
      ref[
        sex == "men" &
          rn_clb < 0.3,
        dif := dif - 1L
      ]
      ref[
        sex == "women" &
          qimd != "1 most deprived" &
          qimd != "5 least deprived" &
          rn_clb < 0.8,
        dif := dif - 1L
      ]
      ref[
        qimd == "1 most deprived" &
          rn_clb < 0.5,
        dif := dif + 1L
      ]
      ref[
        qimd == "5 least deprived" &
          rn_clb < 0.6,
        dif := dif - 1L
      ] # - reduces
      ref[, rn_clb := NULL]
      absorb_dt(sp$pop, ref)

      # when not enough active smokers convert those ex smokers with the shortest quit years
      tt <-
        sp$pop[
          smok_status == 3L,
          .(year, age, sex, qimd, pid, smok_quit_yrs, dif)
        ]
      setnafill(tt, "c", 0, cols = "dif")
      tt[dif < 0, dif := 0L]
      setkey(tt, year, age, sex, qimd, smok_quit_yrs)
      pid_to_conv <-
        tt[
          dif > 0,
          .(pid = head(pid, max(dif))),
          keyby = .(year, age, sex, qimd)
        ]
      sp$pop[
        pid_to_conv,
        on = .(year, pid),
        `:=`(
          smok_status = 4L,
          smok_quit_yrs = 0L,
          smok_dur = smok_dur + 1L
        )
      ] # TODO fix smoking duration

      # when too many never smokers convert to smok status 3
      tt <-
        sp$pop[
          smok_status == 4L,
          .(year, age, sex, qimd, pid, rank_smok_clb, dif)
        ]
      setnafill(tt, "c", 0, cols = "dif")
      tt[dif > 0, dif := 0L]
      tt[, dif := -dif]
      setkey(tt, year, age, sex, qimd, rank_smok_clb)
      # Ensure there are enough people to sample from
      ttt <-
        tt[,
          .(lpid = length(pid), mdif = max(dif)),
          by = .(year, age, sex, qimd)
        ][mdif > lpid, ]
      tt[ttt, on = .NATURAL, dif := i.lpid]
      pid_to_conv <-
        tt[
          dif > 0,
          .(pid = head(pid, max(dif))),
          keyby = .(year, age, sex, qimd)
        ]
      sp$pop[
        pid_to_conv,
        on = .(year, pid),
        `:=`(smok_status = 3L, smok_quit_yrs = 1L)
      ]
      sp$pop[, dif := NULL]
      sp$pop[, rank_smok_clb := NULL]
      rm(tt, ttt, obs, ref, pid_to_conv, cell_rn)
    }

    # Assign smok_cig
    set(sp$pop, NULL, "smok_cig", 0L)

    # smok_cig_curr for current smokers (smok_status == 4)
    idx <- sp$pop[, which(smok_status == 4L)]

    self$design$exposures$smok_cig_curr$generate(sp$pop, self$design, idx)

    # smok_cig_ex for ex-smokers at initialization (pid_mrk & smok_status == 3)

    idx <- sp$pop[, which((pid_mrk) & smok_status == 3L)]

    self$design$exposures$smok_cig_ex$generate(sp$pop, self$design, idx)

    simsmok_cig(sp$pop) # carry forward smok_cig if smok_status == 3

    sp$pop[smok_cig == 0L & smok_status != 1L, smok_cig := 1L]

    if (self$design$sim_prm$simsmok_calibration) {
      simsmok_postcalibration(sp$pop)
    } # need to be post cig simulation

    sp$pop[, smok_status := factor(smok_status)]

    colnam <- c("prb_smok_incid", "prb_smok_cess", "smok_status_ref")
    if (!self$design$sim_prm$keep_simulants_rn) {
      colnam <- c(
        colnam,
        "rankstat_simsmok"
      )
    }
    sp$pop[, (colnam) := NULL]


    # Generate ETS (BI) ----

    # Note at the moment this is independent of smoking prevalence TODO
    # calculate how many each smoker pollutes by year, SHA (not qimd) to
    # be used in scenarios. Ideally correct for mortality
    self$design$exposures$ets$generate(sp$pop, self$design)
    # NOTE for line above invert = true in sim_design yaml so 1 - mu to be
    # equivalent to qbinom(rank_ets, 1, mu). Otherwise rank_ets < mu is
    # equivalent to qbinom(rank_ets, 1, mu, lower.tail = FALSE)). The two
    # have the same prevalence, but correlations are captured correctly only
    # with the 1-mu variant.

    # View(sp$pop[, prop_if(ets == 1)/prop_if(smok_status == "4"), keyby = .(year, sha)])

    # Generate alcohol (ZINBI) ----
    self$design$exposures$alcohol$generate(sp$pop, self$design)

    # Generate BMI (BCPEo) ----
    self$design$exposures$bmi$generate(sp$pop, self$design)

    # Generate SBP (BCPEo) ----
    self$design$exposures$sbp$generate(sp$pop, self$design)

    # Generate BP medication (BI) -----
    # Temporarily round sbp for join, then restore
    sp$pop[, `:=`(
      sbp_acc = sbp,
      sbp = as.integer(round(clamp(sbp, 110, 200), -1))
    )]
    self$design$exposures$bp_med$generate(sp$pop, self$design)
    sp$pop[, `:=`(sbp = sbp_acc, sbp_acc = NULL)]

    # TODO calculate probability of dgn HTN

    # Generate tchol (BCT) ----
    self$design$exposures$tchol$generate(sp$pop, self$design)

    # Generate HDL (to tchol ratio) (GB1) ----
    self$design$exposures$hdl_to_tchol$generate(sp$pop, self$design)


    # NOTE this very highly correlated with hdl level (~0.76) and
    #  highly to tchol (~-0.47). The latter is captured by the correlated RNs

    # Generate statins medication (BI) -----
    # Temporarily round tchol for join, then restore
    sp$pop[, `:=`(tchol_acc = tchol, tchol = round(clamp(tchol, 2, 12), 0))]
    self$design$exposures$statin_px$generate(sp$pop, self$design)
    sp$pop[, `:=`(tchol = tchol_acc, tchol_acc = NULL)]

    sp$pop[,
      statin_adherence := qBE(
        rankstat_statin_adherence,
        self$design$sim_prm$statin_adherence,
        0.2
      )
    ]
    sp$pop[,
      bpmed_adherence := qBE(
        rankstat_bpmed_adherence,
        self$design$sim_prm$bpmed_adherence,
        0.2
      )
    ]
    if (!self$design$sim_prm$keep_simulants_rn) {
      sp$pop[,
        c("rankstat_bpmed_adherence", "rankstat_statin_adherence") := NULL
      ]
    }

    if ("age100" %in% names(sp$pop)) {
      sp$pop[, age := NULL]
      setnames(sp$pop, "age100", "age")
    }

    setnames(sp$pop, xps_nam, c(paste0(xps_nam, "_curr_xps")))
  }
)


IMPACTncd$run(1:n_runs, multicore = TRUE, "future_full_scn")


IMPACTncd$export_summaries(multicore = TRUE)

IMPACTncd$export_tables(multicore = TRUE)
