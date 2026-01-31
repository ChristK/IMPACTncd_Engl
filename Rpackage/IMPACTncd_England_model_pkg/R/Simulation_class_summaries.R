## IMPACTncdEngland is an implementation of the IMPACTncd framework, developed by Chris
## Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz.
##
## Copyright (C) 2018-2026 University of Liverpool, Chris Kypridemos
##
## IMPACTncdEngland is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version. This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details. You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/> or write
## to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
## Boston, MA 02110-1301 USA.

# -----------------------------------------------------------------------------
# Simulation class summary export methods
# This file adds summary export methods to the Simulation class using $set()
# -----------------------------------------------------------------------------


# export_summaries ----
# Exports simulation summaries from lifecourse files.
# See main class documentation in Simulation_class.R for details.
Simulation$set("public", "export_summaries", function(
    multicore = TRUE,
    type = c(
      "le",
      "hle",
      "dis_char",
      "prvl",
      "incd",
      "dis_mrtl",
      "mrtl",
      "all_cause_mrtl_by_dis",
      "cms",
      "qalys",
      "costs"
    ),
    single_year_of_age = FALSE
) {
  if (multicore) {
    arrow::set_cpu_count(1L)
    data.table::setDTthreads(threads = 1L, restore_after_fork = NULL)
    fst::threads_fst(nr_of_threads = 1L, reset_after_fork = NULL)
  } else {
    arrow::set_cpu_count(self$design$sim_prm$clusternumber_export)
    data.table::setDTthreads(
      threads = self$design$sim_prm$clusternumber_export,
      restore_after_fork = NULL
    )
    fst::threads_fst(
      nr_of_threads = self$design$sim_prm$clusternumber_export,
      reset_after_fork = NULL
    )
  }

  lc <- open_dataset(private$output_dir("lifecourse"))

  con <- NULL
  mc_set <- NULL

  tryCatch(
    {
      con <- dbConnect(duckdb::duckdb(), ":memory:", read_only = TRUE)
      if (multicore) {
        private$execute_sql(con, "PRAGMA threads=1")
      } else {
        private$execute_sql(
          con,
          paste0(
            "PRAGMA threads=",
            self$design$sim_prm$clusternumber_export
          )
        )
      }

      duckdb::duckdb_register_arrow(con, "lc_table", lc)
      mc_set <- private$query_sql(
        con,
        "SELECT DISTINCT mc FROM lc_table"
      )$mc
      scn_set <- private$query_sql(
        con,
        "SELECT DISTINCT scenario FROM lc_table"
      )$scenario
    },
    error = function(e) {
      stop(sprintf(
        "Failed to initialize DuckDB connection or query mc_set or query scn_set: %s",
        e$message
      ))
    },
    finally = {
      if (!is.null(con)) {
        try(dbDisconnect(con, shutdown = TRUE), silent = TRUE)
      }
    }
  )

  if (is.null(mc_set) || length(mc_set) == 0) {
    stop("No Monte Carlo iterations found in lifecourse data")
  }
  if (is.null(scn_set) || length(scn_set) == 0) {
    stop("No scenarios found in lifecourse data")
  }

  scn_set <- c(paste0("primary", scn_set), paste0("secondary", scn_set))
  rng_seeds <- sapply(mc_set, function(i) digest2int(scn_set, i))
  colnames(rng_seeds) <- paste0("mc=", mc_set)
  rownames(rng_seeds) <- gsub("primary|secondary", "", scn_set)
  duplicates <- duplicated(as.vector(rng_seeds))
  duplicates <- matrix(
    duplicates,
    nrow = nrow(rng_seeds),
    ncol = ncol(rng_seeds)
  )

  if (any(duplicates)) {
    dup_coords <- which(duplicates, arr.ind = TRUE)
    dup_row_names <- rownames(rng_seeds)[dup_coords[, 1]]
    warning(
      "RNG seeds are not unique for each scenario/mc combination. This may lead to correlated results. ",
      "Try to rename scenario ",
      paste(unique(dup_row_names), collapse = ", "),
      " and rerun the simulation to avoid unintended correlations."
    )
  }

  if ("le" %in% type) {
    file_pth <- private$output_dir("summaries/le_scaled_up")
  } else if ("hle" %in% type) {
    file_pth <- private$output_dir("summaries/hle_1st_cond_scaled_up")
  } else if ("cms" %in% type) {
    file_pth <- private$output_dir("summaries/cms_count_scaled_up")
  } else if ("mrtl" %in% type) {
    file_pth <- private$output_dir("summaries/mrtl_scaled_up")
  } else if ("dis_mrtl" %in% type) {
    file_pth <- private$output_dir("summaries/dis_mrtl_scaled_up")
  } else if ("dis_char" %in% type) {
    file_pth <- private$output_dir("summaries/dis_characteristics_scaled_up")
  } else if ("incd" %in% type) {
    file_pth <- private$output_dir("summaries/incd_scaled_up")
  } else if ("prvl" %in% type) {
    file_pth <- private$output_dir("summaries/prvl_scaled_up")
  } else if ("all_cause_mrtl_by_dis" %in% type) {
    file_pth <- private$output_dir("summaries/all_cause_mrtl_by_dis_scaled_up")
  } else if ("qalys" %in% type) {
    file_pth <- private$output_dir("summaries/qalys_scaled_up")
  } else if ("costs" %in% type) {
    file_pth <- private$output_dir("summaries/costs_scaled_up")
  } else {
    stop("Unknown type of summary")
  }

  if (file.exists(file_pth) && length(list.files(file_pth)) > 0) {
    con2 <- NULL
    mc_toexclude <- NULL

    tryCatch(
      {
        con2 <- dbConnect(duckdb::duckdb(), ":memory:", read_only = TRUE)
        if (multicore) {
          private$execute_sql(con2, "PRAGMA threads=1")
        } else {
          private$execute_sql(
            con2,
            paste0(
              "PRAGMA threads=",
              self$design$sim_prm$clusternumber_export
            )
          )
        }
        duckdb::duckdb_register_arrow(
          con2,
          "tbl",
          open_dataset(file_pth, format = "parquet")
        )
        mc_toexclude <- private$query_sql(
          con2,
          "SELECT DISTINCT mc FROM tbl"
        )$mc
      },
      error = function(e) {
        warning(sprintf("Failed to check existing files: %s", e$message))
        mc_toexclude <- NULL
      },
      finally = {
        if (!is.null(con2)) {
          try(dbDisconnect(con2, shutdown = TRUE), silent = TRUE)
        }
      }
    )

    if (!is.null(mc_toexclude)) {
      mc_set <- mc_set[!mc_set %in% mc_toexclude]
    }
  }

  if (length(mc_set) == 0) {
    if (self$design$sim_prm$logs) {
      message(
        "All required Monte Carlo iterations already processed for type: ",
        paste(type, collapse = ", "),
        ". Skipping summary export."
      )
    }
    return(invisible(self))
  }

  if (multicore) {
    if (self$design$sim_prm$logs) {
      private$time_mark("Start exporting summaries")
    }

    arrow::set_cpu_count(1L)
    data.table::setDTthreads(threads = 1L, restore_after_fork = NULL)
    fst::threads_fst(nr_of_threads = 1L, reset_after_fork = NULL)

    if (.Platform$OS.type == "windows") {
      cl <- makeClusterPSOCK(
        self$design$sim_prm$clusternumber_export,
        dryrun = FALSE,
        quiet = !self$design$sim_prm$logs,
        rscript_startup = quote(local({
          library(CKutils)
          library(IMPACTncdEngland)
          library(R6)
          library(arrow)
          library(duckdb)
          library(data.table)
        })),
        rscript_args = c("--no-init-file", "--no-site-file", "--no-environ"),
        setup_strategy = "parallel"
      )

      on.exit(if (exists("cl")) stopCluster(cl), add = TRUE)

      parLapplyLB(
        cl = cl,
        X = seq_along(mc_set),
        fun = function(i) {
          private$export_summaries_hlpr(
            mcaggr = i,
            type = type,
            single_year_of_age = single_year_of_age,
            implicit_parallelism = FALSE
          )
          NULL
        }
      )
    } else {
      registerDoParallel(self$design$sim_prm$clusternumber_export)
      xps_dt <- foreach(
        i = seq_along(mc_set),
        .inorder = TRUE,
        .options.multicore = list(preschedule = FALSE),
        .verbose = self$design$sim_prm$logs,
        .packages = c(
          "R6",
          "CKutils",
          "IMPACTncdEngland",
          "arrow",
          "duckdb",
          "data.table"
        ),
        .export = NULL,
        .noexport = NULL
      ) %dopar%
        {
          private$export_summaries_hlpr(
            mcaggr = i,
            type = type,
            single_year_of_age = single_year_of_age,
            implicit_parallelism = FALSE
          )
          NULL
        }
    }

    if (self$design$sim_prm$logs) {
      private$time_mark("End of exporting summaries")
    }
  } else {
    if (self$design$sim_prm$logs) {
      private$time_mark("Start of single-core summaries export")
    }

    lapply(seq_along(mc_set), function(i) {
      private$export_summaries_hlpr(
        mcaggr = i,
        type = type,
        single_year_of_age = single_year_of_age,
        implicit_parallelism = TRUE
      )
      NULL
    })

    if (self$design$sim_prm$logs) {
      private$time_mark("End of single-core summaries export")
    }
  }

  if (sink.number(type = "message") > 0L) {
    sink(type = "message")
  }
  while (sink.number(type = "output") > 0L) {
    sink(type = "output")
  }

  invisible(self)
})


# export_summaries_hlpr ----
Simulation$set("private", "export_summaries_hlpr", function(
    mcaggr,
    type = c(
      "le",
      "hle",
      "dis_char",
      "prvl",
      "incd",
      "mrtl",
      "dis_mrtl",
      "all_cause_mrtl_by_dis",
      "cms",
      "qalys",
      "costs"
    ),
    single_year_of_age = FALSE,
    implicit_parallelism
) {
  # Staggered worker startup to reduce memory pressure during parallel runs.
  # Only apply to the first batch of workers (iterations 1 to clusternumber_export).
  # After the first batch, workers finish at different times and naturally
  # pick up new iterations in a staggered manner.
  n_cores <- self$design$sim_prm$clusternumber_export
  if (n_cores > 1L && mcaggr <= n_cores) {
    stagger_delay_sec <- 2L
    worker_position <- mcaggr - 1L
    if (worker_position > 0L) {
      Sys.sleep(worker_position * stagger_delay_sec)
    }
  }

  if (self$design$sim_prm$logs) {
    private$time_mark(paste0("Start mc iteration (summary export) ", mcaggr))
    log_file <- private$output_dir(paste0("logs/log", mcaggr, ".txt"))
    log_con <- file(log_file, open = "at")
    sink(log_con, type = "output", split = FALSE)
    sink(log_con, type = "message")
  }

  if (self$design$sim_prm$logs) {
    message("Exporting summaries...")
  }

  if (implicit_parallelism) {
    arrow::set_cpu_count(self$design$sim_prm$clusternumber_export)
    data.table::setDTthreads(
      threads = self$design$sim_prm$clusternumber_export,
      restore_after_fork = NULL
    )
    fst::threads_fst(
      nr_of_threads = self$design$sim_prm$clusternumber_export,
      reset_after_fork = NULL
    )
  } else {
    arrow::set_cpu_count(1L)
    data.table::setDTthreads(threads = 1L, restore_after_fork = NULL)
    fst::threads_fst(nr_of_threads = 1L, reset_after_fork = NULL)
  }

  lc <- open_dataset(private$output_dir(file.path(
    "lifecourse",
    paste0("mc=", mcaggr)
  )))

  duckdb_con <- dbConnect(duckdb::duckdb(), ":memory:", read_only = FALSE)
  on.exit(dbDisconnect(duckdb_con, shutdown = FALSE), add = TRUE)
  if (implicit_parallelism) {
    private$execute_sql(
      duckdb_con,
      paste0("PRAGMA threads=", self$design$sim_prm$clusternumber_export)
    )
  } else {
    private$execute_sql(duckdb_con, "PRAGMA threads=1")
  }
  duckdb::duckdb_register_arrow(duckdb_con, "lc_table_raw", lc)

  private$execute_sql(
    duckdb_con,
    sprintf(
      "CREATE VIEW lc_table AS SELECT *, %d::INTEGER AS mc FROM lc_table_raw",
      mcaggr
    )
  )

  strata <- c("mc", self$design$sim_prm$strata_for_output)
  strata_noagegrp <- c(
    "mc",
    setdiff(self$design$sim_prm$strata_for_output, c("agegrp"))
  )
  strata_age <- c(strata_noagegrp, "age")

  if (single_year_of_age) {
    strata <- strata_age
  }

  ext <- "parquet"

  if ("le" %in% type) {
    private$export_le_summaries(duckdb_con, mcaggr, strata_noagegrp, ext)
  }
  if ("hle" %in% type) {
    private$export_hle_summaries(duckdb_con, mcaggr, strata_noagegrp, ext)
  }
  if ("dis_char" %in% type) {
    private$export_dis_char_summaries(duckdb_con, mcaggr, strata_noagegrp, ext)
  }
  if ("prvl" %in% type) {
    private$export_prvl_summaries(duckdb_con, mcaggr, strata, ext)
  }
  if ("incd" %in% type) {
    private$export_incd_summaries(duckdb_con, mcaggr, strata, ext)
  }
  if ("mrtl" %in% type) {
    private$export_mrtl_summaries(duckdb_con, mcaggr, strata, ext)
  }
  if ("dis_mrtl" %in% type) {
    private$export_dis_mrtl_summaries(duckdb_con, mcaggr, strata, ext)
  }
  if ("all_cause_mrtl_by_dis" %in% type) {
    private$export_all_cause_mrtl_by_dis_summaries(duckdb_con, mcaggr, strata, ext)
  }
  if ("cms" %in% type) {
    private$export_cms_summaries(duckdb_con, mcaggr, strata, strata_age, ext)
  }
  if ("qalys" %in% type) {
    private$export_qalys_summaries(duckdb_con, mcaggr, strata, ext)
  }
  if ("costs" %in% type) {
    private$export_costs_summaries(duckdb_con, mcaggr, strata, ext)
  }

  if (self$design$sim_prm$logs) {
    sink(type = "message")
    sink(type = "output")
    if (exists("log_con") && inherits(log_con, "connection") && isOpen(log_con)) {
      close(log_con)
    }
  }

  return(invisible(self))
})


# read_summary_dataset ----
Simulation$set("private", "read_summary_dataset", function(summary_type, standardization = "scaled_up") {
  fpth <- private$output_dir(
    paste0("summaries/", summary_type, "_", standardization)
  )
  if (!dir.exists(fpth)) {
    if (self$design$sim_prm$logs) {
      message(fpth, " doesn't exist, skipping...")
    }
    return(NULL)
  }
  tt <- CKutils::read_parquet_dt(fpth)
  return(tt)
})


# calc_QALYs ----
# Calculate QALYs based on Janssen & Szende 2014 (Table 3.7) and
# Sullivan et al. 2011 (Supplementary Tables 4 & 5)
# https://link.springer.com/book/10.1007/978-94-007-7596-1
# https://doi.org/10.1177/0272989X11401031
Simulation$set("private", "calc_QALYs", function(
    duckdb_con,
    mcaggr,
    input_table_name,
    output_view_name
) {
  # Population utility norms by age (from Janssen & Szende 2014)
  # NOTE: Age ranges based on original R implementation using (age-1) indexing
  # Ages 2-17: 1.0, 18-19: 0.94, 20-26: 0.922, 27-36: 0.914,
  # 37-46: 0.888, 47-56: 0.854, 57-66: 0.814, 67-76: 0.775, 77+: 0.706
  utility_pop_norm_expr <- "
    CASE
      WHEN age <= 17 THEN 1.0
      WHEN age BETWEEN 18 AND 19 THEN 0.94
      WHEN age BETWEEN 20 AND 26 THEN 0.922
      WHEN age BETWEEN 27 AND 36 THEN 0.914
      WHEN age BETWEEN 37 AND 46 THEN 0.888
      WHEN age BETWEEN 47 AND 56 THEN 0.854
      WHEN age BETWEEN 57 AND 66 THEN 0.814
      WHEN age BETWEEN 67 AND 76 THEN 0.775
      ELSE 0.706
    END
  "

  # Income utility adjustment (Sullivan et al. 2011)
  # near poor: 0.0150413, low income: 0.0380083, middle: 0.0396568, high: 0.0408501
  utility_income_expr <- "
    CASE income
      WHEN '1 Lowest' THEN 0.0
      WHEN '2' THEN 0.0150413
      WHEN '3' THEN 0.0380083
      WHEN '4' THEN 0.0396568
      WHEN '5 Highest' THEN 0.0408501
      ELSE 0.0
    END
  "

  # Education utility adjustment (Sullivan et al. 2011)
  # high_scl: 0.0028418, other_deg: 0.0056836, bachelor: 0.0060444, ma_phd: 0.0084228
  utility_education_expr <- "
    CASE education
      WHEN 'NVQ4/NVQ5/Degree or equiv' THEN 0.0068372
      WHEN 'Higher ed below degree' THEN 0.0056836
      WHEN 'NVQ3/GCE A Level equiv' THEN 0.0028418
      WHEN 'NVQ2/GCE O Level equiv' THEN 0.0028418
      WHEN 'NVQ1/CSE other grade equiv' THEN 0.0028418
      WHEN 'Foreign/other' THEN 0.0056836
      WHEN 'No qualification' THEN 0.0
      ELSE 0.0
    END
  "

  # Number of chronic conditions adjustment (Sullivan et al. 2011)
  # Based on cms_count clamped to 0-10
  utility_ncc_expr <- "
    CASE LEAST(GREATEST(CAST(cms_count AS INTEGER), 0), 10)
      WHEN 0 THEN 0.0
      WHEN 1 THEN 0.0
      WHEN 2 THEN -0.0528484
      WHEN 3 THEN -0.0415352
      WHEN 4 THEN -0.0202969
      WHEN 5 THEN 0.0083033
      WHEN 6 THEN 0.0408673
      WHEN 7 THEN 0.0668729
      WHEN 8 THEN 0.1158895
      WHEN 9 THEN 0.1344392
      WHEN 10 THEN 0.183614
      ELSE 0.0
    END
  "

  # Sex adjustment: +0.0010046 for men
  utility_sex_expr <- "
    CASE WHEN sex = 'men' THEN 0.0010046 ELSE 0.0 END
  "

  # Ethnicity adjustment (Sullivan et al. 2011)
  # black: -0.0000943, indian: -0.0014681, others: -0.0001887
  utility_ethnicity_expr <- "
    CASE
      WHEN ethnicity = 'white' THEN 0.0
      WHEN ethnicity IN ('black caribbean', 'black african') THEN -0.0000943
      WHEN ethnicity = 'indian' THEN -0.0014681
      WHEN ethnicity IN ('pakistani', 'bangladeshi', 'other asian', 'chinese', 'other') THEN -0.0001887
      ELSE 0.0
    END
  "

  # Disease-specific utility decrements (Sullivan et al. 2011)
  disease_decrements_expr <- "
    - CASE WHEN t2dm_prvl > 0 THEN 0.0714 ELSE 0.0 END
    - CASE WHEN ra_prvl > 0 THEN 0.1568 ELSE 0.0 END
    - CASE WHEN ckd_prvl > 0 THEN 0.1104 ELSE 0.0 END
    - CASE WHEN htn_prvl > 0 THEN 0.0460 ELSE 0.0 END
    - CASE WHEN asthma_prvl > 0 THEN 0.0463 ELSE 0.0 END
    - CASE WHEN helo_prvl > 0 THEN 0.0218 ELSE 0.0 END
    - CASE WHEN alcpr_prvl > 0 THEN 0.0315 ELSE 0.0 END
    - CASE WHEN prostate_ca_prvl > 0 THEN 0.0494 ELSE 0.0 END
    - CASE WHEN ibs_prvl > 0 THEN 0.0727 ELSE 0.0 END
    - CASE WHEN copd_prvl > 0 THEN 0.10291 ELSE 0.0 END
    - CASE WHEN t1dm_prvl > 0 THEN 0.0714 ELSE 0.0 END
    - CASE WHEN ctd_prvl > 0 THEN 0.0833 ELSE 0.0 END
    - CASE WHEN dementia_prvl > 0 THEN 0.2166 ELSE 0.0 END
    - CASE WHEN colorectal_ca_prvl > 0 THEN 0.0674 ELSE 0.0 END
    - CASE WHEN breast_ca_prvl > 0 THEN 0.0194 ELSE 0.0 END
    - CASE WHEN lung_ca_prvl > 0 THEN 0.1192 ELSE 0.0 END
    - CASE WHEN chd_prvl > 0 THEN 0.06713 ELSE 0.0 END
    - CASE WHEN other_ca_prvl > 0 THEN 0.02356 ELSE 0.0 END
    - CASE WHEN af_prvl > 0 THEN 0.0384 ELSE 0.0 END
    - CASE WHEN hf_prvl > 0 THEN 0.1167 ELSE 0.0 END
    - CASE WHEN pain_prvl > 0 THEN 0.34137 ELSE 0.0 END
    - CASE WHEN stroke_prvl > 0 THEN 0.10582 ELSE 0.0 END
    - CASE WHEN epilepsy_prvl > 0 THEN 0.0399 ELSE 0.0 END
    - CASE WHEN andep_prvl > 0 THEN 0.10509 ELSE 0.0 END
    - CASE WHEN psychosis_prvl > 0 THEN 0.11415 ELSE 0.0 END
    - CASE WHEN constipation_prvl > 0 THEN 0.104 ELSE 0.0 END
    - CASE WHEN obesity_prvl > 0 THEN 0.0709 ELSE 0.0 END
  "

  # Combine all components into the full EQ5D expression
  # The raw EQ5D value before death adjustment and clamping
  raw_eq5d_expr <- sprintf(
    "(%s) + (%s) + (%s) + (%s) + (%s) + (%s) %s",
    utility_pop_norm_expr,
    utility_income_expr,
    utility_education_expr,
    utility_ncc_expr,
    utility_sex_expr,
    utility_ethnicity_expr,
    disease_decrements_expr
  )

  # Apply half EQ5D for year of death and clamp to [0, 1]
  # Divide by 2 if dead (all_cause_mrtl > 0), divide by 1 if alive
  eq5d5l_expr <- sprintf(
    "LEAST(GREATEST((%s) / (CASE WHEN all_cause_mrtl > 0 THEN 2.0 ELSE 1.0 END), 0.0), 1.0)",
    raw_eq5d_expr
  )

  create_view_sql <- sprintf(
    "
    CREATE OR REPLACE TEMP VIEW %s AS
    SELECT
      *,
      (%s) AS EQ5D5L
    FROM %s
    WHERE mc = %d;
  ",
    output_view_name,
    eq5d5l_expr,
    input_table_name,
    mcaggr
  )

  private$execute_sql(duckdb_con, create_view_sql)

  NULL
})


# calc_costs ----
# Comprehensive cost and economic output calculation including healthcare, social care,
# informal care, and economic output (productivity value). Based on methodology from
# Sullivan et al. 2011, Janssen & Szende 2014, and UK-specific cost data.
Simulation$set("private", "calc_costs", function(
    duckdb_con,
    mcaggr,
    input_table_name,
    output_view_name
) {
  # Get scenario names
  scnams <- gsub(
    "^scenario=",
    "",
    list.dirs(
      private$output_dir(file.path("lifecourse", paste0("mc=", mcaggr))),
      full.names = FALSE,
      recursive = FALSE
    )
  )

  # ============================================================================
  # STEP 1: Create lookup views for parameters (eliminates repetitive assignments)
  # ============================================================================

  # --- Healthcare costs per disease (2026 values, per case per year) ---
  private$execute_sql(
    duckdb_con,
    "CREATE OR REPLACE TEMP VIEW healthcare_costs_view AS
    SELECT * FROM (VALUES
      ('no_disease', 1254.0),
      ('htn', 124.0), ('af', 1680.0), ('t2dm', 1116.0), ('t1dm', 4756.0),
      ('chd', 5444.0), ('stroke', 8673.0), ('copd', 2344.0), ('ckd', 1221.0),
      ('dementia', 2985.0), ('breast_ca', 2956.0), ('lung_ca', 2778.0),
      ('colorectal_ca', 3755.0), ('ra', 1192.0), ('asthma', 282.0),
      ('helo', 953.0), ('alcpr', 1074.0), ('prostate_ca', 3139.0),
      ('ibs', 514.0), ('ctd', 2002.0), ('other_ca', 3139.0), ('hf', 8233.0),
      ('pain', 795.0), ('epilepsy', 530.0), ('andep', 1348.0),
      ('psychosis', 3314.0), ('constipation', 479.0), ('obesity', 1836.0)
    ) AS t(disease, cost)",
    "healthcare_costs_view"
  )

  # --- Social care disease add-ons ---
  private$execute_sql(
    duckdb_con,
    "CREATE OR REPLACE TEMP VIEW socialcare_addons_view AS
    SELECT * FROM (VALUES
      ('cancer', 1419.082),
      ('chd', 1868.744),
      ('dementia', 12412.006),
      ('stroke', 8417.945)
    ) AS t(disease, addon_cost)",
    "socialcare_addons_view"
  )

  # --- Social care age-based costs ---
  private$execute_sql(
    duckdb_con,
    "CREATE OR REPLACE TEMP VIEW socialcare_age_costs_view AS
    SELECT age, costs FROM (VALUES
      (18, 43.77339), (19, 100.24136), (20, 151.81292), (21, 208.28088), (22, 264.74884),
      (23, 264.74884), (24, 264.74884), (25, 257.53188), (26, 257.53188), (27, 257.53188),
      (28, 257.53188), (29, 257.53188), (30, 247.51041), (31, 247.51041), (32, 247.51041),
      (33, 247.51041), (34, 247.51041), (35, 233.63486), (36, 233.63486), (37, 233.63486),
      (38, 233.63486), (39, 233.63486), (40, 211.38868), (41, 211.38868), (42, 211.38868),
      (43, 211.38868), (44, 211.38868), (45, 167.03307), (46, 218.86003), (47, 270.68700),
      (48, 322.51396), (49, 374.34093), (50, 359.68800), (51, 359.68800), (52, 359.68800),
      (53, 359.68800), (54, 359.68800), (55, 269.06953), (56, 269.06953), (57, 269.06953),
      (58, 269.06953), (59, 269.06953), (60, 138.55678), (61, 138.55678), (62, 138.55678),
      (63, 138.55678), (64, 138.55678), (65, 533.29864), (66, 478.54140), (67, 423.78417),
      (68, 369.02693), (69, 314.26969), (70, 151.53999), (71, 297.48722), (72, 443.43444),
      (73, 589.38167), (74, 735.32890), (75, 434.09010), (76, 639.63391), (77, 845.17772),
      (78, 1050.72153), (79, 1256.26534), (80, 1182.72718), (81, 1873.89682), (82, 2565.06646),
      (83, 3256.23609), (84, 3947.40573), (85, 4683.05748), (86, 6384.48774), (87, 8085.91799),
      (88, 9787.34825), (89, 11488.77851), (90, 10609.59492), (91, 10609.59492), (92, 10609.59492),
      (93, 10609.59492), (94, 10609.59492), (95, 9815.27655), (96, 9815.27655), (97, 9815.27655),
      (98, 9815.27655), (99, 9815.27655), (100, 9815.17399)
    ) AS t(age, costs)",
    "socialcare_age_costs_view"
  )

  # --- Employment rates by age and sex (combined lookup table) ---
  private$execute_sql(
    duckdb_con,
    "CREATE OR REPLACE TEMP VIEW employment_params_view AS
    SELECT age_low, age_high, sex, employment_rate, paid_hrs_week, hrs_pay FROM (VALUES
      -- Men
      (18, 21, 'men', 61.1, 27.2, 12.5),
      (22, 24, 'men', 61.1, 36.5, 18.1),
      (25, 29, 'men', 88.3, 36.5, 18.1),
      (30, 34, 'men', 88.3, 37.6, 23.6),
      (35, 39, 'men', 90.1, 37.6, 23.6),
      (40, 49, 'men', 90.1, 37.3, 26.7),
      (50, 59, 'men', 75.3, 37.2, 26.1),
      (60, 64, 'men', 75.3, 33.5, 22.5),
      (65, 70, 'men', 15.6, 33.5, 22.5),
      (71, 99, 'men', 0.0, 0.0, 0.0),
      -- Women
      (18, 21, 'women', 60.0, 22.1, 12.5),
      (22, 24, 'women', 60.0, 33.1, 17.2),
      (25, 29, 'women', 80.5, 33.1, 17.2),
      (30, 34, 'women', 80.5, 31.4, 20.8),
      (35, 39, 'women', 81.6, 31.4, 20.8),
      (40, 49, 'women', 81.6, 30.7, 22.1),
      (50, 59, 'women', 68.0, 30.2, 20.8),
      (60, 64, 'women', 68.0, 25.5, 18.0),
      (65, 70, 'women', 10.3, 25.5, 18.0),
      (71, 99, 'women', 0.0, 0.0, 0.0)
    ) AS t(age_low, age_high, sex, employment_rate, paid_hrs_week, hrs_pay)",
    "employment_params_view"
  )

  # --- Informal care regression coefficients ---
  # Source: https://doi.org/10.1007/s10198-015-0718-5 (Supplement)
  private$execute_sql(
    duckdb_con,
    "CREATE OR REPLACE TEMP VIEW informal_care_params_view AS
    SELECT * FROM (VALUES
      -- First regression (intensity)
      ('cons1', 2.6540856), ('sex1', -0.0226854), ('age1', 0.0194834),
      ('age_sq1', -0.0001203), ('utility1', -0.8583409), ('comorbidity1', 0.1478759),
      -- ICD chapter coefficients
      ('C', 0.0404685), ('E', -0.0209114), ('F', -0.2746472), ('G', 0.0016516),
      ('H', -0.2065952), ('I', -0.0260221), ('J', -0.1084329), ('K', -0.220285),
      ('M', -0.0010566), ('N', -0.2357236),
      -- Second regression (probability)
      ('cons2', -3.342591), ('sex2', -0.56294), ('age2', 0.0482827),
      ('age_sq2', -0.0004012), ('utility2', 4.122554), ('comorbidity2', -0.3932173),
      -- Care cost parameters (2026 values)
      ('hours_per_day', 3.0), ('cost_per_hour', 10.4)
    ) AS t(param_name, param_value)",
    "informal_care_params_view"
  )

  # ============================================================================
  # STEP 2: Create comprehensive cost calculation view
  # ============================================================================

  # --- Healthcare cost expression (sum of disease-specific costs) ---
  healthcare_cost_expr <- "
    1254.0  -- base cost (no disease)
    + 124.0 * CASE WHEN htn_prvl > 0 THEN 1 ELSE 0 END
    + 1680.0 * CASE WHEN af_prvl > 0 THEN 1 ELSE 0 END
    + 1116.0 * CASE WHEN t2dm_prvl > 0 THEN 1 ELSE 0 END
    + 4756.0 * CASE WHEN t1dm_prvl > 0 THEN 1 ELSE 0 END
    + 5444.0 * CASE WHEN chd_prvl > 0 THEN 1 ELSE 0 END
    + 8673.0 * CASE WHEN stroke_prvl > 0 THEN 1 ELSE 0 END
    + 2344.0 * CASE WHEN copd_prvl > 0 THEN 1 ELSE 0 END
    + 1221.0 * CASE WHEN ckd_prvl > 0 THEN 1 ELSE 0 END
    + 2985.0 * CASE WHEN dementia_prvl > 0 THEN 1 ELSE 0 END
    + 2956.0 * CASE WHEN breast_ca_prvl > 0 THEN 1 ELSE 0 END
    + 2778.0 * CASE WHEN lung_ca_prvl > 0 THEN 1 ELSE 0 END
    + 3755.0 * CASE WHEN colorectal_ca_prvl > 0 THEN 1 ELSE 0 END
    + 1192.0 * CASE WHEN ra_prvl > 0 THEN 1 ELSE 0 END
    + 282.0 * CASE WHEN asthma_prvl > 0 THEN 1 ELSE 0 END
    + 953.0 * CASE WHEN helo_prvl > 0 THEN 1 ELSE 0 END
    + 1074.0 * CASE WHEN alcpr_prvl > 0 THEN 1 ELSE 0 END
    + 3139.0 * CASE WHEN prostate_ca_prvl > 0 THEN 1 ELSE 0 END
    + 514.0 * CASE WHEN ibs_prvl > 0 THEN 1 ELSE 0 END
    + 2002.0 * CASE WHEN ctd_prvl > 0 THEN 1 ELSE 0 END
    + 3139.0 * CASE WHEN other_ca_prvl > 0 THEN 1 ELSE 0 END
    + 8233.0 * CASE WHEN hf_prvl > 0 THEN 1 ELSE 0 END
    + 795.0 * CASE WHEN pain_prvl > 0 THEN 1 ELSE 0 END
    + 530.0 * CASE WHEN epilepsy_prvl > 0 THEN 1 ELSE 0 END
    + 1348.0 * CASE WHEN andep_prvl > 0 THEN 1 ELSE 0 END
    + 3314.0 * CASE WHEN psychosis_prvl > 0 THEN 1 ELSE 0 END
    + 479.0 * CASE WHEN constipation_prvl > 0 THEN 1 ELSE 0 END
    + 1836.0 * CASE WHEN obesity_prvl > 0 THEN 1 ELSE 0 END
  "

  # --- EQ5D calculation (reuses logic from calc_QALYs for consistency) ---
  # Population utility norms by age (Janssen & Szende 2014)
  # NOTE: Age ranges based on original R implementation using (age-1) indexing
  # Using lc.age to avoid ambiguity with socialcare_age_costs_view.age in JOIN
  utility_pop_norm_expr <- "
    CASE
      WHEN lc.age <= 17 THEN 1.0
      WHEN lc.age BETWEEN 18 AND 19 THEN 0.94
      WHEN lc.age BETWEEN 20 AND 26 THEN 0.922
      WHEN lc.age BETWEEN 27 AND 36 THEN 0.914
      WHEN lc.age BETWEEN 37 AND 46 THEN 0.888
      WHEN lc.age BETWEEN 47 AND 56 THEN 0.854
      WHEN lc.age BETWEEN 57 AND 66 THEN 0.814
      WHEN lc.age BETWEEN 67 AND 76 THEN 0.775
      ELSE 0.706
    END
  "

  # Income utility adjustment (Sullivan et al. 2011)
  utility_income_expr <- "
    CASE income
      WHEN '1 Lowest' THEN 0.0
      WHEN '2' THEN 0.0150413
      WHEN '3' THEN 0.0380083
      WHEN '4' THEN 0.0396568
      WHEN '5 Highest' THEN 0.0408501
      ELSE 0.0
    END
  "

  # Education utility adjustment
  utility_education_expr <- "
    CASE education
      WHEN 'NVQ4/NVQ5/Degree or equiv' THEN 0.0068372
      WHEN 'Higher ed below degree' THEN 0.0056836
      WHEN 'NVQ3/GCE A Level equiv' THEN 0.0028418
      WHEN 'NVQ2/GCE O Level equiv' THEN 0.0028418
      WHEN 'NVQ1/CSE other grade equiv' THEN 0.0028418
      WHEN 'Foreign/other' THEN 0.0056836
      WHEN 'No qualification' THEN 0.0
      ELSE 0.0
    END
  "

  # Number of chronic conditions adjustment
  utility_ncc_expr <- "
    CASE LEAST(GREATEST(CAST(cms_count AS INTEGER), 0), 10)
      WHEN 0 THEN 0.0 WHEN 1 THEN 0.0 WHEN 2 THEN -0.0528484
      WHEN 3 THEN -0.0415352 WHEN 4 THEN -0.0202969 WHEN 5 THEN 0.0083033
      WHEN 6 THEN 0.0408673 WHEN 7 THEN 0.0668729 WHEN 8 THEN 0.1158895
      WHEN 9 THEN 0.1344392 WHEN 10 THEN 0.183614 ELSE 0.0
    END
  "

  # Disease-specific utility decrements
  disease_decrements_expr <- "
    - 0.0714 * CASE WHEN t2dm_prvl > 0 THEN 1 ELSE 0 END
    - 0.1568 * CASE WHEN ra_prvl > 0 THEN 1 ELSE 0 END
    - 0.1104 * CASE WHEN ckd_prvl > 0 THEN 1 ELSE 0 END
    - 0.0460 * CASE WHEN htn_prvl > 0 THEN 1 ELSE 0 END
    - 0.0463 * CASE WHEN asthma_prvl > 0 THEN 1 ELSE 0 END
    - 0.0218 * CASE WHEN helo_prvl > 0 THEN 1 ELSE 0 END
    - 0.0315 * CASE WHEN alcpr_prvl > 0 THEN 1 ELSE 0 END
    - 0.0494 * CASE WHEN prostate_ca_prvl > 0 THEN 1 ELSE 0 END
    - 0.0727 * CASE WHEN ibs_prvl > 0 THEN 1 ELSE 0 END
    - 0.10291 * CASE WHEN copd_prvl > 0 THEN 1 ELSE 0 END
    - 0.0714 * CASE WHEN t1dm_prvl > 0 THEN 1 ELSE 0 END
    - 0.0833 * CASE WHEN ctd_prvl > 0 THEN 1 ELSE 0 END
    - 0.2166 * CASE WHEN dementia_prvl > 0 THEN 1 ELSE 0 END
    - 0.0674 * CASE WHEN colorectal_ca_prvl > 0 THEN 1 ELSE 0 END
    - 0.0194 * CASE WHEN breast_ca_prvl > 0 THEN 1 ELSE 0 END
    - 0.1192 * CASE WHEN lung_ca_prvl > 0 THEN 1 ELSE 0 END
    - 0.06713 * CASE WHEN chd_prvl > 0 THEN 1 ELSE 0 END
    - 0.02356 * CASE WHEN other_ca_prvl > 0 THEN 1 ELSE 0 END
    - 0.0384 * CASE WHEN af_prvl > 0 THEN 1 ELSE 0 END
    - 0.1167 * CASE WHEN hf_prvl > 0 THEN 1 ELSE 0 END
    - 0.34137 * CASE WHEN pain_prvl > 0 THEN 1 ELSE 0 END
    - 0.10582 * CASE WHEN stroke_prvl > 0 THEN 1 ELSE 0 END
    - 0.0399 * CASE WHEN epilepsy_prvl > 0 THEN 1 ELSE 0 END
    - 0.10509 * CASE WHEN andep_prvl > 0 THEN 1 ELSE 0 END
    - 0.11415 * CASE WHEN psychosis_prvl > 0 THEN 1 ELSE 0 END
    - 0.104 * CASE WHEN constipation_prvl > 0 THEN 1 ELSE 0 END
    - 0.0709 * CASE WHEN obesity_prvl > 0 THEN 1 ELSE 0 END
  "

  # Combined EQ5D expression (clamped to [0,1])
  eq5d_raw_expr <- sprintf(
    "(%s) + (%s) + (%s) + (%s) + 0.0010046 * CASE WHEN sex = 'men' THEN 1 ELSE 0 END
     + CASE
         WHEN ethnicity = 'white' THEN 0.0
         WHEN ethnicity IN ('black caribbean', 'black african') THEN -0.0000943
         WHEN ethnicity = 'indian' THEN -0.0014681
         ELSE -0.0001887
       END
     %s",
    utility_pop_norm_expr,
    utility_income_expr,
    utility_education_expr,
    utility_ncc_expr,
    disease_decrements_expr
  )

  eq5d_expr <- sprintf(
    "LEAST(GREATEST((%s) / CASE WHEN all_cause_mrtl > 0 THEN 2.0 ELSE 1.0 END, 0.0), 1.0)",
    eq5d_raw_expr
  )

  # --- Comorbidity status for informal care regression ---
  comorbidity_status_expr <- "
    (CASE WHEN prostate_ca_prvl > 0 OR breast_ca_prvl > 0 OR lung_ca_prvl > 0 OR colorectal_ca_prvl > 0 OR other_ca_prvl > 0 THEN 1 ELSE 0 END)
    + (CASE WHEN obesity_prvl > 0 OR t1dm_prvl > 0 OR t2dm_prvl > 0 THEN 1 ELSE 0 END)
    + (CASE WHEN dementia_prvl > 0 OR andep_prvl > 0 OR psychosis_prvl > 0 OR alcpr_prvl > 0 THEN 1 ELSE 0 END)
    + (CASE WHEN epilepsy_prvl > 0 THEN 1 ELSE 0 END)
    + (CASE WHEN helo_prvl > 0 THEN 1 ELSE 0 END)
    + (CASE WHEN htn_prvl > 0 OR chd_prvl > 0 OR af_prvl > 0 OR hf_prvl > 0 OR stroke_prvl > 0 THEN 1 ELSE 0 END)
    + (CASE WHEN copd_prvl > 0 OR asthma_prvl > 0 THEN 1 ELSE 0 END)
    + (CASE WHEN ibs_prvl > 0 OR constipation_prvl > 0 THEN 1 ELSE 0 END)
    + (CASE WHEN ctd_prvl > 0 OR pain_prvl > 0 OR ra_prvl > 0 THEN 1 ELSE 0 END)
    + (CASE WHEN ckd_prvl > 0 THEN 1 ELSE 0 END)
  "

  # --- ICD chapter coefficients for informal care ---
  icd_chapter_expr <- "
    0.0404685 * CASE WHEN prostate_ca_prvl > 0 OR breast_ca_prvl > 0 OR lung_ca_prvl > 0 OR colorectal_ca_prvl > 0 OR other_ca_prvl > 0 THEN 1 ELSE 0 END
    + (-0.0209114) * CASE WHEN obesity_prvl > 0 OR t1dm_prvl > 0 OR t2dm_prvl > 0 THEN 1 ELSE 0 END
    + (-0.2746472) * CASE WHEN dementia_prvl > 0 OR andep_prvl > 0 OR psychosis_prvl > 0 OR alcpr_prvl > 0 THEN 1 ELSE 0 END
    + 0.0016516 * CASE WHEN epilepsy_prvl > 0 THEN 1 ELSE 0 END
    + (-0.2065952) * CASE WHEN helo_prvl > 0 THEN 1 ELSE 0 END
    + (-0.0260221) * CASE WHEN htn_prvl > 0 OR chd_prvl > 0 OR af_prvl > 0 OR hf_prvl > 0 OR stroke_prvl > 0 THEN 1 ELSE 0 END
    + (-0.1084329) * CASE WHEN copd_prvl > 0 OR asthma_prvl > 0 THEN 1 ELSE 0 END
    + (-0.220285) * CASE WHEN ibs_prvl > 0 OR constipation_prvl > 0 THEN 1 ELSE 0 END
    + (-0.0010566) * CASE WHEN ctd_prvl > 0 OR pain_prvl > 0 OR ra_prvl > 0 THEN 1 ELSE 0 END
    + (-0.2357236) * CASE WHEN ckd_prvl > 0 THEN 1 ELSE 0 END
  "

  # --- Employment and productivity parameters by age/sex ---
  # Using lc.age to avoid ambiguity with socialcare_age_costs_view.age in JOIN
  employment_rate_expr <- "
    CASE
      WHEN sex = 'men' AND lc.age BETWEEN 18 AND 24 THEN 61.1
      WHEN sex = 'men' AND lc.age BETWEEN 25 AND 34 THEN 88.3
      WHEN sex = 'men' AND lc.age BETWEEN 35 AND 49 THEN 90.1
      WHEN sex = 'men' AND lc.age BETWEEN 50 AND 64 THEN 75.3
      WHEN sex = 'men' AND lc.age BETWEEN 65 AND 70 THEN 15.6
      WHEN sex = 'women' AND lc.age BETWEEN 18 AND 24 THEN 60.0
      WHEN sex = 'women' AND lc.age BETWEEN 25 AND 34 THEN 80.5
      WHEN sex = 'women' AND lc.age BETWEEN 35 AND 49 THEN 81.6
      WHEN sex = 'women' AND lc.age BETWEEN 50 AND 64 THEN 68.0
      WHEN sex = 'women' AND lc.age BETWEEN 65 AND 70 THEN 10.3
      ELSE 0.0
    END
  "

  paid_hrs_week_expr <- "
    CASE
      WHEN sex = 'men' AND lc.age BETWEEN 18 AND 21 THEN 27.2
      WHEN sex = 'men' AND lc.age BETWEEN 22 AND 29 THEN 36.5
      WHEN sex = 'men' AND lc.age BETWEEN 30 AND 39 THEN 37.6
      WHEN sex = 'men' AND lc.age BETWEEN 40 AND 49 THEN 37.3
      WHEN sex = 'men' AND lc.age BETWEEN 50 AND 59 THEN 37.2
      WHEN sex = 'men' AND lc.age >= 60 THEN 33.5
      WHEN sex = 'women' AND lc.age BETWEEN 18 AND 21 THEN 22.1
      WHEN sex = 'women' AND lc.age BETWEEN 22 AND 29 THEN 33.1
      WHEN sex = 'women' AND lc.age BETWEEN 30 AND 39 THEN 31.4
      WHEN sex = 'women' AND lc.age BETWEEN 40 AND 49 THEN 30.7
      WHEN sex = 'women' AND lc.age BETWEEN 50 AND 59 THEN 30.2
      WHEN sex = 'women' AND lc.age >= 60 THEN 25.5
      ELSE 0.0
    END
  "

  hrs_pay_expr <- "
    CASE
      WHEN sex = 'men' AND lc.age BETWEEN 18 AND 21 THEN 12.5
      WHEN sex = 'men' AND lc.age BETWEEN 22 AND 29 THEN 18.1
      WHEN sex = 'men' AND lc.age BETWEEN 30 AND 39 THEN 23.6
      WHEN sex = 'men' AND lc.age BETWEEN 40 AND 49 THEN 26.7
      WHEN sex = 'men' AND lc.age BETWEEN 50 AND 59 THEN 26.1
      WHEN sex = 'men' AND lc.age >= 60 THEN 22.5
      WHEN sex = 'women' AND lc.age BETWEEN 18 AND 21 THEN 12.5
      WHEN sex = 'women' AND lc.age BETWEEN 22 AND 29 THEN 17.2
      WHEN sex = 'women' AND lc.age BETWEEN 30 AND 39 THEN 20.8
      WHEN sex = 'women' AND lc.age BETWEEN 40 AND 49 THEN 22.1
      WHEN sex = 'women' AND lc.age BETWEEN 50 AND 59 THEN 20.8
      WHEN sex = 'women' AND lc.age >= 60 THEN 18.0
      ELSE 0.0
    END
  "

  # Unpaid hours per year (capped at age 70)
  unpaid_hrs_yr_expr <- "
    CASE
      WHEN sex = 'men' THEN 12.0 * (27.92 + 1.79 * LEAST(lc.age, 70))
      WHEN sex = 'women' THEN 12.0 * (50.03 + 2.3 * LEAST(lc.age, 70))
      ELSE 0.0
    END
  "

  # ============================================================================
  # STEP 3: Create final cost output view for each scenario
  # ============================================================================

  final_view_sql <- sprintf(
    "
    CREATE OR REPLACE TEMP VIEW %s AS
    WITH base_data AS (
      SELECT lc.*,
        -- EQ5D utility value
        (%s) AS eq5d,

        -- Comorbidity status for informal care
        (%s) AS comorbidity_status,

        -- Employment parameters
        (%s) AS employment_rate,
        (%s) AS paid_hrs_week,
        (%s) AS hrs_pay,
        (%s) AS unpaid_hrs_yr,

        -- Social care age-based cost from lookup table
        COALESCE(sc.costs, 0.0) AS socialcare_age_cost

      FROM %s lc
      LEFT JOIN socialcare_age_costs_view sc ON lc.age = sc.age
      WHERE lc.mc = %d AND lc.scenario = %s
    ),
    costs_calculated AS (
      SELECT *,

        -- Healthcare cost (halved for year of death)
        (%s) / CASE WHEN all_cause_mrtl > 0 THEN 2.0 ELSE 1.0 END AS healthcare_cost,

        -- Social care cost (age-based from FST + disease add-ons from FST)
        (socialcare_age_cost
        + 1419.082 * CASE WHEN prostate_ca_prvl > 0 OR breast_ca_prvl > 0 OR lung_ca_prvl > 0 OR colorectal_ca_prvl > 0 OR other_ca_prvl > 0 THEN 1 ELSE 0 END
        + 1868.744 * CASE WHEN chd_prvl > 0 THEN 1 ELSE 0 END
        + 12412.006 * CASE WHEN dementia_prvl > 0 THEN 1 ELSE 0 END
        + 8417.945 * CASE WHEN stroke_prvl > 0 THEN 1 ELSE 0 END
        ) / CASE WHEN all_cause_mrtl > 0 THEN 2.0 ELSE 1.0 END AS socialcare_cost,

        -- Informal care cost (two-stage regression model)
        -- First regression (intensity): cons1 + sex1*female + age1*age + age_sq1*age^2 + utility1*eq5d + comorbidity1*(comorbidity>1) + ICD_chapters
        (
          (1.0 - (EXP(
            -3.342591 + (-0.56294) * CASE WHEN sex = 'women' THEN 1 ELSE 0 END
            + 0.0482827 * age + (-0.0004012) * age * age
            + 4.122554 * eq5d + (-0.3932173) * CASE WHEN comorbidity_status > 1 THEN 1 ELSE 0 END
          ) / (1.0 + EXP(
            -3.342591 + (-0.56294) * CASE WHEN sex = 'women' THEN 1 ELSE 0 END
            + 0.0482827 * age + (-0.0004012) * age * age
            + 4.122554 * eq5d + (-0.3932173) * CASE WHEN comorbidity_status > 1 THEN 1 ELSE 0 END
          ))))
          * EXP(
            2.6540856 + (-0.0226854) * CASE WHEN sex = 'women' THEN 1 ELSE 0 END
            + 0.0194834 * age + (-0.0001203) * age * age
            + (-0.8583409) * eq5d + 0.1478759 * CASE WHEN comorbidity_status > 1 THEN 1 ELSE 0 END
            + (%s)
          )
          / 42.0 * 365.0 * 3.0 * 10.4
        ) / CASE WHEN all_cause_mrtl > 0 THEN 2.0 ELSE 1.0 END AS informalcare_cost,

        -- Productivity cost calculation
        -- paid_yr = (employment_rate/100) * (paid_hrs_week/7) * hrs_pay * 365 * 1.07
        -- unpaid_yr = unpaid_hrs_yr * 13.8
        -- total_yr = paid_yr + unpaid_yr
        -- Relative productivity based on eq5d vs full health
        (
          CASE WHEN age < 18 THEN 0.0 ELSE
            -- Full health productivity (eq5d = 1)
            (
              -- paid_yr at full health
              (employment_rate / 100.0) * (paid_hrs_week / 7.0) * hrs_pay * 365.0 * 1.07
              -- unpaid_yr
              + unpaid_hrs_yr * 13.8
            )
            -- Relative productivity factor based on eq5d
            * (
              -- produc_eq5d / produc_full
              (EXP(
                2.95 * (age / 10.0) + (-0.35) * POWER(age / 10.0, 2)
                + 1.37 * ((-1.0443 * (age / 10.0) + 25.918 * eq5d + 31.0231) / 10.0)
                + (-0.09) * POWER((-1.0443 * (age / 10.0) + 25.918 * eq5d + 31.0231) / 10.0, 2)
                + 1.19 * ((1.0383 * (age / 10.0) + 5.0122 * eq5d + 32.5459) / 10.0)
                + (-0.09) * POWER((1.0383 * (age / 10.0) + 5.0122 * eq5d + 32.5459) / 10.0, 2)
                - 13.2
              ) / (1.0 + EXP(
                2.95 * (age / 10.0) + (-0.35) * POWER(age / 10.0, 2)
                + 1.37 * ((-1.0443 * (age / 10.0) + 25.918 * eq5d + 31.0231) / 10.0)
                + (-0.09) * POWER((-1.0443 * (age / 10.0) + 25.918 * eq5d + 31.0231) / 10.0, 2)
                + 1.19 * ((1.0383 * (age / 10.0) + 5.0122 * eq5d + 32.5459) / 10.0)
                + (-0.09) * POWER((1.0383 * (age / 10.0) + 5.0122 * eq5d + 32.5459) / 10.0, 2)
                - 13.2
              )))
              /
              NULLIF(
                (EXP(
                  2.95 * (age / 10.0) + (-0.35) * POWER(age / 10.0, 2)
                  + 1.37 * ((-1.0443 * (age / 10.0) + 25.918 * 1.0 + 31.0231) / 10.0)
                  + (-0.09) * POWER((-1.0443 * (age / 10.0) + 25.918 * 1.0 + 31.0231) / 10.0, 2)
                  + 1.19 * ((1.0383 * (age / 10.0) + 5.0122 * 1.0 + 32.5459) / 10.0)
                  + (-0.09) * POWER((1.0383 * (age / 10.0) + 5.0122 * 1.0 + 32.5459) / 10.0, 2)
                  - 13.2
                ) / (1.0 + EXP(
                  2.95 * (age / 10.0) + (-0.35) * POWER(age / 10.0, 2)
                  + 1.37 * ((-1.0443 * (age / 10.0) + 25.918 * 1.0 + 31.0231) / 10.0)
                  + (-0.09) * POWER((-1.0443 * (age / 10.0) + 25.918 * 1.0 + 31.0231) / 10.0, 2)
                  + 1.19 * ((1.0383 * (age / 10.0) + 5.0122 * 1.0 + 32.5459) / 10.0)
                  + (-0.09) * POWER((1.0383 * (age / 10.0) + 5.0122 * 1.0 + 32.5459) / 10.0, 2)
                  - 13.2
                ))),
                0.0
              )
            )
          END
        ) / CASE WHEN all_cause_mrtl > 0 THEN 2.0 ELSE 1.0 END AS economic_output

      FROM base_data
    )
    SELECT
      mc, scenario, year, agegrp, sex, dimd, wt, wt_esp,

      -- Individual cost components
      healthcare_cost,
      socialcare_cost,
      informalcare_cost,
      economic_output,

      -- Aggregated costs (subtracting economic_output as it represents value produced, not cost)
      (socialcare_cost + informalcare_cost - economic_output) AS indirect_cost,
      (healthcare_cost + socialcare_cost + informalcare_cost - economic_output) AS total_cost

    FROM costs_calculated;
    ",
    paste0(output_view_name, "_", scnams, "_view"),
    eq5d_expr,
    comorbidity_status_expr,
    employment_rate_expr,
    paid_hrs_week_expr,
    hrs_pay_expr,
    unpaid_hrs_yr_expr,
    input_table_name,
    mcaggr,
    paste0("'", scnams, "'"),
    healthcare_cost_expr,
    icd_chapter_expr
  )

  # Execute view creation for each scenario
  sapply(final_view_sql, function(sql) {
    private$execute_sql(
      duckdb_con,
      sql,
      paste0("Cost calculation view: ", output_view_name)
    )
  })

  return(invisible(NULL))
})


# export_le_summaries ----
Simulation$set("private", "export_le_summaries", function(
    duckdb_con,
    mcaggr,
    strata_noagegrp,
    ext
) {
  lapply(
    paste0(rep(c("le", "le60"), each = 2), "_", c("scaled_up", "esp")),
    function(subdir) {
      private$create_new_folder(private$output_dir(paste0(
        "summaries/",
        subdir
      )))
    }
  )

  group_by_cols <- paste(strata_noagegrp, collapse = ", ")

  # Define configurations for LE calculations
  le_configs <- list(
    list(
      prefix = "le",
      age_filter = "",
      weight_col = "wt",
      suffix = "_scaled_up"
    ),
    list(
      prefix = "le",
      age_filter = "",
      weight_col = "wt_esp",
      suffix = "_esp"
    ),
    list(
      prefix = "le60",
      age_filter = "AND age > 60",
      weight_col = "wt",
      suffix = "_scaled_up"
    ),
    list(
      prefix = "le60",
      age_filter = "AND age > 60",
      weight_col = "wt_esp",
      suffix = "_esp"
    )
  )

  for (config in le_configs) {
    # Skip LE60 calculation if age range doesn't cover 60
    if (
      startsWith(config$prefix, "le60") &&
        !(self$design$sim_prm$ageL < 60L &&
          self$design$sim_prm$ageH > 60L)
    ) {
      next
    }

    query <- sprintf(
      "SELECT %s, SUM(%s) AS popsize, SUM(age * %s) / NULLIF(SUM(%s), 0) AS LE
      FROM lc_table
      WHERE mc = %d AND all_cause_mrtl > 0 %s
      GROUP BY %s
      ORDER BY %s",
      group_by_cols,
      config$weight_col,
      config$weight_col,
      config$weight_col,
      mcaggr,
      config$age_filter,
      group_by_cols,
      group_by_cols
    )

    output_path <- private$output_dir(
      paste0(
        "summaries/",
        config$prefix,
        config$suffix,
        "/",
        mcaggr,
        "_",
        config$prefix,
        config$suffix,
        ".",
        ext
      )
    )

    # Execute query and write result with retry logic
    private$execute_db_diskwrite_with_retry(
      duckdb_con,
      query,
      output_path
    )
  }
  NULL
})


# export_hle_summaries ----
Simulation$set("private", "export_hle_summaries", function(
    duckdb_con,
    mcaggr,
    strata_noagegrp,
    ext
) {
  lapply(
    paste0(
      rep(c("hle_1st_cond", "hle_cmsmm1.5"), each = 2),
      "_",
      c("scaled_up", "esp")
    ),
    function(subdir) {
      private$create_new_folder(private$output_dir(paste0(
        "summaries/",
        subdir
      )))
    }
  )

  group_by_cols <- paste(strata_noagegrp, collapse = ", ")

  # Define configurations for HLE calculations
  hle_configs <- list(
    list(
      prefix = "hle_1st_cond",
      condition = "cms_count = 1",
      weight_col = "wt",
      suffix = "_scaled_up"
    ),
    list(
      prefix = "hle_1st_cond",
      condition = "cms_count = 1",
      weight_col = "wt_esp",
      suffix = "_esp"
    ),
    list(
      prefix = "hle_cmsmm1.5",
      condition = '"cmsmm1.5_prvl" = 1',
      weight_col = "wt",
      suffix = "_scaled_up"
    ),
    list(
      prefix = "hle_cmsmm1.5",
      condition = '"cmsmm1.5_prvl" = 1',
      weight_col = "wt_esp",
      suffix = "_esp"
    )
  )

  for (config in hle_configs) {
    query <- sprintf(
      "SELECT %s, SUM(%s) AS popsize, SUM(age * %s) / NULLIF(SUM(%s), 0) AS HLE
         FROM lc_table
         WHERE mc = %d AND %s
         GROUP BY %s
         ORDER BY %s",
      group_by_cols,
      config$weight_col,
      config$weight_col,
      config$weight_col,
      mcaggr,
      config$condition,
      group_by_cols,
      group_by_cols
    )

    output_path <- private$output_dir(
      paste0(
        "summaries/",
        config$prefix,
        config$suffix,
        "/",
        mcaggr,
        "_",
        config$prefix,
        config$suffix,
        ".",
        ext
      )
    )
    private$execute_db_diskwrite_with_retry(
      duckdb_con,
      query,
      output_path
    )

    NULL
  }
})


# export_dis_char_summaries ----
Simulation$set("private", "export_dis_char_summaries", function(
    duckdb_con,
    mcaggr,
    strata_noagegrp,
    ext
) {
  # Create output directories
  lapply(
    paste0(
      rep(c("dis_characteristics"), each = 2),
      "_",
      c("scaled_up", "esp")
    ),
    function(subdir) {
      private$create_new_folder(private$output_dir(paste0(
        "summaries/",
        subdir
      )))
    }
  )

  # Get disease prevalence columns from DuckDB schema
  lc_table_name <- "lc_table"
  all_cols <- dbListFields(duckdb_con, lc_table_name)
  nm <- grep("_prvl$", all_cols, value = TRUE)

  # Construct the SQL query dynamically
  # Ensure agegrp is excluded from grouping to aggregate over it
  select_cols_no_agegrp <- paste(
    setdiff(strata_noagegrp, "agegrp"),
    collapse = ", "
  )

  if (length(nm) > 0) {
    # Memory-efficient approach: Process diseases one at a time
    disease_results_scaled_up <- list()
    disease_results_esp <- list()

    # Get disease names without _prvl suffix
    disease_names <- gsub("_prvl$", "", nm)

    # Process each disease individually to avoid memory issues
    for (i in seq_along(nm)) {
      disease_col <- nm[i]
      disease_name <- disease_names[i]
      quoted_disease_col <- paste0('"', disease_col, '"')

      # Query for current disease (scaled_up version)
      disease_query_scaled_up <- sprintf(
        "SELECT %s, '%s' AS disease,
         SUM(wt) AS cases,
         SUM(CASE WHEN %s = 1 THEN age * wt ELSE 0 END) / NULLIF(SUM(CASE WHEN %s = 1 THEN wt ELSE 0 END), 0) AS mean_age_incd,
         SUM(age * wt) / NULLIF(SUM(wt), 0) AS mean_age_prvl,
         SUM(%s * wt) / NULLIF(SUM(wt), 0) AS mean_duration,
         SUM(cms_score * wt) / NULLIF(SUM(wt), 0) AS mean_cms_score,
         SUM(cms_count * wt) / NULLIF(SUM(wt), 0) AS mean_cms_count
         FROM %s
         WHERE mc = %d AND %s > 0
         GROUP BY %s",
        select_cols_no_agegrp,
        disease_name,
        quoted_disease_col,
        quoted_disease_col,
        quoted_disease_col,
        lc_table_name,
        mcaggr,
        quoted_disease_col,
        select_cols_no_agegrp
      )

      # Execute and store result for scaled_up
      result_scaled_up <- private$query_sql(
        duckdb_con,
        disease_query_scaled_up,
        paste("Disease characteristics for", disease_name, "(scaled_up)")
      )
      if (nrow(result_scaled_up) > 0) {
        disease_results_scaled_up[[i]] <- result_scaled_up
      }

      # Query for current disease (ESP version)
      disease_query_esp <- gsub("wt", "wt_esp", disease_query_scaled_up)

      # Execute and store result for ESP
      result_esp <- private$query_sql(
        duckdb_con,
        disease_query_esp,
        paste("Disease characteristics for", disease_name, "(ESP)")
      )
      if (nrow(result_esp) > 0) {
        disease_results_esp[[i]] <- result_esp
      }
    }

    # Combine and process scaled_up results
    if (length(disease_results_scaled_up) > 0) {
      # Remove NULL entries
      disease_results_scaled_up <- disease_results_scaled_up[
        !vapply(disease_results_scaled_up, is.null, FUN.VALUE = logical(1))
      ]

      if (length(disease_results_scaled_up) > 0) {
        combined_scaled_up <- rbindlist(disease_results_scaled_up)
        setDT(combined_scaled_up)

        # Get strata columns without agegrp for formula
        strata_no_agegrp <- setdiff(strata_noagegrp, "agegrp")

        # Create pivot transformation using data.table dcast
        pivot_result_scaled_up <- dcast(
          combined_scaled_up,
          formula = as.formula(paste(
            paste(strata_no_agegrp, collapse = " + "),
            "~ disease"
          )),
          value.var = c(
            "cases",
            "mean_age_incd",
            "mean_age_prvl",
            "mean_duration",
            "mean_cms_score",
            "mean_cms_count"
          ),
          fill = 0
        )

        # Order by strata columns (excluding agegrp)
        setkeyv(pivot_result_scaled_up, strata_no_agegrp)

        # Write scaled_up version
        output_path <- private$output_dir(
          paste0(
            "summaries/dis_characteristics_scaled_up/",
            mcaggr,
            "_dis_characteristics_scaled_up.",
            ext
          )
        )
        arrow::write_parquet(pivot_result_scaled_up, output_path)
      }
    }

    # Combine and process ESP results
    if (length(disease_results_esp) > 0) {
      # Remove NULL entries
      disease_results_esp <- disease_results_esp[
        !vapply(disease_results_esp, is.null, FUN.VALUE = logical(1))
      ]

      if (length(disease_results_esp) > 0) {
        combined_esp <- rbindlist(disease_results_esp)
        setDT(combined_esp)

        # Get strata columns without agegrp for formula
        strata_no_agegrp <- setdiff(strata_noagegrp, "agegrp")

        # Create pivot transformation using data.table dcast
        pivot_result_esp <- dcast(
          combined_esp,
          formula = as.formula(paste(
            paste(strata_no_agegrp, collapse = " + "),
            "~ disease"
          )),
          value.var = c(
            "cases",
            "mean_age_incd",
            "mean_age_prvl",
            "mean_duration",
            "mean_cms_score",
            "mean_cms_count"
          ),
          fill = 0
        )

        # Order by strata columns (excluding agegrp)
        setkeyv(pivot_result_esp, strata_no_agegrp)

        # Write ESP version
        output_path_esp <- private$output_dir(
          paste0(
            "summaries/dis_characteristics_esp/",
            mcaggr,
            "_dis_characteristics_esp.",
            ext
          )
        )
        arrow::write_parquet(pivot_result_esp, output_path_esp)
      }
    }
  }

  NULL
})


# export_prvl_summaries ----
Simulation$set("private", "export_prvl_summaries", function(
    duckdb_con,
    mcaggr,
    strata,
    ext
) {
  lapply(
    paste0(rep(c("prvl"), each = 2), "_", c("scaled_up", "esp")),
    function(subdir) {
      private$create_new_folder(private$output_dir(paste0(
        "summaries/",
        subdir
      )))
    }
  )

  # Get disease prevalence columns from DuckDB schema
  lc_table_name <- "lc_table"
  all_cols <- dbListFields(duckdb_con, lc_table_name)
  nm_prvl <- grep("_prvl$", all_cols, value = TRUE)

  # Construct the SQL query dynamically
  select_cols <- paste(strata, collapse = ", ")
  sum_cases_cols <- paste(
    sprintf(
      'SUM(CASE WHEN "%s" > 0 THEN wt ELSE 0 END) AS "%s"',
      nm_prvl,
      nm_prvl
    ),
    collapse = ", "
  )

  sql_query <- sprintf(
    "SELECT %s, SUM(wt) AS popsize, %s
            FROM %s
            WHERE mc = %d
            GROUP BY %s
            ORDER BY %s",
    select_cols,
    sum_cases_cols,
    lc_table_name,
    mcaggr,
    select_cols,
    select_cols
  )

  output_path <- private$output_dir(
    paste0("summaries/prvl_scaled_up/", mcaggr, "_prvl_scale_up.", ext)
  )

  private$execute_db_diskwrite_with_retry(
    duckdb_con,
    sql_query,
    output_path
  )

  # esp
  sql_query_esp <- gsub("wt", "wt_esp", sql_query)
  output_path_esp <- private$output_dir(
    paste0("summaries/prvl_esp/", mcaggr, "_prvl_esp.", ext)
  )

  private$execute_db_diskwrite_with_retry(
    duckdb_con,
    sql_query_esp,
    output_path_esp
  )

  NULL
})


# export_incd_summaries ----
Simulation$set("private", "export_incd_summaries", function(
    duckdb_con,
    mcaggr,
    strata,
    ext
) {
  lapply(
    paste0(rep(c("incd"), each = 2), "_", c("scaled_up", "esp")),
    function(subdir) {
      private$create_new_folder(private$output_dir(paste0(
        "summaries/",
        subdir
      )))
    }
  )

  # Get disease prevalence columns from DuckDB schema
  lc_table_name <- "lc_table"
  all_cols <- dbListFields(duckdb_con, lc_table_name)
  nm_prvl <- grep("_prvl$", all_cols, value = TRUE)

  # Construct the SQL query dynamically
  select_cols <- paste(strata, collapse = ", ")
  sum_cases_cols <- paste(
    sprintf(
      'SUM(CASE WHEN "%s" = 1 THEN wt ELSE 0 END) AS "%s"',
      nm_prvl,
      gsub("_prvl$", "_incd", nm_prvl)
    ),
    collapse = ", "
  )

  sql_query <- sprintf(
    "SELECT %s, SUM(wt) AS popsize, %s
            FROM %s
            WHERE mc = %d
            GROUP BY %s
            ORDER BY %s",
    select_cols,
    sum_cases_cols,
    lc_table_name,
    mcaggr,
    select_cols,
    select_cols
  )

  output_path <- private$output_dir(
    paste0("summaries/incd_scaled_up/", mcaggr, "_incd_scale_up.", ext)
  )

  private$execute_db_diskwrite_with_retry(
    duckdb_con,
    sql_query,
    output_path
  )

  # esp
  sql_query_esp <- gsub("wt", "wt_esp", sql_query)
  output_path_esp <- private$output_dir(
    paste0("summaries/incd_esp/", mcaggr, "_incd_esp.", ext)
  )

  private$execute_db_diskwrite_with_retry(
    duckdb_con,
    sql_query_esp,
    output_path_esp
  )

  NULL
})


# export_mrtl_summaries ----
Simulation$set("private", "export_mrtl_summaries", function(
    duckdb_con,
    mcaggr,
    strata,
    ext
) {
  lapply(
    paste0(rep(c("mrtl"), each = 2), "_", c("scaled_up", "esp")),
    function(subdir) {
      private$create_new_folder(private$output_dir(paste0(
        "summaries/",
        subdir
      )))
    }
  )

  lc_table_name <- "lc_table"
  select_cols_mrtl <- paste(strata, collapse = ", ")
  sql_query <- sprintf(
    "SELECT %s,
              SUM(wt) AS popsize,
              SUM(CASE WHEN all_cause_mrtl > 0 THEN wt ELSE 0 END) AS all_cause_mrtl
       FROM %s
       WHERE mc = %d
       GROUP BY %s
       ORDER BY %s",
    select_cols_mrtl,
    lc_table_name,
    mcaggr,
    select_cols_mrtl,
    select_cols_mrtl
  )
  output_path <- private$output_dir(
    paste0("summaries/mrtl_scaled_up/", mcaggr, "_mrtl_scale_up.", ext)
  )

  private$execute_db_diskwrite_with_retry(
    duckdb_con,
    sql_query,
    output_path
  )

  # esp
  sql_query_esp <- gsub("wt", "wt_esp", sql_query)
  output_path_esp <- private$output_dir(
    paste0("summaries/mrtl_esp/", mcaggr, "_mrtl_esp.", ext)
  )

  private$execute_db_diskwrite_with_retry(
    duckdb_con,
    sql_query_esp,
    output_path_esp
  )

  NULL
})


# export_dis_mrtl_summaries ----
Simulation$set("private", "export_dis_mrtl_summaries", function(
    duckdb_con,
    mcaggr,
    strata,
    ext
) {
  lapply(
    paste0(rep(c("dis_mrtl"), each = 2), "_", c("scaled_up", "esp")),
    function(subdir) {
      private$create_new_folder(private$output_dir(paste0(
        "summaries/",
        subdir
      )))
    }
  )

  lc_table_name <- "lc_table"
  quoted_strata <- paste0('"', strata, '"')
  strata_cols_sql <- paste(quoted_strata, collapse = ", ")

  # All pivoted column names, including potentially 'alive_deaths'
  pivoted_death_col_names_all <- paste0(
    names(private$death_codes),
    "_deaths"
  )
  quoted_pivoted_death_col_names_all <- paste0(
    '"',
    pivoted_death_col_names_all,
    '"'
  )

  # Correctly format the IN clause with quoted aliases for PIVOT (includes 'alive_deaths')
  death_codes_pivot_sql <- paste0(
    "'",
    private$death_codes,
    "' AS ",
    quoted_pivoted_death_col_names_all,
    collapse = ", "
  )

  # Pivoted column names EXCLUDING 'alive_deaths' for the final SELECT list
  pivoted_death_col_names_final <- setdiff(
    pivoted_death_col_names_all,
    "alive_deaths"
  )
  quoted_pivoted_death_col_names_final <- paste0(
    '"',
    pivoted_death_col_names_final,
    '"'
  )

  # Create COALESCE expressions for the final SELECT list (excludes 'alive_deaths')
  coalesce_select_sql_final <- paste0(
    "COALESCE(t.",
    quoted_pivoted_death_col_names_final,
    ", 0) AS ",
    quoted_pivoted_death_col_names_final,
    collapse = ", "
  )

  # Sum of ALL coalesced columns for popsize (includes 'alive_deaths')
  death_cols_sum_sql <- paste0(
    "COALESCE(t.",
    quoted_pivoted_death_col_names_all,
    ", 0)",
    collapse = " + "
  )

  # Select strata columns prefixed with t. and quoted
  select_strata_sql <- paste0("t.", quoted_strata, collapse = ", ")

  # Construct the PIVOT SQL query (inner query) - uses all death codes
  sql_query <- sprintf(
    "WITH AggregatedDeaths AS (
     SELECT
       %s,
       all_cause_mrtl,
       SUM(wt) AS deaths
     FROM %s
     WHERE mc = %d
     GROUP BY %s, all_cause_mrtl
     )
     PIVOT AggregatedDeaths
     ON all_cause_mrtl IN (%s)
     USING SUM(deaths)
     GROUP BY %s
     ",
    strata_cols_sql,
    lc_table_name,
    mcaggr,
    strata_cols_sql,
    death_codes_pivot_sql,
    strata_cols_sql
  )

  # Construct the final SQL query with COALESCE and popsize calculation
  sql_query_final <- sprintf(
    "SELECT
     %s,
     %s,
     %s AS popsize
     FROM (%s) AS t
     ORDER BY %s",
    select_strata_sql,
    coalesce_select_sql_final,
    death_cols_sum_sql,
    sql_query,
    strata_cols_sql
  )

  output_path <- private$output_dir(
    paste0(
      "summaries/dis_mrtl_scaled_up/",
      mcaggr,
      "_dis_mrtl_scale_up.",
      ext
    )
  )

  private$execute_db_diskwrite_with_retry(
    duckdb_con,
    sql_query_final,
    output_path
  )

  # esp
  sql_query_esp <- gsub("wt", "wt_esp", sql_query_final)
  output_path_esp <- private$output_dir(
    paste0("summaries/dis_mrtl_esp/", mcaggr, "_dis_mrtl_esp.", ext)
  )

  private$execute_db_diskwrite_with_retry(
    duckdb_con,
    sql_query_esp,
    output_path_esp
  )

  NULL
})


# export_all_cause_mrtl_by_dis_summaries ----
Simulation$set("private", "export_all_cause_mrtl_by_dis_summaries", function(
    duckdb_con,
    mcaggr,
    strata,
    ext
) {
  lapply(
    paste0(
      rep(c("all_cause_mrtl_by_dis"), each = 2),
      "_",
      c("scaled_up", "esp")
    ),
    function(subdir) {
      private$create_new_folder(private$output_dir(paste0(
        "summaries/",
        subdir
      )))
    }
  )

  lc_table_name <- "lc_table"
  schema <- private$query_sql(
    duckdb_con,
    sprintf("DESCRIBE %s;", lc_table_name)
  )
  prvl_cols <- schema$column_name[endsWith(schema$column_name, "_prvl")]
  disease_names <- gsub("_prvl$", "", prvl_cols)

  # Quote strata columns for safety
  quoted_strata <- paste0('"', strata, '"')
  strata_cols_sql <- paste(quoted_strata, collapse = ", ")

  # Construct CASE WHEN statements for each disease's cases and deaths
  case_statements <- vapply(disease_names, function(dis) {
    prvl_col <- paste0('"', dis, '_prvl"')
    sprintf(
      'SUM(CASE WHEN %s > 0 THEN wt ELSE 0 END) AS "cases_%s"',
      prvl_col,
      dis
    )
  }, FUN.VALUE = character(1))

  death_statements <- vapply(disease_names, function(dis) {
    prvl_col <- paste0('"', dis, '_prvl"')
    sprintf(
      'SUM(CASE WHEN %s > 0 AND all_cause_mrtl > 0 THEN wt ELSE 0 END) AS "deaths_%s"',
      prvl_col,
      dis
    )
  }, FUN.VALUE = character(1))

  # Combine all select parts
  select_parts_sql <- paste(
    c(strata_cols_sql, case_statements, death_statements),
    collapse = ",\n  "
  )

  # Construct the full SQL query
  sql_query_dis_char <- sprintf(
    "
  SELECT
    %s
  FROM %s
  WHERE mc = %d
  GROUP BY %s
  ORDER BY %s
  ",
    select_parts_sql,
    lc_table_name,
    mcaggr,
    strata_cols_sql,
    strata_cols_sql
  )

  output_path <- private$output_dir(
    paste0(
      "summaries/all_cause_mrtl_by_dis_scaled_up/",
      mcaggr,
      "_all_cause_mrtl_by_dis_scale_up.",
      ext
    )
  )

  private$execute_db_diskwrite_with_retry(
    duckdb_con,
    sql_query_dis_char,
    output_path
  )

  # esp
  sql_query_esp <- gsub("wt", "wt_esp", sql_query_dis_char)
  output_path_esp <- private$output_dir(
    paste0(
      "summaries/all_cause_mrtl_by_dis_esp/",
      mcaggr,
      "_all_cause_mrtl_by_dis_esp.",
      ext
    )
  )

  private$execute_db_diskwrite_with_retry(
    duckdb_con,
    sql_query_esp,
    output_path_esp
  )

  NULL
})


# export_cms_summaries ----
Simulation$set("private", "export_cms_summaries", function(
    duckdb_con,
    mcaggr,
    strata,
    strata_age,
    ext
) {
  lapply(
    paste0(
      rep(c("cms_score", "cms_score_by_age", "cms_count"), each = 2),
      "_",
      c("scaled_up", "esp")
    ),
    function(subdir) {
      private$create_new_folder(private$output_dir(paste0(
        "summaries/",
        subdir
      )))
    }
  )

  # Define configurations for CMS calculations
  strata_sql <- paste(strata, collapse = ", ")
  strata_age_sql <- paste(strata_age, collapse = ", ")

  cms_configs <- list(
    list(
      metric = "cms_score",
      group_cols_sql = strata_sql,
      group_cols_r = strata,
      weight_col = "wt",
      suffix = "_scaled_up",
      file_group_suffix = ""
    ),
    list(
      metric = "cms_score",
      group_cols_sql = strata_age_sql,
      group_cols_r = strata_age,
      weight_col = "wt",
      suffix = "_by_age_scaled_up",
      file_group_suffix = "_by_age"
    ),
    list(
      metric = "cms_score",
      group_cols_sql = strata_sql,
      group_cols_r = strata,
      weight_col = "wt_esp",
      suffix = "_esp",
      file_group_suffix = ""
    ),
    list(
      metric = "cms_score",
      group_cols_sql = strata_age_sql,
      group_cols_r = strata_age,
      weight_col = "wt_esp",
      suffix = "_by_age_esp",
      file_group_suffix = "_by_age"
    ),
    list(
      metric = "cms_count",
      group_cols_sql = strata_sql,
      group_cols_r = strata,
      weight_col = "wt",
      suffix = "_scaled_up",
      file_group_suffix = ""
    ),
    list(
      metric = "cms_count",
      group_cols_sql = strata_sql,
      group_cols_r = strata,
      weight_col = "wt_esp",
      suffix = "_esp",
      file_group_suffix = ""
    )
  )

  lc_table_name <- "lc_table"

  for (config in cms_configs) {
    query <- sprintf(
      "SELECT %s, SUM(%s) AS popsize, SUM(%s * %s) / SUM(%s) AS %s
       FROM %s
       WHERE mc = %d
       GROUP BY %s
       ORDER BY %s",
      config$group_cols_sql,
      config$weight_col,
      config$metric,
      config$weight_col,
      config$weight_col,
      config$metric,
      lc_table_name,
      mcaggr,
      config$group_cols_sql,
      config$group_cols_sql
    )

    output_path <- private$output_dir(
      paste0(
        "summaries/",
        config$metric,
        config$file_group_suffix,
        ifelse(config$weight_col == "wt_esp", "_esp", "_scaled_up"),
        "/",
        mcaggr,
        "_",
        config$metric,
        config$file_group_suffix,
        ifelse(config$weight_col == "wt_esp", "_esp", "_scaled_up"),
        ".",
        ext
      )
    )

    # Ensure output directory for the specific file exists
    private$create_new_folder(dirname(output_path))

    private$execute_db_diskwrite_with_retry(
      duckdb_con,
      query,
      output_path
    )
  }

  NULL
})


# export_qalys_summaries ----
Simulation$set("private", "export_qalys_summaries", function(
    duckdb_con,
    mcaggr,
    strata,
    ext
) {
  lapply(
    paste0(rep(c("qalys"), each = 2), "_", c("scaled_up", "esp")),
    function(subdir_suffix) {
      private$create_new_folder(private$output_dir(paste0(
        "summaries/",
        subdir_suffix
      )))
    }
  )

  lc_table_name <- "lc_table"

  # Define the name for the temporary view that calc_QALYs will create
  qaly_view_name <- "lc_with_qalys_view"

  # Call calc_QALYs to create/replace the temporary view with EQ5D5L column.
  private$calc_QALYs(
    duckdb_con = duckdb_con,
    mcaggr = mcaggr,
    input_table_name = lc_table_name,
    output_view_name = qaly_view_name
  )

  # Prepare strata columns for SQL query (quoted)
  quoted_strata_cols_sql <- paste(
    sprintf('"%s"', strata),
    collapse = ", "
  )

  # Define QALY metrics for SELECT statement
  qaly_metrics_select_wt <- 'SUM("EQ5D5L" * wt) AS "EQ5D5L"'
  qaly_metrics_select_wt_esp <- 'SUM("EQ5D5L" * wt_esp) AS "EQ5D5L"'

  # --- Scaled-up QALYs ---
  query_scaled_up <- sprintf(
    "SELECT %s, SUM(wt) AS popsize, %s
       FROM %s
       GROUP BY %s
       ORDER BY %s",
    quoted_strata_cols_sql,
    qaly_metrics_select_wt,
    qaly_view_name,
    quoted_strata_cols_sql,
    quoted_strata_cols_sql
  )
  output_path_scaled_up <- private$output_dir(
    paste0("summaries/qalys_scaled_up/", mcaggr, "_qalys_scaled_up.", ext)
  )

  private$execute_db_diskwrite_with_retry(
    duckdb_con,
    query_scaled_up,
    output_path_scaled_up
  )

  # --- ESP QALYs ---
  query_esp <- sprintf(
    "SELECT %s, SUM(wt_esp) AS popsize, %s
       FROM %s
       GROUP BY %s
       ORDER BY %s",
    quoted_strata_cols_sql,
    qaly_metrics_select_wt_esp,
    qaly_view_name,
    quoted_strata_cols_sql,
    quoted_strata_cols_sql
  )
  output_path_esp <- private$output_dir(
    paste0("summaries/qalys_esp/", mcaggr, "_qalys_esp.", ext)
  )

  private$execute_db_diskwrite_with_retry(
    duckdb_con,
    query_esp,
    output_path_esp
  )

  # drop the temporary view if it's no longer needed for this mcaggr
  private$execute_sql(
    duckdb_con,
    sprintf("DROP VIEW IF EXISTS %s;", qaly_view_name)
  )

  NULL
})


# export_costs_summaries ----
Simulation$set("private", "export_costs_summaries", function(
    duckdb_con,
    mcaggr,
    strata,
    ext
) {
  # Create output directories for scaled-up and ESP-weighted summaries
  lapply(
    paste0(rep(c("costs"), each = 2), "_", c("scaled_up", "esp")),
    function(subdir_suffix) {
      private$create_new_folder(private$output_dir(paste0(
        "summaries/",
        subdir_suffix
      )))
    }
  )

  # get scenario names
  scnams <- gsub(
    "^scenario=",
    "",
    list.dirs(
      private$output_dir(file.path("lifecourse", paste0("mc=", mcaggr))),
      full.names = FALSE,
      recursive = FALSE
    )
  )

  lc_table_name <- "lc_table"

  # Define the name for the temporary view that calc_costs will create
  costs_view_name <- "lc_with_costs"

  costs_scn_views <- paste0(costs_view_name, "_", scnams, "_view")

  # Call calc_costs to create/replace the temporary view with cost columns.
  private$calc_costs(
    duckdb_con = duckdb_con,
    mcaggr = mcaggr,
    input_table_name = lc_table_name,
    output_view_name = costs_view_name
  )

  # Prepare strata columns for SQL query (quoted)
  quoted_strata_cols_sql <- paste(
    sprintf('"%s"', strata),
    collapse = ", "
  )

  # Define cost metrics for SELECT statement (matching calc_costs output columns)
  cost_metrics_select_wt_esp <- paste(
    'SUM(healthcare_cost * wt_esp) AS healthcare_cost',
    'SUM(socialcare_cost * wt_esp) AS socialcare_cost',
    'SUM(informalcare_cost * wt_esp) AS informalcare_cost',
    'SUM(economic_output * wt_esp) AS economic_output',
    'SUM(indirect_cost * wt_esp) AS indirect_cost',
    'SUM(total_cost * wt_esp) AS total_cost',
    sep = ", "
  )

  # Memory-efficient approach: Process scenarios one at a time
  esp_results <- list()
  scaled_up_results <- list()

  # Process each scenario individually to avoid loading all scenarios in RAM
  for (i in seq_along(scnams)) {
    scnam <- scnams[i]
    view_name <- costs_scn_views[i]

    # ESP query for current scenario
    query_esp_scenario <- sprintf(
      "SELECT %s, SUM(wt_esp) AS popsize, %s
       FROM %s
       GROUP BY %s",
      quoted_strata_cols_sql,
      cost_metrics_select_wt_esp,
      view_name,
      quoted_strata_cols_sql
    )

    # Execute and store result for ESP
    esp_results[[i]] <- as.data.table(private$query_sql(
      duckdb_con,
      query_esp_scenario,
      paste("ESP costs for scenario", scnam)
    ))

    # Scaled-up query for current scenario (replace wt_esp with wt)
    cost_metrics_select_wt <- gsub(
      "wt_esp",
      "wt",
      cost_metrics_select_wt_esp
    )
    query_scaled_up_scenario <- sprintf(
      "SELECT %s, SUM(wt) AS popsize, %s
       FROM %s
       GROUP BY %s",
      quoted_strata_cols_sql,
      cost_metrics_select_wt,
      view_name,
      quoted_strata_cols_sql
    )

    # Execute and store result for scaled-up
    scaled_up_results[[i]] <- as.data.table(private$query_sql(
      duckdb_con,
      query_scaled_up_scenario,
      paste("Scaled-up costs for scenario", scnam)
    ))
  }

  # Combine all ESP results and write to disk
  if (length(esp_results) > 0) {
    combined_esp <- rbindlist(esp_results, use.names = TRUE, fill = TRUE)
    rm(esp_results)
    # Order by strata columns
    setkeyv(combined_esp, strata)

    output_path_esp <- private$output_dir(
      paste0("summaries/costs_esp/", mcaggr, "_costs_esp.", ext)
    )

    arrow::write_parquet(combined_esp, output_path_esp)
  }

  # Combine all scaled-up results and write to disk
  if (length(scaled_up_results) > 0) {
    combined_scaled_up <- rbindlist(
      scaled_up_results,
      use.names = TRUE,
      fill = TRUE
    )
    rm(scaled_up_results)
    # Order by strata columns
    setkeyv(combined_scaled_up, strata)

    output_path_scaled_up <- private$output_dir(
      paste0(
        "summaries/costs_scaled_up/",
        mcaggr,
        "_costs_scaled_up.",
        ext
      )
    )

    arrow::write_parquet(combined_scaled_up, output_path_scaled_up)
  }

  # Drop the temporary views if they're no longer needed for this mcaggr
  sapply(costs_scn_views, function(view_name) {
    private$execute_sql(
      duckdb_con,
      sprintf("DROP VIEW IF EXISTS %s;", view_name)
    )
  })

  NULL
})
