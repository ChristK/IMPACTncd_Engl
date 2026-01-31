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
# Calculate QALYs based on Shiroiwa 2021 paper
Simulation$set("private", "calc_QALYs", function(
    duckdb_con,
    mcaggr,
    input_table_name,
    output_view_name,
    include_non_significant = FALSE
) {
  eq5d5l_expr <- "
    0.989
    + CASE agegrp
        WHEN '20-24' THEN -0.018 WHEN '25-29' THEN -0.018
        WHEN '30-34' THEN -0.019 WHEN '35-39' THEN -0.019
        WHEN '40-44' THEN -0.018 WHEN '45-49' THEN -0.018
        WHEN '50-54' THEN -0.028 WHEN '55-59' THEN -0.028
        WHEN '60-64' THEN -0.021 WHEN '65-69' THEN -0.021
        WHEN '70-74' THEN -0.057 WHEN '75-79' THEN -0.057
        WHEN '80-84' THEN -0.129 WHEN '85-89' THEN -0.129
        WHEN '90-94' THEN -0.129 WHEN '95-99' THEN -0.129
        ELSE 0.0
      END
    + CASE WHEN sex = 'women' THEN -0.011 ELSE 0.0 END
    + CASE WHEN chd_prvl = 0 THEN 0.0 ELSE -0.073 END
    + CASE WHEN stroke_prvl = 0 THEN 0.0 ELSE -0.265 END
    + CASE WHEN t2dm_prvl = 0 THEN 0.0 ELSE -0.046 END
  "

  hui3_expr <- "
    0.897
    + CASE agegrp
        WHEN '20-24' THEN -0.023 WHEN '25-29' THEN -0.023
        WHEN '30-34' THEN -0.018 WHEN '35-39' THEN -0.018
        WHEN '40-44' THEN -0.004 WHEN '45-49' THEN -0.004
        WHEN '50-54' THEN -0.021 WHEN '55-59' THEN -0.021
        WHEN '60-64' THEN -0.013 WHEN '65-69' THEN -0.013
        WHEN '70-74' THEN -0.042 WHEN '75-79' THEN -0.042
        WHEN '80-84' THEN -0.145 WHEN '85-89' THEN -0.145
        WHEN '90-94' THEN -0.145 WHEN '95-99' THEN -0.145
        ELSE 0.0
      END
    + CASE WHEN sex = 'women' THEN 0.011 ELSE 0.0 END
    + CASE WHEN chd_prvl = 0 THEN 0.0 ELSE -0.081 END
    + CASE WHEN stroke_prvl = 0 THEN 0.0 ELSE -0.293 END
    + CASE WHEN t2dm_prvl = 0 THEN 0.0 ELSE -0.055 END
  "

  if (!include_non_significant) {
    eq5d5l_expr <- paste0(
      eq5d5l_expr,
      " + CASE WHEN htn_prvl = 0 THEN 0.0 ELSE -0.005 END",
      " + CASE WHEN obesity_prvl = 0 THEN 0.0 ELSE -0.034 END"
    )
    hui3_expr <- paste0(
      hui3_expr,
      " + CASE WHEN htn_prvl = 0 THEN 0.0 ELSE -0.006 END",
      " + CASE WHEN obesity_prvl = 0 THEN 0.0 ELSE 0.019 END"
    )
  }

  create_view_sql <- sprintf(
    "
    CREATE OR REPLACE TEMP VIEW %s AS
    SELECT
      *,
      (%s) AS EQ5D5L,
      (%s) AS HUI3
    FROM %s
    WHERE mc = %d;
  ",
    output_view_name,
    eq5d5l_expr,
    hui3_expr,
    input_table_name,
    mcaggr
  )

  private$execute_sql(duckdb_con, create_view_sql)

  NULL
})


# calc_costs ----
Simulation$set("private", "calc_costs", function(
    duckdb_con,
    mcaggr,
    input_table_name,
    output_view_name
) {
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

  # --- Inflation Factors ---
  prod_informal_inflation_factor <- 1.025
  direct_costs_inflation_factor <- 99.6 / 99.7

  # --- Step 1: Create baseline aggregation views using SQL only ---
  base_agg_sql <- "
    CREATE OR REPLACE TEMP VIEW %s AS
    SELECT agegrp, sex, ROUND(SUM(CASE WHEN %s THEN wt ELSE 0 END)) AS V1
    FROM %s WHERE year = %d AND scenario = 'sc0' AND mc = %d GROUP BY agegrp, sex
    "

  # Create all baseline aggregation views
  aggregation_configs <- list(
    list("chd_prvl_2016_agg_view", "chd_dgns > 0", 2016),
    list("chd_prvl_2019_agg_view", "chd_dgns > 0", 2019),
    list("stroke_prvl_2016_agg_view", "stroke_dgns > 0", 2016),
    list("stroke_prvl_2019_agg_view", "stroke_dgns > 0", 2019),
    list("chd_mrtl_2016_initial_view", "all_cause_mrtl = 2", 2016),
    list("stroke_mrtl_2016_initial_view", "all_cause_mrtl = 3", 2016)
  )

  for (config in aggregation_configs) {
    private$execute_sql(
      duckdb_con,
      sprintf(
        base_agg_sql,
        config[[1]],
        config[[2]],
        input_table_name,
        config[[3]],
        mcaggr
      ),
      config[[1]]
    )
  }

  # --- Step 2: Memory-efficient mortality data handling ---
  # Load only essential columns and filter immediately

  # Load observed population (minimal columns)
  obs_pop_2016 <- read_fst(
    "inputs/pop_estimates_lsoa/national_pop_est.fst",
    columns = c("year", "age", "sex", "pops"),
    as.data.table = TRUE
  )[year == 2016L]

  # Process CHD mortality efficiently (parquet dataset)
  chd_ftlt_2016 <- read_parquet_dt(
    "inputs/disease_burden/chd_ftlt",
    cols = c("year", "age", "sex", "mu2"),
    filter = arrow_in("year", 2016L)
  )

  chd_joined <- chd_ftlt_2016[
    obs_pop_2016,
    on = c("age", "sex"),
    nomatch = 0L
  ][, `:=`(
    deaths_calc = mu2 * pops,
    agegrp = fcase(
      age %between% c(30, 34),
      "30-34",
      age %between% c(35, 39),
      "35-39",
      age %between% c(40, 44),
      "40-44",
      age %between% c(45, 49),
      "45-49",
      age %between% c(50, 54),
      "50-54",
      age %between% c(55, 59),
      "55-59",
      age %between% c(60, 64),
      "60-64",
      age %between% c(65, 69),
      "65-69",
      age %between% c(70, 74),
      "70-74",
      age %between% c(75, 79),
      "75-79",
      age %between% c(80, 84),
      "80-84",
      age %between% c(85, 89),
      "85-89",
      age %between% c(90, 94),
      "90-94",
      age >= 95,
      "95-99",
      default = NA_character_
    )
  )][
    !is.na(agegrp),
    .(calculated_deaths = round(sum(deaths_calc))),
    keyby = .(agegrp, sex)
  ]

  # Register minimal table
  dbWriteTable(
    duckdb_con,
    "chd_ftlt_ext_2016_table",
    chd_joined,
    overwrite = TRUE
  )
  rm(chd_joined) # Immediate cleanup

  # Process stroke mortality efficiently (parquet dataset)
  stroke_ftlt_2016 <- read_parquet_dt(
    "inputs/disease_burden/stroke_ftlt",
    cols = c("year", "age", "sex", "mu2"),
    filter = arrow_in("year", 2016L)
  )

  stroke_joined <- stroke_ftlt_2016[
    obs_pop_2016,
    on = c("age", "sex"),
    nomatch = 0L
  ][, `:=`(
    deaths_calc = mu2 * pops,
    agegrp = fcase(
      age %between% c(30, 34),
      "30-34",
      age %between% c(35, 39),
      "35-39",
      age %between% c(40, 44),
      "40-44",
      age %between% c(45, 49),
      "45-49",
      age %between% c(50, 54),
      "50-54",
      age %between% c(55, 59),
      "55-59",
      age %between% c(60, 64),
      "60-64",
      age %between% c(65, 69),
      "65-69",
      age %between% c(70, 74),
      "70-74",
      age %between% c(75, 79),
      "75-79",
      age %between% c(80, 84),
      "80-84",
      age %between% c(85, 89),
      "85-89",
      age %between% c(90, 94),
      "90-94",
      age >= 95,
      "95-99",
      default = NA_character_
    )
  )][
    !is.na(agegrp),
    .(calculated_deaths = round(sum(deaths_calc))),
    keyby = .(agegrp, sex)
  ]

  dbWriteTable(
    duckdb_con,
    "stroke_ftlt_ext_2016_table",
    as.data.frame(stroke_joined),
    overwrite = TRUE
  )
  rm(stroke_joined) # Immediate cleanup

  # Cleanup large intermediate objects immediately
  rm(obs_pop_2016, chd_ftlt_2016, stroke_ftlt_2016)

  # Update mortality views
  private$execute_sql(
    duckdb_con,
    "
    CREATE OR REPLACE TEMP VIEW chd_mrtl_2016_agg_view AS
    SELECT i.agegrp, i.sex,
           CASE WHEN i.V1 = 0 THEN COALESCE(f.calculated_deaths, i.V1) ELSE i.V1 END AS V1
    FROM chd_mrtl_2016_initial_view i
    LEFT JOIN chd_ftlt_ext_2016_table f ON i.agegrp = f.agegrp AND i.sex = f.sex
  ",
    "chd_mrtl_2016_agg_view"
  )

  private$execute_sql(
    duckdb_con,
    "
    CREATE OR REPLACE TEMP VIEW stroke_mrtl_2016_agg_view AS
    SELECT i.agegrp, i.sex,
           CASE WHEN i.V1 = 0 THEN COALESCE(f.calculated_deaths, i.V1) ELSE i.V1 END AS V1
    FROM stroke_mrtl_2016_initial_view i
    LEFT JOIN stroke_ftlt_ext_2016_table f ON i.agegrp = f.agegrp AND i.sex = f.sex
  ",
    "stroke_mrtl_2016_agg_view"
  )

  # --- Step 3: Memory-efficient cost parameter calculation ---
  # Create parameter tables directly in SQL to avoid R object creation

  # Employee parameters - create as SQL view to avoid R data.table
  private$execute_sql(
    duckdb_con,
    "
    CREATE OR REPLACE TEMP VIEW employee_params_view AS
    SELECT agegrp, sex, CAST(employees AS DOUBLE) AS employees FROM (VALUES
      ('30-34', 'men', 1683780), ('35-39', 'men', 1829610), ('40-44', 'men', 2174550), ('45-49', 'men', 2057710),
      ('50-54', 'men', 1702470), ('55-59', 'men', 1425510), ('60-64', 'men', 963430), ('65-69', 'men', 369640),
      ('70-74', 'men', 106850), ('75-79', 'men', 0), ('80-84', 'men', 0), ('85-89', 'men', 0),
      ('90-94', 'men', 0), ('95-99', 'men', 0),
      ('30-34', 'women', 919700), ('35-39', 'women', 894770), ('40-44', 'women', 1049490), ('45-49', 'women', 1037140),
      ('50-54', 'women', 854970), ('55-59', 'women', 685040), ('60-64', 'women', 376370), ('65-69', 'women', 132470),
      ('70-74', 'women', 44050), ('75-79', 'women', 0), ('80-84', 'women', 0), ('85-89', 'women', 0),
      ('90-94', 'women', 0), ('95-99', 'women', 0)
    ) AS t(agegrp, sex, employees)
  ",
    "employee_params_view"
  )

  # CHD informal care parameters
  private$execute_sql(
    duckdb_con,
    "
    CREATE OR REPLACE TEMP VIEW chd_infm_care_view AS
    SELECT agegrp, sex, CAST(infm_care_hrs AS DOUBLE) AS infm_care_hrs FROM (VALUES
      ('30-34', 'men', 0.030), ('35-39', 'men', 0.030), ('40-44', 'men', 0.030), ('45-49', 'men', 0.030),
      ('50-54', 'men', 0.030), ('55-59', 'men', 0.030), ('60-64', 'men', 0.030), ('65-69', 'men', 0.200),
      ('70-74', 'men', 0.200), ('75-79', 'men', 0.200), ('80-84', 'men', 0), ('85-89', 'men', 0),
      ('90-94', 'men', 0), ('95-99', 'men', 0),
      ('30-34', 'women', 0.030), ('35-39', 'women', 0.030), ('40-44', 'women', 0.030), ('45-49', 'women', 0.030),
      ('50-54', 'women', 0.030), ('55-59', 'women', 0.030), ('60-64', 'women', 0.030), ('65-69', 'women', 0.200),
      ('70-74', 'women', 0.200), ('75-79', 'women', 0.200), ('80-84', 'women', 0), ('85-89', 'women', 0),
      ('90-94', 'women', 0), ('95-99', 'women', 0)
    ) AS t(agegrp, sex, infm_care_hrs)
  ",
    "chd_infm_care_view"
  )

  # Stroke informal care parameters
  private$execute_sql(
    duckdb_con,
    "
    CREATE OR REPLACE TEMP VIEW stroke_infm_care_view AS
    SELECT agegrp, sex, CAST(infm_care_hrs AS DOUBLE) AS infm_care_hrs FROM (VALUES
      ('30-34', 'men', 5.20), ('35-39', 'men', 5.20), ('40-44', 'men', 5.20), ('45-49', 'men', 5.20),
      ('50-54', 'men', 5.20), ('55-59', 'men', 5.20), ('60-64', 'men', 5.20), ('65-69', 'men', 5.03),
      ('70-74', 'men', 5.03), ('75-79', 'men', 5.03), ('80-84', 'men', 9.23), ('85-89', 'men', 9.23),
      ('90-94', 'men', 9.23), ('95-99', 'men', 9.23),
      ('30-34', 'women', 5.20), ('35-39', 'women', 5.20), ('40-44', 'women', 5.20), ('45-49', 'women', 5.20),
      ('50-54', 'women', 5.20), ('55-59', 'women', 5.20), ('60-64', 'women', 5.20), ('65-69', 'women', 5.03),
      ('70-74', 'women', 5.03), ('75-79', 'women', 5.03), ('80-84', 'women', 9.23), ('85-89', 'women', 9.23),
      ('90-94', 'women', 9.23), ('95-99', 'women', 9.23)
    ) AS t(agegrp, sex, infm_care_hrs)
  ",
    "stroke_infm_care_view"
  )

  # Direct cost parameters
  private$execute_sql(
    duckdb_con,
    "
    CREATE OR REPLACE TEMP VIEW chd_direct_tcost_view AS
    SELECT agegrp2, sex, CAST(tcost_val AS DOUBLE) AS tcost_val FROM (VALUES
      ('30-44', 'men', 10300000000.0), ('45-64', 'men', 121000000000.0), ('65-69', 'men', 72300000000.0),
      ('70-74', 'men', 90100000000.0), ('75-99', 'men', 197000000000.0),
      ('30-44', 'women', 2500000000.0), ('45-64', 'women', 22500000000.0), ('65-69', 'women', 18900000000.0),
      ('70-74', 'women', 30300000000.0), ('75-99', 'women', 132800000000.0)
    ) AS t(agegrp2, sex, tcost_val)
  ",
    "chd_direct_tcost_view"
  )

  private$execute_sql(
    duckdb_con,
    "
    CREATE OR REPLACE TEMP VIEW stroke_direct_tcost_view AS
    SELECT agegrp2, sex, CAST(tcost_val AS DOUBLE) AS tcost_val FROM (VALUES
      ('30-44', 'men', 24900000000.0), ('45-64', 'men', 186400000000.0), ('65-69', 'men', 109000000000.0),
      ('70-74', 'men', 144100000000.0), ('75-99', 'men', 465600000000.0),
      ('30-44', 'women', 18000000000.0), ('45-64', 'women', 106800000000.0), ('65-69', 'women', 60900000000.0),
      ('70-74', 'women', 95700000000.0), ('75-99', 'women', 606800000000.0)
    ) AS t(agegrp2, sex, tcost_val)
  ",
    "stroke_direct_tcost_view"
  )

  # Optimised cost parameter calculation - single SQL statement approach
  cost_param_sql <- "
    CREATE OR REPLACE TEMP VIEW %s AS
    WITH joined_data AS (
      SELECT p.agegrp, p.sex, p.%s AS factor_col, agg.V1
      FROM %s p
      JOIN %s agg ON p.agegrp = agg.agegrp AND p.sex = agg.sex
    ),
    weighted_data AS (
      SELECT agegrp, sex, (factor_col * V1) AS weighted_factor
      FROM joined_data
    ),
    total_weighted_sum AS (
      SELECT SUM(weighted_factor) AS total_wt_sum FROM weighted_data
    )
    SELECT wd.agegrp, wd.sex,
           (%.2f * wd.weighted_factor / NULLIF(tws.total_wt_sum, 0)) * %.6f / NULLIF(jd.V1, 0) AS cost_param
    FROM weighted_data wd
    CROSS JOIN total_weighted_sum tws
    JOIN joined_data jd ON wd.agegrp = jd.agegrp AND wd.sex = jd.sex
  "

  # Create all cost parameter views efficiently
  cost_configs <- list(
    list(
      "chd_prvl_prdv_cost_param_view",
      "employees",
      "employee_params_view",
      "chd_prvl_2016_agg_view",
      141000000000.00
    ),
    list(
      "stroke_prvl_prdv_cost_param_view",
      "employees",
      "employee_params_view",
      "stroke_prvl_2016_agg_view",
      322000000000.00
    ),
    list(
      "chd_mrtl_prdv_cost_param_view",
      "employees",
      "employee_params_view",
      "chd_mrtl_2016_agg_view",
      2257000000000.00
    ),
    list(
      "stroke_mrtl_prdv_cost_param_view",
      "employees",
      "employee_params_view",
      "stroke_mrtl_2016_agg_view",
      1352000000000.00
    ),
    list(
      "chd_informal_cost_param_view",
      "infm_care_hrs",
      "chd_infm_care_view",
      "chd_prvl_2016_agg_view",
      291000000000.00
    ),
    list(
      "stroke_informal_cost_param_view",
      "infm_care_hrs",
      "stroke_infm_care_view",
      "stroke_prvl_2016_agg_view",
      1651000000000.00
    )
  )

  for (config in cost_configs) {
    private$execute_sql(
      duckdb_con,
      sprintf(
        cost_param_sql,
        config[[1]],
        config[[2]],
        config[[3]],
        config[[4]],
        config[[5]],
        prod_informal_inflation_factor
      ),
      config[[1]]
    )
  }

  # Direct cost parameters with optimised SQL
  direct_cost_sql <- "
    CREATE OR REPLACE TEMP VIEW %s AS
    WITH lc_with_agegrp2 AS (
      SELECT agegrp, sex, V1,
        CASE
          WHEN agegrp IN ('30-34', '35-39', '40-44') THEN '30-44'
          WHEN agegrp IN ('45-49', '50-54', '55-59', '60-64') THEN '45-64'
          WHEN agegrp = '65-69' THEN '65-69'
          WHEN agegrp = '70-74' THEN '70-74'
          ELSE '75-99'
        END AS agegrp2
      FROM %s
    ),
    agg_by_agegrp2 AS (
      SELECT agegrp2, sex, SUM(V1) AS V1_sum
      FROM lc_with_agegrp2
      GROUP BY agegrp2, sex
    )
    SELECT orig.agegrp, orig.sex,
           (tc.tcost_val * %.6f / NULLIF(agg.V1_sum, 0)) AS cost_param
    FROM %s orig
    JOIN lc_with_agegrp2 lwa ON orig.agegrp = lwa.agegrp AND orig.sex = lwa.sex
    JOIN agg_by_agegrp2 agg ON lwa.agegrp2 = agg.agegrp2 AND lwa.sex = agg.sex
    JOIN %s tc ON agg.agegrp2 = tc.agegrp2 AND agg.sex = tc.sex
  "

  private$execute_sql(
    duckdb_con,
    sprintf(
      direct_cost_sql,
      "chd_direct_cost_param_view",
      "chd_prvl_2019_agg_view",
      direct_costs_inflation_factor,
      "chd_prvl_2019_agg_view",
      "chd_direct_tcost_view"
    ),
    "chd_direct_cost_param_view"
  )

  private$execute_sql(
    duckdb_con,
    sprintf(
      direct_cost_sql,
      "stroke_direct_cost_param_view",
      "stroke_prvl_2019_agg_view",
      direct_costs_inflation_factor,
      "stroke_prvl_2019_agg_view",
      "stroke_direct_tcost_view"
    ),
    "stroke_direct_cost_param_view"
  )

  # --- Step 4: Create Final Output View with All Cost Columns ---
  # One view per scenario named paste0(output_view_name, "_", scnams, "_view")
  final_view_creation_sql <- sprintf(
    "
    CREATE OR REPLACE TEMP VIEW %s AS
    WITH base_filtered AS (
      SELECT mc, scenario, year, agegrp, sex, dimd, chd_dgns, all_cause_mrtl, stroke_dgns, wt, wt_esp
      FROM %s
      WHERE mc = %d AND scenario = %s
      ),
      chd_costs AS (
        SELECT agegrp, sex,
          COALESCE(cppc.cost_param, 0) AS chd_prvl_prdv,
          COALESCE(cpmc.cost_param, 0) AS chd_mrtl_prdv,
          COALESCE(cic.cost_param, 0) AS chd_informal,
          COALESCE(cdc.cost_param, 0) AS chd_direct
        FROM chd_prvl_prdv_cost_param_view cppc
        LEFT JOIN chd_mrtl_prdv_cost_param_view cpmc USING(agegrp, sex)
        LEFT JOIN chd_informal_cost_param_view cic USING(agegrp, sex)
        LEFT JOIN chd_direct_cost_param_view cdc USING(agegrp, sex)
      ),
      stroke_costs AS (
        SELECT agegrp, sex,
          COALESCE(sppc.cost_param, 0) AS stroke_prvl_prdv,
          COALESCE(spmc.cost_param, 0) AS stroke_mrtl_prdv,
          COALESCE(sic.cost_param, 0) AS stroke_informal,
          COALESCE(sdc.cost_param, 0) AS stroke_direct
        FROM stroke_prvl_prdv_cost_param_view sppc
        LEFT JOIN stroke_mrtl_prdv_cost_param_view spmc USING(agegrp, sex)
        LEFT JOIN stroke_informal_cost_param_view sic USING(agegrp, sex)
        LEFT JOIN stroke_direct_cost_param_view sdc USING(agegrp, sex)
      ),
      basic_costs AS (
        SELECT
          m.mc, m.scenario, m.year, m.agegrp, m.sex, m.dimd,
          m.wt, m.wt_esp,

          -- CHD basic cost components
          CASE WHEN m.chd_dgns > 0 THEN cc.chd_prvl_prdv ELSE 0 END AS chd_prvl_prdv_costs,
          CASE WHEN m.all_cause_mrtl = 2 THEN cc.chd_mrtl_prdv ELSE 0 END AS chd_mrtl_prdv_costs,
          CASE WHEN m.chd_dgns > 0 THEN cc.chd_informal ELSE 0 END AS chd_informal_costs,
          CASE WHEN m.chd_dgns > 0 THEN cc.chd_direct ELSE 0 END AS chd_direct_costs,

          -- Stroke basic cost components
          CASE WHEN m.stroke_dgns > 0 THEN sc.stroke_prvl_prdv ELSE 0 END AS stroke_prvl_prdv_costs,
          CASE WHEN m.all_cause_mrtl = 3 THEN sc.stroke_mrtl_prdv ELSE 0 END AS stroke_mrtl_prdv_costs,
          CASE WHEN m.stroke_dgns > 0 THEN sc.stroke_informal ELSE 0 END AS stroke_informal_costs,
          CASE WHEN m.stroke_dgns > 0 THEN sc.stroke_direct ELSE 0 END AS stroke_direct_costs

        FROM base_filtered m
        LEFT JOIN chd_costs cc ON m.agegrp = cc.agegrp AND m.sex = cc.sex
        LEFT JOIN stroke_costs sc ON m.agegrp = sc.agegrp AND m.sex = sc.sex
      )
      SELECT
        mc, scenario, year, agegrp, sex, dimd, wt, wt_esp,

        -- Basic cost components (already calculated)
        chd_prvl_prdv_costs,
        chd_mrtl_prdv_costs,
        chd_informal_costs,
        chd_direct_costs,
        stroke_prvl_prdv_costs,
        stroke_mrtl_prdv_costs,
        stroke_informal_costs,
        stroke_direct_costs,

        -- Aggregated productivity costs
        (chd_prvl_prdv_costs + chd_mrtl_prdv_costs) AS chd_productivity_costs,
        (stroke_prvl_prdv_costs + stroke_mrtl_prdv_costs) AS stroke_productivity_costs,
        (chd_prvl_prdv_costs + chd_mrtl_prdv_costs + stroke_prvl_prdv_costs + stroke_mrtl_prdv_costs) AS cvd_productivity_costs,

        -- Aggregated informal costs
        (chd_informal_costs + stroke_informal_costs) AS cvd_informal_costs,

        -- Aggregated direct costs
        (chd_direct_costs + stroke_direct_costs) AS cvd_direct_costs,

        -- Aggregated indirect costs (productivity + informal)
        (chd_prvl_prdv_costs + chd_mrtl_prdv_costs + chd_informal_costs) AS chd_indirect_costs,
        (stroke_prvl_prdv_costs + stroke_mrtl_prdv_costs + stroke_informal_costs) AS stroke_indirect_costs,
        (chd_prvl_prdv_costs + chd_mrtl_prdv_costs + chd_informal_costs + stroke_prvl_prdv_costs + stroke_mrtl_prdv_costs + stroke_informal_costs) AS cvd_indirect_costs,

        -- Total costs (indirect + direct)
        (chd_prvl_prdv_costs + chd_mrtl_prdv_costs + chd_informal_costs + chd_direct_costs) AS chd_total_costs,
        (stroke_prvl_prdv_costs + stroke_mrtl_prdv_costs + stroke_informal_costs + stroke_direct_costs) AS stroke_total_costs,
        (chd_prvl_prdv_costs + chd_mrtl_prdv_costs + chd_informal_costs + chd_direct_costs + stroke_prvl_prdv_costs + stroke_mrtl_prdv_costs + stroke_informal_costs + stroke_direct_costs) AS cvd_total_costs

      FROM basic_costs;
  ",
    paste0(output_view_name, "_", scnams, "_view"),
    input_table_name,
    mcaggr,
    paste0("'", scnams, "'")
  )

  sapply(final_view_creation_sql, function(sql) {
    private$execute_sql(
      duckdb_con,
      sql,
      paste0("Final cost views per scenario:", output_view_name)
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

  # Call calc_QALYs to create/replace the temporary view with EQ5D5L and HUI3 columns.
  private$calc_QALYs(
    duckdb_con = duckdb_con,
    mcaggr = mcaggr,
    input_table_name = lc_table_name,
    output_view_name = qaly_view_name,
    include_non_significant = FALSE
  )

  # Prepare strata columns for SQL query (quoted)
  quoted_strata_cols_sql <- paste(
    sprintf('"%s"', strata),
    collapse = ", "
  )

  # Define QALY metrics for SELECT statement
  qaly_metrics_select_wt <- 'SUM("EQ5D5L" * wt) AS "EQ5D5L", SUM("HUI3" * wt) AS "HUI3"'
  qaly_metrics_select_wt_esp <- 'SUM("EQ5D5L" * wt_esp) AS "EQ5D5L", SUM("HUI3" * wt_esp) AS "HUI3"'

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

  # Define cost metrics for SELECT statement
  cost_metrics_select_wt_esp <- paste(
    'SUM(chd_direct_costs * wt_esp) AS chd_direct_costs',
    'SUM(stroke_direct_costs * wt_esp) AS stroke_direct_costs',
    'SUM(cvd_direct_costs * wt_esp) AS cvd_direct_costs',
    'SUM(chd_productivity_costs * wt_esp) AS chd_productivity_costs',
    'SUM(stroke_productivity_costs * wt_esp) AS stroke_productivity_costs',
    'SUM(cvd_productivity_costs * wt_esp) AS cvd_productivity_costs',
    'SUM(chd_informal_costs * wt_esp) AS chd_informal_costs',
    'SUM(stroke_informal_costs * wt_esp) AS stroke_informal_costs',
    'SUM(cvd_informal_costs * wt_esp) AS "cvd_informal_costs"',
    'SUM(chd_indirect_costs * wt_esp) AS chd_indirect_costs',
    'SUM(stroke_indirect_costs * wt_esp) AS stroke_indirect_costs',
    'SUM(cvd_indirect_costs * wt_esp) AS cvd_indirect_costs',
    'SUM(chd_total_costs * wt_esp) AS chd_total_costs',
    'SUM(stroke_total_costs * wt_esp) AS stroke_total_costs',
    'SUM(cvd_total_costs * wt_esp) AS cvd_total_costs',
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
