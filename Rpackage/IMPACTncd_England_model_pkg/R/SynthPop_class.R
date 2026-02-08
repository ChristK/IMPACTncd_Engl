## IMPACTncdEngland is an implementation of the IMPACTncd framework, developed by
## Chris Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz.
##
## Copyright (C) 2018-2026 University of Liverpool, Chris Kypridemos
##
## IMPACTncdEngland is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by the Free
## Software Foundation; either version 3 of the License, or (at your option) any
## later version. This program is distributed in the hope that it will be
## useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
## Public License for more details. You should have received a copy of the GNU
## General Public License along with this program; if not, see
## <http://www.gnu.org/licenses/> or write to the Free Software Foundation,
## Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.



# From
# https://stackoverflow.com/questions/33424233/how-do-i-tell-an-r6-class-what-to-do-with-square-brackets
# Allows data.table syntax to the R6class object directly. Assumes it has a
# field 'pop' that is a data.table
#' @export
`[.SynthPop` <- function(x, ...) x$pop[...]

#' @title SynthPop Class
#' @name SynthPop
#'
#' @description
#' An R6 class representing a synthetic population in the IMPACTncd simulation framework.
#' This class encapsulates synthetic simulants, manages file-based persistence,
#' simulates exposures, and adjusts synthetic population weights to match national or local projections.
#'
#' @details
#' The `SynthPop` class stores and operates on a synthetic population used for simulation of public health interventions.
#' It ensures consistency between simulations by persisting data to disk and using consistent checksums and filenames.
#' The synthetic population is composed of individuals (simulants) with simulated life-course data and exposures,
#' generated using a design object and stratified sampling.
#'
#' Key features include:
#' \itemize{
#'   \item Generation of synthetic demographic and exposure data.
#'   \item Calibration of weights to match national/local population projections.
#'   \item Modular and parallel generation of population files.
#'   \item Tools for risk storage, deletion, and validation of synthpop files.
#' }
#'
#' @section Public Fields:
#' \describe{
#'   \item{\code{mc}}{Monte Carlo iteration number.}
#'   \item{\code{mc_aggr}}{Monte Carlo aggregation group.}
#'   \item{\code{metadata}}{YAML metadata associated with the synthpop.}
#'   \item{\code{pop}}{A \code{data.table} containing the synthetic population life-course data.}
#' }
#'
#' @section Public Methods:
#' \describe{
#'   \item{\code{initialize(mc_, design_)}}{Initializes or loads the synthpop.}
#'   \item{\code{update_design(design_)}}{Updates the internal design object.}
#'   \item{\code{update_pop_weights(scenario_nam)}}{Adjusts weights for baseline or scenario populations.}
#'   \item{\code{delete_synthpop(mc_)}}{Deletes synthpop files from disk.}
#'   \item{\code{delete_incomplete_synthpop()}}{Removes orphaned or incomplete synthpop files.}
#'   \item{\code{check_integridy()}}{Checks and optionally removes malformed synthpop files.}
#'   \item{\code{count_synthpop()}}{Reports file counts and size.}
#'   \item{\code{get_checksum()}}{Returns the checksum of the current synthpop.}
#'   \item{\code{get_filename()}}{Prints synthpop filenames.}
#'   \item{\code{get_design()}}{Returns the internal design object.}
#'   \item{\code{get_dir()}}{Prints the synthpop storage directory.}
#'   \item{\code{gen_synthpop_demog()}}{Generates baseline sociodemographic structure.}
#'   \item{\code{write_synthpop(mc_)}}{Generates and writes synthpop files in parallel.}
#'   \item{\code{get_risks(disease_nam)}}{Returns disease-specific risk tables.}
#'   \item{\code{store_risks(disease_nam)}}{Stores and removes risk-related columns from `pop`.}
#'   \item{\code{print()}}{Prints summary of synthpop metadata.}
#' }
#'
#' @section Private Methods:
#' Many helper methods support internal logic for data integrity, file generation, risk storage,
#' LSOA and LAD resolution, checksum creation, and exposure generation.
#'
#' @seealso \code{\link[R6]{R6Class}}, \code{\link[yaml]{read_yaml}}, \code{\link[fst]{read_fst}}, \code{\link[data.table]{data.table}}#'
#'
#' @export
SynthPop <-
  R6::R6Class(
    classname = "SynthPop",

    # public ------------------------------------------------------------------
    public = list(
      #' @field mc The Monte Carlo iteration of the synthetic population. Every
      #'   integer generates a unique synthetic population.
      mc = NA,

      #' @field mc_aggr The Monte Carlo iteration of the synthetic population to
      #'   be used when multiple synthetic populations getting aggregated. It
      #'   ensures correct seeds for the RNGs during the simulation for the RRs
      #'   and the lags.
      mc_aggr = NA,

      #' @field metadata Metadata of the synthpop.
      metadata = NA,

      #' @field pop The data.table that contains the life-course of simulants.
      #'   If the file exists, it is loaded from disk. If it doesn't, it is
      #'   first generated, then saved to disk, and then loaded from disk.
      pop = NA,

      #' Initialize SynthPop Object
      #' @description
      #' Constructs a `SynthPop` object for a specified Monte Carlo iteration using a given `Design` object.
      #' If the corresponding synthpop files already exist on disk, they are loaded.
      #' If not, the synthetic population is generated from scratch, saved, and then loaded.
      #'
      #' This method ensures reproducibility and safe concurrent execution via deterministic filenames and file checks.
      #'
      #' @param mc_ Integer. The Monte Carlo iteration number (≥ 0). If set to 0, initializes an empty synthpop object.
      #' @param design_ A valid \code{\link{Design}} object containing simulation parameters.
      #'
      #' @return An invisible `SynthPop` object.
      #'
      #' @examples
      #' design <- Design$new("inputs/sim_design.yaml")
      #' sp <- SynthPop$new(1, design)
      #' sp$print()
      # initialize ----
      initialize = function(mc_, design_) {
        stopifnot(length(mc_) == 1L, is.numeric(mc_), ceiling(mc_) >= 0L)
        stopifnot("Design" %in% class(design_))

        mc_ <- as.integer(ceiling(mc_))
        # Create synthpop_dir if it doesn't exists
        # NOTE code below is duplicated in Simulation class. This is intentional
        if (!dir.exists(design_$sim_prm$synthpop_dir)) {
          dir.create(design_$sim_prm$synthpop_dir, recursive = TRUE)
          message(paste0(
            "Folder ",
            design_$sim_prm$synthpop_dir,
            " was created"
          ))
        }

        # get unique lsoas
        # lsoas <- private$get_unique_LSOAs(design_)

        private$checksum <- private$gen_checksum(design_)

        self$mc <- mc_
        self$mc_aggr <-
          as.integer(ceiling(mc_ / design_$sim_prm$num_chunks))

        private$design <- design_
        private$synthpop_dir <- design_$sim_prm$synthpop_dir
        private$ladcontribution <- private$LAD_contribution()

        if (mc_ > 0) {
          private$filename <- private$gen_synthpop_filename(
            mc_,
            private$checksum,
            design_
          )
          # logic for the synthpop load
          files_exist <- vapply(
            private$filename,
            file.exists,
            FUN.VALUE = logical(1)
          )
          if (all(!files_exist)) {
            # No files exist. Create the synthpop and store the file on disk (no
            # parallelism)
            private$gen_synthpop(mc_, private$filename, design_)
          } else if (
            file.exists(private$filename$metafile) &&
              !all(files_exist)
          ) {
            # Metafile exists but not all three files. It means that most likely
            # a generate_synthpop() is still running. So the function waits
            # until the file is created before it proceeds to load it. Note that
            # if this is not the case then the loop is infinite!!!
            while (
              !all(vapply(
                private$filename,
                file.exists,
                FUN.VALUE = logical(1)
              ))
            ) {
              # Sys.sleep(5)
              # if (design_$sim_prm$logs) {
              #   message(
              #     "Metafile exists without a synthpop file. Check for synthpop after 5 sec."
              #   )}
              self$delete_incomplete_synthpop()
              private$gen_synthpop(mc_, private$filename, design_)
            }
            # Ensure the file write is complete (size stable)
            if (design_$sim_prm$logs) {
              message("Synthpop file found.")
            }

            sz1 <- file.size(private$filename$synthpop)
            Sys.sleep(3)
            sz2 <- file.size(private$filename$synthpop)
            while (sz1 != sz2) {
              if (design_$sim_prm$logs) {
                message("Synthpop file size increases.")
              }
              sz1 <- file.size(private$filename$synthpop)
              Sys.sleep(3)
              sz2 <- file.size(private$filename$synthpop)
            }
            if (design_$sim_prm$logs) {
              message("Synthpop file stabilised.")
            }
          } else if (
            !file.exists(private$filename$metafile) &&
              !all(files_exist)
          ) {
            # Metafile doesn't exist but some other files exist. In this case
            # delete everything and start from scratch
            self$delete_incomplete_synthpop()
            private$gen_synthpop(mc_, private$filename, design_)
          }
          # No need to provision for case when all file present. The following
          # lines handle this case anyway

          self$pop <- private$get_synthpop()
          self$metadata <- yaml::read_yaml(private$filename$metafile)

          if (design_$sim_prm$logs) self$print()
        }
        invisible(self)
      },

      # update_design ----
      #' @description
      #' Updates the Design object that is stored in the SynthPop object.
      #' @param design_ A design object with the simulation parameters.
      #' @return The invisible self for chaining.

      update_design = function(design_) {
        if (!inherits(design_, "Design")) {
          stop("Argument design_ needs to be a Design object.")
        }

        private$design <- design_
        invisible(self)
      },

      # update_pop_weights ----
      #' @description
      #' Updates the wt to reflect population change.
      #' @param scenario_nam If "sc0" (the baseline scenario) update weights to
      #'   scale to ONS projections. Else copy weights from the baseline
      #'   scenario to the current scenario for the common person-years.
      #' @return The invisible self for chaining.
      update_pop_weights = function(scenario_nam = "sc0") {
        strata <- c("year", "age", "sex")
        if (private$design$sim_prm$calibrate_to_pop_projections_by_LAD) {
          strata <- c(strata, "LAD17CD")
        }

        if (scenario_nam == "sc0") {
          # baseline # & !"wt" %in% names(self$pop)
          tt <- private$get_pop_size(
            private$design$sim_prm$calibrate_to_pop_projections_by_LAD
          )
          # fix for ICBs that may have only a part of LADs. Note that the fix
          # assumes that the portion of the LAD that is part of the ICB has the
          # same age/sex distribution as the rest of the LAD.
          if (private$design$sim_prm$calibrate_to_pop_projections_by_LAD) {
            tt[private$ladcontribution, on = "LAD17CD", pops := pops * spTOons]
          }
          self$pop[all_cause_mrtl >= 0, wt := .N, by = strata] # Long deads are NAs and those died in the year need to be counted
          absorb_dt(self$pop, tt)
          self$pop[,
            wt := pops / (wt * private$design$sim_prm$num_chunks)
          ]
          self$pop[is.na(all_cause_mrtl), wt := 0]
          self$pop[, pops := NULL]

          # NOTE that the sum of the weights in the above may be smaller that
          # the total population/private$design$sim_prm$num_chunks
          # if n is small because not all year/age/sex/(LAD17CD) are represented
          # in the sample. The recent change to stratified sampling by age sex
          # should reduce the chance of this issue but is almost certainly
          # happening when design$sim_prm$calibrate_to_pop_projections_by_LAD = TRUE.
          # cbind(
          #   tt[age >= private$design$sim_prm$ageL & year >= private$design$sim_prm$init_year, sum(pops), keyby = year],
          #   self$pop[age >= private$design$sim_prm$ageL & year >= private$design$sim_prm$init_year, sum(wt) * 2, keyby = year]
          # )

          # Fix for missing pop segments when locality is England. Calibrate to
          # latest ONS projections if possible.
          if (private$design$sim_prm$locality == "England") {
            refpop <- private$get_pop_size(FALSE)[
              age >= private$design$sim_prm$ageL &
                year >= private$design$sim_prm$init_year
            ] # Use national projections
            ttt <- self$pop[
              age >= private$design$sim_prm$ageL &
                year >= private$design$sim_prm$init_year,
              sum(wt, na.rm = TRUE),
              keyby = .(year, age, sex)
            ]
            refpop[
              ttt,
              on = c("year", "age", "sex"),
              correction := pops /
                (i.V1 * private$design$sim_prm$num_chunks)
            ]
            self$pop[
              refpop,
              on = c("year", "age", "sex"),
              wt := i.correction * wt
            ]
          } else if (
            # For LADs and regions use the LAD projections
            # TODO cover cases when locality contains both LADS and Regions
            all(
              private$design$sim_prm$locality %in%
                unique(
                  read_fst(
                    "./inputs/pop_estimates_lsoa/lsoa_to_locality_indx.fst",
                    columns = "LAD17NM",
                    as.data.table = TRUE
                  )$LAD17NM
                )
            ) ||
              all(
                private$design$sim_prm$locality %in%
                  unique(
                    read_fst(
                      "./inputs/pop_estimates_lsoa/lsoa_to_locality_indx.fst",
                      columns = "RGN11NM",
                      as.data.table = TRUE
                    )$RGN11NM
                  )
              )
          ) {
            refpop <- private$get_pop_size(TRUE)[
              age >= private$design$sim_prm$ageL &
                year >= private$design$sim_prm$init_year
            ] # Use LAD projections
            refpop <- refpop[, .(pops = sum(pops)), keyby = .(year, age, sex)]
            ttt <- self$pop[
              age >= private$design$sim_prm$ageL &
                year >= private$design$sim_prm$init_year,
              sum(wt, na.rm = TRUE),
              keyby = .(year, age, sex)
            ]
            refpop[
              ttt,
              on = c("year", "age", "sex"),
              correction := pops /
                (i.V1 * private$design$sim_prm$num_chunks)
            ]
            self$pop[
              refpop,
              on = c("year", "age", "sex"),
              wt := i.correction * wt
            ]
          } else {
            # for ICBs
            refpop <- private$get_pop_size(TRUE)[
              age >= private$design$sim_prm$ageL &
                year >= private$design$sim_prm$init_year
            ] # Use LAD projections
            refpop[
              private$ladcontribution,
              on = "LAD17CD",
              pops := pops * spTOons
            ]
            refpop <- refpop[, .(pops = sum(pops)), keyby = .(year, age, sex)]
            ttt <- self$pop[
              age >= private$design$sim_prm$ageL &
                year >= private$design$sim_prm$init_year,
              sum(wt, na.rm = TRUE),
              keyby = .(year, age, sex)
            ]
            refpop[
              ttt,
              on = c("year", "age", "sex"),
              correction := pops /
                (i.V1 * private$design$sim_prm$num_chunks)
            ]
            self$pop[
              refpop,
              on = c("year", "age", "sex"),
              wt := i.correction * wt
            ]
          }

          # NOTE no correction is needed for ICBs
        } else if (scenario_nam != "sc0") {
          # & !"wt" %in% names(self$pop)
          # For policy scenarios - read from partitioned parquet format
          fnam <- file.path(
            private$design$sim_prm$output_dir,
            "lifecourse",
            paste0("mc=", self$mc_aggr),
            "scenario=sc0"
          )
          t0 <- read_parquet_dt(
            fnam,
            cols = c("pid", "year", "wt", "wt_esp")
          )
          # For some reason pid and year get read incorrectly as character sometimes
          # TODO check if this is still necessary
          if (!is.integer(t0$pid)) {
            t0[, pid := as.integer(pid)]
          }
          if (!is.integer(self$pop$pid)) {
            self$pop[, pid := as.integer(pid)]
          }
          if (!is.integer(t0$year)) {
            t0[, year := as.integer(year)]
          }
          if (!is.integer(self$pop$year)) {
            self$pop[, year := as.integer(year)]
          }

          setkeyv(t0, c("pid", "year"))
          self$pop[
            t0,
            on = c("pid", "year"),
            `:=`(wt = i.wt, wt_esp = i.wt_esp)
          ]

          # New way of calculating policy scenario population weights
          setkeyv(self$pop, c("pid", "year"))
          setnafill(self$pop, type = "locf", cols = c("wt", "wt_esp"))
          self$pop[is.na(all_cause_mrtl), `:=`(wt = 0, wt_esp = 0)]
        } else {
          stop(
            "The baseline scenario need to be named 'sc0' and simulated first, before any policy scenarios."
          ) # TODO more informative message
        }

        invisible(self)
      },

      # delete_synthpop ----
      #' @description
      #' Delete (all) synthpop files in the synthpop directory.
      #' @param mc_ If `mc_ = NULL`, delete all files in the synthpop directory.
      #'   If `mc_` is an integer vector delete the specific synthpop files
      #'   including the metadata and index files.
      #' @param check_checksum If  `TRUE` only delete files with the same
      #'   checksum as the synthpop. Only relevant when `mc_ = NULL`.
      #' @param invert If `TRUE` (default is `FALSE`) keeps files with the same
      #'   checksum as the synthpop and deletes all other synthpops. Only
      #'   relevant when `mc_ = NULL` and `check_checksum = TRUE`.
      #' @return The invisible `SynthPop` object.
      delete_synthpop = function(mc_, check_checksum = TRUE, invert = FALSE) {
        if (missing(mc_)) {
          stop("Use mc_ = NULL if you want to delete all synthpop files.")
        }
        if (is.null(mc_)) {
          if (check_checksum) {
            fl <- list.files(
              private$synthpop_dir,
              pattern = paste0("^synthpop_", private$checksum),
              full.names = TRUE,
              recursive = TRUE
            )
            if (invert) {
              fl2 <- list.files(
                private$synthpop_dir,
                pattern = "^synthpop_",
                full.names = TRUE,
                recursive = TRUE
              )
              fl <- setdiff(fl2, fl)
            }
          } else {
            fl <- list.files(
              private$synthpop_dir,
              pattern = "^synthpop_",
              full.names = TRUE,
              recursive = TRUE
            )
          }
          file.remove(fl)
        } else if (
          length(mc_) == 1L &&
            is.numeric(mc_) &&
            ceiling(mc_) > 0L
        ) {
          fl <- unlist(
            private$gen_synthpop_filename(mc_, private$checksum, private$design)
          )
          file.remove(fl)
        } else if (
          length(mc_) > 1L &&
            all(is.numeric(mc_)) &&
            all(ceiling(mc_) > 0L)
        ) {
          fl <-
            lapply(
              mc_,
              private$gen_synthpop_filename,
              private$checksum,
              private$design
            )
          fl <- unlist(fl)
          file.remove(fl)
        } else {
          message("mc_ need to be NULL or numeric. Nothing was deleted.")
        }

        return(invisible(self))
      },

      # delete_incomplete_synthpop ----
      #' @description
      #' Check that every synthpop file has a metafile and an index file. Delete
      #' any orphan files.
      #' @param check_checksum If  `TRUE` only delete incomplete group files
      #'   with the same checksum as the synthpop.
      #' @return The invisible `SynthPop` object.
      delete_incomplete_synthpop = function(check_checksum = TRUE) {
        if (check_checksum) {
          f1 <- paste0("^synthpop_", private$checksum, ".*\\.fst$")
          f2 <- paste0("^synthpop_", private$checksum, ".*_meta\\.yaml$")
        } else {
          f1 <- "^synthpop_.*\\.fst$"
          f2 <- "^synthpop_.*_meta\\.yaml$"
        }

        files <-
          list.files(private$synthpop_dir, f1)
        # remove indx files
        files <- sub("\\.fst$", "", files)
        metafiles <-
          list.files(private$synthpop_dir, f2)
        metafiles <- sub("_meta\\.yaml$", "", metafiles)

        to_remove <- setdiff(metafiles, files)
        if (length(to_remove) > 0) {
          to_remove <- paste0(to_remove, "_meta.yaml")
          file.remove(file.path(private$synthpop_dir, to_remove))
        }

        to_remove <- setdiff(files, metafiles)
        if (length(to_remove) > 0) {
          to_remove2 <- paste0(to_remove, ".fst")
          file.remove(file.path(private$synthpop_dir, to_remove2))
        }

        return(invisible(self))
      },

      # check_integridy ----
      #' @description
      #' Check the integrity of (and optionally delete) .fst files by checking
      #' their metadata are readable.
      #' @param remove_malformed If `TRUE`, delete all malformed .fst files and
      #'   their associated files.
      #' @param check_checksum If  `TRUE` only check files with the same
      #'   checksum as the synthpop.
      #' @return The invisible `SynthPop` object.
      check_integrity = function(
        remove_malformed = FALSE,
        check_checksum = TRUE
      ) {
        if (check_checksum) {
          pat <- paste0("^synthpop_", private$checksum, ".*\\.fst$")
        } else {
          pat <- "^synthpop_.*\\.fst$"
        }

        files <-
          list.files(private$synthpop_dir, pat, full.names = TRUE)
        if (length(files) > 0L) {
          malformed <- vapply(
            files,
            function(x) {
              out <- try(metadata_fst(x), silent = TRUE)
              inherits(out, "try-error")
            },
            FUN.VALUE = logical(1),
            USE.NAMES = FALSE
          )

          des <- sum(malformed)

          if (remove_malformed) {
            if (des == 0L) {
              message(paste0(des, " malformed fst file(s)"))
            } else {
              # des != 0L
              message(paste0(des, " malformed fst file(s)..."))
              to_remove <- files[malformed]

              # then remove other files
              to_remove <- gsub(".fst$", "", to_remove)

              # _meta.yaml
              tr <- paste0(to_remove, "_meta.yaml")
              file.remove(tr[file.exists(tr)])
              # .fst
              tr <- paste0(to_remove, ".fst")
              file.remove(tr[file.exists(tr)])

              message("...now deleted!")
            }
          } else {
            # remove_malformed = FALSE
            message(paste0(des, " malformed fst file(s)"))
          }
        } else {
          # if length(files) == 0
          message("no .fst files found.")
        }
        return(invisible(self))
      },

      # count_synthpop ----
      #' @description
      #' Count the synthpop files in a directory. It includes files without
      #' metafiles and index files.
      #' @return The invisible `SynthPop` object.
      count_synthpop = function() {
        out <- list()
        # folder size
        files <-
          list.files(private$synthpop_dir, full.names = TRUE)
        if (length(files) > 0L) {
          vect_size <- vapply(files, file.size, FUN.VALUE = numeric(1))
          out$`synthpop folder size (Gb)` <-
            signif(sum(vect_size) / (1024^3), 4) # Gb

          # synthpops with same checksum
          files <- list.files(
            private$synthpop_dir,
            paste0("^synthpop_", private$checksum, ".*\\.fst$")
          )

          out$`synthpop meta files with same checksum` <-
            length(list.files(
              private$synthpop_dir,
              paste0("^synthpop_", private$checksum, ".*_meta\\.yaml$")
            ))

          # synthpops with any checksum
          files <-
            list.files(private$synthpop_dir, "^synthpop_.*\\.fst$")

          out$`synthpop meta files with any checksum` <-
            length(list.files(
              private$synthpop_dir,
              "^synthpop_.*_meta\\.yaml$"
            ))

          cat(paste0(names(out), ": ", out, "\n"))
        } else {
          # if length(files) == 0L
          cat("no files found.")
        }
        return(invisible(self))
      },

      # get_checksum ----
      #' @description
      #' Get the synthpop file paths.
      #' @param x One of "all", "synthpop" or "metafile". Can be abbreviated.
      #' @return The invisible `SynthPop` object.
      get_checksum = function() {
        out <- private$checksum
        names(out) <- "Checksum"
        cat(paste0(names(out), ": ", out))
        invisible(self)
      },

      # get_filename ----
      #' @description
      #' Get the synthpop file paths.
      #' @param x One of "all", "synthpop" or "metafile". Can be abbreviated.
      #' @return The invisible `SynthPop` object.
      get_filename = function(x = c("all", "synthpop", "metafile")) {
        if (self$mc == 0L) {
          print("Not relevant because mc = 0L")
        } else {
          x <- match.arg(x)
          switch(
            x,
            all = print(private$filename),
            synthpop = print(private$filename[["synthpop"]]),
            metafile = print(private$filename[["metafile"]])
          )
        }
        invisible(self)
      },

      # get_design ----
      #' @description
      #' Get the synthpop design.
      #' @return The invisible `SynthPop` object.
      get_design = function() {
        # print(private$design)
        # invisible(self)
        private$design
      },

      # get_dir ----
      #' @description
      #' Get the synthpop dir.
      #' @return The invisible `SynthPop` object.
      get_dir = function() {
        print(private$synthpop_dir)
        invisible(self)
      },

      # gen_synthpop_demog ----
      #' @description
      #' Generate synthpop sociodemographics for a random sample of the population in the initial year.
      #' @param design_ A Design object,
      #' @param month April or July are accepted. Use July for mid-year
      #'   population estimates.
      #' @return An invisible `data.table` with sociodemographic information.
      gen_synthpop_demog = function(design_, month = "July") {
        stopifnot(
          "Argument month need to be April or July" = month %in%
            c("April", "July")
        )
        # Use month = July for mid-year
        lsoas_ <- private$get_unique_LSOAs(design_)
        # load dtb
        if (month == "April") {
          file <- "./inputs/pop_estimates_lsoa/LSOA_1st_April_population_estimates"
        } else {
          file <- "./inputs/pop_estimates_lsoa/LSOA_mid_year_population_estimates"
        }

        dtb <- read_parquet_dt(
          file,
          filter = arrow::Expression$field_ref("year") ==
            design_$sim_prm$init_year &
            arrow_in("LSOA11CD", lsoas_)
        )
        dtb[, `99` := `99` + `100`] # 99 means 99+ to be consistent with the rest of the model
        # delete unwanted ages
        dtb[,
          c(
            paste0(0:(design_$sim_prm$ageL - 1L)),
            c(paste0((design_$sim_prm$ageH + 1L):100))
          ) := NULL
        ]

        dtb <-
          melt(
            dtb,
            grep("^[0-9]", names(dtb), value = TRUE, invert = TRUE),
            variable.name = "age",
            value.name = "population_size",
            variable.factor = FALSE
          )
        dtb[, age := as.integer(age)]
        dtb[, year := as.integer(year)]
        dtb[,
          population_proportion := population_size / sum(population_size),
          by = .(age, sex)
        ] # by age & sex for stratified sampling
        dtb[, population_size := NULL]
        # load ethnicity proportions by lsoa
        file <- "./inputs/pop_estimates_lsoa/ethn2011_pct.fst"
        dt_meta <- metadata_fst(file)
        if (!identical("LSOA11CD", dt_meta$keys[1])) {
          stop("Ethnicity file need to be keyed by LSOA")
        }
        file_indx <- read_fst(
          file,
          as.data.table = TRUE,
          columns = "LSOA11CD"
        )[, .(from = min(.I), to = max(.I)), keyby = "LSOA11CD"][
          LSOA11CD %in% lsoas_,
          .("from" = min(from), "to" = max(to))
        ]
        ethn <- read_fst(
          file,
          from = file_indx$from,
          to = file_indx$to,
          as.data.table = TRUE
        )[LSOA11CD %in% lsoas_]

        absorb_dt(dtb, ethn)
        .ethn_nam <- c(
          "white",
          "indian",
          "pakistani",
          "bangladeshi",
          "other asian",
          "black caribbean",
          "black african",
          "chinese",
          "other"
        )

        for (j in .ethn_nam) {
          set(dtb, NULL, j, dtb[, get(j) * population_proportion])
        }
        dtb[, population_proportion := NULL]
        dtb <-
          melt(
            dtb,
            measure.vars = .ethn_nam,
            variable.name = "ethnicity",
            variable.factor = TRUE,
            value.name = "prbl"
          )

        # The following is a new addition to the model on 18/06/24. I need to
        # ensure that all age/sex combinations are represented in the sample.
        # Therefore instead of simple random sampling that was the default
        # till now, I implemented stratified random sampling by age and sex. The
        # people that are sampled could be slightly more than n because I
        # truncate upwards.
        strata_size <- ceiling(
          design_$sim_prm$n /
            ((design_$sim_prm$ageH - design_$sim_prm$ageL + 1L) * 2)
        ) # times two for sex

        # I do not explicitly set.seed because I do so in the gen_synthpop()
        dtinit <- dtb[
          dtb[, sample(.I, strata_size, TRUE, prbl), keyby = .(age, sex)]$V1
        ]

        # Generate the cohorts of 30 year old to enter every year
        # as sim progress these will become 30 yo
        # no population growth here as I will calibrate to pop
        # projections and it fluctuates at +-2% anyways (for age == 30).

        if (design_$sim_prm$logs) {
          message(
            "Generating the cohorts of ",
            design_$sim_prm$ageL,
            " year old"
          )
        }
        # Assuming that the new cohorts of 30 year olds entering the model
        # have the same age-sex-dimd-ethnicity distribution as in the initial
        # year. This is not a huge problem if the init year is post 2010
        # because the distribution is fairly stable since then with perhaps an
        # increase in deprivation. There are very strong trends before 2010.
        dtb <- dtb[age == design_$sim_prm$ageL]
        dtfut <- dtb[
          dtb[,
            sample(
              .I,
              strata_size * design_$sim_prm$sim_horizon_max,
              TRUE,
              prbl
            ),
            keyby = sex
          ]$V1
        ]
        dtfut[,
          age := age - rep(design_$sim_prm$sim_horizon_max:1, strata_size)
        ]

        dtb <- rbind(dtfut, dtinit)
        dtb[, prbl := NULL]

        indx_hlp <-
          read_fst(
            "./inputs/pop_estimates_lsoa/lsoa_to_locality_indx.fst",
            as.data.table = TRUE
          )

        dtb[
          indx_hlp,
          on = "LSOA11CD",
          `:=`(
            # tds = i.tds,
            # tds_quintile = i.tds_quintile,
            # imd = i.imd,
            qimd = i.qimd,
            dimd = i.dimd,
            sha = i.SHA11NM,
            # CCG17CDH = i.CCG17CDH,
            LAD17CD = i.LAD17CD # needed for calibration of weights to pop if design$sim_prm$calibrate_to_pop_projections_by_LAD: yes
          )
        ]
        setkey(dtb, year, age, sex)
        return(invisible(dtb))
      },

      # write_synthpop ----
      #' @description
      #' Generate synthpop files in parallel, using foreach, and writes them to
      #' disk. It skips files that are already on disk.
      #' Note: the backend for foreach needs to be initialised before calling
      #' the function.
      #' @param mc_ An integer vector for the Monte Carlo iteration of the
      #'   synthetic population. Each integer generates a unique synthetic
      #'   population.
      #' @return The invisible `SynthPop` object.
      write_synthpop = function(mc_) {
        stopifnot(all(is.numeric(mc_)), all(ceiling(mc_) > 0L))
        on.exit(self$delete_incomplete_synthpop(), add = TRUE)
        mc_ <- as.integer(ceiling(mc_))

        # get unique lsoas
        lsoas <- private$get_unique_LSOAs(private$design)

        if (.Platform$OS.type == "windows") {
          # TODO update to make compatible with windows
          cl <-
            makeCluster(private$design$sim_prm$clusternumber) # used for clustering. Windows compatible
          registerDoParallel(cl)
        } else {
          registerDoParallel(private$design$sim_prm$clusternumber) # used for forking. Only Linux/OSX compatible
        }

        foreach(
          mc_iter = mc_,
          .inorder = FALSE,
          .options.multicore = list(preschedule = FALSE),
          .verbose = private$design$sim_prm$logs,
          .packages = c(
            "R6",
            "gamlss.dist",
            # For distr in prevalence.R
            "dqrng",
            "qs2",
            "fst",
            "arrow",
            "CKutils",
            "IMPACTncdEngland",
            "data.table"
          ),
          .export = NULL,
          .noexport = NULL # c("time_mark")
        ) %dopar%
          {
            # Staggered worker startup to reduce memory pressure during parallel runs.
            # Only apply to the first batch of workers (iterations 1 to clusternumber).
            # After the first batch, workers finish at different times and naturally
            # pick up new iterations in a staggered manner.
            n_cores <- private$design$sim_prm$clusternumber
            if (n_cores > 1L && mc_iter <= n_cores) {
              stagger_delay_sec <- 2L # seconds between worker startups
              worker_position <- mc_iter - 1L
              if (worker_position > 0L) {
                Sys.sleep(worker_position * stagger_delay_sec)
              }
            }

            arrow::set_cpu_count(1L)
            data.table::setDTthreads(
              threads = 1L,
              restore_after_fork = NULL
            )
            fst::threads_fst(
              nr_of_threads = 1L,
              reset_after_fork = NULL
            )
            filename <- private$gen_synthpop_filename(
              mc_iter,
              private$checksum,
              private$design
            )

            # logic for the synthpop load
            files_exist <- vapply(
              filename,
              file.exists,
              FUN.VALUE = logical(1)
            )
            if (all(!files_exist)) {
              # No files exist. Create the synthpop and store
              # the file on disk
              private$gen_synthpop(mc_iter, filename, private$design)
            } else if (
              file.exists(filename$metafile) &&
                !all(files_exist)
            ) {
              # Metafile exists but not all three files. It means
              # that most likely a generate_synthpop() is still running. So the
              # function waits until the file is created before it proceeds to
              # load it. Note that if this is not the case then the loop is
              # infinite!!!
              while (
                !all(vapply(
                  filename,
                  file.exists,
                  FUN.VALUE = logical(1)
                ))
              ) {
                Sys.sleep(5)
              }

              # Ensure the file write is complete (size stable)
              sz1 <- file.size(filename$synthpop)
              Sys.sleep(3)
              sz2 <- file.size(filename$synthpop)
              while (sz1 != sz2) {
                sz1 <- file.size(filename$synthpop)
                Sys.sleep(3)
                sz2 <- file.size(filename$synthpop)
              }
            } else if (
              !file.exists(filename$metafile) &&
                !all(files_exist)
            ) {
              # Metafile doesn't exist but some other files exist. In this case
              # delete everything and start from scratch
              self$delete_incomplete_synthpop()
              private$gen_synthpop(mc_iter, filename, private$design)
            }
            # No need to provision for case when all files present.

            return(NULL)
          }
        if (exists("cl")) {
          stopCluster(cl)
        }

        invisible(self)
      },

      # get_risks  ----
      #' @description Get the risks for all individuals in a synthetic
      #'   population for a disease.
      #' @param disease_nam The disease that the risks will be returned.
      #' @return A data.table with columns for pid, year, and all associated
      #'   risks if disease_nam is specified. Else a list of data.tables for all
      #'   diseases.
      get_risks = function(disease_nam) {
        if (missing(disease_nam)) {
          return(private$risks)
        } else {
          stopifnot(is.character(disease_nam))
          return(private$risks[[disease_nam]])
        }
      },

      # store_risks ----
      #' @description Stores the disease risks for all individuals in a synthetic
      #'   population in a private list.
      #' @param disease_nam The disease that the risks will be stored.
      #' @return The invisible self for chaining.
      store_risks = function(disease_nam) {
        stopifnot(is.character(disease_nam))

        nam <- grep("_rr$", names(self$pop), value = TRUE)

        private$risks[[disease_nam]] <-
          self$pop[, .SD, .SDcols = c("pid", "year", nam)]

        self$pop[, (nam) := NULL]
        invisible(self)
      },

      # print ----
      #' @description
      #' Prints the synthpop object metadata.
      #' @return The invisible `SynthPop` object.
      print = function() {
        print(c(
          "path" = ifelse(
            self$mc == 0L,
            "Not relevant because mc = 0L",
            private$filename$synthpop
          ),
          "checksum" = private$checksum,
          "mc" = self$mc,
          self$metadata
        ))
        invisible(self)
      }
    ),

    # private -----------------------------------------------------------------
    private = list(
      filename = NA,
      checksum = NA,
      # The design object with the simulation parameters.
      design = NA,
      synthpop_dir = NA,
      risks = list(), # holds the risks for all individuals
      ladcontribution = data.table(), # holds the LAD contribution to the synthpop
      # Special deep copy for data.table. Use POP$clone(deep = TRUE) to
      # dispatch. Otherwise a reference is created
      # deep_clone ----
      # @param name The name associated with the object.
      # @param value The object to be cloned.
      #
      # @return A deep clone of the input object.
      #
      # @details
      # The function uses methods for cloning based on the class of the input object.
      # It performs a deep copy for data.table objects using the `data.table::copy` function,
      # and for R6 objects, it utilizes the `clone` method. For other classes, a shallow copy is returned.
      #
      deep_clone = function(name, value) {
        if ("data.table" %in% class(value)) {
          data.table::copy(value)
        } else if ("R6" %in% class(value)) {
          value$clone()
        } else {
          # For everything else, just return it. This results in a shallow
          # copy of s3.
          value
        }
      },

      # get all unique LSOAs included in locality vector
      # get_unique_LSOAs ----
      # @param design_ The design object containing simulation parameters and locality information.
      #
      # @return A sorted vector of unique LSOA codes.
      #
      # @details
      # The function reads an index file containing LSOA information and extracts unique LSOAs based on the
      # specified locality criteria in the design object. If "England" is included in the locality criteria,
      # it returns all LSOAs at the national level; otherwise, it filters LSOAs based on the specified locality criteria.
      #
      get_unique_LSOAs = function(design_) {
        indx_hlp <-
          read_fst(
            "./inputs/pop_estimates_lsoa/lsoa_to_locality_indx.fst",
            as.data.table = TRUE,
            columns = c("LSOA11CD", "LAD17NM", "RGN11NM", "ICB22NM")
          )

        if ("England" %in% design_$sim_prm$locality) {
          lsoas <- indx_hlp[, unique(LSOA11CD)] # national
        } else {
          lsoas <-
            indx_hlp[
              LSOA11CD %in%
                design_$sim_prm$locality |
                LAD17NM %in% design_$sim_prm$locality |
                RGN11NM %in% design_$sim_prm$locality |
                ICB22NM %in% design_$sim_prm$locality,
              unique(LSOA11CD)
            ]
        }
        return(sort(lsoas))
      },

      # get all unique LADs included in locality vector.
      # get_unique_LADs ----
      # @param design_ The design object containing simulation parameters and locality information.
      #
      # @return A sorted vector of unique LAD codes.
      #
      # @details
      # The function reads an index file containing LAD information and extracts unique LADs based on the
      # specified locality criteria in the design object. If "England" is included in the locality criteria,
      # it returns all LADs at the national level; otherwise, it filters LADs based on the specified locality criteria.
      #
      get_unique_LADs = function(design_) {
        indx_hlp <-
          read_fst(
            "./inputs/pop_estimates_lsoa/lsoa_to_locality_indx.fst",
            as.data.table = TRUE,
            columns = c("LSOA11CD", "LAD17NM", "RGN11NM", "ICB22NM", "LAD17CD")
          )

        if ("England" %in% design_$sim_prm$locality) {
          lads <- indx_hlp[, unique(LAD17CD)] # national
        } else {
          lads <-
            indx_hlp[
              LSOA11CD %in%
                design_$sim_prm$locality |
                LAD17NM %in% design_$sim_prm$locality |
                RGN11NM %in% design_$sim_prm$locality |
                ICB22NM %in% design_$sim_prm$locality,
              unique(LAD17CD)
            ]
        }
        return(sort(lads))
      },

      # get a smaller design list only with characteristics that are important
      # for synthpop creation and define the uniqueness of the object. I.e. if
      # these parameters are different the synthpop has to have different
      # filename and vice-versa
      # get_unique_characteristics ----
      get_unique_characteristics = function(design_) {
        design_$sim_prm[c(
          "n",
          "sim_horizon_max",
          "init_year_long",
          "maxlag",
          "smoking_relapse_limit",
          "ageL",
          "ageH",
          "jumpiness",
          "simsmok_calibration",
          "statin_adherence",
          "bpmed_adherence"
        )]
      },

      # gen synthpop unique checksum for the given set of inputs
      # gen_checksum ----
      # @param design_ The design object containing simulation parameters.
      #
      # @return A data frame containing unique simulation characteristics.
      #
      # @details
      # The function extracts specific simulation characteristics from the design parameters, including:
      # - n: Number of simulations
      # - sim_horizon_max: Maximum simulation horizon
      # - init_year_long: Initial year for long simulations
      # - maxlag: Maximum lag for time-dependent covariates
      # - smoking_relapse_limit: Smoking relapse limit
      # - ageL: Lower age limit
      # - ageH: Upper age limit
      # - jumpiness: Jumpiness parameter
      # - simsmok_calibration: Calibration parameter for simsmok model
      # - statin_adherence: Adherence parameter for statin medication
      # - bpmed_adherence: Adherence parameter for blood pressure medication
      #
      gen_checksum = function(design_) {
        # get a md5 checksum based on function arguments
        # First get function call arguments
        fcall <- private$get_unique_characteristics(design_)

        lsoas_ <- private$get_unique_LSOAs(design_)

        locality_years_age_id <-
          digest(
            paste(lsoas_, fcall, sep = ",", collapse = ","),
            serialize = FALSE
          )
        return(locality_years_age_id)
      },

      # gen_synthpop_filename ----
      # for the given set of inputs
      # @param mc_ The Monte Carlo iteration.
      # @param checksum_ The checksum used for file identification.
      # @param design_ The design object containing simulation parameters and file directory information.
      #
      # @return A list containing the paths to the synthpop data file ("synthpop") and its metafile ("metafile").
      #
      # @details
      # The function constructs filenames based on the specified parameters including the Monte Carlo iteration,
      # checksum, and the directory specified in the design object. The resulting paths are normalized using `normalizePath`.
      #
      gen_synthpop_filename = function(mc_, checksum_, design_) {
        return(
          list(
            "synthpop" = normalizePath(
              paste0(
                design_$sim_prm$synthpop_dir,
                "/synthpop_",
                checksum_,
                "_",
                mc_,
                ".fst"
              ),
              mustWork = FALSE
            ),
            "metafile" = normalizePath(
              paste0(
                design_$sim_prm$synthpop_dir,
                "/synthpop_",
                checksum_,
                "_",
                mc_,
                "_meta.yaml"
              ),
              mustWork = FALSE
            )
          )
        )
      },

      # del_incomplete ----
      # @param filename_ A list containing the paths to the synthpop data file ("synthpop") and its metafile ("metafile").
      #
      # @details
      # The function checks if the metafile exists and the synthpop file does not exist. If these conditions are met,
      # it deletes both files using \code{vapply(filename_, file.remove, logical(1))}.
      #
      del_incomplete = function(filename_) {
        if (
          file.exists(filename_$metafile) &&
            (!file.exists(filename_$synthpop))
        ) {
          suppressWarnings(vapply(
            filename_,
            file.remove,
            FUN.VALUE = logical(1)
          ))
        }
      },

      # gen_synthpop ----
      # @details
      # Generates a complete synthetic population for a single Monte Carlo iteration,
      # including demographics, correlated exposure trajectories, and risk factor
      # distributions. The function creates a lifecourse dataset spanning from
      # sim_horizon_min - maxlag to sim_horizon_max, applying random walks to ranks
      # and generating exposures from statistical distributions (ZINBI, BCPEo, DEL, etc.).
      # The result is written directly to disk as .fst files with associated YAML metadata.
      # returns NULL. Writes synthpop on disk
      gen_synthpop = function(mc_, filename_, design_) {
        # increase design_$sim_prm$jumpiness for more erratic jumps in
        # trajectories

        # In Shiny app this function runs as a future. It is not
        # straightforward to check whether the future has been resolved or
        # not. To circumvent the problem I will save the metafile here (almost
        # function beginning) and the synthpop file at the end. So if both
        # files exist the function has finished. If only metafile exists the
        # function probably still runs.

        # Save synthpop metadata
        if (!file.exists(filename_$metafile)) {
          yaml::write_yaml(
            private$get_unique_characteristics(design_),
            filename_$metafile
          )
        }
        # NOTE In shiny app if 2 users click the  button at the same time, 2
        # functions will run almost concurrently with potential race condition

        # To avoid edge cases when the function stopped prematurely and a
        # metafile was created while the file was not. On.exit ensures that
        # either both files exist or none.

        on.exit(private$del_incomplete(filename_), add = TRUE)

        dqRNGkind("pcg64")
        SEED <-
          2121870L # sample(1e7, 1) # Hard-coded for reproducibility
        set.seed(SEED + mc_)
        dqset.seed(SEED, mc_)

        # Generate synthpops with sociodemographic and exposures information.

        dtb <- self$gen_synthpop_demog(design_, month = "April")

        # NOTE!! from now on year in the short form i.e. 13 not 2013
        dtb[, `:=`(pid = .I)]
        new_n <- nrow(dtb)

        # Generate correlated ranks for the individuals ----
        if (design_$sim_prm$logs) {
          message("Generate correlated ranks for the individuals")
        }

        cm_mean <- as.matrix(
          read_parquet_dt(
            "./inputs/exposure_distributions/exposure_corr_mean/part-0.parquet"
          ),
          rownames = "rn"
        )

        tr <- which(
          colnames(cm_mean) %in%
            c("af_r", "ckd_r", "famcvd_r", "dm_r", "dm_dgn_r")
        )
        cm_mean <- cm_mean[-tr, -tr]

        rank_mtx <- generate_corr_unifs(new_n, cm_mean)
        if (design_$sim_prm$logs) {
          message("generate correlated uniforms")
        }

        # Restrict the range of some RNs to avoid unrealistic exposures
        # This scaling does not affect correlations
        # /0.999 because I multiplied all the columns below
        # rank_mtx <- rank_mtx * 0.999
        # rank_mtx[, "frtpor_r"] <- rank_mtx[, "frtpor_r"] * 0.99 / 0.999
        # rank_mtx[, "vegpor_r"] <- rank_mtx[, "vegpor_r"] * 0.93 / 0.999
        # rank_mtx[, "smok_cig_ex_r"] <-
        #   rank_mtx[, "smok_cig_ex_r"] * 0.99 / 0.999
        # rank_mtx[, "totalwu_r"] <- rank_mtx[, "totalwu_r"] * 0.99 / 0.999
        # rank_mtx[, "smok_quit_yrs_r"] <-
        #   rank_mtx[, "smok_quit_yrs_r"] * 0.99 / 0.999
        # rank_mtx[, "smok_dur_ex_r"] <-
        #   rank_mtx[, "smok_dur_ex_r"] * 0.99 / 0.999
        # rank_mtx[, "smok_dur_curr_r"] <-
        #   rank_mtx[, "smok_dur_curr_r"] * 0.88 / 0.999
        # rank_mtx[, "bmi_r"] <-
        #   rank_mtx[, "bmi_r"] * 0.95 / 0.999
        # sum((cor(rank_mtx) - cm_mean) ^ 2)
        if (design_$sim_prm$logs) {
          message("correlated ranks matrix to data.table")
        }

        rank_mtx <- data.table(rank_mtx)

        # NOTE rankstat_* is unaffected by the RW. Stay constant through the lifecourse
        dtb[,
          c(
            "rank_education",
            "rank_income",
            "rank_pa",
            "rank_fruit",
            "rank_veg",
            "rankstat_smok",
            "rankstat_smok_quit_yrs",
            "rankstat_smok_dur_ex",
            "rankstat_smok_dur_curr",
            "rankstat_smok_cig_ex",
            "rankstat_smok_cig_curr",
            "rank_ets",
            "rank_alcohol",
            "rank_bmi",
            "rank_sbp",
            "rank_bpmed",
            "rank_tchol",
            "rank_hdl",
            "rank_statin_px"
          ) := rank_mtx
        ]

        rm(rank_mtx)

        # add non-correlated RNs
        rank_cols <-
          c(
            # "rankstat_ncc",
            "rankstat_simsmok",
            "rankstat_pa_dur",
            "rankstat_pa_met",
            "rankstat_bpmed_adherence",
            "rankstat_statin_adherence"
          )

        for (nam in rank_cols) {
          set(dtb, NULL, nam, dqrunif(new_n))
        } # NOTE do not replace with generate_rns function.

        # Generate education (exception as it remains stable through lifecourse) ----
        if (design_$sim_prm$logs) {
          message("Generate education")
        }
        if (max(dtb$age) > 90L) {
          dtb[, age100 := age]
          dtb[age > 90L, age := 90L]
        }

        tbl <- design_$exposures$education$get_table()

        nam <- intersect(names(dtb), names(tbl))
        # logic necessary for new cohorts entering the simulation that currently age < 30
        # These will have the same distribution as if 30 years old
        tt <- tbl[age == min(age)]
        tt <- clone_dt(tt, design_$sim_prm$sim_horizon_max) # TODO adding design_$sim_prm$sim_horizon_max + design_$sim_prm$maxlag for longer projections
        tt[, age := age - .id] # as the sim progress these will become 30 yo
        # increase population by 0.5% every year
        tt[, .id := NULL]
        tbl <- rbind(tt, tbl)
        dtb[
          tbl,
          education := (rank_education > ed1) +
            (rank_education > ed2) +
            (rank_education > ed3) +
            (rank_education > ed4) +
            (rank_education > ed5) +
            (rank_education > ed6) +
            1L,
          on = nam
        ]
        dtb[,
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
        if (!design_$sim_prm$keep_simulants_rn) {
          dtb[, rank_education := NULL]
        }

        if ("age100" %in% names(dtb)) {
          dtb[, age := NULL]
          setnames(dtb, "age100", "age")
        }

        # Project forward for simulation and back project for lags  ----
        if (design_$sim_prm$logs) {
          message("Project forward and back project")
        }

        dtb <-
          clone_dt(
            dtb,
            design_$sim_prm$sim_horizon_max +
              design_$sim_prm$maxlag +
              1L
          )

        dtb[
          .id <= design_$sim_prm$maxlag,
          `:=`(age = age - .id, year = year - .id)
        ]
        dtb[
          .id > design_$sim_prm$maxlag,
          `:=`(
            age = age + .id - design_$sim_prm$maxlag - 1L,
            year = year + .id - design_$sim_prm$maxlag - 1L
          )
        ]
        # dtb <-
        #   dtb[between(age, design_$sim_prm$ageL - design_$sim_prm$maxlag, design_$sim_prm$ageH)]
        # delete unnecessary ages
        del_dt_rows(
          dtb,
          !between(
            dtb$age,
            design_$sim_prm$ageL - design_$sim_prm$maxlag,
            design_$sim_prm$ageH
          ),
          environment()
        )

        dtb[, `:=`(.id = NULL)]

        if (max(dtb$age) > 90L) {
          dtb[, age100 := age]
          dtb[age > 90L, age := 90L]
        }

        # to_agegrp(dtb, 20L, 85L, "age", "agegrp20", to_factor = TRUE)
        # to_agegrp(dtb, 10L, 85L, "age", "agegrp10", to_factor = TRUE)
        # to_agegrp(dtb,  5L, 85L, "age", "agegrp5" , to_factor = TRUE)

        # Simulate exposures -----

        # Random walk for ranks ----
        if (design_$sim_prm$logs) {
          message("Random walk for ranks")
        }

        setkeyv(dtb, c("pid", "year"))
        setindexv(dtb, c("year", "age", "sex", "sha", "qimd", "ethnicity"))

        dtb[, pid_mrk := mk_new_simulant_markers(pid)]

        dtb[,
          lapply(
            .SD,
            fscramble_trajectories,
            pid_mrk,
            design_$sim_prm$jumpiness
          ),
          .SDcols = patterns("^rank_")
        ]
        # ggplot2::qplot(year, rank_income, data = dtb[pid %in% sample(1e5, 1)], ylim = c(0,1))

        # Generate income ----
        design_$exposures$income$generate(dtb, design_)
        if (!design_$sim_prm$keep_simulants_rn) {
          dtb[, rank_income := NULL]
        }

        # Generate active days ----
        if (design_$sim_prm$logs) {
          message("Generate active days")
        }
        design_$exposures$active_days$generate(dtb, design_)
        if (!design_$sim_prm$keep_simulants_rn) {
          dtb[, rank_pa := NULL]
        }

        # dtb[, active_days := factor(active_days, levels = 0:7)]

        # convert active days to MET mins/Week I have number of days each
        # week individuals make at least moderate activity for at least 30
        # min. From
        # http://healthsurvey.hscic.gov.uk/media/63730/HSE16-Adult-phy-act.pdf
        # Moderate physical activity (MPA) includes activities with
        # estimated intensity levels of 3-6 METs; vigorous physical
        # activities (VPA) are those with estimated intensity levels of 6
        # METs or higher.

        # For MET distribution I will use hist(3 + rbinom(1e5, 8, 3/11))
        # arbitrarily Note that GBD has MET RR between 0 and 4200. Shouldn't
        # depend on age too much as MET are relative to the basic metabolism

        # For the total min that someone is active I will use the distribution hist(30 +
        # rexp(1e5, 1/7))
        dtb[,
          met := as.integer(floor(
            active_days *
              (3L + qbinom(rankstat_pa_met, 8, 3 / 11)) *
              (30 + qexp(rankstat_pa_dur, 1 / 7)) /
              100
          ))
        ] # TODO make data driven
        if (!design_$sim_prm$keep_simulants_rn) {
          dtb[, c("rankstat_pa_met", "rankstat_pa_dur") := NULL]
        }

        # Generate fruit consumption (ZISICHEL) ----
        design_$exposures$fruit$generate(dtb, design_)
        if (!design_$sim_prm$keep_simulants_rn) {
          dtb[, c("rankstat_pa_met", "rankstat_pa_dur") := NULL]
        }

        # Generate veg consumption (DEL) ----
        design_$exposures$veg$generate(dtb, design_)
        if (!design_$sim_prm$keep_simulants_rn) {
          dtb[, rank_veg := NULL]
        }
        # Smoking simulation ----
        if (design_$sim_prm$logs) {
          message("Smoking simulation")
        }

        # Assign smok_status when pid_mrk == true (the first year an individual enters the simulation (with lags))
        design_$exposures$smok_status$generate(dtb, design_)
        dtb[, smok_status_ref := smok_status]

        # Initialize smoking duration variables
        set(dtb, NULL, "smok_quit_yrs", 0L)
        set(dtb, NULL, "smok_dur", 0L)

        # Assign smok_quit_yrs and smok_dur for ex-smokers (pid_mrk & smok_status %in% 2:3)
        # Only generate values for rows that need them for efficiency
        idx <- dtb[, which((pid_mrk) & smok_status %in% 2:3)]
        design_$exposures$smok_quit_yrs$generate(dtb, design_, idx)
        design_$exposures$smok_dur_ex$generate(dtb, design_, idx)

        # Assign smok_dur for current smokers (pid_mrk & smok_status == 4)
        idx <- dtb[, which((pid_mrk) & smok_status == 4L)]
        design_$exposures$smok_dur_curr$generate(dtb, design_, idx)

        # Clean up rank variables if not keeping
        if (!design_$sim_prm$keep_simulants_rn) {
          for (rn in c(
            "rankstat_smok_quit_yrs",
            "rankstat_smok_dur_ex",
            "rankstat_smok_dur_curr"
          )) {
            if (rn %in% names(dtb)) dtb[, (rn) := NULL]
          }
        }

        # Ensure smoking histories start from age 12
        dtb[age - smok_quit_yrs < 12L, smok_quit_yrs := age - 12L]
        dtb[age - smok_dur < 12L, smok_dur := age - 12L]
        dtb[
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
        design_$exposures$smok_incid$generate(dtb, design_)

        # Assign smok_cessation probabilities
        design_$exposures$smok_cess$generate(dtb, design_)

        # Handle smok_relapse probabilities
        tbl <-
          read_parquet_dt("./inputs/exposure_distributions/smok_relapse")
        tbl <-
          dcast(tbl, sex + qimd ~ smok_quit_yrs, value.var = "pr")
        nam <- tbl[, paste0(sex, " ", qimd)]
        tbl <-
          as.matrix(tbl[, mget(paste0(1:15))], rownames = nam)

        simsmok(dtb, tbl, design_$sim_prm$smoking_relapse_limit)
        # dtb[!(pid_mrk), table(smok_status)]
        # dtb[pid == 1, plot(year, smok_status, ylim = c(0, 4))]
        # dtb[pid == 10, .(age, smok_status, smok_quit_yrs, smok_dur)]
        # dtb[, sum(smok_status == 4)/.N, keyby = year]

        if (design_$sim_prm$simsmok_calibration) {
          # calculate dif between ref (multinom) and simsmok
          # I will further calibrate to better match HSE
          resample <-
            function(x, ...) {
              x[sample.int(length(x), ...)]
            }
          obs <-
            dtb[smok_status == 1L, .(nsa = .N), keyby = .(year, age, sex, qimd)] # ? add sha/ethn
          ref <-
            dtb[
              smok_status_ref == 1L,
              .(nsr = .N),
              keyby = .(year, age, sex, qimd)
            ] # ? add sha/ethn
          absorb_dt(ref, obs)
          setnafill(ref, "c", 0L, cols = "nsa")
          ref[, `:=`(dif = nsr - nsa, nsr = NULL, nsa = NULL)]
          # Further calibrate to match better with HSE
          ref[
            sex == "men" &
              rbinom(.N, 1L, 0.9) == 1L,
            dif := dif - 2L
          ]
          ref[
            sex == "men" &
              rbinom(.N, 1L, 1 * clamp((year - min(year)) / 10, 0, 1)) == 1L,
            dif := dif + 2L
          ]
          ref[
            sex == "women" &
              rbinom(.N, 1L, 0.2) == 1L,
            dif := dif - 1L
          ]
          ref[
            age < 49 &
              rbinom(.N, 1L, 0.5) == 1L,
            dif := dif + 1L
          ]
          ref[
            age < 49 &
              sex == "men" &
              qimd == "3" &
              rbinom(.N, 1L, 0.5) == 1L,
            dif := dif - 1L
          ]
          ref[
            age < 49 &
              sex == "women" &
              qimd %in% c("4", "5 least deprived") &
              rbinom(.N, 1L, 0.4) == 1L,
            dif := dif - 1L
          ]
          ref[
            between(age, 50, 69) &
              rbinom(.N, 1L, 0.2) == 1L,
            dif := dif + 1L
          ]
          ref[
            between(age, 50, 69) &
              sex == "women" &
              qimd == "5 least deprived" &
              rbinom(.N, 1L, 0.4) == 1L,
            dif := dif + 1L
          ]
          ref[
            between(age, 70, 89) &
              rbinom(.N, 1L, 0.4) == 1L,
            dif := dif + 1L
          ]
          ref[
            between(age, 70, 89) &
              sex == "men" &
              qimd %in% c("1 most deprived", "2") &
              rbinom(.N, 1L, 0.4) == 1L,
            dif := dif - 1L
          ]
          ref[
            between(age, 70, 89) &
              sex == "women" &
              qimd %in% c("4d", "2") &
              rbinom(.N, 1L, 0.4) == 1L,
            dif := dif + 1L
          ]
          absorb_dt(dtb, ref)

          # when not enough never smokers convert those ex smokers with the longer quit years
          tt <-
            dtb[
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
          dtb[
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
            dtb[smok_status %in% 1, .(year, age, sex, qimd, pid, dif)]
          setnafill(tt, "c", 0, cols = "dif")
          tt[dif > 0, dif := 0L]
          tt[, dif := -dif]
          setkey(tt, year, age, sex, qimd)
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
              .(pid = resample(pid, max(dif))),
              keyby = .(year, age, sex, qimd)
            ]
          dtb[
            pid_to_conv,
            on = .(year, pid),
            `:=`(
              smok_status = 2L,
              smok_quit_yrs = dtb[smok_status == 2L, median(smok_quit_yrs)],
              smok_dur = dtb[smok_status == 2L, median(smok_dur)],
              smok_cig = 1L
            )
          ]
          dtb[, dif := NULL]

          # Same logic for active smokers
          # I will further calibrate to better match HSE
          # calculate dif between ref (multinom) and simsmok
          obs <-
            dtb[smok_status == 4L, .(nsa = .N), keyby = .(year, age, sex, qimd)] # ? add sha/ethn
          ref <-
            dtb[
              smok_status_ref == 4L,
              .(nsr = .N),
              keyby = .(year, age, sex, qimd)
            ] # ? add sha/ethn
          absorb_dt(ref, obs)
          setnafill(ref, "c", 0, cols = "nsa")
          ref[, `:=`(dif = nsr - nsa, nsr = NULL, nsa = NULL)]
          # Further calibrate to match better with HSE
          ref[
            sex == "men" &
              rbinom(.N, 1L, 0.3) == 1L,
            dif := dif - 1L
          ]
          ref[
            sex == "women" &
              qimd != "1 most deprived" &
              qimd != "5 least deprived" &
              rbinom(.N, 1L, 0.8) == 1L,
            dif := dif - 1L
          ]
          ref[
            qimd == "1 most deprived" &
              rbinom(.N, 1L, 0.5) == 1L,
            dif := dif + 1L
          ]
          ref[
            qimd == "5 least deprived" &
              rbinom(.N, 1L, 0.6) == 1L,
            dif := dif - 1L
          ] # - reduces
          absorb_dt(dtb, ref)

          # when not enough active smokers convert those ex smokers with the shortest quit years
          tt <-
            dtb[
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
          dtb[
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
            dtb[smok_status == 4L, .(year, age, sex, qimd, pid, dif)]
          setnafill(tt, "c", 0, cols = "dif")
          tt[dif > 0, dif := 0L]
          tt[, dif := -dif]
          setkey(tt, year, age, sex, qimd)
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
              .(pid = resample(pid, max(dif))),
              keyby = .(year, age, sex, qimd)
            ]
          dtb[
            pid_to_conv,
            on = .(year, pid),
            `:=`(smok_status = 3L, smok_quit_yrs = 1L)
          ]
          dtb[, dif := NULL]

          rm(tt, ttt, obs, ref, pid_to_conv)
        }

        # Assign smok_cig
        set(dtb, NULL, "smok_cig", 0L)

        # smok_cig_curr for current smokers (smok_status == 4)
        idx <- dtb[, which(smok_status == 4L)]
        design_$exposures$smok_cig_curr$generate(dtb, design_, idx)

        # smok_cig_ex for ex-smokers at initialization (pid_mrk & smok_status == 3)
        idx <- dtb[, which((pid_mrk) & smok_status == 3L)]
        design_$exposures$smok_cig_ex$generate(dtb, design_, idx)

        simsmok_cig(dtb) # carry forward smok_cig if smok_status == 3
        dtb[smok_cig == 0L & smok_status > 1L, smok_cig := 1L]

        if (design_$sim_prm$simsmok_calibration) {
          simsmok_postcalibration(dtb)
        } # need to be post cig simulation

        dtb[, smok_status := factor(smok_status)]

        colnam <- c("prb_smok_incid", "prb_smok_cess", "smok_status_ref")
        if (!design_$sim_prm$keep_simulants_rn) {
          colnam <- c(
            colnam,
            "rankstat_simsmok"
          )
        }
        dtb[, (colnam) := NULL]

        # Generate ETS (BI) ----

        # Note at the moment this is independent of smoking prevalence TODO
        # calculate how many each smoker pollutes by year, SHA (not qimd) to
        # be used in scenarios. Ideally correct for mortality
        design_$exposures$ets$generate(dtb, design_)
        if (!design_$sim_prm$keep_simulants_rn) {
          dtb[, rank_ets := NULL]
        }
        # NOTE for line above invert = true in sim_design yaml so 1 - mu to be
        # equivalent to qbinom(rank_ets, 1, mu). Otherwise rank_ets < mu is
        # equivalent to qbinom(rank_ets, 1, mu, lower.tail = FALSE)). The two
        # have the same prevalence, but correlations are captured correctly only
        # with the 1-mu variant.

        # View(dtb[, prop_if(ets == 1)/prop_if(smok_status == "4"), keyby = .(year, sha)])

        # Generate alcohol (ZINBI) ----
        design_$exposures$alcohol$generate(dtb, design_)
        if (!design_$sim_prm$keep_simulants_rn) {
          dtb[, rank_alcohol := NULL]
        }
        # Generate BMI (BCPEo) ----
        design_$exposures$bmi$generate(dtb, design_)
                if (!design_$sim_prm$keep_simulants_rn) {
                  dtb[, rank_bmi := NULL]
                }

        # Generate SBP (BCPEo) ----
        design_$exposures$sbp$generate(dtb, design_)
        if (!design_$sim_prm$keep_simulants_rn) {
          dtb[, rank_sbp := NULL]
        }
        # Generate BP medication (BI) -----
        # Temporarily round sbp for join, then restore
        dtb[, `:=`(
          sbp_acc = sbp,
          sbp = as.integer(round(clamp(sbp, 110, 200), -1))
        )]
        design_$exposures$bp_med$generate(dtb, design_)
        dtb[, `:=`(sbp = sbp_acc, sbp_acc = NULL)]

        # TODO calculate probability of dgn HTN

        # Generate tchol (BCT) ----
        design_$exposures$tchol$generate(dtb, design_)
        if (!design_$sim_prm$keep_simulants_rn) {
          dtb[, rank_tchol := NULL]
        }

        # Generate HDL (to tchol ratio) (GB1) ----
        design_$exposures$hdl_to_tchol$generate(dtb, design_)
        if (!design_$sim_prm$keep_simulants_rn) {
          dtb[, rank_hdl := NULL]
        }
        # NOTE this very highly correlated with hdl level (~0.76) and
        #  highly to tchol (~-0.47). The latter is captured by the correlated RNs

        # Generate statins medication (BI) -----
        # Temporarily round tchol for join, then restore
        dtb[, `:=`(tchol_acc = tchol, tchol = round(clamp(tchol, 2, 12), 0))]
        design_$exposures$statin_px$generate(dtb, design_)
        dtb[, `:=`(tchol = tchol_acc, tchol_acc = NULL)]

        # Estimate number of comorbidities (ncc) calculation ----
        # to be used in QALY
        # if (design_$sim_prm$logs) message("Generate ncc")
        #
        # dtb[, ncc := as.integer(clamp(qbinom(rankstat_ncc, ceiling(age / 8L),
        #                                     fifelse(age < 55, 0.25, 0.40)),
        #                              0, 10))]
        # calibrated to Sullivan et all 2011 (web table 1)
        # to_agegrp(output, 10L, 89L, "age", "agegrp10", to_factor = TRUE)
        # output[, round(mean(ncc), 1), keyby = agegrp10]
        # target by agegrp 1.1  1.6  2.4  3.1  4.0  4.4 from

        # dtb[, `:=` (
        #   pid_mrk = NULL,
        #   # to be recreated when loading synthpop
        #   rankstat_ncc = NULL
        # )]

        dtb[,
          statin_adherence := qBE(
            rankstat_statin_adherence,
            design_$sim_prm$statin_adherence,
            0.2
          )
        ]
        dtb[,
          bpmed_adherence := qBE(
            rankstat_bpmed_adherence,
            design_$sim_prm$bpmed_adherence,
            0.2
          )
        ]
        if (!design_$sim_prm$keep_simulants_rn) {
          dtb[,
            c("rankstat_bpmed_adherence", "rankstat_statin_adherence") := NULL
          ]
        }

        xps_tolag <- c(
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
        xps_nam <- paste0(xps_tolag, "_curr_xps")
        setnames(dtb, xps_tolag, xps_nam)

        # Prune & write synthpop to disk ----
        # del rn as they are reproducible
        # dtb[, ("LSOA11CD") := NULL] # needed for ICB weights

        if ("age100" %in% names(dtb)) {
          dtb[, age := NULL]
          setnames(dtb, "age100", "age")
        }

        setkey(dtb, pid, year) # Just in case
        setcolorder(dtb, c("pid", "year", "age", "sex", "dimd"))
        setindexv(dtb, c("year", "age", "sex", "dimd", "ethnicity"))
        if (design_$sim_prm$logs) {
          message("Writing synthpop to disk")
        }
        write_fst(dtb, filename_$synthpop, 90) # 100 is too slow
        return(invisible(NULL))
      },

      #  LAD_contribution ----
      # used to calculate the LAD contribution to the population size. Useful for ICBs that include portions of some LADs.
      # @return a data.table with cols LAD17CD and spTOons (spTOons := sppop / i.onspop).
      LAD_contribution = function() {
        lsoas_ <- private$get_unique_LSOAs(private$design)
        file <- "./inputs/pop_estimates_lsoa/LSOA_mid_year_population_estimates"
        popLSOA <- read_parquet_dt(
          file,
          filter = arrow::Expression$field_ref("year") ==
            private$design$sim_prm$init_year &
            arrow_in("LSOA11CD", lsoas_)
        )
        popLSOA[, `99` := `99` + `100`]
        popLSOA[,
          c(
            paste0(0:(private$design$sim_prm$ageL - 1L)),
            c(paste0((private$design$sim_prm$ageH + 1L):100))
          ) := NULL
        ]
        popLSOA <- melt(
          popLSOA,
          grep("^[0-9]", names(popLSOA), value = TRUE, invert = TRUE),
          variable.name = "age",
          value.name = "population_size",
          variable.factor = FALSE
        )
        popLSOA[, age := as.integer(age)]
        popLSOA[, year := as.integer(year)]
        indx_hlp <-
          read_fst(
            "./inputs/pop_estimates_lsoa/lsoa_to_locality_indx.fst",
            as.data.table = TRUE
          )
        popLSOA[indx_hlp, on = "LSOA11CD", LAD17CD := i.LAD17CD]

        strata <- c("year", "age", "sex")
        if (private$design$sim_prm$calibrate_to_pop_projections_by_LAD) {
          strata <- c(strata, "LAD17CD")
        }
        popLAD <- private$get_pop_size(
          private$design$sim_prm$calibrate_to_pop_projections_by_LAD
        )

        minage <- max(popLSOA[, min(age)], popLAD[, min(age)])
        minyear <- popLSOA[, min(year)]

        popLSOA <- popLSOA[age >= minage & year == minyear]
        popLAD <- popLAD[age >= minage & year == minyear]
        sppoplad <- popLSOA[, .(sppop = sum(population_size)), keyby = LAD17CD]
        onspoplad <- popLAD[, .(onspop = sum(pops)), keyby = LAD17CD]
        sppoplad[onspoplad, on = "LAD17CD", spTOons := sppop / i.onspop][,
          sppop := NULL
        ]

        return(invisible(sppoplad))
      },

      # Load a synthpop file from disk in full or in chunks.
      #  get_synthpop ----
      # @param exclude_cols A vector specifying columns to be excluded from the synthetic population data.
      #
      # @return An invisible data table containing the processed synthetic population data.
      #
      # @details
      # The function reads synthetic population data stored in an fst file, filters the data based on the specified criteria, ensures uniqueness of individual identifiers (pid), generates population weights, and makes additional adjustments to the data.
      #
      get_synthpop = function(exclude_cols = c()) {
        mm_synthpop <- metadata_fst(private$filename$synthpop)
        mm_synthpop <- setdiff(mm_synthpop$columnNames, exclude_cols)

        # Read synthpop

        dtb <- read_fst(
          private$filename$synthpop,
          columns = mm_synthpop,
          as.data.table = TRUE
        )
        dtb <- dtb[
          between(
            year,
            private$design$sim_prm$init_year - private$design$sim_prm$maxlag,
            private$design$sim_prm$init_year +
              private$design$sim_prm$sim_horizon_fromGUI
          ) &
            between(
              age,
              private$design$sim_prm$ageL - private$design$sim_prm$maxlag,
              private$design$sim_prm$ageH
            )
        ]

        # Ensure pid does not overlap for files from different mc
        new_n <-
          it <- as.integer(ceiling(
            self$mc %% private$design$sim_prm$num_chunks
          ))
        if ((max(dtb$pid) + it * 1e8) >= .Machine$integer.max) {
          stop("pid larger than int32 limit.")
        }
        dtb[, pid := as.integer(pid + it * 1e8)]

        dtb[, pid_mrk := mk_new_simulant_markers(pid)] # TODO Do I need this?

        dtb[,
          smok_packyrs_curr_xps := as.integer(round(
            smok_cig_curr_xps * smok_dur_curr_xps / 20
          ))
        ]
        set(dtb, NULL, "all_cause_mrtl", 0L)
        set(dtb, NULL, "cms_score", 0) # CMS score of diagnosed conditions
        set(dtb, NULL, "cms_count", 0L) # Count of diagnosed CMS conditions

        setkey(dtb, pid, year)
        # dtb[, dead := identify_longdead(all_cause_mrtl, pid_mrk)]
        # dtb[, ncc := clamp(
        #   ncc - (chd_prvl > 0) - (stroke_prvl > 0) -
        #     (poststroke_dementia_prvl > 0) -
        #     (htn_prvl > 0) - (t2dm_prvl > 0) - (af_prvl > 0) -
        #     (copd_prvl > 0) - (lung_ca_prvl > 0) -
        #     (colon_ca_prvl > 0) -
        #     (breast_ca_prvl > 0),
        #   0L,
        #   10L
        # )]
        # to be added back in the qaly fn. Otherwise when I prevent disease
        # the ncc does not decrease.

        invisible(dtb)
      },

      # get_pop_size ----
      # Get Population Size By Group
      #
      # @return An invisible modified data table with population size by group.
      #
      # @details
      # The function reads population estimates and projections data from an fst
      # file, filters the data based on specified criteria, and stratifies the
      # population size by prespicified groups. If the time horizon is longer
      # than the available projections it repeats the last year of the
      # projection. The max age is max age+. I.e 99 is 99+ and includes all ages
      # above 99
      get_pop_size = function(
        return_pop_projections_by_LAD = private$design$sim_prm$calibrate_to_pop_projections_by_LAD
      ) {
        minage <- private$design$sim_prm$ageL - private$design$sim_prm$maxlag
        minyear <- private$design$sim_prm$init_year -
          private$design$sim_prm$maxlag
        maxyear <- private$design$sim_prm$init_year +
          private$design$sim_prm$sim_horizon_fromGUI

        if (return_pop_projections_by_LAD) {
          strata <- c("year", "age", "sex", "LAD17CD")
          lads <- private$get_unique_LADs(private$design)
          tt <- read_fst(
            "./inputs/pop_projections/lad17_proj.fst",
            as.data.table = TRUE
          )[
            LAD17CD %in% lads & age >= minage & between(year, minyear, maxyear),
            .(pops = sum(pops)),
            keyby = strata
          ]
        } else {
          # Use national pop estimates and projections
          strata <- c("year", "age", "sex")
          tt <- rbind(
            read_fst(
              "inputs/pop_estimates_lsoa/national_pop_est.fst",
              as.data.table = TRUE
            ),
            read_fst(
              "inputs/pop_projections/national_proj.fst",
              as.data.table = TRUE
            )
          )[
            between(age, minage, private$design$sim_prm$ageH) &
              between(year, minyear, maxyear),
            .(pops = sum(pops)),
            keyby = strata
          ]
        }
        # include all ages > max age to the max age
        # TODO consider make it a function argument
        tt[
          age > private$design$sim_prm$ageH,
          age := private$design$sim_prm$ageH
        ]
        tt <- tt[, .(pops = sum(pops)), keyby = strata]

        # Repeat last year if necessary
        if (maxyear > tt[, max(year)]) {
          ttt <- tt[year == max(year), ]
          ttt <- clone_dt(ttt, maxyear - tt[, max(year)])
          ttt[, `:=`(year = year + .id, .id = NULL)]
          tt <- rbind(tt, ttt)
          setkeyv(tt, strata)
        }

        # WiP (based on the fact that in minyear the synthpop size is private$design$sim_prm$n by design.
        # Also in init_year the synthpop size is private$design$sim_prm$n by design for ages between ageL an ageH.)
        # tt <- tt[LAD17CD %in% lads &
        #    between(age, private$design$sim_prm$ageL, private$design$sim_prm$ageH) &
        #    between(year, private$design$sim_prm$init_year, maxyear),
        #  .(pops = sum(pops)), keyby = strata] # Note ages < ageL and years before init_year are not included
        # absorb_dt(tt, tt[year == private$design$sim_prm$init_year, .(pops_init_year = sum(pops)), keyby = eval(strata[strata != "year"])])
        # tt[, wt := (tt[year == private$design$sim_prm$init_year, sum(pops)]/ private$design$sim_prm$n) * pops / pops_init_year]

        # The expectation of population size in the sample, 'expected_sample_pops' is
        # tt[, expected_sample_pops := private$design$sim_prm$n * pops / pops_init_year]
        # Note that some pops = 0 when calibrate_to_pop_projections_by_LAD == TRUE
        # ttt <- sp$pop[between(age, private$design$sim_prm$ageL, private$design$sim_prm$ageH) &
        #    between(year, private$design$sim_prm$init_year, maxyear),
        #     .N, keyby = strata]
        # absorb_dt(ttt, tt)
        # ttt[, sum(wt * N), keyby = year]
        # # So the weight 'wt' for the expected population is
        # tt[, wt := pops / expected_sample_pops]
        # tt[, sum(wt), keyby = year]

        # sp$pop[, c("pops", "wt", "annual_pops", "mpops") := NULL]
        # absorb_dt(sp$pop, tt)

        # sp$pop[between(age, private$design$sim_prm$ageL, private$design$sim_prm$ageH) &
        #    between(year, private$design$sim_prm$init_year, maxyear),
        #     .(unique(annual_pops), sum(wt, T), sum(wt_immrtl, T)), keyby = year]

        # WiP2
        # sp$pop[, c("pops", "wt", "annual_pops", "mpops") := NULL]
        # tt[, `:=` (wt = pops / sum(pops), annual_pops = sum(pops)), by = year]
        # absorb_dt(sp$pop, tt)
        # sp$pop[, mpops := (wt * unique(annual_pops) / sum(wt)) / private$design$sim_prm$num_chunks, by = year] # mpops are the new weights
        # sp$pop[, .(unique(annual_pops), sum(mpops), sum(wt_immrtl)), by = year]
        # sp$pop[, summary(mpops/mpops2)]
        # setnames(sp$pop, "mpops", "mpops2")
        # sp$pop[, mpops2 := NULL]
        # ttt <- sp$pop[, .N, keyby = strata]
        # absorb_dt(tt, ttt)
        # setnafill(tt, "c", 0L, cols = "N")
        # tt[, wt_immrtl := pops / (N * private$design$sim_prm$num_chunks)]

        invisible(tt)
      }
    )
  )
