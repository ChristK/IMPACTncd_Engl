## IMPACTncdEngl is an implementation of the IMPACTncd framework, developed by Chris
## Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz.
##
## Copyright (C) 2018-2020 University of Liverpool, Chris Kypridemos
##
## IMPACTncdEngl is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version. This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details. You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/> or write
## to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
## Boston, MA 02110-1301 USA.



# From
# https://stackoverflow.com/questions/33424233/how-do-i-tell-an-r6-class-what-to-do-with-square-brackets
# Allows data.table syntax to the R6class object directly. Assumes it has a
# field 'pop' that is a data.table
#' @export
`[.SynthPop` <- function(x, ...) x$pop[...]

#' R6 Class representing a synthetic population
#'
#' @description
#' A synthpop has a `pop` field that contains the life course of simulants in a
#' `data.table`.
#'
#' @details
#' To be completed...
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
      #'   be used when multiple synthetic populations getting aggregated.It
      #'   ensures correct seeds for the RNGs during the simulation for the RRs
      #'   and the lags.
      mc_aggr = NA,

      #' @field metadata Metadata of the synthpop.
      metadata = NA,

      #' @field pop The data.table that contains the life-course of simulants.
      #'   If the file exists, it is loaded from disk. If it doesn't, it is
      #'   first generated, then saved to disk, and then loaded from disk.
      pop = NA,


      #' @description Create a new SynthPop object.
      #' If a synthpop file in \code{\link[fst]{fst-package}} format already
      #' exists, then the synthpop is loaded from there. Otherwise it is
      #' generated from scratch and then saved as `filename` in
      #' \code{\link[fst]{fst-package}} format. Two additional files are saved
      #' for each 'synthpop'. A metadata file, and an index file.
      #' @param mc_ The Monte Carlo iteration of the synthetic population. Each
      #'   integer generates a unique synthetic population. If `mc = 0` an
      #'   object with an empty synthpop is initiated.
      #' @param design_ A \code{\link[IMPACTncdEngl]{Design}} object.
      #' @param synthpop_dir_ The directory where 'SynthPop' objects are stored.
      #'   The synthpop file in \code{\link[fst]{fst-package}} format. If
      #'   `filename` already exists, then the synthpop is loaded from there.
      #'   Otherwise it is generated from scratch and then saved as `filename`
      #'   in \code{\link[fst]{fst-package}} format. Two additional files are
      #'   saved for each 'synthpop'. A metadata file, and an index file.
      #' @return A new `SynthPop` object.
      #' @examples
      #' design <- Design$new("./validation/design_for_trends_validation.yaml")
      #' POP$write_synthpop(1:6, design)
      #' POP <- SynthPop$new(4L, design)
      #' POP$print()
      #' POP$count_synthpop()
      #'
      #' POP$delete_synthpop(1L)
      #' POP$delete_synthpop(5:6)
      #' POP$get_filename()
      initialize = function(mc_, design_) {
        stopifnot(length(mc_) == 1L, is.numeric(mc_), ceiling(mc_) >= 0L)
        stopifnot("Design" %in% class(design_))

        mc_ <- as.integer(ceiling(mc_))
        # Create synthpop_dir if it doesn't exists
        if (!dir.exists(design_$sim_prm$synthpop_dir)) {
          dir.create(design_$sim_prm$synthpop_dir, recursive = TRUE)
          message(paste0("Directory ", design_$sim_prm$synthpop_dir,
                         " was created"))
        }

        # get unique lsoas
        # lsoas <- private$get_unique_LSOAs(design_)

        private$checksum <- private$gen_checksum(design_)

        self$mc <- mc_
        self$mc_aggr <- as.integer(ceiling(mc_ / design_$sim_prm$n_synthpop_aggregation))

        private$design <- design_
        private$synthpop_dir <- design_$sim_prm$synthpop_dir

        if (mc_ > 0) {
          private$filename <- private$gen_synthpop_filename(mc_,
                                                            private$checksum,
                                                            design_)
          # logic for the synthpop load
          files_exist <- sapply(private$filename, file.exists)
          if (all(!files_exist)) {
            # No files exist. Create the synthpop and store the file on disk (no
            # parallelism)
            private$gen_synthpop(mc_,
                                 private$filename,
                                 design_)

          } else if (file.exists(private$filename$metafile) &&
                     !all(files_exist)) {
            # Metafile exists but not all three files. It means that most likely
            # a generate_synthpop() is still running. So the function waits
            # until the file is created before it proceeds to load it. Note that
            # if this is not the case then the loop is infinite!!!
            while (!all(sapply(private$filename, file.exists))) {
              Sys.sleep(5)
              if (design_$sim_prm$logs) {
                message("Metafile exists without a synthpop file. Check for synthpop after 5 sec.")
              }
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

          } else if (!file.exists(private$filename$metafile) &&
                     !all(files_exist)) {
            # Metafile doesn't exist but some other files exist. In this case
            # delete everything and start from scratch
            self$delete_incomplete_synthpop()
            private$gen_synthpop(mc_,
                                 private$filename,
                                 design_)
          }
          # No need to provision for case when all file present. The following
          # lines handle this case anyway

          self$pop <- private$get_synthpop()
          self$metadata <- yaml::read_yaml(private$filename$metafile)

          if (design_$sim_prm$logs) self$print()
        }
        invisible(self)
      },

      #' @description
      #' Updates the Design object that is stored in the SynthPop object.
      #' @param design_ A design object with the simulation parameters.
      #' @return The invisible self for chaining.

      update_design = function(design_ = design) {
        if (!inherits(design, "Design"))
          stop("Argument design_ needs to be a Design object.")

        private$design <- design
        invisible(self)
      },

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
        if (missing(mc_)) stop("Use mc_ = NULL if you want to delete all synthpop files.")
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

        } else if (length(mc_) == 1L &&
                   is.numeric(mc_) && ceiling(mc_) > 0L) {
          fl <- unlist(
            private$gen_synthpop_filename(mc_,
              private$checksum,
              private$design))
          file.remove(fl)

        } else if (length(mc_) > 1L &&
                   all(is.numeric(mc_)) && all(ceiling(mc_) > 0L)) {
          fl <-
            lapply(mc_,
                   private$gen_synthpop_filename,
                   private$checksum,
                   private$design)
          fl <- unlist(fl)
          file.remove(fl)

        } else
          message("mc_ need to be NULL or numeric. Nothing was deleted.")

        return(invisible(self))
      },

      #' @description
      #' Check that every synthpop file has a metafile and an index file. Delete
      #' any orphan files.
      #' @param check_checksum If  `TRUE` only delete incomplete group files
      #'   with the same checksum as the synthpop.
      #' @return The invisible `SynthPop` object.
      delete_incomplete_synthpop =
        function(check_checksum = TRUE) {
          if (check_checksum) {
            f1 <- paste0("^synthpop_", private$checksum , ".*\\.fst$")
            f2 <- paste0("^synthpop_", private$checksum , ".*_meta\\.yaml$")
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

      #' @description
      #' Check the integrity of (and optionally delete) .fst files by checking
      #' their metadata are readable.
      #' @param remove_malformed If `TRUE`, delete all malformed .fst files and
      #'   their associated files.
      #' @param check_checksum If  `TRUE` only check files with the same
      #'   checksum as the synthpop.
      #' @return The invisible `SynthPop` object.
      check_integridy =
        function(remove_malformed = FALSE,
          check_checksum = TRUE) {
          if (check_checksum) {
            pat <- paste0("^synthpop_", private$checksum , ".*\\.fst$")
          } else {
            pat <- "^synthpop_.*\\.fst$"
          }

          files <-
            list.files(private$synthpop_dir,
              pat,
              full.names = TRUE)
          if (length(files) > 0L) {
            malformed <- sapply(files, function(x) {
              out <- try(metadata_fst(x), silent = TRUE)
              out <- inherits(out, "try-error")
              out
            }, USE.NAMES = FALSE)


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
                tr <-  paste0(to_remove, ".fst")
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




      #' @description
      #' Count the synthpop files in a directory. It includes files without
      #' metafiles and index files.
      #' @return The invisible `SynthPop` object.
      count_synthpop =
        function() {
          out <- list()
          # folder size
          files <-
            list.files(private$synthpop_dir, full.names = TRUE)
          if (length(files) > 0L) {
            vect_size <- sapply(files, file.size)
            out$`synthpop folder size (Gb)` <-
              signif(sum(vect_size) / (1024 ^ 3), 4) # Gb

            # synthpops with same checksum
            files <- list.files(private$synthpop_dir,
              paste0("^synthpop_", private$checksum , ".*\\.fst$"))

            out$`synthpop meta files with same checksum` <-
              length(list.files(
                private$synthpop_dir,
                paste0("^synthpop_", private$checksum , ".*_meta\\.yaml$")
              ))

            # synthpops with any checksum
            files <-
              list.files(private$synthpop_dir, "^synthpop_.*\\.fst$")

            out$`synthpop meta files with any checksum` <-
              length(list.files(private$synthpop_dir, "^synthpop_.*_meta\\.yaml$"))


            cat(paste0(names(out), ": ", out, "\n"))
          } else {
            # if length(files) == 0L
            cat("no files found.")
          }
          return(invisible(self))
        },

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
            all      = print(private$filename),
            synthpop = print(private$filename[["synthpop"]]),
            metafile = print(private$filename[["metafile"]])
          )
        }
        invisible(self)
      },

      #' @description
      #' Get the synthpop design.
      #' @return The invisible `SynthPop` object.
      get_design = function() {
        # print(private$design)
        # invisible(self)
        private$design
      },

      #' @description
      #' Get the synthpop dir.
      #' @return The invisible `SynthPop` object.
      get_dir = function() {
        print(private$synthpop_dir)
        invisible(self)
      },

      #' @description
      #' Generate synthpop sociodemographics, random sample of the population.
      #' @param design_ A Design object,
      #' @param month April or July are accepted. Use July for mid-year
      #'   population estimates.
      #' @return An invisible `data.table` with sociodemographic information.
      gen_synthpop_demog =
        function(design_, month = "April") {
          stopifnot("Argument month need to be April or July" = month %in% c("April", "July"))
          # Use month = July for mid-year
          lsoas_ <- private$get_unique_LSOAs(design_)
          # load dt
          if (month == "April") {
            file <- "./inputs/pop_estimates_lsoa/LSOA_1st_April_population_estimates.fst"
          } else {
            file <- "./inputs/pop_estimates_lsoa/LSOA_mid_year_population_estimates.fst"
          }
          dt_meta <- metadata_fst(file)
          stopifnot("Population size file need to be keyed by year" =
                      identical("year", dt_meta$keys[1]))
          file_indx <- read_fst(file, as.data.table = TRUE, columns = "year"
          )[, .(from = min(.I), to = max(.I)), keyby = "year"][year == design_$sim_prm$init_year]
          dt <-
            read_fst(file, from = file_indx$from, to = file_indx$to,
                     as.data.table = TRUE)[LSOA11CD %in% lsoas_]
          # delete unwanted ages
          dt[, c(paste0(0:(design_$sim_prm$ageL - 1L)), c(paste0((
            design_$sim_prm$ageH + 1L
          ):100))) := NULL]

          dt <-
            melt(
              dt,
              grep("^[0-9]", names(dt), value = TRUE, invert = TRUE),
              variable.name = "age",
              value.name = "population_size",
              variable.factor = FALSE
            )
          dt[, age := as.integer(age)]
          dt[, year := as.integer(year)]
          dt[, population_size := population_size / sum(population_size)]

          # load ethnicity proportions by lsoa
          file <- "./inputs/pop_estimates_lsoa/ethn2011_pct.fst"
          dt_meta <- metadata_fst(file)
          stopifnot("Ethnicity file need to be keyed by LSOA" = identical("LSOA11CD", dt_meta$keys[1]))
          file_indx <- read_fst(file, as.data.table = TRUE, columns = "LSOA11CD"
          )[, .(from = min(.I), to = max(.I)), keyby = "LSOA11CD"][LSOA11CD %in% lsoas_, .("from" = min(from), "to" = max(to))]
          ethn <- read_fst(
            file,
            from = file_indx$from,
            to = file_indx$to,
            as.data.table = TRUE
          )[LSOA11CD %in% lsoas_]

          absorb_dt(dt, ethn)
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

          for (j in .ethn_nam)
            set(dt, NULL, j, dt[, get(j) * population_size])
          dt[, population_size := NULL]
          dt <-
            melt(
              dt,
              measure.vars = .ethn_nam,
              variable.name = "ethnicity",
              variable.factor = TRUE,
              value.name = "prbl"
            )

          # I do not explicitly set.seed because I do so in the gen_synthpop()
          dtinit <- dt[sample(.N, design_$sim_prm$n, TRUE, prbl)]

          # Generate the cohorts of 30 year old to enter every year
          # as sim progress these will become 30 yo
          # no population growth here as I will calibrate to dt
          # projections and it fluctuates at +-2% anyways.

          # tt1 <-
          #   read_fst(
          #     "./inputs/pop_estimates_lsoa/national_mid_year_population_estimates.fst",
          #     as.data.table = TRUE
          #   )[age == design_$sim_prm$ageL & year >= design_$sim_prm$init_year]
          # tt2 <-
          #   read_fst("./inputs/pop_projections/national_proj.fst",
          #            as.data.table = TRUE)[age == design_$sim_prm$ageL &
          # year <= (design_$sim_prm$init_year + design_$sim_prm$sim_horizon_max)]
          # doubleyrs <- intersect(unique(tt1$year), unique(tt2$year))
          # if (length(doubleyrs)) {
          #   tt2 <- tt2[!year %in% doubleyrs]
          # }
          # tt <- rbind(tt1, tt2)
          # setkey(tt, year)
          # tt[, growth := shift(pops)/pops, keyby = sex]

          if (design_$sim_prm$logs)
            message("Generate the cohorts of ", design_$sim_prm$ageL," year old")

          dt <- dt[age == design_$sim_prm$ageL]
          siz <- dtinit[age == design_$sim_prm$ageL, .N]
          dtfut <- dt[sample(.N, siz * design_$sim_prm$sim_horizon_max, TRUE, prbl)]
          dtfut[, age := age - rep(1:design_$sim_prm$sim_horizon_max, siz)]

          dt <- rbind(dtfut, dtinit)
          dt[, prbl := NULL]

          indx_hlp <-
            read_fst("./inputs/pop_estimates_lsoa/lsoa_to_locality_indx.fst",
                     as.data.table = TRUE)

          dt[indx_hlp, on = "LSOA11CD", `:=` (
            # tds = i.tds,
            # tds_quintile = i.tds_quintile,
            # imd = i.imd,
            qimd = i.qimd,
            dimd = i.dimd,
            sha = i.SHA11NM,
            CCG17CDH = CCG17CDH
          )]

          return(invisible(dt))
        },

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

        if (Sys.info()["sysname"] == "Windows") {
          cl <-
            makeCluster(private$design$sim_prm$clusternumber) # used for clustering. Windows compatible
          registerDoParallel(cl)
        } else {
          registerDoParallel(private$design$sim_prm$clusternumber) # used for forking. Only Linux/OSX compatible
        }

        foreach(
          mc_iter = mc_,
          .inorder = FALSE,
          .verbose = private$design$sim_prm$logs,
          .packages = c(
            "R6",
            "gamlss.dist",
            # For distr in prevalence.R
            "dqrng",
            "qs",
            "fst",
            "CKutils",
            "IMPACTncdEngl",
            "data.table"
          ),
          .export = NULL,
          .noexport = NULL # c("time_mark")
        ) %dopar% {
          data.table::setDTthreads(private$design$sim_prm$n_cpus)
          fst::threads_fst(private$design$sim_prm$n_cpus)
          filename <-
            private$gen_synthpop_filename(mc_iter,
                                          private$checksum,
                                          private$design)

          # logic for the synthpop load
          files_exist <- sapply(filename, file.exists)
          if (all(!files_exist)) {
            # No files exist. Create the synthpop and store
            # the file on disk
            private$gen_synthpop(mc_iter,
                                 filename,
                                 private$design)

          } else if (file.exists(filename$metafile) &&
                     !all(files_exist)) {
            # Metafile exists but not all three files. It means
            # that most likely a generate_synthpop() is still running. So the
            # function waits until the file is created before it proceeds to
            # load it. Note that if this is not the case then the loop is
            # infinite!!!
            while (!all(sapply(filename, file.exists)))
              Sys.sleep(5)

            # Ensure the file write is complete (size stable)
            sz1 <- file.size(filename$synthpop)
            Sys.sleep(3)
            sz2 <- file.size(filename$synthpop)
            while (sz1 != sz2) {
              sz1 <- file.size(filename$synthpop)
              Sys.sleep(3)
              sz2 <- file.size(filename$synthpop)
            }

          } else if (!file.exists(filename$metafile) &&
                     !all(files_exist)) {
            # Metafile doesn't exist but some other files exist. In this case
            # delete everything and start from scratch
            self$delete_incomplete_synthpop()
            private$gen_synthpop(mc_iter,
                                 filename,
                                 private$design)
          }
          # No need to provision for case when all files present.

          return(NULL)
        }
        if (exists("cl"))
          stopCluster(cl)

        invisible(self)
      },

      #' @description
      #' Prints the synthpop object metadata.
      #' @return The invisible `SynthPop` object.
      print = function() {
        print(c(
          "path" = ifelse(self$mc == 0L,
                          "Not relevant because mc = 0L",
                          private$filename$synthpop),
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

      # Special deep copy for data.table. Use POP$clone(deep = TRUE) to
      # dispatch. Otherwise a reference is created
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
      get_unique_LSOAs = function(design_) {
        indx_hlp <-
          read_fst("./inputs/pop_estimates_lsoa/lsoa_to_locality_indx.fst",
                   as.data.table = TRUE, columns = c("LSOA11CD", "LAD17NM", "RGN11NM"))

        if ("England" %in% design_$sim_prm$locality) {
          lsoas <- indx_hlp[, unique(LSOA11CD)] # national
        } else {
          lsoas <-
            indx_hlp[LAD17NM %in% design_$sim_prm$locality |
                       RGN11NM %in% design_$sim_prm$locality, unique(LSOA11CD)]
        }
        return(sort(lsoas))
      },

      # get all unique LADs included in locality vector.
      get_unique_LADs = function(design_) {
        indx_hlp <-
          read_fst("./inputs/pop_estimates_lsoa/lsoa_to_locality_indx.fst",
                   as.data.table = TRUE, columns = c("LAD17CD", "LAD17NM", "RGN11NM"))

        if ("England" %in% design_$sim_prm$locality) {
          lads <- indx_hlp[, unique(LAD17CD)] # national
        } else {
          lads <-
            indx_hlp[LAD17NM %in% design_$sim_prm$locality |
                       RGN11NM %in% design_$sim_prm$locality, unique(LAD17CD)]
        }
        return(sort(lads))
      },

      # get a smaller design list only with characteristics that are important
      # for synthpop creation and define the uniqueness of the object. I.e. if
      # these parameters are different the synthpop has to have different
      # filename and vice-versa
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
      gen_checksum =
        function(design_) {
          # get a md5 checksum based on function arguments
          # First get function call arguments
          fcall <- private$get_unique_characteristics(design_)

          lsoas_ <- private$get_unique_LSOAs(design_)

          locality_years_age_id <-
            digest::digest(paste(lsoas_, fcall, sep = ",", collapse = ","),
                           serialize = FALSE)
          return(locality_years_age_id)
        },

      # gen synthpop filename for the given set of inputs
      gen_synthpop_filename =
        function(mc_,
                 checksum_,
                 design_) {
          return(
            list(
              "synthpop" = normalizePath(
                paste0(design_$sim_prm$synthpop_dir,
                       "/synthpop_",
                       checksum_,
                       "_",
                       mc_,
                       ".fst"),
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

      del_incomplete = function(filename_) {
        if (file.exists(filename_$metafile) &&
            (!file.exists(filename_$synthpop)
            )) {
          suppressWarnings(sapply(filename_, file.remove))
        }
      },

      gen_synthpop = # returns NULL. Writes synthpop on disk
        function(mc_,
                 filename_,
                 design_) {
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
            yaml::write_yaml(private$get_unique_characteristics(design_),
                             filename_$metafile)
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

            dt <- self$gen_synthpop_demog(design_, month = "April")

            # NOTE!! from now on year in the short form i.e. 13 not 2013
            dt[, `:=`(pid  = .I)]
            new_n <- nrow(dt)


            # Generate correlated ranks for the individuals ----
            if (design_$sim_prm$logs) message("Generate correlated ranks for the
                                              individuals")

            cm_mean <- as.matrix(
              read_fst(
                "./inputs/exposure_distributions/exposure_corr_mean.fst",
                as.data.table = TRUE
              ),
              rownames = "rn"
            )

            tr <- which(colnames(cm_mean) %in%
                          c("af_r", "ckd_r", "famcvd_r", "dm_r", "dm_dgn_r"))
            cm_mean <- cm_mean[-tr, -tr]

            rank_mtx <- generate_corr_unifs(new_n, cm_mean)
            if (design_$sim_prm$logs) message("generate correlated uniforms")

            # Restrict the range of some RNs to avoid unrealistic exposures
            # This scaling does not affect correlations
            # /0.999 because I multiplied all the columns below
            rank_mtx <- rank_mtx * 0.999
            rank_mtx[, "frtpor_r"] <- rank_mtx[, "frtpor_r"] * 0.99 / 0.999
            rank_mtx[, "vegpor_r"] <- rank_mtx[, "vegpor_r"] * 0.93 / 0.999
            rank_mtx[, "smok_cig_ex_r"] <-
              rank_mtx[, "smok_cig_ex_r"] * 0.99 / 0.999
            rank_mtx[, "totalwu_r"] <- rank_mtx[, "totalwu_r"] * 0.99 / 0.999
            rank_mtx[, "smok_quit_yrs_r"] <-
              rank_mtx[, "smok_quit_yrs_r"] * 0.99 / 0.999
            rank_mtx[, "smok_dur_ex_r"] <-
              rank_mtx[, "smok_dur_ex_r"] * 0.99 / 0.999
            rank_mtx[, "smok_dur_curr_r"] <-
              rank_mtx[, "smok_dur_curr_r"] * 0.88 / 0.999
            rank_mtx[, "bmi_r"] <-
              rank_mtx[, "bmi_r"] * 0.95 / 0.999
            # sum((cor(rank_mtx) - cm_mean) ^ 2)
            if (design_$sim_prm$logs) message("correlated ranks matrix to data.table")

            rank_mtx <- data.table(rank_mtx)

            # NOTE rankstat_* is unaffected by the RW. Stay constant through the lifecourse
            dt[, c(
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
            ) := rank_mtx]

            rm(rank_mtx)

            # add non-correlated RNs
            rank_cols <-
              c(
                "rankstat_ncc",
                "rankstat_ca_history",
                "rankstat_famlungca",
                "rankstat_pa_dur",
                "rankstat_pa_met",
                "rankstat_bpmed_adherence",
                "rankstat_statin_adherence"
              )


            for (nam in rank_cols)
              set(dt, NULL, nam, dqrunif(new_n)) # NOTE do not replace with generate_rns function.

            # Generate education (exception as it remains stable through lifecourse) ----
            if (design_$sim_prm$logs) message("Generate education")

            tbl <-
              read_fst("./inputs/exposure_distributions/education_table.fst",
                       as.data.table = TRUE)
            nam <- intersect(names(dt), names(tbl))
            # logic necessary for new cohorts entering the simulation that currently age < 30
            # These will have the same distribution as if 30 years old
            tt <- tbl[age == min(age)]
            tt <- clone_dt(tt, design_$sim_prm$sim_horizon_max)
            tt[, age := age - .id] # as the sim progress these will become 30 yo
            # increase population by 0.5% every year
            tt[, .id := NULL]
            tbl <- rbind(tt, tbl)
            dt[tbl, education := (rank_education > ed1) + (rank_education > ed2) +
                 (rank_education > ed3) + (rank_education > ed4) +
                 (rank_education > ed5) + (rank_education > ed6) + 1L,
               on = nam]
            dt[, education := factor(
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
            )]
            dt[, rank_education := NULL]

            # Project forward for simulation and back project for lags  ----
            if (design_$sim_prm$logs) message("Project forward and back project")

            dt <-
              clone_dt(dt,
                       design_$sim_prm$sim_horizon_max +
                         design_$sim_prm$maxlag + 1L)

            dt[.id <= design_$sim_prm$maxlag, `:=` (age  = age  - .id,
                                                    year = year - .id)]
            dt[.id > design_$sim_prm$maxlag, `:=` (
              age  = age  + .id - design_$sim_prm$maxlag - 1L,
              year = year + .id - design_$sim_prm$maxlag - 1L
            )]
            # dt <-
            #   dt[between(age, design_$sim_prm$ageL - design_$sim_prm$maxlag, design_$sim_prm$ageH)]
            # delete unnecessary ages
            del_dt_rows(
              dt,
              !between(
                dt$age,
                design_$sim_prm$ageL - design_$sim_prm$maxlag,
                design_$sim_prm$ageH
              ),
              environment()
            )

            dt[, `:=` (.id = NULL)]

            if (max(dt$age) > 90L) {
              dt[, age100 := age]
              dt[age > 90L, age := 90L]
            }

            # to_agegrp(dt, 20L, 85L, "age", "agegrp20", to_factor = TRUE)
            # to_agegrp(dt, 10L, 85L, "age", "agegrp10", to_factor = TRUE)
            # to_agegrp(dt,  5L, 85L, "age", "agegrp5" , to_factor = TRUE)

            # Simulate exposures -----

            # Random walk for ranks ----
            if (design_$sim_prm$logs) message("Random walk for ranks")

            setkeyv(dt, c("pid", "year"))
            setindexv(dt, c("year", "age", "sex", "sha", "qimd", "ethnicity"))

            dt[, pid_mrk := mk_new_simulant_markers(pid)]

            dt[, lapply(.SD,
                        fscramble_trajectories,
                        pid_mrk,
                        design_$sim_prm$jumpiness),
               .SDcols = patterns("^rank_")]
            # ggplot2::qplot(year, rank_income, data = dt[pid %in% sample(1e5, 1)], ylim = c(0,1))

            # Generate income ----
            if (design_$sim_prm$logs) message("Generate income")

            tbl <-
              read_fst("./inputs/exposure_distributions/income_table.fst", as.data.table = TRUE)
            nam <- intersect(names(dt), names(tbl))
            dt[tbl, income := (rank_income > inc1) + (rank_income > inc2) +
                 (rank_income > inc3) + (rank_income > inc4) + 1L,
               on = nam]
            dt[, income := factor(
              income,
              levels = 1:5,
              labels = c("1 Highest", "2", "3", "4", "5 Lowest")
            )]
            dt[, rank_income := NULL]

            # Generate active days ----
            if (design_$sim_prm$logs) message("Generate active days")

            tbl <-
              read_fst("./inputs/exposure_distributions/active_days_table.fst",
                       as.data.table = TRUE)
            nam <- intersect(names(dt), names(tbl))
            dt[tbl, active_days := (rank_pa > pa0) + (rank_pa > pa1) + (rank_pa > pa2) +
                 (rank_pa > pa3) + (rank_pa > pa4) + (rank_pa > pa5) + (rank_pa > pa6),
               on = nam]
            dt[, rank_pa := NULL]
            # dt[, active_days := factor(active_days, levels = 0:7)]

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
            dt[, met := as.integer(floor(active_days * (3L + qbinom(rankstat_pa_met, 8, 3/11)) *
                 (30 + qexp(rankstat_pa_dur, 1/7)) / 100))] # TODO make data driven
            dt[, c("rankstat_pa_met", "rankstat_pa_dur") := NULL]

            # Generate fruit consumption (ZISICHEL) ----
            if (design_$sim_prm$logs) message("Generate fruit consumption")

            tbl <-
              read_fst("./inputs/exposure_distributions/frtpor_table.fst", as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(dt), names(tbl)))
            absorb_dt(dt, tbl)
            dt[, fruit :=
                 my_qZISICHEL(rank_fruit,
                              mu, sigma, nu, tau, n_cpu = design_$sim_prm$n_cpu) * 80L]  # g/d
            dt[, (col_nam) := NULL]
            dt[, rank_fruit := NULL]

            # Generate veg consumption (DEL) ----
            if (design_$sim_prm$logs) message("Generate veg consumption")

            tbl <-
              read_fst("./inputs/exposure_distributions/vegpor_table.fst", as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(dt), names(tbl)))
            absorb_dt(dt, tbl)
            dt[, veg :=
                 my_qDEL(rank_veg, mu, sigma, nu, n_cpu = design_$sim_prm$n_cpu) * 80L]  # g/d
            dt[, (col_nam) := NULL]
            dt[, rank_veg := NULL]

            # Smoking simulation ----
            if (design_$sim_prm$logs) message("Smoking simulation")

            # Assign smok_status when pid_mrk == true (the first year an individual enters the simulation (with lags))
            tbl <-
              read_fst("./inputs/exposure_distributions/smok_status_table.fst",
                       as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(dt), names(tbl)))
            absorb_dt(dt, tbl)
            dt[, smok_status_ref := my_qMN4(rankstat_smok, mu, sigma, nu)] # for calibration
            dt[(pid_mrk), smok_status := smok_status_ref]
            dt[, (col_nam) := NULL]
            dt[, rankstat_smok := dqrunif(.N)] # this is now used for simsmoke. There shouldn't be correlated any more (colname hardcoded in C++ code).

            # Assign smok_quit_yrs when pid_mrk == true (the first year an
            # individual enters the simulation) I could use these estimates for
            # calibration but I need to calculate mortality first
            tbl <-
              read_fst("./inputs/exposure_distributions/smok_quit_yrs_table.fst",
                       as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(dt), names(tbl)))
            absorb_dt(dt, tbl)
            set(dt, NULL, "smok_quit_yrs", 0L)
            dt[(pid_mrk) &
                 smok_status %in% 2:3,
               smok_quit_yrs := my_qDPO(rankstat_smok_quit_yrs, mu, sigma)]
            dt[, rankstat_smok_quit_yrs := NULL]
            dt[, (col_nam) := NULL]

            # Assign smok_dur_ex when pid_mrk == true (the first year an individual enters the simulation)
            tbl <-
              read_fst("./inputs/exposure_distributions/smok_dur_ex_table.fst",
                       as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(dt), names(tbl)))
            absorb_dt(dt, tbl)
            set(dt, NULL, "smok_dur", 0L)
            dt[(pid_mrk) &
                 smok_status %in% 2:3, smok_dur := my_qDPO(rankstat_smok_dur_ex, mu, sigma)]
            dt[, rankstat_smok_dur_ex := NULL]
            dt[, (col_nam) := NULL]

            # Assign smok_dur_curr when pid_mrk == true (the first year an individual enters the simulation)
            tbl <-
              read_fst("./inputs/exposure_distributions/smok_dur_curr_table.fst",
                       as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(dt), names(tbl)))
            absorb_dt(dt, tbl)
            dt[(pid_mrk) &
                 smok_status == 4, smok_dur := as.integer(round(qNBI(rankstat_smok_dur_curr, mu, sigma)))]
            dt[, rankstat_smok_dur_curr := NULL]
            dt[, (col_nam) := NULL]

            # Ensure smoking histories start from age 12
            dt[age - smok_quit_yrs < 12L, smok_quit_yrs := age - 12L]
            dt[age - smok_dur < 12L, smok_dur := age - 12L]
            dt[age - smok_dur - smok_quit_yrs < 12L ,
               `:=`(smok_dur = as.integer(smok_dur / ((
                 smok_dur + smok_quit_yrs
               ) / (age - 12L))),
               smok_quit_yrs = as.integer(smok_quit_yrs / ((
                 smok_dur + smok_quit_yrs
               ) / (age - 12L))))]

            # Assign smok_incid probabilities
            tbl <-
              read_fst("./inputs/exposure_distributions/smok_incid_table.fst",
                       as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(dt), names(tbl)))
            absorb_dt(dt, tbl)
            setnames(dt, "mu", "prb_smok_incid")

            # Assign smok_cessation probabilities
            tbl <-
              read_fst("./inputs/exposure_distributions/smok_cess_table.fst",
                       as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(dt), names(tbl)))
            absorb_dt(dt, tbl)
            setnames(dt, "mu", "prb_smok_cess")

            # Handle smok_relapse probabilities
            tbl <-
              read_fst("./inputs/exposure_distributions/smok_relapse_table.fst",
                       as.data.table = TRUE)
            tbl <-
              dcast(tbl, sex + qimd ~ smok_quit_yrs, value.var = "pr")
            nam <- tbl[, paste0(sex, " ", qimd)]
            tbl <-
              as.matrix(tbl[, mget(paste0(1:15))], rownames = nam)

            simsmok(dt, tbl, design_$sim_prm$smoking_relapse_limit)
            # dt[!(pid_mrk), table(smok_status)]
            # dt[pid == 1, plot(year, smok_status, ylim = c(0, 4))]
            # dt[pid == 10, .(age, smok_status, smok_quit_yrs, smok_dur)]
            # dt[, sum(smok_status == 4)/.N, keyby = year]

            if (design_$sim_prm$simsmok_calibration) {
              # calculate dif between ref (multinom) and simsmok
              # I will further calibrate to better match HSE
              resample <-
                function(x, ...)
                  x[sample.int(length(x), ...)]
              obs <-
                dt[smok_status == 1L, .(nsa = .N), keyby = .(year, age, sex, qimd)] # ? add sha/ethn
              ref <-
                dt[smok_status_ref == 1L, .(nsr = .N), keyby = .(year, age, sex, qimd)] # ? add sha/ethn
              absorb_dt(ref, obs)
              setnafill(ref, "c", 0L, cols = "nsa")
              ref[, `:=`(dif = nsr - nsa,
                         nsr = NULL,
                         nsa = NULL)]
              # Further calibrate to match better with HSE
              ref[sex == "men" &
                    rbinom(.N, 1L, 0.9) == 1L, dif := dif - 2L]
              ref[sex == "men" &
                    rbinom(.N, 1L, 1 * clamp((year - min(year)) / 10, 0, 1)) == 1L, dif := dif + 2L]
              ref[sex == "women" &
                    rbinom(.N, 1L, 0.2) == 1L, dif := dif - 1L]
              ref[age < 49 &
                    rbinom(.N, 1L, 0.5) == 1L, dif := dif + 1L]
              ref[age < 49 &
                    sex == "men" &
                    qimd == "3" &
                    rbinom(.N, 1L, 0.5) == 1L, dif := dif - 1L]
              ref[age < 49 &
                    sex == "women" &
                    qimd %in% c("4", "5 least deprived") &
                    rbinom(.N, 1L, 0.4) == 1L, dif := dif - 1L]
              ref[between(age, 50, 69) &
                    rbinom(.N, 1L, 0.2) == 1L, dif := dif + 1L]
              ref[between(age, 50, 69) &
                    sex == "women" &
                    qimd == "5 least deprived" &
                    rbinom(.N, 1L, 0.4) == 1L, dif := dif + 1L]
              ref[between(age, 70, 89) &
                    rbinom(.N, 1L, 0.4) == 1L, dif := dif + 1L]
              ref[between(age, 70, 89) &
                    sex == "men" &
                    qimd %in% c("1 most deprived", "2") &
                    rbinom(.N, 1L, 0.4) == 1L, dif := dif - 1L]
              ref[between(age, 70, 89) &
                    sex == "women" &
                    qimd %in% c("4d", "2") &
                    rbinom(.N, 1L, 0.4) == 1L, dif := dif + 1L]
              absorb_dt(dt, ref)

              # when not enough never smokers convert those ex smokers with the longer quit years
              tt <-
                dt[smok_status %in% 2:3, .(year, age, sex, qimd, pid, smok_quit_yrs, dif)]
              setnafill(tt, "c", 0, cols = "dif")
              tt[dif < 0, dif := 0L]
              setkey(tt, year, age, sex, qimd, smok_quit_yrs)
              pid_to_conv <-
                tt[dif > 0, .(pid = tail(pid, max(dif))), keyby = .(year, age, sex, qimd)]
              dt[pid_to_conv, on = .(year, pid), `:=`(
                smok_status = 1L,
                smok_quit_yrs = 0L,
                smok_dur = 0L,
                smok_cig = 0L
              )]

              # when too many never smokers convert to smok status 2 (occasional)
              tt <-
                dt[smok_status %in% 1, .(year, age, sex, qimd, pid, dif)]
              setnafill(tt, "c", 0, cols = "dif")
              tt[dif > 0, dif := 0L]
              tt[, dif := -dif]
              setkey(tt, year, age, sex, qimd)
              # Ensure there are enough people to sample from
              ttt <-
                tt[, .(lpid = length(pid), mdif =  max(dif)), by = .(year, age, sex, qimd)][mdif >
                                                                                              lpid, ]
              tt[ttt, on = .NATURAL, dif := i.lpid]
              pid_to_conv <-
                tt[dif > 0, .(pid = resample(pid, max(dif))), keyby = .(year, age, sex, qimd)]
              dt[pid_to_conv, on = .(year, pid), `:=`(
                smok_status = 2L,
                smok_quit_yrs = dt[smok_status == 2L, median(smok_quit_yrs)],
                smok_dur = dt[smok_status == 2L, median(smok_dur)],
                smok_cig = 1L
              )]
              dt[, dif := NULL]

              # Same logic for active smokers
              # I will further calibrate to better match HSE
              # calculate dif between ref (multinom) and simsmok
              obs <-
                dt[smok_status == 4L, .(nsa = .N), keyby = .(year, age, sex, qimd)] # ? add sha/ethn
              ref <-
                dt[smok_status_ref == 4L, .(nsr = .N), keyby = .(year, age, sex, qimd)] # ? add sha/ethn
              absorb_dt(ref, obs)
              setnafill(ref, "c", 0, cols = "nsa")
              ref[, `:=`(dif = nsr - nsa,
                         nsr = NULL,
                         nsa = NULL)]
              # Further calibrate to match better with HSE
              ref[sex == "men" &
                    rbinom(.N, 1L, 0.3) == 1L, dif := dif - 1L]
              ref[sex == "women" &
                    qimd != "1 most deprived" &
                    qimd != "5 least deprived" &
                    rbinom(.N, 1L, 0.8) == 1L, dif := dif - 1L]
              ref[qimd == "1 most deprived" &
                    rbinom(.N, 1L, 0.5) == 1L, dif := dif + 1L]
              ref[qimd == "5 least deprived" &
                    rbinom(.N, 1L, 0.6) == 1L, dif := dif - 1L] # - reduces
              absorb_dt(dt, ref)

              # when not enough active smokers convert those ex smokers with the shortest quit years
              tt <-
                dt[smok_status == 3L, .(year, age, sex, qimd, pid, smok_quit_yrs, dif)]
              setnafill(tt, "c", 0, cols = "dif")
              tt[dif < 0, dif := 0L]
              setkey(tt, year, age, sex, qimd, smok_quit_yrs)
              pid_to_conv <-
                tt[dif > 0, .(pid = head(pid, max(dif))), keyby = .(year, age, sex, qimd)]
              dt[pid_to_conv, on = .(year, pid),
                 `:=`(
                   smok_status = 4L,
                   smok_quit_yrs = 0L,
                   smok_dur = smok_dur + 1
                 )] # TODO fix smoking duration

              # when too many never smokers convert to smok status 3
              tt <-
                dt[smok_status == 4L, .(year, age, sex, qimd, pid, dif)]
              setnafill(tt, "c", 0, cols = "dif")
              tt[dif > 0, dif := 0L]
              tt[, dif := -dif]
              setkey(tt, year, age, sex, qimd)
              # Ensure there are enough people to sample from
              ttt <-
                tt[, .(lpid = length(pid), mdif =  max(dif)), by = .(year, age, sex, qimd)][mdif >
                                                                                              lpid,]
              tt[ttt, on = .NATURAL, dif := i.lpid]
              pid_to_conv <-
                tt[dif > 0, .(pid = resample(pid, max(dif))), keyby = .(year, age, sex, qimd)]
              dt[pid_to_conv, on = .(year, pid), `:=`(smok_status = 3L, smok_quit_yrs = 1L)]
              dt[, dif := NULL]

              rm(tt, ttt, obs, ref, pid_to_conv)
            }

            # Assign smok_cig_curr when pid_mrk == true (the first year an individual enters the simulation)
            set(dt, NULL, "smok_cig", 0L)

            tbl <-
              read_fst("./inputs/exposure_distributions/smok_cig_curr_table.fst",
                       as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(dt), names(tbl)))
            absorb_dt(dt, tbl)
            dt[smok_status == 4L,
               smok_cig := qZINBI(rankstat_smok_cig_curr, mu, sigma, nu)]
            dt[, (col_nam) := NULL]

            # Assign smok_cig_ex when pid_mrk == true (the first year an individual enters the simulation)
            # dt[smok_status == 2, smok_cig := 1L]


            tbl <-
              read_fst("./inputs/exposure_distributions/smok_cig_ex_table.fst",
                       as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(dt), names(tbl)))
            absorb_dt(dt, tbl)
            dt[(pid_mrk) &
                 smok_status == 3L,
               smok_cig := my_qZABNB(rankstat_smok_cig_ex,
                                     mu,
                                     sigma,
                                     nu,
                                     tau,
                                     n_cpu = design_$sim_prm$n_cpu)]
            dt[, (col_nam) := NULL]

            simsmok_cig(dt) # carry forward smok_cig if smok_status == 3
            dt[smok_cig == 0L & smok_status > 1L, smok_cig := 1L]

            if (design_$sim_prm$simsmok_calibration)
              simsmok_postcalibration(dt) # need to be post cig simulation

            dt[, smok_status := factor(smok_status)]

            dt[, c(
              "rankstat_smok",
              "rankstat_smok_cig_curr",
              "rankstat_smok_cig_ex",
              "prb_smok_incid",
              "prb_smok_cess",
              "smok_status_ref"
            ) := NULL]

            # Generate ETS (BI) ----
            if (design_$sim_prm$logs) message("Generate ETS")

            # Note at the moment this is independent of smoking prevalence TODO
            # calculate how many each smoker pollutes by year, SHA (not qimd) to
            # be used in scenarios. Ideally correct for mortality
            tbl <-
              read_fst("./inputs/exposure_distributions/ets_table.fst", as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(dt), names(tbl)))
            absorb_dt(dt, tbl)
            dt[, ets := as.integer(rank_ets < mu)]
            dt[, rank_ets := NULL]
            dt[, (col_nam) := NULL]
            # View(dt[, prop_if(ets == 1)/prop_if(smok_status == "4"), keyby = .(year, sha)])

            # Generate alcohol (ZINBI) ----
            if (design_$sim_prm$logs) message("Generate alcohol")

            tbl <-
              read_fst("./inputs/exposure_distributions/alcohol_table.fst", as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(dt), names(tbl)))
            absorb_dt(dt, tbl)
            dt[, alcohol := as.integer(qZINBI(rank_alcohol, mu, sigma, nu))]
            dt[, rank_alcohol := NULL]
            dt[, (col_nam) := NULL]

            # Generate BMI (BCPEo) ----
            if (design_$sim_prm$logs) message("Generate BMI")

            tbl <-
              read_fst("./inputs/exposure_distributions/bmi_table.fst", as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(dt), names(tbl)))
            absorb_dt(dt, tbl)
            dt[, bmi := my_qBCPEo(rank_bmi, mu, sigma, nu, tau, n_cpu = design_$sim_prm$n_cpu)]
            dt[, rank_bmi := NULL]
            dt[, (col_nam) := NULL]

            # Generate SBP (BCPEo) ----
            if (design_$sim_prm$logs) message("Generate SBP")

            tbl <-
              read_fst("./inputs/exposure_distributions/sbp_table.fst", as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(dt), names(tbl)))
            absorb_dt(dt, tbl)
            dt[, sbp := my_qBCPEo(rank_sbp, mu, sigma, nu, tau, n_cpu = design_$sim_prm$n_cpu)]
            dt[, rank_sbp := NULL]
            dt[, (col_nam) := NULL]

            # Generate BP medication (BI) -----
            if (design_$sim_prm$logs) message("Generate BP medication")

            dt[, `:=` (sbp_acc = sbp,
                       sbp = round(clamp(sbp, 110, 200), -1))]
            tbl <-
              read_fst("./inputs/exposure_distributions/bp_med_table.fst", as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(dt), names(tbl)))
            absorb_dt(dt, tbl)
            dt[, bpmed := as.integer(rank_bpmed < mu)]
            dt[, rank_bpmed := NULL]
            dt[, (col_nam) := NULL]
            dt[, `:=` (sbp = sbp_acc,
                       sbp_acc = NULL)]

            # TODO calculate probability of dgn HTN

            # Generate tchol (BCT) ----
            if (design_$sim_prm$logs) message("Generate tchol")

            tbl <-
              read_fst("./inputs/exposure_distributions/tchol_table.fst", as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(dt), names(tbl)))
            absorb_dt(dt, tbl)
            dt[, tchol := my_qBCT(rank_tchol, mu, sigma, nu, tau, n_cpu = design_$sim_prm$n_cpu)]
            dt[, rank_tchol := NULL]
            dt[, (col_nam) := NULL]

            # Generate HDL (to tchol ratio) (GB1) ----
            if (design_$sim_prm$logs) message("Generate HDL (to tchol ratio)")

            # NOTE this very highly correlated with hdl level (~0.76) and
            #  highly to tchol (~-0.47). The latter is captured by the correlated RNs
            tbl <-
              read_fst("./inputs/exposure_distributions/hdl_to_tchol_table.fst",
                       as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(dt), names(tbl)))
            absorb_dt(dt, tbl)
            dt[, tchol_hdl_ratio := 1 / qGB1(rank_hdl, mu, sigma, nu, tau)]
            dt[, rank_hdl := NULL]
            dt[, (col_nam) := NULL]

            # Generate statins medication (BI) -----
            if (design_$sim_prm$logs) message("Generate statins medication")

            dt[, `:=` (tchol_acc = tchol,
                       tchol = round(clamp(tchol, 2, 12), 0))]
            tbl <-
              read_fst("./inputs/exposure_distributions/statin_px_table.fst",
                       as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(dt), names(tbl)))
            absorb_dt(dt, tbl)
            dt[, statin_px := as.integer(rank_statin_px < mu)]
            dt[, rank_statin_px := NULL]
            dt[, (col_nam) := NULL]
            dt[, `:=` (tchol = tchol_acc,
                       tchol_acc = NULL)]


            # Estimate number of comorbidities (ncc) calculation ----
            # to be used in QALY
            if (design_$sim_prm$logs) message("Generate ncc")

            dt[, ncc := as.integer(clamp(qbinom(rankstat_ncc, ceiling(age / 8L),
                                                fifelse(age < 55, 0.25, 0.40)),
                                         0, 10))]
            # calibrated to Sullivan et all 2011 (web table 1)
            # to_agegrp(output, 10L, 89L, "age", "agegrp10", to_factor = TRUE)
            # output[, round(mean(ncc), 1), keyby = agegrp10]
            # target by agegrp 1.1  1.6  2.4  3.1  4.0  4.4 from

            dt[, `:=` (
              pid_mrk = NULL,
              # to be recreated when loading synthpop
              rankstat_ncc = NULL
            )]

            # Estimate history of cancer
            tbl <-
              read_fst("./inputs/exposure_distributions/history_of_cancer.fst",
                       as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(dt), names(tbl)))
            absorb_dt(dt, tbl)
            dt[, history_of_ca := as.integer(rankstat_ca_history < mu)]
            dt[, c(col_nam, "rankstat_ca_history") := NULL]

            # Estimate family history of lung ca (crude approx)
            # The prob of having lung ca is around 50/1e5.
            # So the probability of having at least
            # one of 3 family members with lung ca is 1 - (1-50/1e5)^3.
            tt <- 1 - (1 - 50 / 1e5) ^ 3
            dt[, fam_lung_ca := as.integer(rankstat_famlungca < tt)]
            dt[, ("rankstat_famlungca") := NULL]


            dt[, statin_adherence := qBE(rankstat_statin_adherence, design_$sim_prm$statin_adherence, 0.2)]
            dt[, bpmed_adherence := qBE(rankstat_bpmed_adherence, design_$sim_prm$bpmed_adherence, 0.2)]
            dt[, c("rankstat_bpmed_adherence", "rankstat_statin_adherence") := NULL]

            exps_tolag <- c(
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
            exps_nam <-  paste0(exps_tolag, "_curr_xps")
            setnames(dt, exps_tolag, exps_nam)


            # Prune & write synthpop to disk ----
            # del rn as they are reproducible
            nam <- c("LSOA11CD",
                     # "LAD11CD",
                     # "LAD11NM",
                     # "tds_quintile",
                     # "imd",
                     # "sha", # NEEDED for social scenarios
                     # "pid_mrk",
                     "CCG17CDH"
                     )
            dt[, (nam) := NULL]

            if ("age100" %in% names(dt)) {
              dt[, age := NULL]
              setnames(dt, "age100", "age")
            }

            setkey(dt, pid, year) # Just in case
            setcolorder(dt, c("pid", "year", "age", "sex", "dimd"))
            setindexv(dt, c("year", "age", "sex", "dimd", "ethnicity"))
            if (design_$sim_prm$logs) message("Writing synthpop to disk")
            write_fst(dt,
                      filename_$synthpop,
                      90) # 100 is too slow
          return(invisible(NULL))
        },


      # Load a synthpop file from disk in full or in chunks.
      get_synthpop =
        function(exclude_cols = c()) {
          mm_synthpop <- metadata_fst(private$filename$synthpop)
          mm_synthpop <- setdiff(mm_synthpop$columnNames, exclude_cols)

          # Read synthpop

          dt <- read_fst(private$filename$synthpop,
                     columns = mm_synthpop,
                     as.data.table = TRUE)
          dt <- dt[between(
            year,
            private$design$sim_prm$init_year - private$design$sim_prm$maxlag,
            private$design$sim_prm$init_year + private$design$sim_prm$sim_horizon_fromGUI
          ) &
            between(age,
                    private$design$sim_prm$ageL - private$design$sim_prm$maxlag,
                    private$design$sim_prm$ageH)]

          dt[, pid_mrk := mk_new_simulant_markers(pid)]
          # Above necessary because of pruning  and potential merging above

          # Ensure pid does not overlap for files from different mc
          new_n <- uniqueN(dt$pid)
          it <- as.integer(ceiling(self$mc %% private$design$sim_prm$n_synthpop_aggregation))
          it[it == 0L] <- private$design$sim_prm$n_synthpop_aggregation
          it <- it - 1L
          if (max(dt$pid + (private$design$sim_prm$n_synthpop_aggregation - 1) * new_n) < .Machine$integer.max) {
            dt[, pid := as.integer(pid + it * new_n)]
          } else stop("pid larger than int32 limit.")
          # generate population weights
          private$gen_pop_weights(dt, private$design) # TODO replace?
          dt[, smok_packyrs_curr_xps := as.integer(round(smok_cig_curr_xps * smok_dur_curr_xps / 20))]
          set(dt, NULL, "all_cause_mrtl", 0L)

          # dt[, dead := identify_longdead(all_cause_mrtl, pid_mrk)]
          # dt[, ncc := clamp(
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


          invisible(dt)
        },

      # Calculate weights so that their sum is the population of the area based
      # on ONS. It takes into account synthpop aggregation. So you need to sum
      # all the synthpops belong to the same aggregation to reach the total pop.
      # NOTE Wt are still incomplete because they assume everyone remains alive.
      # So baseline population underestimated as clearly some die
      gen_pop_weights = function(dt, design) {
        tt <-
          read_fst("./inputs/pop_projections/lad17_proj.fst", as.data.table = TRUE)
        lads <- private$get_unique_LADs(design)
        tt <- tt[LAD17CD %in% lads &
                   between(age, min(dt$age), max(dt$age)) &
                   between(year, min(dt$year), max(dt$year)),
                 .(pops = sum(pops)), keyby = .(year, age, sex)]
        dt[, wt := .N, by = .(year, age, sex)]
        absorb_dt(dt, tt)
        dt[, wt := pops / (wt * design$sim_prm$n_synthpop_aggregation)]
        dt[, pops := NULL]

        invisible(dt)
      }
    )
  )
