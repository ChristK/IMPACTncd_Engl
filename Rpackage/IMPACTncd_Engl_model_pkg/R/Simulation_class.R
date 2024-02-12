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
# field 'output' that is a data.table

#' @export
`[.Simulation` <- function(x, ...) x$output[...]

#' R6 Class representing a simulation environment
#' @description A simulation environment.
#' @details To be completed...
#' @export
Simulation <-
  R6::R6Class(
    classname = "Simulation",

    # public ------------------------------------------------------------------
    public = list(
      #' @field design A Design object.
      design = NA,

      #' @field diseases A list of Disease objects.
      diseases = NA,

      #' @field RR A list of RR for the simulated exposures.
      RR = NA,

      #' @field scenarios A list of scenario objects.
      scenarios = NA,

      # initialise ----
      #' @description Create a new simulation object.
      #' @param sim_prm Either a path to a yaml file or a Design object.
      #' @return A new `Simulation` object.
      initialize = function(sim_prm) {
        if (is.character(sim_prm))
          self$design <- Design$new(sim_prm)
        else if (inherits(sim_prm, "Design"))
          self$design <- sim_prm$clone(deep = TRUE)
        else
          stop("sim_prm need to be a path to an appropriate yaml file or a Design object")

        data.table::setDTthreads(threads = self$design$sim_prm$clusternumber,
                                 restore_after_fork = NULL)
        fst::threads_fst(
          nr_of_threads = self$design$sim_prm$clusternumber,
          reset_after_fork = NULL
        )

        deployArgs <- list(
          sTag = self$design$sim_prm$sTag,
          bOverwriteFilesOnDeploy = self$design$sim_prm$bOverwriteFilesOnDeploy,
          sDeployToRootDirPath = self$design$sim_prm$sDeployToRootDirPath,
          sUploadSrcDirPath = self$design$sim_prm$sUploadSrcDirPath,
          sToken = self$design$sim_prm$sToken
        )
        self$DeployGitHubAssets(deployArgs)
        message("Loading exposures.")
        # RR Create a named list of Exposure objects for the files in
        # ./inputs/RR
        fl <- list.files(path = "./inputs/RR", pattern = ".csvy$", full.names = TRUE)
        # RR <- future_lapply(fl, Exposure$new, future.seed = 950480304L)
        self$RR <- lapply(fl, Exposure$new, design = self$design)
        names(self$RR) <- sapply(self$RR, function(x) x$get_name())
        # invisible(future_lapply(RR, function(x) {
        #   x$gen_stochastic_effect(design, overwrite = FALSE, smooth = FALSE)
        # }, future.seed = 627524136L))
        invisible(lapply(self$RR, function(x) {
          x$gen_stochastic_effect(self$design, overwrite = FALSE, smooth = FALSE)
        }))
        # NOTE smooth cannot be exported to Design for now, because the first
        # time this parameter changes we need logic to overwrite unsmoothed
        # files
        rm(fl)

        # Generate diseases
        message("Loading diseases.")
        self$diseases <- lapply(self$design$sim_prm$diseases, function(x) {
          x[["design_"]] <- self$design
          x[["RR"]] <- self$RR
          do.call(Disease$new, x)
        })
        names(self$diseases) <- sapply(self$design$sim_prm$diseases, `[[`, "name")

        message("Generating microsimulation structure.")
        # Generate the graph with the causality structure
        ds <- unlist(strsplit(names(self$RR), "~"))
        ds[grep("^smok_", ds)] <- "smoking"
        ds <- gsub("_prvl$", "", ds)

        ds1 <- ds[as.logical(seq_along(ds) %% 2)]
        ds2 <- ds[!as.logical(seq_along(ds) %% 2)]
        ds <- unique(data.table(ds1, ds2))

        private$causality_structure <- make_graph(unlist(transpose(ds)),
                                                  directed = TRUE)

        # European standardised population 2013 (esp) weights
        tt <- data.table(agegrp = agegrp_name(0, 99),
                         wt_esp  = c(1000, 4000, 5500, 5500, 5500, 6000, 6000, 6500,
                                     7000, 7000, 7000, 7000, 6500, 6000, 5500, 5000,
                                     4000, 2500, 1500, 800, 200))
        esp <- CJ(agegrp = agegrp_name(0, 99),
                  sex = c("men", "women"),
                  dimd = c("1 most deprived", as.character(2:9), "10 least deprived")
        )

        private$esp_weights <- copy(absorb_dt(esp, tt))

        private$death_codes <- unlist(lapply(self$diseases, function(x)
          x$meta$mortality$code))
        private$death_codes[["alive"]] <- 0L

        invisible(self)
      },

      # UploadGitHubAssets ----
      # Upload GitHub assets
      #' @description Scan a specified directory for large asset files and upload these to GitHub.
      #' @param sToken GitHub personal access token (default: NULL).
      #' @param sRepo GitHub repository (default: NULL).
      #' @param sTag GitHub tag or release (default: NULL).
      #' @param sDeployToRootDirPath Path to the root directory for deployment (default: NULL).
      #' @param bOverwriteFilesOnDeploy Logical indicating whether to overwrite existing files on deployment (default: FALSE).
      #' @param sUploadSrcDirPath Path to the source directory for upload (default: NULL).
      #' @param iTestWithFirstNAssets Number of assets to upload as a test
      #'
      #' @return Invisible self for chaining
      #' @export
      UploadGitHubAssets = function(sToken = NULL, iTestWithFirstNAssets = 0,
                                     sTag = NULL, sUploadSrcDirPath = NULL,
                                     sDeployToRootDirPath = NULL,
                                     bOverwriteFilesOnDeploy) { # consider as a Simulation class method
        # set the GitHub repo
        sRepo <- "ChristK/IMPACTncd_Engl"
        # get GitHub asset route data
        GetGitHubAssetRouteInfo(sRepo, sTag, sUploadSrcDirPath, sDeployToRootDirPath,
                                bOverwriteFilesOnDeploy, sToken, iTestWithFirstNAssets)
        sUploadSrcDirPath <- TrimSlashes(sUploadSrcDirPath, bRidStartSlash = FALSE)
        # pb_new_release(sRepo, sTag) # Only need to run the first time a github repo is created
        # pb_release_delete(sRepo, sTag)
        # find local assets, excluding secure_data, "*tmp.qs", and "*parf/PARF_*.qs" files.
        lsLocalAssetPathNames <- list.files(sUploadSrcDirPath, pattern = ".fst$|.xls$|.xlsx$|.qs$",
                                            all.files = TRUE, full.names = TRUE, recursive = TRUE)
        lsLocalAssetPathNames <- grep("secure_data", lsLocalAssetPathNames, value = TRUE, invert = TRUE)
        # NOTE currently no necessary files ar .qs. For future proof I add them above
        # and exclude the below
        lsLocalAssetPathNames <- grep("tmp.qs$", lsLocalAssetPathNames, value = TRUE, invert = TRUE)
        lsLocalAssetPathNames <- grep("/simulation/parf/PARF_.*\\.qs$", lsLocalAssetPathNames, value = TRUE, invert = TRUE)
        lsLocalAssetPathNames <- grep("/DELETEME/", lsLocalAssetPathNames, value = TRUE, invert = TRUE)

        # write table with original and *unique sanitised* file names and directories
        lsSanitisedFileNames <- gsub("[^[:alnum:]&&^\\.]", ".", basename(lsLocalAssetPathNames)) # replace all not (alphanumeric or '.') with .
        lsSanitisedFileNames <- sub("^[.]", "default.", lsSanitisedFileNames)
        dtOriginalAndSanitisedFilePathNames <- data.table(
          orig_file = basename(lsLocalAssetPathNames),
          sanit_file = lsSanitisedFileNames, abs_dir = dirname(lsLocalAssetPathNames),
          rel_dir = gsub(sUploadSrcDirPath, "", dirname(lsLocalAssetPathNames)), key = "orig_file"
        ) # rel_dir to keep initial slash
        # print(nrow(dtOriginalAndSanitisedFilePathNames))
        # return(0)
        if (any(duplicated(dtOriginalAndSanitisedFilePathNames$sanit_file))) stop("Duplicated filenames found")
        if (iTestWithFirstNAssets != 0) {
          dtOriginalAndSanitisedFilePathNames <- dtOriginalAndSanitisedFilePathNames[1:iTestWithFirstNAssets, ]
        }
        fwrite(dtOriginalAndSanitisedFilePathNames, "./auxil/filindx.csv")

        lsSrcFilePath <- file.path(dtOriginalAndSanitisedFilePathNames$abs_dir, dtOriginalAndSanitisedFilePathNames$orig_file)
        lsHttpResponses <- piggyback::pb_upload(lsSrcFilePath,
                                                repo = sRepo, tag = sTag, name = dtOriginalAndSanitisedFilePathNames$sanit_file,
                                                .token = sToken, use_timestamps = FALSE, overwrite = TRUE
        ) # piggyback:: prefix necessary to use R.utils::reassignInPackage() injected code modification

        StopOnHttpFailure(lsHttpResponses, TRUE)
        invisible(self)
      },

      # DeployGitHubAssets ----
      # Download GitHub Assets
      #' @description This function downloads assets from a GitHub repository and organizes them
      #' in a specified directory structure.
      #'
      #' @param sToken GitHub personal access token (default: NULL).
      #' @param sRepo GitHub repository (default: NULL).
      #' @param sTag GitHub tag or release (default: NULL).
      #' @param sDeployToRootDirPath Path to the root directory for deployment (default: NULL).
      #' @param bOverwriteFilesOnDeploy Logical indicating whether to overwrite existing files on deployment (default: FALSE).
      #' @param sUploadSrcDirPath Path to the source directory for upload (default: NULL).
      #' @param args a list of the above specified parameters
      #'
      #' @return Invisible self for chaining
      #' @export
      DeployGitHubAssets = function(args = NULL) {
        sTag <- ifelse(is.null(args$sTag), NULL, args$sTag)
        bOverwriteFilesOnDeploy <- ifelse(is.null(args$bOverwriteFilesOnDeploy), NULL, args$bOverwriteFilesOnDeploy)
        sDeployToRootDirPath <- ifelse(is.null(args$sDeployToRootDirPath), NULL, args$sDeployToRootDirPath)
        sUploadSrcDirPath <- ifelse(is.null(args$sUploadSrcDirPath), NULL, args$sUploadSrcDirPath)
        sToken <- ifelse(is.null(args$sToken), gh::gh_token(), args$sToken)
        sRepo <- "ChristK/IMPACTncd_Engl"

        GetGitHubAssetRouteInfo(sRepo = sRepo, sTag = sTag, iTestWithFirstNAssets,
                                sDeployToRootDirPath = sDeployToRootDirPath,
                                sUploadSrcDirPath = sUploadSrcDirPath,
                                bOverwriteFilesOnDeploy = bOverwriteFilesOnDeploy,
                                sToken = sToken)

        sDeployToRootDirPath <- TrimSlashes(sDeployToRootDirPath, bRidStartSlash = FALSE)

        sanitisedToOriginalFilePaths <- fread(file.path(sDeployToRootDirPath, "/auxil/filindx.csv"), key = "orig_file")

        if (iTestWithFirstNAssets != 0) {
          sanitisedToOriginalFilePaths <- sanitisedToOriginalFilePaths[1:iTestWithFirstNAssets, ]
        }

        subDirectoryPaths <- sapply(sanitisedToOriginalFilePaths$rel_dir, TrimSlashes)

        sapply(file.path(sDeployToRootDirPath, unique(subDirectoryPaths)), dir.create, showWarnings = FALSE, recursive = TRUE)

        all_files <- file.path(sDeployToRootDirPath, subDirectoryPaths, sanitisedToOriginalFilePaths$orig_file)

        if (!all(sapply(all_files, file.exists))) {
          missing_files <- sanitisedToOriginalFilePaths$sanit_file[!file.exists(all_files)]
          missing_orig_files <- sanitisedToOriginalFilePaths$orig_file[which(sanitisedToOriginalFilePaths$sanit_file %in% missing_files)]
          missing_sub_dirs <- subDirectoryPaths[which(sanitisedToOriginalFilePaths$sanit_file %in% missing_files)]

          lsHttpResponses <- piggyback::pb_download(
            file = missing_files,
            dest = file.path(sDeployToRootDirPath, missing_sub_dirs, missing_orig_files),
            repo = sRepo, tag = sTag, overwrite = bOverwriteFilesOnDeploy,
            use_timestamps = FALSE, .token = sToken)

          StopOnHttpFailure(lsHttpResponses, FALSE)
        } else {
          if (bOverwriteFilesOnDeploy) {
            lsHttpResponses <- piggyback::pb_download(
              file = sanitisedToOriginalFilePaths$sanit_file,
              dest = file.path(sDeployToRootDirPath, subDirectoryPaths, sanitisedToOriginalFilePaths$orig_file),
              repo = sRepo, tag = sTag, overwrite = bOverwriteFilesOnDeploy,
              use_timestamps = FALSE, .token = sToken)

            StopOnHttpFailure(lsHttpResponses, FALSE)
          } else {
            print("All assets are already downloaded.")
          }
        }
      },

      # run ----
      #' @description Runs a simulation
      #' @param mc A positive sequential integer vector with the Monte Carlo
      #'   iterations of synthetic population to simulate, or a scalar.
      #' @param multicore If TRUE run the simulation in parallel.
      #' @param scenario_nam A string for the scenario name (i.e. sc1)
      #' @return The invisible self for chaining.
      run = function(mc, multicore = TRUE, scenario_nam) {

        if (!is.integer(mc)) stop("mc need to be an integer")
        if (any(mc <= 0)) stop("mc need to be positive integer")

        # check if sequential vector. Necessary if
        # design$sim_prm$n_synthpop_aggregation > 1
        if (anyNA(mc) || any(is.infinite(mc)) || length(mc) < 1L ||
            (length(mc) > 1L && diff(mc[1:2]) == 0) ||
            (length(mc) > 1L && diff(range(diff(mc))) > sqrt(.Machine$double.eps)))
          stop("mc need to be a sequential integer vector, or a scalar")
        # NOTE mc is in fact mc_aggr. mc_ is the mc of the synthpop
        mc_sp <-
          (
            min(mc) * self$design$sim_prm$n_synthpop_aggregation -
              self$design$sim_prm$n_synthpop_aggregation + 1L
          ):(max(mc) * self$design$sim_prm$n_synthpop_aggregation)

        # Create folders if don't exist
        # TODO write hlp function and use lapply
        if (file.exists(self$design$sim_prm$output_dir) && file.access(self$design$sim_prm$output_dir, mode = 2) == -1L)
          stop("You don't have write access to the output folder. Please change the permissions or the path.")

        message("Creating output subfolders.")
        private$create_new_folder(self$design$sim_prm$output_dir, self$design$sim_prm$logs)
        private$create_new_folder(private$output_dir("summaries/"), self$design$sim_prm$logs)
        private$create_new_folder(private$output_dir("tables/"), self$design$sim_prm$logs)
        private$create_new_folder(private$output_dir("plots/"), self$design$sim_prm$logs)
        private$create_new_folder(private$output_dir("lifecourse/"), self$design$sim_prm$logs)
        if (self$design$sim_prm$export_PARF) {
          private$create_new_folder(private$output_dir("parf/"), self$design$sim_prm$logs)
        }
        if (self$design$sim_prm$export_xps) {
          private$create_new_folder(private$output_dir("xps/"), self$design$sim_prm$logs)
        }
        if (self$design$sim_prm$logs) {
          private$create_new_folder(private$output_dir("logs/"), self$design$sim_prm$logs)
        }

        if (file.exists(self$design$sim_prm$synthpop_dir) && file.access(self$design$sim_prm$synthpop_dir, mode = 2) == -1L)
          stop("You don't have write access to the synthpop folder. Please change the permissions or the path.")

        # NOTE code below is duplicated in Synthpop class. This is intentional
        private$create_new_folder(self$design$sim_prm$synthpop_dir, self$design$sim_prm$logs)


        if (any(file.exists( # TODO fix when lifecourse is not saved
          file.path(
            self$design$sim_prm$output_dir,
            "lifecourse",
            paste0(mc, "_lifecourse.cs")
          )
        ))) {
          # stop("Results from a previous simulation exists in the output
          #      folder. Please remove them before run a new one.")
          message(
            "Results from a previous simulation exists in the output folder. Please remove them if this was unintentional."
          )
        }



        # Generate PARF files if they don't exist. Note that generation is
        # multicore
        lapply(self$diseases, function(x) {
          x$gen_parf_files(self$design, self$diseases)
        })

        if (multicore) {

          if (Sys.info()["sysname"] == "Windows") {
            cl <-
              makeCluster(self$design$sim_prm$clusternumber) # used for clustering. Windows compatible
            registerDoParallel(cl)
          } else {
            registerDoParallel(self$design$sim_prm$clusternumber) # used for forking. Only Linux/OSX compatible
          }

          if (self$design$sim_prm$logs)
            private$time_mark("Start of parallelisation")

          xps_dt <- foreach(
            mc_iter = mc_sp,
            .inorder = FALSE,
            .options.multicore = list(preschedule = FALSE),
            .verbose = self$design$sim_prm$logs,
            .packages = c(
              "R6",
              "gamlss.dist",
              "dqrng",
              "CKutils",
              "IMPACTncdEngl",
              "fst",
              "data.table"
            ),
            .export = NULL,
            .noexport = NULL # c("time_mark")
          ) %dopar% {

            private$run_sim(mc_ = mc_iter, scenario_nam)

          }

          if (exists("cl")) stopCluster(cl)

          if (self$design$sim_prm$logs) private$time_mark("End of parallelisation")


        } else {
          if (self$design$sim_prm$logs)
            private$time_mark("Start of single-core run")

          lapply(mc_sp, private$run_sim, scenario_nam)

          if (self$design$sim_prm$logs)
            private$time_mark("End of single-core run")

        }

        if (self$design$sim_prm$avoid_appending_csv) {
          # collect the lifecourse fragmented file
          private$collect_files("lifecourse", "_lifecourse.csv$", to_mc_aggr = TRUE)

          if (self$design$sim_prm$export_xps) {
            private$collect_files("xps", "_xps20.csv$", to_mc_aggr = FALSE)
            private$collect_files("xps", "_xps_esp.csv$", to_mc_aggr = FALSE)
          }

          if (self$design$sim_prm$logs)
            private$time_mark("End of collecting mc lifecourse files")
        }

        while (sink.number() > 0L) sink()

        invisible(self)
      },

      # export_summaries ----

      #' @description Process the lifecourse files
      #' @param multicore If TRUE run the simulation in parallel.
      #' @param type The type of summary to extract.
      #' @return The invisible self for chaining.
      export_summaries = function(multicore = TRUE,
                                  type = c("le", "hle", "dis_char", "prvl",
                                           "incd", "dis_mrtl", "mrtl",
                                           "allcause_mrtl_by_dis", "cms")) {

        fl <- list.files(private$output_dir("lifecourse"), full.names = TRUE)

        # logic to avoid inappropriate dual processing of already processed mcs
        # TODO take into account scenarios
        if ("le" %in% type) file_pth <- private$output_dir("summaries/le_scaled_up.csv.gz") else
          if ("hle" %in% type) file_pth <- private$output_dir("summaries/hle_1st_cond_scaled_up.csv.gz") else
            if ("cms" %in% type) file_pth <- private$output_dir("summaries/cms_count_scaled_up.csv.gz") else
              if ("mrtl" %in% type) file_pth <- private$output_dir("summaries/mrtl_scaled_up.csv.gz") else
                if ("dis_mrtl" %in% type) file_pth <- private$output_dir("summaries/dis_mrtl_scaled_up.csv.gz") else
                  if ("dis_char" %in% type) file_pth <- private$output_dir("summaries/dis_characteristics_scaled_up.csv.gz") else
                    if ("incd" %in% type) file_pth <- private$output_dir("summaries/incd_scaled_up.csv.gz") else
                      if ("prvl" %in% type) file_pth <- private$output_dir("summaries/prvl_scaled_up.csv.gz") else
                        if ("allcause_mrtl_by_dis" %in% type) file_pth <- private$output_dir("summaries/all_cause_mrtl_by_dis_scaled_up.csv.gz")


                        if (file.exists(file_pth)) {
                          tt <- unique(fread(file_pth, select = "mc")$mc)
                          for (i in seq_along(tt)) {
                            fl <- grep(paste0("/", tt[[i]], "_lifecourse.csv.gz$"), fl,
                                       value = TRUE, invert = TRUE)
                          }
                        }
                        # end of logic

                        if (multicore) {

                          if (Sys.info()["sysname"] == "Windows") {
                            cl <-
                              makeCluster(self$design$sim_prm$clusternumber_export) # used for clustering. Windows compatible
                            registerDoParallel(cl)
                          } else {
                            registerDoParallel(self$design$sim_prm$clusternumber_export) # used for forking. Only Linux/OSX compatible
                          }

                          if (self$design$sim_prm$logs)
                            private$time_mark("Start exporting summaries")

                          void <- foreach(
                            i = seq_along(fl),
                            .inorder = TRUE,
                            .options.multicore = list(preschedule = FALSE),
                            .verbose = self$design$sim_prm$logs,
                            .packages = c(
                              "R6",
                              "CKutils",
                              "IMPACTncdEngl",
                              "data.table"
                            ),
                            .export = NULL,
                            .noexport = NULL # c("time_mark")
                          ) %dopar% {

                            lc <-   fread(fl[i], stringsAsFactors = TRUE, key = c("scenario", "pid", "year"))
                            private$export_summaries_hlpr(lc, type = type)
                            NULL
                          }

                          if (exists("cl")) stopCluster(cl)

                          if (self$design$sim_prm$logs)
                            private$time_mark("End of exporting summuries")


                        } else {
                          if (self$design$sim_prm$logs)
                            private$time_mark("Start of single-core run")

                          lapply(seq_along(fl), function(i) {
                            lc <-   fread(fl[i], stringsAsFactors = TRUE, key = c("scenario", "pid", "year"))
                            private$export_summaries_hlpr(lc, type = type)
                            NULL
                          })

                          if (self$design$sim_prm$logs)
                            private$time_mark("End of single-core run")

                        }

                        if (self$design$sim_prm$avoid_appending_csv) {
                          # collect the summary fragmentrd file
                          if ("le" %in% type) {
                            private$collect_files("summaries", "_le_scaled_up.csv$", to_mc_aggr = FALSE)
                            private$collect_files("summaries", "_le_esp.csv$", to_mc_aggr = FALSE)
                            private$collect_files("summaries", "_le60_scaled_up.csv$", to_mc_aggr = FALSE)
                            private$collect_files("summaries", "_le60_esp.csv$", to_mc_aggr = FALSE)
                          }
                          if ("hle" %in% type) {
                            private$collect_files("summaries", "_hle_1st_cond_scaled_up.csv$", to_mc_aggr = FALSE)
                            private$collect_files("summaries", "_hle_1st_cond_esp.csv$", to_mc_aggr = FALSE)
                            private$collect_files("summaries", "_hle_cmsmm1.5_scaled_up.csv$", to_mc_aggr = FALSE)
                            private$collect_files("summaries", "_hle_cmsmm1.5_esp.csv$", to_mc_aggr = FALSE)
                          }
                          if ("cms" %in% type) {
                            private$collect_files("summaries", "_cms_score_scaled_up.csv$", to_mc_aggr = FALSE)
                            private$collect_files("summaries", "_cms_score_esp.csv$", to_mc_aggr = FALSE)
                            private$collect_files("summaries", "_cms_score_by_age_scaled_up.csv$", to_mc_aggr = FALSE)
                            private$collect_files("summaries", "_cms_score_by_age_esp.csv$", to_mc_aggr = FALSE)
                            private$collect_files("summaries", "_cms_count_scaled_up.csv$", to_mc_aggr = FALSE)
                            private$collect_files("summaries", "_cms_count_esp.csv$", to_mc_aggr = FALSE)
                          }
                          if ("mrtl" %in% type) {
                            private$collect_files("summaries", "_mrtl_scaled_up.csv$", to_mc_aggr = FALSE)
                            private$collect_files("summaries", "_mrtl_esp.csv$", to_mc_aggr = FALSE)
                          }
                          if ("dis_mrtl" %in% type) {
                            private$collect_files("summaries", "_dis_mrtl_scaled_up.csv$", to_mc_aggr = FALSE)
                            private$collect_files("summaries", "_dis_mrtl_esp.csv$", to_mc_aggr = FALSE)
                          }
                          if ("dis_char" %in% type) {
                            private$collect_files("summaries", "_dis_characteristics_scaled_up.csv$", to_mc_aggr = FALSE)
                            private$collect_files("summaries", "_dis_characteristics_esp.csv$", to_mc_aggr = FALSE)
                          }
                          if ("incd" %in% type) {
                            private$collect_files("summaries", "_incd_scaled_up.csv$", to_mc_aggr = FALSE)
                            private$collect_files("summaries", "_incd_esp.csv$", to_mc_aggr = FALSE)
                          }
                          if ("prvl" %in% type) {
                            private$collect_files("summaries", "_prvl_scaled_up.csv$", to_mc_aggr = FALSE)
                            private$collect_files("summaries", "_prvl_esp.csv$", to_mc_aggr = FALSE)
                          }
                          if ("allcause_mrtl_by_dis" %in% type) {
                            private$collect_files("summaries", "_all_cause_mrtl_by_dis_scaled_up.csv$", to_mc_aggr = FALSE)
                            private$collect_files("summaries", "_all_cause_mrtl_by_dis_esp.csv$", to_mc_aggr = FALSE)
                          }

                          if (self$design$sim_prm$logs)
                            private$time_mark("End of collecting mc_aggr summary files")
                        }

                        while (sink.number() > 0L) sink()

                        invisible(self)
      },

      # get_causal_structure ----

      #' @description Returns the causality matrix and optionally plots the
      #'   causality structure.
      #' @param processed If `TRUE` generates the causality matrix from the
      #'   graph.
      #' @param print_plot If `TRUE` prints the causal structure graph.
      #' @param focus If missing the whole causal structure is returned.
      #'  Otherwise, if a named node only the subgraph of the 1st order
      #'  neighbours that point to the given vertrice is returned.
      #' @return The processed causality matrix if `processed = TRUE` or the
      #'   graph otherwise.
      get_causal_structure = function(processed = TRUE, print_plot = FALSE, focus = FALSE) {
        if (missing(focus)) {
          graph <- private$causality_structure
        } else {
          if (length(focus) > 1L) stop("focus need to be scalar string.")
          if (!focus %in% self$get_node_names()) stop("focus need to be a node name. Use get_node_names() to get the list of eligible values.")
          graph <- make_ego_graph(private$causality_structure, order = 1,  nodes = focus, mode = "in")[[1]]
        }
        if (print_plot) {
          print(
            plot(
              graph,
              vertex.shape = "none",
              edge.arrow.size = .3,
              vertex.label.font = 2,
              vertex.label.color = "gray40",
              edge.arrow.width = .5,
              vertex.label.cex = .7,
              edge.color = "gray85",
              layout = layout_components
            )
          )
        }

        if (processed) {
          graph <- as.matrix(as_adjacency_matrix(graph))
          n <- sapply(self$diseases, `[[`, "name")
          graph <- graph[rowSums(graph) > 0, colnames(graph) %in% n]
        }

        return(graph)
      },

      # get_node_names ----

      #' @description Returns the names of all exposures and diseases.
      #' @return A string vector.
      get_node_names = function() {
        return(V(private$causality_structure)$name)
      },

      # get_causal_path ----

      #' @description Returns the causal paths between an exposure and an outcome (disease).
      #' @param from the beginning of the path (an exposure) as a string. Use `get_node_names` for available nodes.
      #' @param to the end of the path (a disease) as a string. Use `get_node_names` for available nodes.
      #' @param shortest_paths Boolean. If true, only returns the paths with the smallest number of nodes. Else, all possible paths (excluding multiple and loop edges) are returned.
      #' @return A list with all the possible paths between exposure and disease.
      get_causal_path = function(from, to, shortest_paths = FALSE) {
        nm <- V(private$causality_structure)$name
        from <- which(nm == from)
        to <- which(nm == to)
        if (shortest_paths) {
          out <- get.all.shortest.paths(private$causality_structure, from, to, mode = "out")
        } else {
          out <- all_simple_paths(private$causality_structure, from, to, mode = "out")
        }
        return(out)
      },

      # update_design ----

      #' @description Updates the Design object that is stored in the Simulation
      #'   object.
      #' @param new_design A design object with the simulation parameters.
      #' @return The invisible self for chaining.
      update_design = function(new_design) {
        if (!inherits(new_design, "Design"))
          stop("Argument new_design needs to be a Design object.")

        self$design <- new_design

        invisible(self)
      },

      # del_outputs ----

      #' @description Delete all output files.
      #' @return The invisible self for chaining.
      del_outputs = function() {

        if (dir.exists(self$design$sim_prm$output_dir)) {

          # Check for safety that folders /lifecourse, /tables, /plots, and /summaries exist to avoid accidental deletes of other folders
          if (dir.exists(file.path(self$design$sim_prm$output_dir, "lifecourse")) &&
              dir.exists(file.path(self$design$sim_prm$output_dir, "summaries")) &&
              dir.exists(file.path(self$design$sim_prm$output_dir, "tables")) &&
              dir.exists(file.path(self$design$sim_prm$output_dir, "plots"))) {

            fl <- list.files(self$design$sim_prm$output_dir,
                             full.names = TRUE,
                             recursive = TRUE
            )
            file.remove(fl)

            if (length(fl) > 0 && self$design$sim_prm$logs) {
              message("Output files deleted.")
            }
          } else {
            message("Output folder doesn't contain the expected subfolders. Please check the output folder path.")
          }
        } else { # If output folder doesn't exist
          message("Output folder doesn't exist.")
        }

        invisible(self)
      },

      # del_logs ----

      #' @description Delete log files.
      #' @return The invisible self for chaining.
      del_logs = function() {

        fl <- list.files(private$output_dir("logs/"), full.names = TRUE)

        file.remove(fl)

        if (length(fl) > 0 && self$design$sim_prm$logs)
          message("Log files deleted.")

        invisible(self)
      },

      # get_esp ----

      #' @description Get the European Standardised Population 2013 by sex and
      #'   dimd.
      #' @return A data.table with the European Standardised Population 2013.
      get_esp = function() {
        private$esp_weights
      },

      # get_mm_weights ----

      #' @description Get the disease multimorbidity weights (i.e. Cambridge
      #'   Morbidity Score weights).
      #' @return A named vector with disease weights.
      get_mm_weights = function() {
        unlist(sapply(self$diseases, function(x) x$meta$diagnosis$mm_wt))
      },

      # allow_universal_output_folder_access ----

      #' @description Make output folder available to all users (Linux specific).
      #' @return The invisible self for chaining.
      allow_universal_output_folder_access = function() {
        if (Sys.info()["sysname"] == "Linux") system2("chmod", paste0("ugo+rwx ", self$design$sim_prm$output_dir)) else
          message("This function is only available in Linux.")
        invisible(self)
      },

      # update_output_path ----

      #' @description Updates the output path.
      #' @param new_path A string with the new output path (absolute).
      #' @return The invisible self for chaining.
      update_output_path = function(new_path) {
        if (!is.character(new_path)) stop("new_path needs to be a string.")
        self$design$sim_prm$output_dir <- new_path
        invisible(self)
      },

      # update_synthpop_path ----

      #' @description Updates the synthpop path.
      #' @param new_path A string with the new synthpop path (absolute).
      #' @return The invisible self for chaining.
      update_synthpop_path = function(new_path) {
        if (!is.character(new_path)) stop("new_path needs to be a string.")
        self$design$sim_prm$synthpop_dir <- new_path
        invisible(self)
      },


      # print ----
      # ? ARU check
      #' @description Prints the simulation object metadata.
      #' @return The invisible `SynthPop` object.
      print = function() {
        print(c(
          "TODO..."
        ))
        invisible(self)
      }
    ),



    # private -----------------------------------------------------------------
    private = list(
      synthpop_dir = NA,
      causality_structure = NA,
      death_codes = NA,
      # diseasenam_hlp = NA,
      esp_weights = data.table(),

      # run_sim ----
      # ? ARU check
      # Runs the simulation in one core. mc is scalar
      # @param mc_ Monte Carlo iteration number.
      # @param scenario_nam Name of the scenario. Defaults to "sc0" if not provided.
      #
      # @details
      # This function performs the following steps:
      # - Initializes a synthetic population.
      # - Applies mortality calibration to the population.
      # - Applies primary and secondary prevention scenarios.
      # - Initializes the scenario and runs the simulation using simcpp.
      # - Updates population weights based on the scenario.
      # - Prunes the population based on specified criteria.
      # - Applies weights.
      # - Exports exposures if export is enabled.
      # - Writes the lifecourse data to a CSV or compressed CSV file.
      #
      run_sim = function(mc_, scenario_nam = "") {
        if (!nzchar(scenario_nam)) scenario_nam <- "sc0"

        if (self$design$sim_prm$logs) {
          private$time_mark(paste0("Start mc iteration ", mc_))
          sink(
            file = private$output_dir(paste0("logs/log", mc_, ".txt")),
            append = TRUE,
            type = "output",
            split = FALSE
          )
        }

        sp <- SynthPop$new(mc_, self$design)
        e <- read_fst("./inputs/mortality/mrtl_clb.fst", as.data.table = TRUE) # mortality calibration
        lookup_dt(sp$pop, e, check_lookup_tbl_validity = self$design$sim_prm$logs)
        setnafill(sp$pop, "const", 1, cols = "mrtl_clbr")
        rm(e)

        # From Karl, somehow scenario_fn() makes init_prvl different. The
        # following code solves the problem.
        # TODO: investigate the root cause


        lapply(self$diseases, function(x) {
          if (self$design$sim_prm$logs) print(x$name)
          x$
            gen_parf(sp, self$design, self$diseases)$
            set_init_prvl(sp = sp, design_ = self$design)
        })

        scenario_fn_primary_prevention(sp) # apply primary pevention scenario

        lapply(self$diseases, function(x) {
          x$set_rr(sp, self$design)$
            set_incd_prb(sp, self$design)$
            set_dgns_prb(sp, self$design)$
            set_mrtl_prb(sp, self$design)
        })

        scenario_fn_secondary_prevention(sp) # apply secondary pevention scenario

        # ds <- copy(self$diseases) # Necessary for parallelisation
        # lapply(self$diseases, function(x) {
        #   if (self$design$sim_prm$logs) print(x$name)
        #   x$gen_parf(sp, self$design, self$diseases)$
        #     set_init_prvl(sp, self$design)$
        #     set_rr(sp, self$design)$
        #     set_incd_prb(sp, self$design)$
        #     set_dgns_prb(sp, self$design)$
        #     set_mrtl_prb(sp, self$design)
        # })

        l <- private$mk_scenario_init(sp, scenario_nam)
        if (!identical(key(sp$pop), c("pid", "year"))) stop("synthpop key is not as expected")
        simcpp(sp$pop, l, sp$mc)
        # it doesn't matter if mc or mc_aggr is used in the above, because it is
        # only used for the RNG stream and the pid are different in each mc_aggr
        # pop

        sp$update_pop_weights(scenario_nam)

        # Prune pop (NOTE that assignment in the function env makes this
        # data.table local)
        sp$pop <- sp$pop[all_cause_mrtl >= 0L &
                           year >= self$design$sim_prm$init_year &
                           between(age, self$design$sim_prm$ageL, self$design$sim_prm$ageH), ]
        setkey(sp$pop, pid, year)
        sp$pop[, pid_mrk := mk_new_simulant_markers(pid)]

        # apply ESP weights
        to_agegrp(sp$pop, 5, 99)
        absorb_dt(sp$pop, private$esp_weights)
        sp$pop[, wt_esp := wt_esp * unique(wt_esp) / sum(wt_esp),
               by = .(year, agegrp, sex, dimd)] # NOTE keyby changes the key

        if (self$design$sim_prm$export_xps) {
          if (self$design$sim_prm$logs) message("Exporting exposures...")
          private$export_xps(sp, scenario_nam)
        }

        nam <- c(self$design$sim_prm$cols_for_output,
                 grep("^cms_|_prvl$|_dgns$|_mrtl$", names(sp$pop), value = TRUE))
        nam <- grep("^prb_", nam, value = TRUE, invert = TRUE) # exclude prb_ ... _dgns
        sp$pop[, setdiff(names(sp$pop), nam) := NULL]
        sp$pop[, mc := sp$mc_aggr]


        # TODO add logic for the years of having MM. Currently 1 is not the real
        # incidence. It is still prevalence
        sp$pop[, `:=` (
          cms1st_cont_prvl   = carry_forward_incr(as.integer(cms_count == 1),
                                                  pid_mrk, TRUE, 1L, byref = TRUE),
          cmsmm0_prvl   = carry_forward_incr(as.integer(cms_score > 0),
                                             pid_mrk, TRUE, 1L, byref = TRUE),
          cmsmm1_prvl   = carry_forward_incr(as.integer(cms_score > 1),
                                             pid_mrk, TRUE, 1L, byref = TRUE),
          cmsmm1.5_prvl = carry_forward_incr(as.integer(cms_score > 1.5),
                                             pid_mrk, TRUE, 1L, byref = TRUE),
          cmsmm2_prvl   = carry_forward_incr(as.integer(cms_score > 2),
                                             pid_mrk, TRUE, 1L, byref = TRUE)
        )]

        sp$pop[, scenario := scenario_nam]

        setkeyv(sp$pop, c("pid", "year"))

        # Write lifecourse
        if (self$design$sim_prm$logs) message("Exporting lifecourse...")

        if (self$design$sim_prm$avoid_appending_csv) {
          fnam <- private$output_dir(paste0(
            "lifecourse/", sp$mc_aggr, "_", sp$mc, "_lifecourse.csv"
          ))
        } else {
          fnam <- private$output_dir(paste0(
            "lifecourse/", sp$mc_aggr, "_lifecourse.csv.gz"
          ))
        }
        fwrite_safe(sp$pop, fnam)


        if (self$design$sim_prm$logs) {
          private$time_mark(paste0("End mc iteration ", mc_))
          sink()
        }

        NULL
      },

      # creates the list that is used in c++ side sp is needed for sp$mc_aggr in
      # to_cpp()

      # mk_scenario_init ----
      # functions - ⁠carry_forward_inc;
      # Create initialization parameters for a simulation scenario.
      # @param sp A SynthPop object.
      # @param scenario_name A character string specifying the name of the scenario.
      #
      # @return A list of initialization parameters for the scenario.
      #
      mk_scenario_init = function(sp, scenario_name) {
        # scenario_suffix_for_pop <- paste0("_", scenario_name)

        # TODO the next line counteracts the commented line above. This is
        # intentional until we finalise the scenario mechanism
        scenario_suffix_for_pop <- ""

        list(
          "exposures"          = self$design$sim_prm$exposures,
          "scenarios"          = self$design$sim_prm$scenarios, # to be generated programmatically
          "scenario"           = scenario_name,
          "kismet"             = self$design$sim_prm$kismet, # If TRUE random numbers are the same for each scenario.
          "init_year"          = self$design$sim_prm$init_year,
          "pids"               = "pid",
          "years"              = "year",
          "ages"               = "age",
          "sexs"               = "sex",
          "dimds"              = "dimd",
          "ageL"               = self$design$sim_prm$ageL,
          "all_cause_mrtl"     = paste0("all_cause_mrtl", scenario_suffix_for_pop),
          "cms_score"          = paste0("cms_score", scenario_suffix_for_pop),
          "cms_count"          = paste0("cms_count", scenario_suffix_for_pop),
          # "strata_for_outputs" = c("pid", "year", "age", "sex", "dimd"),
          "diseases"           = lapply(self$diseases, function(x)
            x$to_cpp(sp, self$design, scenario_name, scenario_suffix_for_pop))
        )
      },

      # export_xps ----
      # ? ARU check details
      # Export Exposure Data
      # @param sp A SynthPop object.
      # @param scenario_nam A character string specifying the scenario name.
      # @return NULL
      #
      # @details This function exports exposure data from a SynthPop object to CSV files.
      #
      export_xps = function(sp, scenario_nam) {
        # NOTE no need to check validity of inputs here as it is only used
        # internally

        to_agegrp(sp$pop, grp_width = 20L, max_age = self$design$sim_prm$ageH,
                  min_age = self$design$sim_prm$ageL, age_colname = "age",
                  agegrp_colname = "agegrp20", to_factor = TRUE)

        sp$pop[, smok_never_curr_xps := fifelse(smok_status_curr_xps == "1", 1L, 0L)]
        sp$pop[, smok_active_curr_xps := fifelse(smok_status_curr_xps == "4", 1L, 0L)]

        xps <- grep("_curr_xps$", names(sp$pop), value = TRUE)
        xps <- grep("_prvl_curr_xps$", xps, value = TRUE, invert = TRUE)
        xps <- xps[-which(xps %in% c("smok_status_curr_xps", "met_curr_xps",
                                     "bpmed_curr_xps"))]
        sp$pop[smok_status_curr_xps == "1", `:=` (
          smok_packyrs_curr_xps = NA,
          smok_quit_yrs_curr_xps = NA,
          smok_dur_curr_xps = NA,
          smok_cig_curr_xps = NA
        )]
        sp$pop[smok_status_curr_xps == "4", `:=` (
          smok_quit_yrs_curr_xps = NA)]

        out_xps20 <- groupingsets(
          sp$pop[all_cause_mrtl >= 0L &
                   year >= self$design$sim_prm$init_year &
                   age >= self$design$sim_prm$ageL, ],
          j = lapply(.SD, weighted.mean, wt, na.rm = TRUE),
          by = c("year", "sex", "agegrp20", "qimd"), # "ethnicity", "sha"
          .SDcols = xps,
          sets = list(
            "year",
            c("year", "agegrp20"),
            c("year", "sex"),
            c("year", "qimd"),
            c("year", "agegrp20", "sex"),
            c("year", "sex", "agegrp20", "qimd")

            # c("year", "ethnicity"),
            # c("year", "sha")
          )
        )[, `:=` (year = year + 2000L, mc = sp$mc, scenario = scenario_nam)]
        # TODO above mc could also be mc_aggr. Getting the uncertainty right here is tricky

        for (j in seq_len(ncol(out_xps20)))
          set(out_xps20, which(is.na(out_xps20[[j]])), j, "All")
        setkey(out_xps20, year)
        if (self$design$sim_prm$avoid_appending_csv) {
          fwrite_safe(out_xps20, private$output_dir(paste0("xps/", sp$mc, "_xps20.csv")))
        } else {
          fwrite_safe(out_xps20, private$output_dir("xps/xps20.csv.gz"))
        }

        # TODO link strata in the outputs to the design.yaml
        out_xps5 <- groupingsets(
          sp$pop[all_cause_mrtl >= 0L &
                   year >= self$design$sim_prm$init_year &
                   age >= self$design$sim_prm$ageL, ],
          j = lapply(.SD, weighted.mean, wt_esp, na.rm = TRUE),
          by = c("year", "sex", "qimd"), # "ethnicity", "sha"
          .SDcols = xps,
          sets = list(
            "year",
            c("year", "sex"),
            c("year", "qimd"),
            c("year", "sex", "qimd")
            # c("year", "ethnicity"),
            # c("year", "sha")
          )
        )[, `:=` (year = year + 2000L, mc = sp$mc, scenario = scenario_nam)]
        for (j in seq_len(ncol(out_xps5)))
          set(out_xps5, which(is.na(out_xps5[[j]])), j, "All")
        setkey(out_xps5, year)
        if (self$design$sim_prm$avoid_appending_csv) {
          fwrite_safe(out_xps5, private$output_dir(paste0("xps/", sp$mc, "_xps_esp.csv")))
        } else {
          fwrite_safe(out_xps5, private$output_dir("xps/xps_esp.csv.gz"))
        }

        # Tidy up
        sp$pop[, c(
          "agegrp20",
          "smok_never_curr_xps",
          "smok_active_curr_xps"
        ) := NULL]
        sp$pop[smok_status_curr_xps == "1", `:=` (
          smok_packyrs_curr_xps = 0,
          smok_quit_yrs_curr_xps = 0,
          smok_dur_curr_xps = 0,
          smok_cig_curr_xps = 0
        )]
        sp$pop[smok_status_curr_xps == "4", `:=` (
          smok_quit_yrs_curr_xps = 0)]

        NULL
      },

      # time_mark ----
      # Function for timing log
      # @param text_id A character string representing the text identifier.
      # @return NULL
      #
      # @details
      # The function logs a timestamp, along with the provided text identifier, to a file.
      # The file is specified as "logs/times.txt" within the output directory. The logging
      # is done in append mode, meaning new entries are added to the end of the existing file.
      # The logged information includes the text identifier, timestamp, and a newline character.
      #
      time_mark = function(text_id) {
        sink(
          file = private$output_dir("logs/times.txt"),
          append = TRUE,
          type = "output",
          split = FALSE
        )
        cat(paste0(text_id, " at: ", Sys.time(), "\n"))
        sink()
      },

      # output_dir ----
      output_dir = function(x = "") {
        file.path(self$design$sim_prm$output_dir, x)
      },

      # export_summaries_hlpr ----
      # function to export summaries from lifecourse files.
      # lc is a lifecourse file
      # ? ARU check short forms
      # Export Summaries Helper Function
      # @param lc Data table containing lifecourse information.
      # @param type Character vector specifying the types of summaries to export. Possible values include "le" (life expectancy), "hle" (healthy life expectancy), "dis_char" (disease characteristics), "prvl" (prevalence), "incd" (incidence), "mrtl" (mortality), "dis_mrtl" (disease-specific mortality), "allcause_mrtl_by_dis" (all-cause mortality by disease), and "cms" (CMS - Chronic Multimorbidity Score).
      # @param design An object representing the design of the simulation. It should contain simulation parameters.
      #
      # @details
      # The function exports various summaries based on the specified types. For example, it exports life expectancy, healthy life expectancy, disease characteristics, prevalence, incidence, mortality, disease-specific mortality, all-cause mortality by disease, and Chronic Multimorbidity Score (CMS) information.
      #
      # @return
      # Returns an invisible updated simulation object.
      #
      export_summaries_hlpr = function(lc, type = c("le", "hle", "dis_char",
                                                    "prvl", "incd", "mrtl",  "dis_mrtl",
                                                    "allcause_mrtl_by_dis", "cms")) {
        if (self$design$sim_prm$logs) message("Exporting summaries...")

        strata <- c("mc", self$design$sim_prm$strata_for_output)
        strata_noagegrp <- c("mc",
                             setdiff(self$design$sim_prm$strata_for_output, c("agegrp")))
        strata_age <- c(strata_noagegrp, "age")

        setkeyv(lc, c("scenario", "pid", "year")) # necessary for age_onset

        mcaggr <- ifelse(self$design$sim_prm$avoid_appending_csv, paste0(lc$mc[1], "_"), "")
        ext <- ifelse(self$design$sim_prm$avoid_appending_csv, ".csv", ".csv.gz")

        # Life expectancy ----
        # NOTE for scaled_up LE weights need to apply from the very beginning.
        # Also note that currently this ignores the deaths for people younger
        # than min_age so not a true LE at birth
        if ("le" %in% type) {
          # fwrite_safe(lc[all_cause_mrtl > 0, .("popsize" = (.N), LE = mean(age)),
          #                keyby = strata_noagegrp],
          #             private$output_dir(paste0("summaries/", "le_out.csv.gz"
          #             )))
          fwrite_safe(lc[all_cause_mrtl > 0, .("popsize" = sum(wt), LE = weighted.mean(age, wt)),  keyby = strata_noagegrp],
                      private$output_dir(paste0("summaries/", mcaggr, "le_scaled_up", ext
                      )))
          fwrite_safe(lc[all_cause_mrtl > 0, .("popsize" = sum(wt_esp), LE = weighted.mean(age, wt_esp)),  keyby = strata_noagegrp],
                      private$output_dir(paste0("summaries/", mcaggr, "le_esp", ext
                      )))
          # Life expectancy at 60 ----

          if (self$design$sim_prm$ageL < 60L && self$design$sim_prm$ageH > 60L) {
            # fwrite_safe(lc[all_cause_mrtl > 0 & age > 60, .("popsize" = (.N), LE60 = mean(age)),  keyby = strata_noagegrp],
            #             private$output_dir(paste0("summaries/", "le60_out", ext
            #             )))
            fwrite_safe(lc[all_cause_mrtl > 0 & age > 60, .("popsize" = sum(wt), LE60 = weighted.mean(age, wt)),  keyby = strata_noagegrp],
                        private$output_dir(paste0("summaries/", mcaggr, "le60_scaled_up", ext
                        )))
            fwrite_safe(lc[all_cause_mrtl > 0 & age > 60, .("popsize" = sum(wt_esp), LE60 = weighted.mean(age, wt_esp)),  keyby = strata_noagegrp],
                        private$output_dir(paste0("summaries/", mcaggr, "le60_esp", ext
                        )))
          }
          # Note: for less aggregation use wtd.mean with popsize i.e le_out[,
          # weighted.mean(LE, popsize), keyby = year]
        }

        # Healthy life expectancy ----
        if ("hle" %in% type) {
          # TODO currently some individuals are counted more than once because
          # disease counter and score can be reduced.
          # Ideally only the first reach to the threshold should be counted
          # fwrite_safe(lc[cms_count == 1L, .("popsize" = (.N), HLE = mean(age)),
          #                keyby = strata_noagegrp],
          #             private$output_dir(paste0("summaries/", "hle_1st_cond_out", ext)))
          fwrite_safe(lc[cms_count == 1L,
                         .("popsize" = sum(wt), HLE = weighted.mean(age, wt)),
                         keyby = strata_noagegrp],
                      private$output_dir(paste0(
                        "summaries/", mcaggr, "hle_1st_cond_scaled_up", ext
                      )))
          fwrite_safe(lc[cms_count == 1L,
                         .("popsize" = sum(wt_esp), HLE = weighted.mean(age, wt_esp)),
                         keyby = strata_noagegrp],
                      private$output_dir(paste0("summaries/", mcaggr, "hle_1st_cond_esp", ext
                      )))

          # fwrite_safe(lc[cmsmm1.5_prvl == 1L, .("popsize" = (.N), HLE = mean(age)),
          #                keyby = strata_noagegrp],
          #             private$output_dir(paste0("summaries/", "hle_cmsmm1.5_out", ext)))
          fwrite_safe(lc[cmsmm1.5_prvl == 1L,
                         .("popsize" = sum(wt), HLE = weighted.mean(age, wt)),
                         keyby = strata_noagegrp],
                      private$output_dir(paste0(
                        "summaries/", mcaggr, "hle_cmsmm1.5_scaled_up", ext
                      )))
          fwrite_safe(lc[cmsmm1.5_prvl == 1L,
                         .("popsize" = sum(wt_esp), HLE = weighted.mean(age, wt_esp)),
                         keyby = strata_noagegrp],
                      private$output_dir(paste0("summaries/", mcaggr, "hle_cmsmm1.5_esp", ext
                      )))
        }

        # Disease characteristics----
        if ("dis_char" %in% type) {
          nm <- grep("_prvl$", names(lc), value = TRUE)

          # tt <- rbindlist(lapply(nm, function(x) {
          #   # sr are the rows the 1st episode occurs per pid
          #   # Need to be sorted on year
          #   sr <- lc[get(x) > 0L, .I[match(1L, get(x))], by = .(pid, scenario)]$V1
          #   sr <- sr[!is.na(sr)]
          #   lc[sr, age_onset := age] # age at 1st ever event
          #
          #   lc[get(x) > 0L, .("disease" = gsub("_prvl$", "", x),
          #                     "cases" = (.N),
          #                     "mean_age_incd" = mean(age[get(x) == 1L]),
          #                     "mean_age_1st_onset" = mean(age_onset, na.rm = TRUE),
          #                     "mean_age_prvl" = mean(age),
          #                     "mean_duration" = mean(get(x)), # Note get(x) very slow here. Implementation with .SDcols also slow because of cases
          #                     "mean_cms_score" = mean(cms_score),
          #                     "mean_cms_count" = mean(cms_count)),
          #      keyby = strata_noagegrp]
          #   lc[, age_onset := NULL]
          # }))
          # tt <- rbindlist(lapply(nm, function(x) {
          #   lc[get(x) > 0L, lapply(.SD, mean), .SDcols = c(x, "age", "cms_score", "cms_count"), keyby = strata_noagegrp]
          # })) # Fast but without cases
          # tt <-
          #   dcast(tt, as.formula(paste0(paste(strata_noagegrp, collapse = "+"), "~disease")),
          #         fill = 0L, value.var = c("cases", "mean_duration", "mean_age_incd",
          #                                  "mean_age_prvl", "mean_cms_score",
          #                                  "mean_cms_count"))
          # fwrite_safe(tt,
          #             private$output_dir(paste0("summaries/", "dis_characteristics_out", ext
          #             )))

          tt <- rbindlist(lapply(nm, function(x) {
            # sr are the rows the 1st episode occurs per pid
            # Need to be sorted on year
            sr <- lc[get(x) > 0L, .I[match(1L, get(x))], by = .(pid, scenario)]$V1
            sr <- sr[!is.na(sr)]
            lc[, wt1st := 0]
            lc[sr, `:=` (age_onset = age, wt1st = wt)] # age at 1st ever event

            ans <- lc[get(x) > 0L, .("disease" = gsub("_prvl$", "", x),
                                     "cases" = sum(wt),
                                     "mean_age_incd" = weighted.mean(age[get(x) == 1L],
                                                                     wt[get(x) == 1L]),
                                     "mean_age_1st_onset" = weighted.mean(age_onset, wt1st, na.rm = TRUE),

                                     "mean_age_prvl" = weighted.mean(age, wt),
                                     "mean_duration" = weighted.mean(get(x), wt), # Note get(x) very slow here. Implementation with .SDcols also slow because of cases
                                     "mean_cms_score" = weighted.mean(cms_score, wt),
                                     "mean_cms_count" = weighted.mean(cms_count, wt)),
                      keyby = strata_noagegrp]
            lc[, c("age_onset", "wt1st") := NULL]
            ans
          }))
          tt <-
            dcast(tt, as.formula(paste0(paste(strata_noagegrp, collapse = "+"), "~disease")),
                  fill = 0L, value.var = c("cases", "mean_duration", "mean_age_incd",
                                           "mean_age_1st_onset",
                                           "mean_age_prvl", "mean_cms_score", "mean_cms_count"))
          fwrite_safe(tt,
                      private$output_dir(paste0("summaries/", mcaggr, "dis_characteristics_scaled_up", ext
                      )))

          tt <- rbindlist(lapply(nm, function(x) {
            # sr are the rows the 1st episode occurs per pid
            # Need to be sorted on year
            sr <- lc[get(x) > 0L, .I[match(1L, get(x))], by = .(pid, scenario)]$V1
            sr <- sr[!is.na(sr)]
            lc[, wt1st := 0]
            lc[sr, `:=` (age_onset = age, wt1st = wt_esp)] # age at 1st ever event

            ans <- lc[get(x) > 0L, .("disease" = gsub("_prvl$", "", x),
                                     "cases" = sum(wt_esp),
                                     "mean_age_incd" = weighted.mean(age[get(x) == 1L],
                                                                     wt_esp[get(x) == 1L]),
                                     "mean_age_1st_onset" = weighted.mean(age_onset, wt1st, na.rm = TRUE),
                                     "mean_age_prvl" = weighted.mean(age, wt_esp),
                                     "mean_duration" = weighted.mean(get(x), wt_esp), # Note get(x) very slow here. Implementation with .SDcols also slow because of cases
                                     "mean_cms_score" = weighted.mean(cms_score, wt_esp),
                                     "mean_cms_count" = weighted.mean(cms_count, wt_esp)),
                      keyby = strata_noagegrp]
            lc[, c("age_onset", "wt1st") := NULL]
            ans
          }))
          tt <-
            dcast(tt, as.formula(paste0(paste(strata_noagegrp, collapse = "+"), "~disease")),
                  fill = 0L, value.var = c("cases", "mean_duration", "mean_age_incd",
                                           "mean_age_1st_onset",
                                           "mean_age_prvl", "mean_cms_score",
                                           "mean_cms_count"))
          fwrite_safe(tt,
                      private$output_dir(paste0("summaries/", mcaggr, "dis_characteristics_esp", ext
                      )))
          rm(tt)
        }

        # prvl ----
        if ("prvl" %in% type) {
          # Note for mortality this exports qx directly. mx is defined as the
          # number of deaths during the year divided by the average number alive
          # during the year, i.e. This differs slightly from qx , which is the
          # number of deaths during the year divided by the number alive at the
          # beginning of the year.
          # fwrite_safe(lc[, c("popsize" = (.N),
          #                    lapply(.SD, function(x) sum(x > 0))),
          #                .SDcols = patterns("_prvl$"), keyby = strata],
          #             private$output_dir(paste0("summaries/", "prvl_out", ext
          #             )))
          fwrite_safe(lc[, c("popsize" = sum(wt),
                             lapply(.SD, function(x, wt) sum((x > 0) * wt), wt)),
                         .SDcols = patterns("_prvl$"), keyby = strata],
                      private$output_dir(paste0("summaries/", mcaggr, "prvl_scaled_up", ext
                      )))
          fwrite_safe(lc[, c("popsize" = sum(wt_esp),
                             lapply(.SD, function(x, wt) sum((x > 0) * wt), wt_esp)),
                         .SDcols = patterns("_prvl$"), keyby = strata],
                      private$output_dir(paste0("summaries/", mcaggr, "prvl_esp", ext
                      )))
        }

        # incd ----
        if ("incd" %in% type) {
          # NOTE incd includes prevalent cases in denominator
          # fwrite_safe(lc[, c("popsize" = (.N),
          #                    lapply(.SD, function(x) sum(x == 1))),
          #                .SDcols = patterns("_prvl$"), keyby = strata],
          #             private$output_dir(paste0("summaries/", "incd_out", ext
          #             )))
          incdtbl <- lc[, c("popsize" = sum(wt),
                            lapply(.SD, function(x, wt) sum((x == 1) * wt), wt)),
                        .SDcols = patterns("_prvl$"), keyby = strata]
          nm <- grep("_prvl$", names(incdtbl), value = TRUE)
          setnames(incdtbl, nm, gsub("_prvl$", "_incd", nm))
          fwrite_safe(incdtbl,
                      private$output_dir(paste0("summaries/", mcaggr, "incd_scaled_up", ext
                      )))

          incdtbl <- lc[, c("popsize" = sum(wt_esp),
                            lapply(.SD, function(x, wt) sum((x == 1) * wt), wt_esp)),
                        .SDcols = patterns("_prvl$"), keyby = strata]
          nm <- grep("_prvl$", names(incdtbl), value = TRUE)
          setnames(incdtbl, nm, gsub("_prvl$", "_incd", nm))
          fwrite_safe(incdtbl,
                      private$output_dir(paste0("summaries/", mcaggr, "incd_esp", ext
                      )))

          rm(incdtbl, nm)
        }

        # mrtl ----
        if ("mrtl" %in% type) {
          # fwrite_safe(lc[, .("popsize" = (.N),
          #                    "all_cause_mrtl" = sum(all_cause_mrtl > 0)),
          #                keyby = strata],
          #             private$output_dir(paste0("summaries/", "mrtl_out", ext
          #             )))
          fwrite_safe(lc[, .("popsize" = sum(wt),
                             "all_cause_mrtl" = sum((all_cause_mrtl > 0) * wt)),
                         keyby = strata],
                      private$output_dir(paste0("summaries/", mcaggr, "mrtl_scaled_up", ext
                      )))
          fwrite_safe(lc[, .("popsize" = sum(wt_esp),
                             "all_cause_mrtl" = sum((all_cause_mrtl > 0) * wt_esp)),
                         keyby = strata],
                      private$output_dir(paste0("summaries/", mcaggr, "mrtl_esp", ext
                      )))
        }

        # disease specific mortality ----
        if ("dis_mrtl" %in% type) {
          # dis_mrtl_out <-
          #   dcast(
          #     lc[, .("deaths" = (.N)),
          #        keyby = c(strata, "all_cause_mrtl")],
          #     formula = as.formula(paste0(
          #       paste(strata, collapse = "+"), "~all_cause_mrtl"
          #     )),
          #     fill = 0L,
          #     value.var = "deaths"
          #   )
          #
          # setnames(dis_mrtl_out, as.character(private$death_codes),
          #          names(private$death_codes), skip_absent = TRUE)
          # dis_mrtl_out[, `:=` (
          #   popsize = Reduce(`+`, .SD),
          #   alive = NULL
          # ), .SDcols = !strata]
          # fwrite_safe(dis_mrtl_out,
          #             private$output_dir(paste0("summaries/", "dis_mrtl_out", ext
          #             )))

          dis_mrtl_out <- # scale up
            dcast(
              lc[, .("deaths" = sum(wt)),
                 keyby = c(strata, "all_cause_mrtl")],
              formula = as.formula(paste0(
                paste(strata, collapse = "+"), "~all_cause_mrtl"
              )),
              fill = 0L,
              value.var = "deaths"
            )

          setnames(dis_mrtl_out, as.character(private$death_codes),
                   paste0(names(private$death_codes), "_deaths"), skip_absent = TRUE)
          dis_mrtl_out[, `:=` (
            popsize = Reduce(`+`, .SD), # it includes alive so it is the pop at the start of the year
            alive_deaths = NULL
          ), .SDcols = !strata]
          fwrite_safe(dis_mrtl_out,
                      private$output_dir(paste0("summaries/", mcaggr, "dis_mrtl_scaled_up", ext
                      )))

          dis_mrtl_out <- # scale up esp
            dcast(
              lc[, .("deaths" = sum(wt_esp)),
                 keyby = c(strata, "all_cause_mrtl")],
              formula = as.formula(paste0(
                paste(strata, collapse = "+"), "~all_cause_mrtl"
              )),
              fill = 0L,
              value.var = "deaths"
            )

          setnames(dis_mrtl_out, as.character(private$death_codes),
                   paste0(names(private$death_codes), "_deaths"), skip_absent = TRUE)
          dis_mrtl_out[, `:=` (
            popsize = Reduce(`+`, .SD),
            alive_deaths = NULL
          ), .SDcols = !strata]
          fwrite_safe(dis_mrtl_out,
                      private$output_dir(paste0("summaries/", mcaggr, "dis_mrtl_esp", ext
                      )))
          rm(dis_mrtl_out)
        }

        # All-cause mrtl by disease ----
        if ("allcause_mrtl_by_dis" %in% type) {
          nm <- grep("_prvl$", names(lc), value = TRUE)

          # tt <- lapply(nm, function(x) {
          #   lc[get(x) > 0L, .("disease" = gsub("_prvl$", "", x), "cases" = (.N), "deaths" = sum(all_cause_mrtl > 0)), keyby = strata]
          # })
          # tt <-
          # dcast(rbindlist(tt), as.formula(paste0(paste(strata, collapse = "+"), "~disease")),
          #       fill = 0L, value.var = c("deaths", "cases"))
          # fwrite_safe(tt,
          #             private$output_dir(paste0("summaries/", "all_cause_mrtl_by_dis_out", ext
          #             )))

          tt <- lapply(nm, function(x) {
            lc[get(x) > 0L, .("disease" = gsub("_prvl$", "", x), "cases" = sum(wt), "deaths" = sum(wt * (all_cause_mrtl > 0))), keyby = strata]
          })
          tt <-
            dcast(rbindlist(tt), as.formula(paste0(paste(strata, collapse = "+"), "~disease")),
                  fill = 0L, value.var = c("deaths", "cases"))
          fwrite_safe(tt,
                      private$output_dir(paste0("summaries/", mcaggr, "all_cause_mrtl_by_dis_scaled_up", ext
                      )))

          tt <- lapply(nm, function(x) {
            lc[get(x) > 0L, .("disease" = gsub("_prvl$", "", x), "cases" = sum(wt_esp), "deaths" = sum(wt_esp * (all_cause_mrtl > 0))), keyby = strata]
          })
          tt <-
            dcast(rbindlist(tt), as.formula(paste0(paste(strata, collapse = "+"), "~disease")),
                  fill = 0L, value.var = c("deaths", "cases"))
          fwrite_safe(tt,
                      private$output_dir(paste0("summaries/", mcaggr, "all_cause_mrtl_by_dis_esp", ext
                      )))
          rm(tt)
        }


        # CMS mean ----
        if ("cms" %in% type) {
          # fwrite_safe(lc[, .("popsize" = (.N), cms_score = mean(cms_score)), keyby = strata],
          #             private$output_dir(paste0("summaries/", "cms_score_out", ext
          #             )))
          fwrite_safe(lc[, .("popsize" = sum(wt), cms_score = weighted.mean(cms_score, wt)),  keyby = strata],
                      private$output_dir(paste0("summaries/", mcaggr, "cms_score_scaled_up", ext
                      )))
          fwrite_safe(lc[, .("popsize" = sum(wt), cms_score = weighted.mean(cms_score, wt)),  keyby = strata_age],
                      private$output_dir(paste0("summaries/", mcaggr, "cms_score_by_age_scaled_up", ext
                      )))

          fwrite_safe(lc[, .("popsize" = sum(wt_esp), cms_score = weighted.mean(cms_score, wt_esp)),  keyby = strata],
                      private$output_dir(paste0("summaries/", mcaggr, "cms_score_esp", ext
                      )))

          fwrite_safe(lc[, .("popsize" = sum(wt_esp), cms_score = weighted.mean(cms_score, wt_esp)),  keyby = strata_age],
                      private$output_dir(paste0("summaries/", mcaggr, "cms_score_by_age_esp", ext
                      )))

          # CMS count ----
          # fwrite_safe(lc[, .("popsize" = (.N), cms_count = mean(cms_count)), keyby = strata],
          #             private$output_dir(paste0("summaries/", "cms_count_out", ext
          #             )))
          fwrite_safe(lc[, .("popsize" = sum(wt), cms_count = weighted.mean(cms_count, wt)),  keyby = strata],
                      private$output_dir(paste0("summaries/", mcaggr, "cms_count_scaled_up", ext
                      )))
          fwrite_safe(lc[, .("popsize" = sum(wt_esp), cms_count = weighted.mean(cms_count, wt_esp)),  keyby = strata],
                      private$output_dir(paste0("summaries/", mcaggr, "cms_count_esp", ext
                      )))
        }

        if (!self$design$sim_prm$keep_lifecourse) file.remove(pth)

        return(invisible(self))
      },

      # Special deep copy for data.table. Use POP$clone(deep = TRUE) to
      # dispatch. Otherwise a reference is created
      # @param name A character string representing the name or identifier of the object being cloned.
      # @param value The object to be cloned.
      #
      # @details
      # The function checks the class of the input object and performs the appropriate cloning method. For data tables, it uses data.table::copy, and for R6 objects, it calls the clone() method. For other object types, a shallow copy is returned.
      #
      # @return
      # A deep clone of the input object.
      #
      deep_clone = function(name, value) {
        if ("data.table" %in% class(value)) {
          data.table::copy(value)
        } else if ("R6" %in% class(value)) {
          value$clone()
        } else {
          # For everything else, just return it. This results in a shallow copy
          # of s3.
          value
        }
      },

      # collect_files ----
      # Collect files written by mc_aggr or mc_aggr_mc in a folder and combine
      # them into one file
      # @param folder_name A character string specifying the name of the folder containing the files to be collected.
      # @param pattern A regular expression pattern to filter files for collection. Only files matching this pattern will be considered.
      # @param to_mc_aggr Logical, indicating whether to aggregate files for Monte Carlo iteration (to_mc_aggr = TRUE) or not.
      #
      # @details
      # The function reads each file from the specified folder, aggregates its content, and, if needed, compresses the aggregated result. The aggregation process involves rewriting the files by applying certain string patterns to the filenames.
      #
      # @return
      # This function does not return a value directly. It collects and aggregates files, and compresses them if specified, cleaning up the original files in the process.
      #
      collect_files = function(folder_name, pattern = NULL, to_mc_aggr = FALSE) {
        if (self$design$sim_prm$logs) message("Collecting mc files...")
        if (to_mc_aggr) {
          string1 <- "_[0-9]+_"
          string2 <- "_"
        } else {
          string1 <- "[0-9]+_"
          string2 <- ""
        }
        sapply(
          list.files(path = private$output_dir(folder_name), pattern = pattern, full.names = TRUE),
          function(fnam) {
            fwrite_safe(fread(fnam), file = sub(string1, string2, fnam))
            file.remove(fnam)
          }
        )
        # gzip the .csv files to .csv.gz (faster than using gzip() and same speed/compression as with fst 80. But fst reads faster)
        if (self$design$sim_prm$logs) message("Compressing aggregated files...")
        sapply(
          list.files(path = private$output_dir(folder_name), pattern = sub("^_", "", pattern), full.names = TRUE),
          function(fnam) {
            fwrite_safe(fread(fnam), file = gsub(".csv$", ".csv.gz", fnam))
            file.remove(fnam)
          }
        )
        NULL
      },

      # create_new_folder ----
      # @details
      # The function checks if the specified folder exists. If it does not exist, it creates the directory using the \code{dir.create} function. The \code{recursive} argument is set to TRUE to create parent directories if needed.
      #
      # @return
      # This function does not return a value directly. It creates the specified folder and reports the status if specified.
      #
      # @description Create folder if doesn't exist. Stops on failure.
      # @param sDirPathName String folder path and name.
      # @param bReport Bool report folder creation.
      create_new_folder = function(sDirPathName,bReport) {
        if (!dir.exists(sDirPathName)) {
          bSuccess <- dir.create(sDirPathName, recursive=TRUE)
          if (!bSuccess) stop (paste("Failed creating directory",sDirPathName))
          if (bReport) message(paste0("Folder ",sDirPathName," was created"))
        }
      }
    )
  )

