# changed this to have similarity with gh_deploy.R
if (!require(piggyback)) {
  if (!nzchar(system.file(package = "pak"))) install.packages("pak")
  dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)
  pak::pkg_install("piggyback", lib = Sys.getenv("R_LIBS_USER"))
  library(piggyback)
}
if (!require(data.table)) {
  dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)
  pak::pkg_install("data.table", lib = Sys.getenv("R_LIBS_USER"))
  library(data.table)
}
source("ghAssetUtils.R")

#' @description Scan a specified directory for large asset files and upload these to GitHub.
#' 	Additional data is held in an asset config YAML file. BACKGROUND: Standard GitHub repositories do not allow large files exceeding 100 MB (as of 22-08-05); however, such files ('assets'), each up to 2 GB in size, may be attached separately to the repository.
#' 	May either execute from the console:
#' 	Rscript gh_upload.R [[<GitHubAssetRouteId> [<GitHubToken>]]]
#' 	or alternatively, may execute directly within R or RStudio:
#' 		source("gh_upload.R")
#' 		UploadGitHubAssets(sToken=<GitHubToken>)
#' 	where in the above, <GitHubToken> is a GitHub personal access token (PAT);
#' 		if omitted, seeks token from gh::gh_token().
#' @param sToken GitHub personal access token (default: NULL).
#' @param sRepo GitHub repository (default: NULL).
#' @param sTag GitHub tag or release (default: NULL).
#' @param sDeployToRootDirPath Path to the root directory for deployment (default: NULL).
#' @param bOverwriteFilesOnDeploy Logical indicating whether to overwrite existing files on deployment (default: FALSE).
#' @param sUploadSrcDirPath Path to the source directory for upload (default: NULL).
#'
#' @return None
UploadGitHubAssets <- function(sToken = NULL, iTestWithFirstNAssets = 0,
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
}

if (sys.nframe() == 0) {
    UploadGitHubAssets(sTag = "v0.0.5", bOverwriteFilesOnDeploy = 0,
                       sDeployToRootDirPath = "D:/Dropbox/ILKConsultancy/IMPACTncd_Engl",
                       sUploadSrcDirPath = "D:/Dropbox/ILKConsultancy/IMPACTncd_Engl",
                       sToken = gh_token())
    #        sAssetConfigFilePath = "./auxil/ghAssetConfig.yaml",
    #       sGitHubAssetRouteId = "local_Chris_IMPACTncd_Engl_0.0.4"
    #  )
} # execute with defaults if run from topmost frame (under Rscript)
