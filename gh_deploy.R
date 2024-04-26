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
#' @description Deploy large GitHub asset files into their relevant locations.
#' 	Additional data is held in an asset config YAML file. BACKGROUND: Standard GitHub repositories do not allow large files exceeding 100 MB (as of 22-08-05); however, such files ('assets'), each up to 2 GB in size, may be attached separately to the repository.
#' 	May either execute from the console:
#' 	Rscript gh_deploy.R [<AssetFilePathName> [<GitHubAssetRouteId> [<GitHubToken>]]]
#' 	or alternatively, may execute directly within R or RStudio:
#' 		source("gh_deploy.R")
#' 		DeployGitHubAssets(sAssetConfigFilePath=<AssetFilePathName>,
#' 			sGitHubAssetRouteId=<GitHubAssetRouteId>,sToken=<GitHubToken>)
#' 	where in the above, <AssetFilePathName> is an optional asset config file's path and name;
#' 			if omitted, tries the default ./auxil/ghAssetConfig.yaml.
#' 	<GitHubAssetRouteId> is an asset route's [id] within <AssetFilePathName>;
#' 		if omitted, seeks ID matching Sys.info()[['user']] name. If <AssetFilePathName> also omitted, assumes running from console.
#' 	<GitHubToken> is a GitHub personal access token (PAT);
#' 		if omitted, seeks token from gh::gh_token().
#' @param sAssetConfigFilePath string (optional): path to ghAssetConfig.yaml file - provided for interactive calls (RStudio).
#' @param sGitHubAssetRouteId string (optional): ID designating asset route in asset config.yaml file - provided for interactive calls (RStudio).
#' @param sToken string (optional): GitHub personal access token (PAT).
DeployGitHubAssets <- function(sToken = NULL, sRepo = NULL, sTag = NULL,
                               sDeployToRootDirPath = NULL,
                               bOverwriteFilesOnDeploy,
                               sUploadSrcDirPath = NULL) {
  # Added repo information here as the function here is more accessible than `GetGitHubAssetRouteInfo` and required for pb_download as well
  # set the GitHub repo
  sRepo <- "ChristK/IMPACTncd_Engl"
  # sUploadSrcDirPath only in `GetGitHubAssetRouteInfo` (as this function used for upload and deploy) removed from `DeployGitHubAssets`
  GetGitHubAssetRouteInfo(sRepo = sRepo, sTag = sTag, iTestWithFirstNAssets,
                          sDeployToRootDirPath = sDeployToRootDirPath,
                          sUploadSrcDirPath = sUploadSrcDirPath,
                          bOverwriteFilesOnDeploy = T, sToken = sToken)
  sDeployToRootDirPath <- TrimSlashes(sDeployToRootDirPath, bRidStartSlash = FALSE)
  # get table which maps SANITISED to ORIGINAL filenames
  sanitisedToOriginalFilePaths <- fread(file.path(sDeployToRootDirPath,
                                                  "/simulation/filindx.csv"),
                                        key = "orig_file")
  if (iTestWithFirstNAssets != 0) { # if testing, only use first N assets
    sanitisedToOriginalFilePaths <- sanitisedToOriginalFilePaths[1:iTestWithFirstNAssets, ]
  }
  # create sub-directories below root directory
  subDirectoryPaths <- sapply(sanitisedToOriginalFilePaths$rel_dir, TrimSlashes)
  sapply(file.path(sDeployToRootDirPath, unique(subDirectoryPaths)),
         dir.create, showWarnings = FALSE, recursive = TRUE)
  # download each Github asset into appropriate sub-directory
  # OLD: finalAssetPaths<- file.path(sDeployToRootDirPath,subDirectoryPaths, sanitisedToOriginalFilePaths$orig_file)
  lsHttpResponses <- piggyback::pb_download(
    file = sanitisedToOriginalFilePaths$sanit_file,
    dest = file.path(sDeployToRootDirPath, subDirectoryPaths,
                     sanitisedToOriginalFilePaths$orig_file),
    repo = sRepo, tag = sTag, overwrite = bOverwriteFilesOnDeploy,
    use_timestamps = FALSE, .token = sToken)
  # 1. WARNING! confusing pb_download() behaviour: 'dest=' *directory* required for single asset; destination *filenames* required for multiple assets.
  # 2. piggyback:: prefix necessary to use R.utils::reassignInPackage() injected code modification.
  StopOnHttpFailure(lsHttpResponses, FALSE)
  # OLD: restore original filename
  # file.rename(finalAssetPaths, file.path(sDeployToRootDirPath,subDirectoryPaths,sanitisedToOriginalFilePaths$orig_file))
}

# if (sys.nframe() == 0) DeployGitHubAssets() # execute with defaults if run from topmost frame (under Rscript)
DeployGitHubAssets(sTag = "v0.0.5", bOverwriteFilesOnDeploy = 1,
                   sDeployToRootDirPath = "D:/Dropbox/ILKConsultancy/IMPACTncd_Engl",
                   sUploadSrcDirPath = "D:/Dropbox/ILKConsultancy/IMPACTncd_Engl",
                   sToken = gh_token())
