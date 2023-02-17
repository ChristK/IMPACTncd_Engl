if(!require(piggyback)) {
  install.packages("piggyback")
  library(piggyback)
}
if(!require(data.table)) {
  install.packages("data.table")
  library(data.table)
}
source("ghAssetUtils.R")

#' @description Scan a specified directory for large asset files and upload these to GitHub.
#'	Additional data is held in an asset config YAML file. BACKGROUND: Standard GitHub repositories do not allow large files exceeding 100 MB (as of 22-08-05); however, such files ('assets'), each up to 2 GB in size, may be attached separately to the repository.
#'	May either execute from the console:
#' 	Rscript gh_upload.R [<AssetFilePathName> [<GitHubAssetRouteId> [<GitHubToken>]]]
#'	or alternatively, may execute directly within R or RStudio:
#'		source("gh_upload.R")
#'		UploadGitHubAssets(sAssetConfigFilePath=<AssetFilePathName>,
#'			sGitHubAssetRouteId=<GitHubAssetRouteId>,sToken=<GitHubToken>)
#'	where in the above, <AssetFilePathName> is an optional asset config file's path and name;
#'			if omitted, tries the default ./auxil/ghAssetConfig.yaml.
#' 	<GitHubAssetRouteId> is an asset route's [id] within <AssetFilePathName>;
#' 		if omitted, seeks ID matching Sys.info()[['user']] name. If <AssetFilePathName> also omitted, assumes running from console.
#' 	<GitHubToken> is a GitHub personal access token (PAT);
#' 		if omitted, seeks token from gh::gh_token().
#' @param sAssetConfigFilePath string (optional): path to ghAssetConfig.yaml file - provided for interactive calls (RStudio).
#' @param sGitHubAssetRouteId string (optional): ID designating asset route in asset config.yaml file - provided for interactive calls (RStudio).
#' @param sToken string (optional): GitHub personal access token (PAT).
UploadGitHubAssets<- function(sAssetConfigFilePath=NULL,sGitHubAssetRouteId=NULL,sToken=NULL) { # consider as a Simulation class method
	# get GitHub asset route data
	GetGitHubAssetRouteInfo(sId,sRepo,sTag,sUploadSrcDirPath,sDeployToRootDirPath,bOverwriteFilesOnDeploy,
		iTestWithFirstNAssets,sToken=sToken,sAssetConfigFilePath=sAssetConfigFilePath,sGitHubAssetRouteId=sGitHubAssetRouteId)
	sUploadSrcDirPath<- TrimSlashes(sUploadSrcDirPath,bRidStartSlash=FALSE)

	# pb_new_release(sRepo, sTag) # Only need to run the first time a github repo is created
	# pb_release_delete(sRepo, sTag)

	# find local assets, excluding secure_data, "*tmp.qs", and "*parf/PARF_*.qs" files.
	lsLocalAssetPathNames<- list.files(sUploadSrcDirPath, pattern = ".fst$|.xls$|.xlsx$|.qs$", all.files = TRUE, full.names = TRUE, recursive = TRUE)
	lsLocalAssetPathNames<- grep("secure_data", lsLocalAssetPathNames, value = TRUE, invert = TRUE)
	# NOTE currently no necessary files ar .qs. For future proof I add them above
	# and exclude the below
	lsLocalAssetPathNames<- grep("tmp.qs$", lsLocalAssetPathNames, value = TRUE, invert = TRUE)
	lsLocalAssetPathNames<- grep("/simulation/parf/PARF_.*\\.qs$", lsLocalAssetPathNames, value = TRUE, invert = TRUE)
	lsLocalAssetPathNames<- grep("/DELETEME/",lsLocalAssetPathNames,value=TRUE,invert=TRUE)

	# write table with original and *unique sanitised* file names and directories
	dtOriginalAndSanitisedFilePathNames<- data.table(orig_file = basename(lsLocalAssetPathNames),
		sanit_file = gsub("[^[:alnum:]&&^\\.]", ".", basename(lsLocalAssetPathNames)), # replace all not (alphanumeric or '.') with .
		abs_dir = dirname(lsLocalAssetPathNames),
		rel_dir = gsub(sUploadSrcDirPath, "", dirname(lsLocalAssetPathNames)), key = "orig_file") # rel_dir to keep initial slash
	if (any(duplicated(dtOriginalAndSanitisedFilePathNames$sanit_file))) stop("Duplicated filenames found")
	if(iTestWithFirstNAssets!=0)
		dtOriginalAndSanitisedFilePathNames<- dtOriginalAndSanitisedFilePathNames[1:iTestWithFirstNAssets,]
	fwrite(dtOriginalAndSanitisedFilePathNames, "./auxil/filindx.csv")

	lsSrcFilePath<- file.path(dtOriginalAndSanitisedFilePathNames$abs_dir,dtOriginalAndSanitisedFilePathNames$orig_file)
 lsHttpResponses <- piggyback::pb_upload(lsSrcFilePath,
     repo = sRepo, tag = sTag, name = dtOriginalAndSanitisedFilePathNames$sanit_file,
     .token = sToken, use_timestamps = FALSE, overwrite = TRUE
 ) # piggyback:: prefix necessary to use R.utils::reassignInPackage() injected code modification

	StopOnHttpFailure(lsHttpResponses,TRUE)
}

if (sys.nframe() == 0) {
    UploadGitHubAssets(
        sAssetConfigFilePath = "./auxil/ghAssetConfig.yaml",
        sGitHubAssetRouteId = "local_Chris_IMPACTncd_Engl_0.0.4"
    )
} # execute with defaults if run from topmost frame (under Rscript)

