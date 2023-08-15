# SUMMARY: remove redundant asset files from GitHub: those which are not assets below local source directory specified in [uploadSrcDirectory] in auxil/ghAssetConfig.yaml.

library("piggyback",lib.loc="/home/pp0u8134/R/x86_64-pc-linux-gnu-library/4.2")
library("data.table")
source("ghAssetUtils.R")

# get GitHub asset route data
GetGitHubAssetRouteInfo(sId,sRepo,sTag,sToken,sUploadSrcDirPath,
	sDeployToRootDirPath,bOverwriteFilesOnDeploy,iTestWithFirstNAssets)
sDeployToRootDirPath<- TrimSlashes(sDeployToRootDirPath,bRidStartSlash=FALSE)

# find local assets, excluding secure_data, "*tmp.qs", and "*parf/PARF_*.qs" files.
lsLocalAssetPathNames<- list.files(sUploadSrcDirPath, pattern = ".fst$|.xls$|.xlsx$|.qs$", full.names = TRUE, recursive = TRUE)
lsLocalAssetPathNames<- grep("secure_data", lsLocalAssetPathNames, value = TRUE, invert = TRUE)
# NOTE currently no necessary files ar .qs. For future proof I add them above
# and exclude the below
lsLocalAssetPathNames<- grep("tmp.qs$", lsLocalAssetPathNames, value = TRUE, invert = TRUE)
lsLocalAssetPathNames<- grep("/simulation/parf/PARF_.*\\.qs$", lsLocalAssetPathNames, value = TRUE, invert = TRUE)
lsLocalAssetPathNames<- grep("simulation/DELETEME", lsLocalAssetPathNames, value = TRUE, invert = TRUE)

# get *unique* sanitised file names and directories
dtOriginalAndSanitisedFilePathNames<- data.table(orig_file = basename(lsLocalAssetPathNames),
	sanit_file = gsub("[^[:alnum:]&&^\\.]", ".", basename(lsLocalAssetPathNames)), # replace all not (alphanumeric or '.') with .
	abs_dir = dirname(lsLocalAssetPathNames),
	rel_dir = gsub(sUploadSrcDirPath, "", dirname(lsLocalAssetPathNames)), key = "orig_file") # rel_dir to keep initial slash
if (any(duplicated(dtOriginalAndSanitisedFilePathNames$sanit_file))) stop("Duplicated filenames found")

dfAssetFilesMetaData<- pb_info(sRepo, sTag, sToken)
message("initial length dfAssetFilesMetaData : ",nrow(dfAssetFilesMetaData))
dfAssetFilesMetaData<- dfAssetFilesMetaData[-which(dfAssetFilesMetaData$file_name %in% dtOriginalAndSanitisedFilePathNames$sanit_file),]
message("after removing matches, length dfAssetFilesMetaData : ",nrow(dfAssetFilesMetaData))

# delete unwanted GitHub asset files
pb_delete(dfAssetFilesMetaData$file_name,sRepo,sTag,sToken)



