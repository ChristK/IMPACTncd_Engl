setOptions_for_repo <- function() {
  chooseCRANmirror(ind = 1)
  repos <- getOption("repos")
}
setOptions_for_repo()

# temporary fix: R.utils::reassignInPackage() used to inject our revised functions (see below) into the [piggyback] package, prior to their expected adoption by the piggyback team
if (!require(R.utils)) {
  dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)
  install.packages("R.utils", lib = Sys.getenv("R_LIBS_USER"), repos = "https://cran.rstudio.com/")
  library(R.utils)
}
if (!require(yaml)) {
  dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)
  install.packages("yaml", lib = Sys.getenv("R_LIBS_USER"), repos = "https://cran.rstudio.com/")
  library(yaml)
}

#' @description Remove slashes from start and end of path.
#' @param sPath: string initial path.
#' @param bRidStartSlash: bool remove slash at start of path
#' @return string possibly modified path.
TrimSlashes <- function(sPath, bRidStartSlash = TRUE) {
  if (bRidStartSlash && substr(sPath, 1, 1) == "/") {
    sPath <- substr(sPath, 2, nchar(sPath))
  }
  if (substr(sPath, nchar(sPath), nchar(sPath)) == "/") {
    sPath <- substr(sPath, 1, nchar(sPath) - 1)
  }
  return(sPath)
}

#' @description Helper to create descriptive list of HTTP headers from the given response object.
#' @param httpResponse HTTP response object as provided by [httr] package.
#' @return String giving each header's name and value as: <name1>=<value1>; <name2>=<value2>; ...
HttpHeaderNamesAndValues <- function(httpResponse) {
  sHeaders <- ""
  iIndex <- 1
  while (iIndex < length(httpResponse$headers)) {
    sHeaders <- paste0(sHeaders, attr(httpResponse$headers[iIndex], "name"), "=", httpResponse$headers[[iIndex]], "; ")
    iIndex <- iIndex + 1
  }
  return(sHeaders)
}

#' @description Stop if get HTTP response indicating missing file or unexcepted file size.
#' @param lsHttpResponses List of HTTP response objects, each provided by [httr] package.
#' @param bUploadedFiles boolean, files have been uploaded.
StopOnHttpFailure <- function(lsHttpResponses, bUploadedFiles) {
  for (httpResponse in lsHttpResponses)
  {
    if (is.null(httpResponse)) next # is NULL when logic decides not to proceed with upload

    # test: failure error code
    if (httpResponse$status_code != 200 & httpResponse$status_code != 201) {
      stop(paste("Error: unexpected HTTP", httpResponse$status_code, " status for asset", str(httpResponse$content), "with headers", HttpHeaderNamesAndValues(httpResponse)))
    }

    if (!bUploadedFiles) {
      # test: unexpected file size
      sFileNamePath <- httpResponse$content[1]
      iFileSizeOnDisk <- as.numeric(file.info(sFileNamePath)$size)
      iFileSizeFromGitHub <- as.numeric(httpResponse$headers$`content-length`)
      if (iFileSizeOnDisk != iFileSizeFromGitHub) {
        stop(paste(
          "Error: asset", sFileNamePath,
          "has unexpected file size (", iFileSizeOnDisk, " (disk),",
          iFileSizeFromGitHub, " (Github)) with headers", HttpHeaderNamesAndValues(httpResponse)
        ))
      }
    }
  }
}

#' @description Get GitHub asset route data from YAML AssetConfigFile or command-line variables.
#' Execution options:
#' 	Rscript <source-file> [<AssetFilePathName> [<GitHubAssetRouteId> [<GitHubToken>]]]
#' 		where <source-file> is the R script name,
#' 			<AssetFilePathName> is an optional asset config file's path and name;
#' 				if omitted, tries IMPACTncd_Engl/auxil/ghAssetConfig.yaml.
#' 			<GitHubAssetRouteId> is an asset route's [id] within <AssetFilePathName>;
#' 				if omitted, seeks ID from Sys.info()[['user']] name.
#' 			<GitHubToken> is a GitHub personal access token (PAT);
#' 				if omitted, seeks token from gh::gh_token().
#' or, source('gh_deploy.R') # loads <AssetFilePathName> via IMPACTncd_Engl/auxil/ghAssetConfig.yaml.
#' @param sRepo string (out param): GitHub repository name.
#' @param sTag string (out param): GitHub repository tag.
#' @param iTestWithFirstNAssets int (out param): only download first [iTestWithFirstNAssets] assets (for testing).
#' @param sUploadSrcDirPath string (out param): source directory path to scan for uploading assets to GitHub.
#' @param sDeployToRootDirPath string (out param): deployment directory path for downloading assets from GitHub.
#' @param bOverwriteFilesOnDeploy bool (out param): overwrite files during deployment.
#' @param sToken string (in|out param): GitHub personal access token (PAT).
GetGitHubAssetRouteInfo <- function(sRepo, sTag, sUploadSrcDirPath, sDeployToRootDirPath,
                                    bOverwriteFilesOnDeploy, sToken = NULL,
                                    iTestWithFirstNAssets) {
  lsCommandArgs <- commandArgs(TRUE)
  iNumCmdArgs <- length(lsCommandArgs)
  # get GitHub access token
  if (iNumCmdArgs > 0) {
    sGitHubToken <- lsCommandArgs[1]
    iNumCmdArgs <- iNumCmdArgs - 1
  } else if (is.null(sToken)) {
    sGitHubToken <- as.character(gh::gh_token())
  } # use token set in environmental variable (if given)
  else {
    sGitHubToken <- sToken
  }
    if (iNumCmdArgs > 0) { # quit - as didn't expect additional variables
    stop("Execute from the console:
					Rscript <scriptR> [<GitHubAssetRouteId> [<GitHubToken>]]
				or alternatively, may execute directly within R or RStudio:
					source(\"<scriptR>\")
					UploadGitHubAssets(sToken=<GitHubToken>)
				where in the above, <scriptR> is the R script name, either gh_deploy.R or gh_upload.R,
					<GitHubAssetRouteId> is an asset route's [id];
						if omitted, seeks ID matching Sys.info()[['user']] name.
					<GitHubToken> is a GitHub personal access token (PAT);
						if omitted, seeks token from gh::gh_token().")
  }
  # get desired asset route's properties
  gitHubAssetRouteFinal <- list(
    repo = sRepo,
    tag = sTag,
    personalAccessToken = sGitHubToken,
    uploadSrcDirectory = sUploadSrcDirPath,
    deployToRootDirectory = sDeployToRootDirPath,
    overwriteFilesOnDeploy = bOverwriteFilesOnDeploy
  )
  if (sGitHubToken == "") { # no access token given previously
    sGitHubToken <- gitHubAssetRouteFinal$personalAccessToken
    if (is.null(sGitHubToken)) warning("Failed reading GitHub personal access token (PAT) from GITHUB_PAT environmental variable.
					May set PAT on command-line or in asset config file [personalAccessToken] variable.")
  }
  eval.parent(substitute(sRepo <- gitHubAssetRouteFinal$repo))
  eval.parent(substitute(sTag <- gitHubAssetRouteFinal$tag))
  eval.parent(substitute(sToken <- sGitHubToken))
  eval.parent(substitute(sUploadSrcDirPath <- gitHubAssetRouteFinal$uploadSrcDirectory))
  eval.parent(substitute(sDeployToRootDirPath <- gitHubAssetRouteFinal$deployToRootDirectory))
  eval.parent(substitute(bOverwriteFilesOnDeploy <- if (gitHubAssetRouteFinal$overwriteFilesOnDeploy == "1") TRUE else FALSE))
  eval.parent(substitute(iTestWithFirstNAssets <- if (is.null(gitHubAssetRouteFinal$testWithFirstNAssets)) {
    0
  } else {
    as.integer(gitHubAssetRouteFinal$testWithFirstNAssets)
  }))
  return(0)
  stop(paste0("Failed finding GitHub asset route data"))
}

####################################################################################
# Temporary fix: below revised [piggyback] package functions injected into this package,
#  prior to the expected adoption of these functions by the piggyback team.
#  NOTE: these functions are as-submitted for [piggyback] distribution Github pull-request.
# Functions: pb_upload_liverpool(), pb_upload_file_liverpool(), pb_download_liverpool().
####################################################################################

#' Upload data to an existing release
#'
#' NOTE: you must first create a release if one does not already exists.
#' @param file path to file to be uploaded
#' @param repo Repository name in format "owner/repo". Defaults to `guess_repo()`.
#' @param tag  tag for the GitHub release to which this data should be attached.
#' @param name name for uploaded file. If not provided will use the basename of
#' `file` (i.e. filename without directory)
#' @param overwrite overwrite any existing file with the same name already
#'  attached to the on release? Default behavior is based on timestamps,
#'  only overwriting those files which are older.
#' @param use_timestamps DEPRECATED.
#' @param show_progress logical, show a progress bar be shown for uploading?
#' Defaults to `[interactive()]` - can also set globally with options("piggyback.verbose")
#' @param .token GitHub authentication token, see `[gh::gh_token()]`
#' @param dir directory relative to which file names should be based, defaults to NULL for current working directory.
#' @examples
#' \dontrun{
#' # Needs your real token to run
#'
#' readr::write_tsv(mtcars, "mtcars.tsv.xz")
#' pb_upload("mtcars.tsv.xz", "cboettig/piggyback-tests")
#' }
#' @export
#'
pb_upload_liverpool <- function(file,
                                repo = guess_repo(),
                                tag = "latest",
                                name = NULL,
                                overwrite = "use_timestamps",
                                use_timestamps = NULL,
                                show_progress = getOption("piggyback.verbose", default = interactive()),
                                .token = gh::gh_token(),
                                dir = NULL) {
  stopifnot(
    is.character(repo),
    is.character(tag),
    length(tag) == 1,
    length(repo) == 1
  )

  releases <- pb_releases(repo, .token)

  if (tag == "latest" && length(releases$tag_name) > 0 && !"latest" %in% releases$tag_name) {
    if (getOption("piggyback.verbose", default = interactive())) {
      cli::cli_alert_info("Uploading to latest release: {.val {releases$tag_name[[1]]}}.")
    }
    tag <- releases$tag_name[[1]]
  }

  if (!tag %in% releases$tag_name && !interactive()) {
    cli::cli_abort("Release {.val {tag}} not found in {.val {repo}}. No upload performed.")
  }

  if (!tag %in% releases$tag_name) {
    cli::cli_alert_warning("Release {.val {tag}} not found in {.val {repo}}.")

    run <- utils::menu(
      choices = c("Yes", "No"),
      title = glue::glue("Would you like to create a new release now?")
    )

    if (run == 2) {
      return(invisible(NULL))
    }
    if (run == 1) pb_new_release(repo = repo, tag = tag, .token = .token)
  }

  ## start fresh
  memoise::forget(pb_info)

  out <- lapply(seq_along(file), function(iFileIndex) {
    pb_upload_file(
      file[iFileIndex],
      repo,
      tag,
      if (is.null(name)) NULL else name[iFileIndex],
      overwrite,
      use_timestamps,
      show_progress,
      .token,
      dir
    )
  })

  ## break cache when done
  memoise::forget(pb_info)
  invisible(out)
}

pb_upload_file_liverpool <- function(file,
                                     repo = guess_repo(),
                                     tag = "latest",
                                     name = NULL,
                                     overwrite = "use_timestamps",
                                     use_timestamps = NULL,
                                     show_progress = getOption("piggyback.verbose", default = interactive()),
                                     .token = gh::gh_token(),
                                     dir = NULL) {
  file_path <- do.call(file.path, compact(list(dir, file)))

  ## Uses NULL as default dir, drops it with compact, then
  ## does the file.path call with what's left
  ##
  ## This is better than using "." as default dir because if you try to pass an
  ## absolute path with "." e.g. file.path(".","C:/Users/Tan") it will
  ## return "./C:/Users/Tan" which is not desired.

  if (!file.exists(file_path)) {
    cli::cli_warn("File {.file {file_path}} does not exist.")
    return(NULL)
  }

  if (!is.null(use_timestamps)) {
    cli::cli_warn("{.code use_timestamps} argument is deprecated, please set {.code overwrite = 'use_timestamps'} instead")
  }

  ## Yeah, having two separate arguments was clearly a mistake!
  ## Code has been partially refactored now so that user just
  ## sets `overwrite` and we handle the twisted logic internally here:
  use_timestamps <- switch(as.character(overwrite),
    "TRUE" = FALSE,
    "FALSE" = FALSE,
    "use_timestamps" = TRUE
  )
  overwrite <- switch(as.character(overwrite),
    "TRUE" = TRUE,
    "FALSE" = FALSE,
    "use_timestamps" = TRUE
  )

  progress <- httr::progress("up")
  if (!show_progress) progress <- NULL

  if (is.null(name)) {
    ## name is name on GitHub, technically need not be name of local file
    name <- basename(file_path)
  }

  ## memoised for piggyback_cache_duration
  df <- pb_info(repo, tag, .token)

  i <- which(df$file_name == name)

  if (length(i) > 0) { # File of same name is on GitHub

    if (use_timestamps) {
      local_timestamp <- fs::file_info(file_path)$modification_time

      no_update <- local_timestamp <= df[i, "timestamp"]
      if (no_update) {
        cli::cli_warn("Matching or more recent version of {.file {file_path}} found on GH, not uploading.")
        return(invisible(NULL))
      }
    }

    if (overwrite) {
      ## If we find matching id, Delete file from release.
      gh::gh("DELETE /repos/:owner/:repo/releases/assets/:id",
        owner = df$owner[[1]],
        repo = df$repo[[1]],
        id = df$id[i],
        .token = .token
      )
    } else {
      cli::cli_warn("Skipping upload of {.file {df$file_name[i]}} as file exists on GitHub and {.code overwrite = FALSE}")
      return(invisible(NULL))
    }
  }

  if (show_progress) cli::cli_alert_info("Uploading {.file {name}} ...")

  releases <- pb_releases(repo = repo, .token = .token)
  upload_url <- releases$upload_url[releases$tag_name == tag]

  r <- httr::RETRY(
    verb = "POST",
    url = sub("\\{.+$", "", upload_url), # rid anything from { to end.
    query = list(name = name),
    httr::add_headers(Authorization = paste("token", .token)),
    body = httr::upload_file(file_path),
    progress,
    terminate_on = c(400, 401, 403, 404, 422)
  )

  if (show_progress) httr::warn_for_status(r)

  ## Release info changed, so break cache
  try({
    memoise::forget(pb_info)
  })
  invisible(r)
}

#' Download data from an existing release
#'
#' @param file name or vector of names of files to be downloaded. If `NULL`,
#' all assets attached to the release will be downloaded.
#' @param dest name of vector of names of where file should be downloaded.
#' Can be a directory or a list of filenames the same length as `file`
#' vector. Any directories in the path provided must already exist. WARNING! confusing behaviour: *directory* required for single asset; destination *filenames* required for multiple assets.
#' @param overwrite Should any local files of the same name be overwritten?
#'  default `TRUE`.
#' @param ignore a list of files to ignore (if downloading "all" because
#'  `file=NULL`).
#' @inheritParams pb_upload
#'
#' @export
#' @examples \dontrun{
#' ## Download a specific file.
#' ## (dest can be omitted when run inside and R project)
#' piggyback::pb_download("iris.tsv.gz",
#'   repo = "cboettig/piggyback-tests",
#'   dest = tempdir()
#' )
#' }
#' \dontrun{
#' ## Download all files
#' piggyback::pb_download(
#'   repo = "cboettig/piggyback-tests",
#'   dest = tempdir()
#' )
#' }
pb_download_liverpool <- function(file = NULL,
                                  dest = ".",
                                  repo = guess_repo(),
                                  tag = "latest",
                                  overwrite = TRUE,
                                  ignore = "manifest.json",
                                  use_timestamps = TRUE,
                                  show_progress = getOption("piggyback.verbose", default = interactive()),
                                  .token = gh::gh_token()) {
  progress <- httr::progress("down")

  if (!show_progress) progress <- NULL

  df <- pb_info(repo, tag, .token)

  ## drop failed upload states from list
  df <- df[df$state != "starter", ]

  if (!is.null(file)) {
    i <- which(df$file_name %in% file)
    if (length(i) < 1) {
      cli::cli_warn("file(s) {.file {file}} not found in repo {.val {repo}}")
    }

    df <- df[i, ]
  } else {
    i <- which(df$file_name %in% ignore)
    if (length(i) >= 1) {
      df <- df[-i, ]
    }
    file <- df$file_name
  }


  ## if dest paths are not provided, we will write all files to dest dir
  # User is responsible for making sure dest dir exists!
  if (length(dest) == 1) { # mbirkett: this unexpectedly causes different behaviour for single/multiple assets (i.e. length of [file] list). Destination _filenames_ expected for multiple assets, while single assets require destination _directory_. For latter, perhaps should use destination filename if given, else assuming destination directory. Current behaviour can cause errors and some confusion.
    i <- which(df$file_name %in% file)
    dest <- file.path(dest, df$file_name[i])
  }
  # dest should now be of length df
  if (nrow(df) != length(file)) {
    iMissingFileIndices <- which(!(file %in% df$file_name))
    stop(
      "Error: following files requested for download, yet not available on Github:\n",
      paste(file[iMissingFileIndices], collapse = ";")
    )
  }
  # GitHub and local files listed in different sequence; need to carefully match destination namepath by given filename
  iGitHubFileRow <- 1
  while (iGitHubFileRow <= nrow(df)) {
    iRequestedFileIndex <- which(file %in% df[iGitHubFileRow, c("file_name")])
    df[iGitHubFileRow, c("dest")] <- dest[iRequestedFileIndex]
    iGitHubFileRow <- iGitHubFileRow + 1
  }


  if (use_timestamps) {
    local_timestamp <- fs::file_info(dest)$modification_time
    update <- df$timestamp > local_timestamp
    update[is.na(update)] <- TRUE # we'll download if missing locally
    df <- df[update, ]

    if (dim(df)[[1]] < 1) {
      cli::cli_alert_info("All local files already up-to-date!")
      return(invisible(NULL))
    }
  }

  resp <- lapply(seq_along(df$id), function(i) {
    gh_download_asset(df$owner[[1]],
      df$repo[[1]],
      id = df$id[i],
      destfile = df$dest[i],
      overwrite = overwrite,
      .token = .token,
      progress = progress
    )
  })
  return(invisible(resp))
}

####################################################################################

reassignInPackage("pb_upload", pkgName = "piggyback", pb_upload_liverpool)
reassignInPackage("pb_upload_file", pkgName = "piggyback", pb_upload_file_liverpool)
reassignInPackage("pb_download", pkgName = "piggyback", pb_download_liverpool)
