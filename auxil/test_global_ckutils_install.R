#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# test_global_ckutils_install.R
#
# Tests how global.R installs CKutils: over HTTPS from GitHub, falling back to
# a git clone over SSH for networks that block HTTPS to github.com but allow
# SSH, and finally to an already-installed copy (see install_ckutils() in
# global.R).
#
# global.R is a bootstrap script, not package code, so this cannot live in the
# package's tinytest suite; and it needs the real network and an SSH key
# registered with GitHub, so it cannot run in CI either. Run it by hand after
# changing that block of global.R:
#
#   Rscript auxil/test_global_ckutils_install.R
#
# Requirements: HTTPS access to github.com and CRAN, git on the PATH, `ssh -T
# git@github.com` authenticating, and CKutils already installed (its
# dependencies are borrowed from your library). One scenario compiles CKutils
# from source (~1-2 min). Exit status 0 only if every check passes.
#
# How the conditions are produced:
#   * "HTTPS to GitHub blocked": https_proxy points at a closed local port, so
#     every HTTPS request fails at once with "Connection refused". no_proxy
#     exempts the CRAN mirrors, as a corporate firewall that blocks only
#     GitHub would. SSH does not use these variables, so it still works.
#   * "SSH blocked": GIT_SSH_COMMAND points ssh at the same closed port.
#   * Each scenario runs in a fresh R process (callr) whose ONLY writable
#     library is a throwaway "shadow" library. It holds symlinks to every
#     package in your user library except CKutils, plus a COPY of your
#     CKutils. Your user library is not on the child's library path, so
#     nothing can be installed into it; the final check confirms it is
#     byte-for-byte unchanged.
#   * "CKutils needs updating": the copy's recorded commit (RemoteSha) is
#     replaced with a fake one, so the SSH route must really clone, compile
#     and install. remotes takes the same path for a package that is not
#     installed at all: both a missing and a stale SHA differ from the remote.
#   * "CKutils not installed": R always keeps the system library (.Library)
#     on the path, so a CKutils installed there cannot be hidden from a child.
#     The child instead replaces system.file() with a test double reporting
#     CKutils as absent, which is the one check install_ckutils() makes.
#
# Caveat: unauthenticated GitHub API calls are rate-limited (60/hour). If the
# HTTPS scenarios fail with HTTP 403, wait or set GITHUB_PAT.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(callr))

global_r <- normalizePath("global.R", mustWork = FALSE)
if (!file.exists(global_r)) {
  stop("Run this from the project root (global.R not found).", call. = FALSE)
}
ssh_url <- "git@github.com:ChristK/CKutils.git"
real_ckutils <- find.package("CKutils") # errors if CKutils is not installed
real_lib <- dirname(real_ckutils)
if (real_lib %in% c(.Library, .Library.site)) {
  stop(
    "CKutils was found only in a system library; install it in your user ",
    "library first, which these tests borrow from.",
    call. = FALSE
  )
}

# PASS/FAIL reporter, as in auxil/verify_init_prvl_crn.R.
ok <- TRUE
say <- function(pass, msg) {
  pass <- isTRUE(pass)
  ok <<- ok && pass
  cat(sprintf("[%s] %s\n", if (pass) " OK " else "FAIL", msg))
}

# --- Fixtures --------------------------------------------------------------

work_dir <- tempfile("ckutils_install_test_")
dir.create(work_dir)

# Children inherit this process's environment, so a route set in the tester's
# own environment would leak into every scenario. Unset it here; scenarios
# that need it set it explicitly.
Sys.unsetenv("IMPACTNCD_CKUTILS_INSTALL")

# A shadow library: every package of the real user library as a symlink,
# except CKutils, which is copied (with_ckutils = TRUE) or left out.
make_shadow_lib <- function(name, with_ckutils) {
  lib <- file.path(work_dir, name)
  dir.create(lib)
  pkgs <- setdiff(list.files(real_lib), "CKutils")
  pkgs <- pkgs[!startsWith(pkgs, "00")] # skip stale install locks
  stopifnot(all(file.symlink(file.path(real_lib, pkgs), file.path(lib, pkgs))))
  if (with_ckutils) {
    stopifnot(file.copy(real_ckutils, lib, recursive = TRUE, copy.date = TRUE))
  }
  lib
}

# Makes the copied CKutils in `lib` claim to be built from commit `sha`.
# packageDescription() -- and so remotes -- reads the installed metadata from
# Meta/package.rds rather than from the DESCRIPTION file, so edit both.
set_recorded_sha <- function(lib, sha) {
  desc_file <- file.path(lib, "CKutils", "DESCRIPTION")
  stopifnot(!nzchar(Sys.readlink(desc_file))) # a copy, never the real file
  d <- read.dcf(desc_file)
  d[1, "RemoteSha"] <- sha
  write.dcf(d, desc_file)
  meta_file <- file.path(lib, "CKutils", "Meta", "package.rds")
  meta <- readRDS(meta_file)
  meta$DESCRIPTION[["RemoteSha"]] <- sha
  saveRDS(meta, meta_file)
}

# Fingerprint of the real user library, to prove the tests never wrote to it.
lib_fingerprint <- function() {
  descs <- file.path(real_lib, list.files(real_lib), "DESCRIPTION")
  descs <- descs[file.exists(descs)]
  tools::md5sum(descs)
}
real_lib_before <- lib_fingerprint()

remote_head <- system2("git", c("ls-remote", ssh_url, "HEAD"), stdout = TRUE)
remote_head <- sub("\t.*", "", remote_head[1])
say(
  grepl("^[0-9a-f]{40}$", remote_head),
  paste("git ls-remote over SSH returns the HEAD commit:", remote_head)
)

closed_port <- "http://127.0.0.1:9"
cran_hosts <- "cloud.r-project.org,packagemanager.posit.co"
block_https <- c(
  https_proxy = closed_port,
  HTTPS_PROXY = closed_port,
  no_proxy = cran_hosts,
  NO_PROXY = cran_hosts
)
block_ssh <- c(
  GIT_SSH_COMMAND = "ssh -o BatchMode=yes -o HostName=127.0.0.1 -o Port=9"
)

# --- Child-process runners -------------------------------------------------

# Runs install_ckutils(), as defined in global.R, in a fresh R process. Only
# that definition is evaluated, so the route is chosen exactly as global.R
# would choose it, from IMPACTNCD_CKUTILS_INSTALL in the child's environment.
run_install <- function(
  lib,
  env = character(),
  quiet = TRUE,
  pretend_absent = FALSE
) {
  child <- function(global_r, quiet, pretend_absent) {
    if (pretend_absent) {
      # Test double: install_ckutils() is defined in this environment below,
      # so it finds this system.file() before base's.
      system.file <- function(..., package = "base") {
        if (identical(package, "CKutils")) {
          return("")
        }
        base::system.file(..., package = package)
      }
    }
    exprs <- parse(global_r, keep.source = FALSE)
    is_def <- vapply(
      exprs,
      function(e) {
        is.call(e) &&
          identical(e[[1]], as.name("<-")) &&
          identical(e[[2]], as.name("install_ckutils"))
      },
      logical(1)
    )
    stopifnot(sum(is_def) == 1L)
    eval(exprs[[which(is_def)]])

    msgs <- character()
    warns <- character()
    res <- withCallingHandlers(
      tryCatch(
        list(route = install_ckutils(quiet = quiet), error = NA_character_),
        error = function(e) {
          list(route = NA_character_, error = conditionMessage(e))
        }
      ),
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    c(
      res,
      list(
        messages = paste(msgs, collapse = "\n"),
        warnings = paste(warns, collapse = "\n"),
        libpaths = .libPaths()
      )
    )
  }
  callr::r(
    child,
    args = list(
      global_r = global_r,
      quiet = quiet,
      pretend_absent = pretend_absent
    ),
    libpath = c(lib, .Library),
    env = c(callr::rcmd_safe_env(), R_LIBS_USER = lib, env)
  )
}

installed_desc <- function(lib) file.path(lib, "CKutils", "DESCRIPTION")
desc_mtime <- function(lib) file.mtime(installed_desc(lib))
has <- function(x, pattern) grepl(pattern, x, fixed = TRUE)

# --- Scenarios with CKutils already installed --------------------------------

cat("\n== CKutils installed ==\n")
lib_inst <- make_shadow_lib("installed", with_ckutils = TRUE)
mtime0 <- desc_mtime(lib_inst)

r <- run_install(lib_inst)
say(
  !real_lib %in% r$libpaths && r$libpaths[1] == lib_inst,
  "the child's library path is the shadow library, not the real one"
)
say(
  identical(r$route, "https") && !nzchar(r$messages) && !nzchar(r$warnings),
  sprintf("auto, open network: HTTPS route, silently (route = %s)", r$route)
)

r <- run_install(lib_inst, env = block_https)
say(
  identical(r$route, "ssh"),
  sprintf("auto, HTTPS blocked: falls back to SSH (route = %s)", r$route)
)
say(
  has(r$messages, "over HTTPS failed; trying SSH") &&
    has(r$messages, "IMPACTNCD_CKUTILS_INSTALL=ssh"),
  "auto, HTTPS blocked: says it is falling back, and how to skip HTTPS"
)
say(!nzchar(r$warnings), "auto, HTTPS blocked: no warning")

r <- run_install(
  lib_inst,
  env = c(block_https, IMPACTNCD_CKUTILS_INSTALL = " SSH ")
)
say(
  identical(r$route, "ssh") && !nzchar(r$messages),
  "mode ' SSH ' (case and blanks ignored): SSH only, HTTPS never tried"
)

r <- run_install(
  lib_inst,
  env = c(block_https, IMPACTNCD_CKUTILS_INSTALL = "https")
)
say(
  identical(r$route, "installed") &&
    has(r$warnings, "[https]") &&
    !has(r$warnings, "[ssh]") &&
    has(r$warnings, "using the installed version"),
  "mode 'https', HTTPS blocked: SSH not tried; warns and keeps the installed copy"
)

r <- run_install(
  lib_inst,
  env = c(block_https, block_ssh, IMPACTNCD_CKUTILS_INSTALL = "")
)
say(
  identical(r$route, "installed") &&
    has(r$warnings, "[https]") &&
    has(r$warnings, "[ssh]"),
  "variable set but empty = auto; both blocked: tries both, keeps the installed copy"
)
say(
  has(r$warnings, "Connection refused") &&
    has(r$warnings, "Could not read from remote repository"),
  "both blocked: the SSH failure carries git's own reason, not just its exit code"
)
say(
  identical(desc_mtime(lib_inst), mtime0),
  "no scenario above reinstalled CKutils (its commit was already current)"
)

# --- Scenarios with CKutils not installed -----------------------------------

cat("\n== CKutils not installed (system.file() test double) ==\n")
lib_empty <- make_shadow_lib("empty", with_ckutils = FALSE)

r <- run_install(lib_empty, env = c(IMPACTNCD_CKUTILS_INSTALL = "sftp"))
say(
  is.na(r$route) &&
    has(r$error, "IMPACTNCD_CKUTILS_INSTALL must be") &&
    has(r$error, "\"sftp\""),
  "an invalid mode is an error naming the variable and the bad value"
)

r <- run_install(lib_empty, env = c(block_https, block_ssh), pretend_absent = TRUE)
say(
  is.na(r$route) &&
    has(r$error, "CKutils is not installed") &&
    has(r$error, "[https]") &&
    has(r$error, "[ssh]") &&
    has(r$error, "ssh -T git@github.com"),
  "both blocked, nothing installed: an error listing both failures and the remedy"
)
say(
  !file.exists(installed_desc(lib_empty)),
  "both blocked: nothing was installed"
)

# --- A real install over SSH --------------------------------------------------

cat("\n== Install over SSH: recorded commit stale (compiles CKutils) ==\n")
lib_stale <- make_shadow_lib("stale", with_ckutils = TRUE)
set_recorded_sha(lib_stale, strrep("0", 40))
say(
  identical(packageDescription("CKutils", lib.loc = lib_stale)$RemoteSha, strrep("0", 40)),
  "fixture: R now sees the shadow CKutils as built from a stale commit"
)
mtime_stale <- desc_mtime(lib_stale)

r <- run_install(lib_stale, env = block_https)
say(
  identical(r$route, "ssh") && is.na(r$error),
  sprintf("auto, HTTPS blocked: installs from an SSH clone (route = %s)", r$route)
)
say(
  !identical(desc_mtime(lib_stale), mtime_stale),
  "CKutils was really reinstalled into the shadow library"
)
d <- read.dcf(installed_desc(lib_stale))
say(
  identical(unname(d[1, "RemoteType"]), "xgit") &&
    identical(unname(d[1, "RemoteUrl"]), ssh_url),
  "the install records where it came from (RemoteType xgit, RemoteUrl)"
)
say(
  identical(unname(d[1, "RemoteSha"]), remote_head),
  "the install records the commit it was built from (RemoteSha = remote HEAD)"
)
loaded_from <- callr::r(
  function() {
    library(CKutils)
    normalizePath(getNamespaceInfo(asNamespace("CKutils"), "path"))
  },
  libpath = c(lib_stale, .Library),
  env = c(callr::rcmd_safe_env(), R_LIBS_USER = lib_stale)
)
say(
  identical(loaded_from, normalizePath(file.path(lib_stale, "CKutils"))),
  "the SSH-installed CKutils attaches, from the shadow library"
)

mtime1 <- desc_mtime(lib_stale)
r <- run_install(
  lib_stale,
  env = c(block_https, IMPACTNCD_CKUTILS_INSTALL = "ssh"),
  quiet = FALSE
)
say(
  identical(r$route, "ssh") &&
    has(r$messages, "Skipping install of 'CKutils' over SSH") &&
    identical(desc_mtime(lib_stale), mtime1),
  "SSH again: skipped, because the recorded commit matches (no recompile)"
)

r <- run_install(lib_stale, quiet = FALSE)
say(
  identical(r$route, "https") &&
    has(r$messages, "Skipping install of 'CKutils' from a github remote") &&
    identical(desc_mtime(lib_stale), mtime1),
  "HTTPS after an SSH install: also skipped -- both routes read the same RemoteSha"
)

# --- global.R itself ----------------------------------------------------------

# Evaluates global.R from the top up to and including the library(CKutils)
# block that follows rm(install_ckutils), with HTTPS to GitHub blocked. The
# rest of global.R (dependency installs, building IMPACTncdEngland) is out of
# scope and would rebuild the model package as a side effect.
cat("\n== global.R, HTTPS to GitHub blocked ==\n")
child_global <- function(global_r) {
  exprs <- parse(global_r, keep.source = FALSE)
  at <- which(vapply(
    exprs,
    function(e) identical(e, quote(rm(install_ckutils))),
    logical(1)
  ))
  stopifnot(length(at) == 1L)
  msgs <- character()
  withCallingHandlers(
    for (e in exprs[seq_len(at + 1L)]) eval(e, globalenv()),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  list(
    attached = "package:CKutils" %in% search(),
    leftover = exists("install_ckutils", envir = globalenv(), inherits = FALSE),
    messages = paste(msgs, collapse = "\n")
  )
}
lib_global <- make_shadow_lib("global", with_ckutils = TRUE)
g <- callr::r(
  child_global,
  args = list(global_r = global_r),
  libpath = c(lib_global, .Library),
  env = c(callr::rcmd_safe_env(), R_LIBS_USER = lib_global, block_https),
  wd = dirname(global_r)
)
say(
  isTRUE(g$attached) && has(g$messages, "trying SSH"),
  "global.R falls back to SSH and attaches CKutils"
)
say(
  identical(g$leftover, FALSE),
  "global.R removes install_ckutils() from the global environment"
)

# --- Why the SSH route compares commits itself --------------------------------

# remotes finds an xgit remote's package name with `git archive --remote`,
# which GitHub refuses; the name is then NA, so remotes cannot find the
# installed commit and reinstalls on every run. install_ckutils() therefore
# compares commits itself. Informational only: if this ever reports a name,
# remotes has fixed it and the comparison in global.R could be dropped.
pkg_name <- tryCatch(
  asNamespace("remotes")$remote_package_name(
    asNamespace("remotes")$git_remote(ssh_url, git = "external", credentials = NULL)
  ),
  error = function(e) NA_character_
)
cat(sprintf(
  "\n[INFO] remotes %s reads the package name of an SSH GitHub remote as %s%s\n",
  format(packageVersion("remotes")),
  deparse(pkg_name),
  if (is.na(pkg_name)) " (so it cannot skip an up-to-date install)" else ""
))

# --- Safety -------------------------------------------------------------------

cat("\n")
say(
  identical(lib_fingerprint(), real_lib_before),
  "the real user library was not modified"
)

unlink(work_dir, recursive = TRUE) # removes the symlinks, not their targets
cat(sprintf(
  "\n%s\n",
  if (ok) "ALL CHECKS PASSED" else "*** SOME CHECKS FAILED ***"
))
quit(status = if (ok) 0L else 1L)
