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
# Working-directory guard (GitHub issue #11)
#
# The model resolves ~90 bare "./inputs/..." and "./simulation/..." paths
# against getwd() at CALL time, while Disease and ExposureEffect fix their
# cache locations when they are CONSTRUCTED. A working directory that changes
# in between -- typically RStudio evaluating a vignette's chunks in the
# vignette's own folder -- used to make gen_parf_files() look for the PARF
# registry in the wrong folder, judge the cache invalid and delete the real
# PARF cache, before failing on a relative read. From a sibling checkout with
# the same layout the relative reads would instead SUCCEED and a run would
# silently mix two trees.
#
# Two defences, both keyed on ONE root that Design$initialize() captures as
# getwd():
#   * anchor: Disease and ExposureEffect build their registry paths (and the
#     base to_rel_path() strips) from that root instead of calling getwd();
#   * guard:  methods that touch ./inputs or ./simulation call .assert_wd()
#     first and stop if the working directory has changed.
# The root is deliberately NOT the YAML's folder (design_dir): that is
# scenarios/<X>/ for scenario drivers and a symlinked staging mount in Docker.
# -----------------------------------------------------------------------------

# .canon_dir ----
# Canonical form of a directory path for comparison: "/" separators (getwd()
# uses "/", normalizePath() defaults to "\\" on Windows), symlinks resolved for
# existing paths, a trailing "/" dropped (but "/" and "C:/" kept), and lower
# case on Windows, where drive-letter case can differ.
# @param p Character vector.
# @return Character vector.
.canon_dir <- function(p) {
  p <- normalizePath(p, winslash = "/", mustWork = FALSE)
  p <- sub("([^:/])/+$", "\\1", p) # strip trailing "/", keep "/" and "C:/"
  if (.Platform$OS.type == "windows") tolower(p) else p
}

# .same_dir ----
# TRUE if two paths name the same directory (see .canon_dir()).
# @param a,b Character scalars.
# @return Logical scalar.
.same_dir <- function(a, b) {
  identical(.canon_dir(a), .canon_dir(b))
}

# .is_abs_path ----
# TRUE for an absolute path: "/...", "\\...", "~...", or a Windows drive
# ("C:" / "C:/" / "C:\\"), which also covers UNC paths ("\\\\server\\share").
# @param p Character vector.
# @return Logical vector.
.is_abs_path <- function(p) {
  grepl("^([/\\\\~]|[A-Za-z]:)", p)
}

# .rel_to_root ----
# `path` relative to `root` where it lies under it, else unchanged; the form
# the git-tracked PARF/RR generation records store, so they are portable
# across users and machines. On Windows both sides are normalised to "/" and
# compared case-insensitively first: normalizePath() of an RR file gives "\\"
# separators, which never matched the "/" of getwd(), so those records used to
# keep absolute machine paths. On other platforms this is exactly the old
# to_rel_path() (both sides are already physical paths with "/").
# @param path Character vector (NULL/empty is returned as is).
# @param root Character scalar.
# @return Character vector.
.rel_to_root <- function(path, root) {
  if (!length(path)) return(path)
  if (.Platform$OS.type == "windows") {
    p <- normalizePath(path, winslash = "/", mustWork = FALSE)
    pre <- paste0(sub("/+$", "", normalizePath(root, winslash = "/", mustWork = FALSE)), "/")
    hit <- startsWith(tolower(p), tolower(pre))
  } else {
    p <- path
    pre <- paste0(root, "/")
    hit <- startsWith(p, pre)
  }
  ifelse(hit, substring(p, nchar(pre) + 1L), path)
}

# .assert_wd ----
# Stop unless the working directory is still `root`.
# @param root Character scalar: the root captured by Design$initialize(), or
#   NULL for an object that predates root tracking (then nothing is checked).
# @param what Character scalar naming the operation, for the message.
# @param cause "entry" (a method was called from another folder) or
#   "scenario" (a scenario function changed the folder during run()): the two
#   have different causes and different remedies.
# @return TRUE, invisibly, if the check passes.
.assert_wd <- function(root, what = "This operation",
                       cause = c("entry", "scenario")) {
  cause <- match.arg(cause)
  if (is.null(root) || !length(root) || is.na(root)) return(invisible(TRUE))
  wd <- getwd()
  if (.same_dir(wd, root)) return(invisible(TRUE))
  msg <- if (cause == "entry") {
    sprintf(
      paste0(
        "%s: working directory changed since this model was set up ",
        "(was '%s', now '%s'). IMPACTncd resolves ./inputs and ./simulation ",
        "against the working directory, so continuing could read or delete ",
        "the wrong files. Run setwd(\"%s\") at the console and try again. ",
        "In RStudio this usually means the code was run inside an .Rmd ",
        "chunk: chunks are evaluated in the document's folder and a setwd() ",
        "inside one is undone when it ends, so run the code from the console ",
        "or a script instead (or set Tools > Global Options > R Markdown > ",
        "'Evaluate chunks in directory' to 'Project')."
      ),
      what, root, wd, root
    )
  } else {
    sprintf(
      paste0(
        "%s changed the working directory (it was '%s' and is now '%s'). ",
        "IMPACTncd resolves ./inputs and ./simulation against the working ",
        "directory, so the rest of the run would read the wrong files. A ",
        "scenario function must leave the working directory as it found it: ",
        "if yours calls setwd(), restore it before returning, e.g. ",
        "old <- setwd(dir); on.exit(setwd(old), add = TRUE)."
      ),
      what, root, wd
    )
  }
  stop(msg, call. = FALSE)
}
