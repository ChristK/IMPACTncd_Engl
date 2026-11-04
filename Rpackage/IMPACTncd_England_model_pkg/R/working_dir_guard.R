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

# .same_dir ----
# TRUE if two paths name the same directory. Normalises separators (getwd()
# uses "/", normalizePath() defaults to "\\" on Windows), resolves symlinks
# for existing paths, drops a trailing "/" (but keeps "/" and "C:/"), and
# compares case-insensitively on Windows, where drive-letter case can differ.
# @param a,b Character scalars.
# @return Logical scalar.
.same_dir <- function(a, b) {
  canon <- function(p) {
    p <- normalizePath(p, winslash = "/", mustWork = FALSE)
    p <- sub("([^:/])/+$", "\\1", p) # strip a trailing "/", keep "/" and "C:/"
    if (.Platform$OS.type == "windows") tolower(p) else p
  }
  identical(canon(a), canon(b))
}

# .assert_wd ----
# Stop unless the working directory is still `root`.
# @param root Character scalar: the root captured by Design$initialize(), or
#   NULL for an object that predates root tracking (then nothing is checked).
# @param what Character scalar naming the operation, for the message.
# @return TRUE, invisibly, if the check passes.
.assert_wd <- function(root, what = "This operation") {
  if (is.null(root) || !length(root) || is.na(root)) return(invisible(TRUE))
  wd <- getwd()
  if (!.same_dir(wd, root)) {
    stop(
      sprintf(
        paste0(
          "%s: working directory changed since this model was set up ",
          "(was '%s', now '%s'). IMPACTncd resolves ./inputs and ./simulation ",
          "against the working directory, so continuing could read or delete ",
          "the wrong files. Run setwd(\"%s\") and retry. (In RStudio this ",
          "usually means code was run from inside an .Rmd file, whose chunks ",
          "are evaluated in the document's folder.)"
        ),
        what, root, wd, root
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}
