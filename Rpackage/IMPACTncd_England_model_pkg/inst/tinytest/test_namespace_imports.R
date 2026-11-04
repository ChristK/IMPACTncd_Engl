# Tests: every function of an Imports package that the package calls
# unqualified is imported (originally DBI/duckdb; widened on review).
#
# Regression (found 2026-11-04 while testing the rewritten how_to_test_run.Rmd):
# `export_summaries(multicore = FALSE)` failed in a clean session with
#   could not find function "dbListFields"
# because five call sites in Simulation_class_summaries.R call dbListFields()
# unqualified while the NAMESPACE imported only dbConnect, dbDisconnect,
# dbExecute and dbGetQuery from DBI. An unimported name falls through the
# namespace and its imports to the SEARCH PATH, so it resolved only where
# DBI happened to be attached. The multicore path masked the bug:
# foreach(.packages = c(..., "duckdb", ...)) attaches duckdb, which Depends on
# DBI; the single-core path -- the one the test-run vignettes recommend --
# attaches nothing.
#
# The general check walks every function in the namespace, INCLUDING the
# methods of R6 class generators (public, private, active, and those added
# with $set()), and fails if any of them calls a function exported by DBI or
# duckdb that is neither imported nor defined in the package. It inspects
# code, so it does not depend on what the test session has attached.
#
# Uses codetools::findGlobals(); codetools is a recommended package shipped
# with R, but the file exits cleanly if it is missing.

if (!requireNamespace("codetools", quietly = TRUE)) {
  exit_file("codetools not available")
}

ns <- asNamespace("IMPACTncdEngland")
imports_env <- parent.env(ns)

# ---- Direct pin: the import that was missing --------------------------------
expect_true(
  exists("dbListFields", envir = imports_env, inherits = FALSE),
  info = "DBI::dbListFields must be imported (used unqualified in export_summaries helpers)"
)

# ---- Collect every function the package defines -----------------------------
pkg_functions <- list()
for (nm in ls(ns, all.names = TRUE)) {
  obj <- get(nm, envir = ns)
  if (is.function(obj)) {
    pkg_functions[[nm]] <- obj
  } else if (inherits(obj, "R6ClassGenerator")) {
    for (slot in c("public_methods", "private_methods", "active")) {
      for (m in names(obj[[slot]])) {
        if (is.function(obj[[slot]][[m]])) {
          pkg_functions[[paste0(nm, "$", m)]] <- obj[[slot]][[m]]
        }
      }
    }
  }
}

# Sanity: the walk reached the R6 methods, including one added via $set()
# from a satellite file (export_summaries, Simulation_class_summaries.R).
expect_true("Simulation$export_summaries" %in% names(pkg_functions))
expect_true(length(pkg_functions) > 100L)

# ---- General check: every package in Imports --------------------------------
# Widened on review (2026-11-04) from DBI/duckdb to every package the
# DESCRIPTION Imports: the same latent bug existed for igraph
# (get.all.shortest.paths() in get_causal_path()) and parallel (makeCluster()
# in SynthPop$write_synthpop()). Names that also exist in base R's default
# packages are skipped, so a deliberate base call is never flagged.
imports_field <- packageDescription("IMPACTncdEngland")$Imports
pkgs <- trimws(sub("\\(.*", "", strsplit(imports_field, ",")[[1]]))
pkgs <- setdiff(pkgs[nzchar(pkgs)], c("R", "Rcpp", "BH"))
pkg_exports <- unique(unlist(lapply(pkgs, function(p) {
  if (requireNamespace(p, quietly = TRUE)) getNamespaceExports(p)
})))
base_names <- unique(unlist(lapply(
  c("base", "stats", "utils", "methods", "graphics", "grDevices"),
  function(p) ls(asNamespace(p), all.names = TRUE)
)))
imported <- ls(imports_env, all.names = TRUE)

# Known, pre-existing and out of scope here: validate() and
# calibrate_incd_ftlt() are developer tools that pipe arrow queries through
# dplyr verbs (collect(), select(); dplyr is not in Imports at all), so they
# already need dplyr attached; their `%>%` is igraph's re-export. Listed so a
# NEW occurrence anywhere else still fails.
allowed <- list(
  `Simulation$validate` = "%>%",
  `Simulation$calibrate_incd_ftlt` = "%>%"
)

offenders <- character(0)
for (nm in names(pkg_functions)) {
  # Qualified calls (pkg::f) show up as the global `::`, not as `f`, and
  # $-method calls are not globals at all, so only bare calls are reported.
  g <- codetools::findGlobals(pkg_functions[[nm]], merge = FALSE)$functions
  bad <- setdiff(intersect(g, pkg_exports), c(imported, base_names, allowed[[nm]]))
  bad <- bad[!vapply(bad, exists, logical(1), envir = ns, inherits = FALSE)]
  if (length(bad)) {
    offenders <- c(offenders, paste0(nm, ": ", paste(sort(bad), collapse = ", ")))
  }
}

expect_identical(
  offenders, character(0),
  info = paste(
    "Unqualified calls to functions of an Imports package that are not",
    "imported (add them to @importFrom in R/IMPACTncdEngland.R and",
    "re-document, or qualify them):",
    paste(offenders, collapse = "; ")
  )
)
