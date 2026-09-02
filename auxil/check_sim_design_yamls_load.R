## Check that every sim_design*.yaml in the repo still constructs a valid
## Design, and report the values it resolves for the OPTIONAL parameters.
##
## Why this exists: optional `sim_prm` keys are defaulted inside
## `Design$initialize()` rather than listed in `required_params`, precisely so
## that a YAML written before a key existed keeps loading. That contract is easy
## to break -- adding a key to `required_params`, or validating a default
## against something a given design cannot supply -- and the failure is not
## visible until someone runs an older scenario. This walks every shipped YAML
## and proves it.
##
## Designs are built from a CLEAN working directory so `has_input_data()` is
## FALSE and the constructor stays in skeleton mode: no data load, so this runs
## in seconds regardless of where it is launched from.
##
## Usage: Rscript auxil/check_sim_design_yamls_load.R

suppressMessages({
  library(yaml)
  library(IMPACTncdEngland)
})

repo <- normalizePath(".")
yamls <- setdiff(
  list.files(repo, pattern = "^sim_design.*\\.yaml$", recursive = TRUE,
             full.names = TRUE),
  character(0)
)
yamls <- yamls[!grepl("/\\.history/|/claude_process/|/graphify-out/", yamls)]
yamls <- sort(yamls)

clean_wd <- file.path(tempdir(), "yaml_check_clean_wd")
dir.create(clean_wd, showWarnings = FALSE, recursive = TRUE)
old_wd <- setwd(clean_wd)
on.exit(setwd(old_wd), add = TRUE)

# The optional keys to report. Extend when a new one is added.
OPTIONAL <- c("export_xps_quantiles", "xps_quantile_probs",
              "xps_quantile_strata")

fmt <- function(v) {
  if (is.null(v)) return("<NULL>")
  if (is.list(v)) {
    return(paste0("[", paste(vapply(v, function(e)
      paste0("(", paste(e, collapse = ","), ")"), character(1)),
      collapse = " "), "]"))
  }
  paste(v, collapse = ",")
}

cat(sprintf("checking %d sim_design YAMLs\n\n", length(yamls)))
fails <- character(0)
for (y in yamls) {
  rel <- sub(paste0(repo, "/"), "", y, fixed = TRUE)
  dl <- tryCatch(read_yaml(y), error = function(e) e)
  if (inherits(dl, "error")) {
    cat(sprintf("  FAIL  %-56s unreadable: %s\n", rel, conditionMessage(dl)))
    fails <- c(fails, rel)
    next
  }
  present <- OPTIONAL[OPTIONAL %in% names(dl)]
  d <- tryCatch(suppressMessages(Design$new(dl)), error = function(e) e)
  if (inherits(d, "error")) {
    cat(sprintf("  FAIL  %-56s %s\n", rel, conditionMessage(d)))
    fails <- c(fails, rel)
    next
  }
  cat(sprintf("  ok    %-56s %s\n", rel,
              if (length(present)) paste0("declares: ",
                                          paste(present, collapse = ", "))
              else "declares none of the optional keys"))
  for (k in OPTIONAL) {
    cat(sprintf("          %-24s %s%s\n", k, fmt(d$sim_prm[[k]]),
                if (k %in% present) "  (from the YAML)" else "  (default)"))
  }
}

cat("\n")
if (length(fails)) {
  cat(sprintf("%d of %d FAILED: %s\n", length(fails), length(yamls),
              paste(fails, collapse = ", ")))
  quit(status = 1L)
}
cat(sprintf("all %d YAMLs load; every optional key resolves.\n", length(yamls)))
