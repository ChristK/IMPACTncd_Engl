# Tests: with `logs: yes`, an error inside run() / export_summaries() /
# export_tables() is released from the console log before R prints it.
#
# Regression (2026-11-04, review of the issue #11 guard): start_console_log()
# diverts the message stream to <output_dir>/logs/console.txt. R prints an
# uncaught error to whichever message sink is active when stop() is called
# -- BEFORE on.exit() runs -- so the error went only to console.txt and the
# terminal showed just "Execution halted". The methods' bodies are now
# wrapped in a calling handler (private$release_console_log_on_error()) that
# writes "Error: <message>" to console.txt and releases the redirect first.
#
# The terminal print itself is R's own error printer and cannot be captured
# here; what is asserted is the mechanism: the redirect is released, the
# caller's sinks are intact, and console.txt holds the "Error:" line, which
# only the handler writes (on.exit() alone would leave console.txt without it).

library(tinytest)
library(IMPACTncdEngland)

old_wd <- getwd()
root <- normalizePath(tempfile("conlog_err_"), mustWork = FALSE)
dir.create(root)
setwd(root)
out <- file.path(root, "out")

dl <- list(
  locality = "England", clusternumber = 1, clusternumber_export = 1,
  logs = TRUE, cols_for_output = c("year", "age", "sex"),
  strata_for_output = c("year", "age", "sex"),
  exposures_for_output = list("bmi"),
  exposure_definitions = list(list(
    name = "bmi", file_name = "bmi", var_name = "bmi",
    rank_var = "rank_bmi", distribution = "continuous"
  )),
  n = 10, init_year_long = 2010, sim_horizon_max = 2030, ageL = 35, ageH = 80,
  diseases = list(chd = list(name = "chd", meta = list(
    incidence = list(influenced_by_disease_name = character(0))
  ))),
  maxlag = 5, smoking_relapse_limit = 3, stochastic = TRUE, kismet = 42,
  jumpiness = 0.2, statin_adherence = 0.85, bpmed_adherence = 0.9,
  export_xps = FALSE, export_PARF = FALSE, simsmok_calibration = FALSE,
  output_dir = out, synthpop_dir = file.path(root, "sp"),
  validation = FALSE, iteration_n_max = 2, num_chunks = 1
)
design <- suppressMessages(suppressWarnings(Design$new(dl)))
sim <- suppressMessages(suppressWarnings(Simulation$new(design)))

out_depth <- sink.number(type = "output")
msg_con <- sink.number(type = "message")

# The lifecourse folder is empty, so export_summaries() stops inside its
# wrapped body ("No Monte Carlo iterations found", or a DuckDB error on the
# empty dataset), with the console log active.
err <- tryCatch({ sim$export_summaries(multicore = FALSE); "<no error>" },
                error = conditionMessage)
expect_false(identical(err, "<no error>"), info = "the fixture does make it fail")
expect_identical(sink.number(type = "output"), out_depth,
                 info = "no output sink left behind")
expect_identical(sink.number(type = "message"), msg_con,
                 info = "the message stream is back on the terminal")

con_txt <- file.path(out, "logs", "console.txt")
log_lines <- if (file.exists(con_txt)) readLines(con_txt, warn = FALSE) else character(0)
expect_true(any(grepl("export_summaries()", log_lines, fixed = TRUE)),
            info = "the console log was active for the call")
expect_true(any(startsWith(log_lines, "Error: ")),
            info = "the handler recorded the error in console.txt before releasing it")

# A caller's own sink survives an error too.
tf <- tempfile()
sink(tf)
invisible(tryCatch(sim$export_summaries(multicore = FALSE), error = function(e) NULL))
depth_inside <- sink.number(type = "output")
sink()
expect_identical(depth_inside, out_depth + 1L,
                 info = "the caller's sink is still there after the error")

setwd(old_wd)
unlink(root, recursive = TRUE)
