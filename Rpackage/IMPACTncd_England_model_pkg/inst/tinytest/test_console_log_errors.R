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
# What is asserted is the state AT THE MOMENT R WOULD PRINT THE ERROR: an
# outer calling handler runs after the package's inner one and before any
# unwinding (on.exit), which is exactly when R's default printer runs. So a
# handler that kept the error in console.txt fails here, and one that gets it
# to the terminal passes, whatever it writes to console.txt.
#
# run() is wrapped the same way but not exercised here: on a skeleton object
# it stops on `data_initialized` BEFORE its console log starts, and a real run
# needs the model data.

library(tinytest)
library(IMPACTncdEngland)

old_wd <- getwd()
# Simulation$new() and export_summaries() change the session's thread
# settings: save them now, restore at the end.
state0 <- list(dt = data.table::getDTthreads(), fst = fst::threads_fst(),
               arrow = arrow::cpu_count())

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
con_txt <- file.path(out, "logs", "console.txt")

msg_con0 <- sink.number(type = "message")
out_depth0 <- sink.number(type = "output")

# Call `expr`, and record the message connection when the error reaches an
# outer calling handler, i.e. when R's printer would run.
state_at_print <- function(expr) {
  at <- NA_integer_
  err <- tryCatch(
    withCallingHandlers({ force(expr); "<no error>" },
                        error = function(e) at <<- sink.number(type = "message")),
    error = conditionMessage
  )
  list(err = err, msg_con_at_print = at)
}

# export_summaries(): the lifecourse folder is empty, so it stops inside its
# wrapped body, with the console log active.
r <- state_at_print(sim$export_summaries(multicore = FALSE))
expect_false(identical(r$err, "<no error>"), info = "export_summaries(): the fixture makes it fail")
expect_identical(r$msg_con_at_print, msg_con0,
                 info = "export_summaries(): the message stream is on the terminal when the error is printed")
expect_identical(sink.number(type = "output"), out_depth0, info = "no output sink left behind")
log_lines <- if (file.exists(con_txt)) readLines(con_txt, warn = FALSE) else character(0)
expect_true(any(grepl("export_summaries()", log_lines, fixed = TRUE)),
            info = "the console log was active for the call")
expect_true(any(startsWith(log_lines, "Error: ")),
            info = "the error is also recorded in console.txt")

# export_tables(): an invalid equity_ridit_reference fails its match.arg(),
# the first statement inside the wrapped body, with the console log active.
r <- state_at_print(sim$export_tables(multicore = FALSE, equity_ridit_reference = "bogus"))
expect_false(identical(r$err, "<no error>"), info = "export_tables(): the fixture makes it fail")
expect_identical(r$msg_con_at_print, msg_con0,
                 info = "export_tables(): the message stream is on the terminal when the error is printed")

# A caller's own sink survives an error too.
tf <- tempfile()
sink(tf)
invisible(tryCatch(sim$export_summaries(multicore = FALSE), error = function(e) NULL))
depth_inside <- sink.number(type = "output")
sink()
unlink(tf)
expect_identical(depth_inside, out_depth0 + 1L,
                 info = "the caller's sink is still there after the error")

setwd(old_wd)
unlink(root, recursive = TRUE)
data.table::setDTthreads(state0$dt)
fst::threads_fst(state0$fst)
arrow::set_cpu_count(state0$arrow)
