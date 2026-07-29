# Tests for private$start_console_log() / private$stop_console_log(), which
# redirect R console output (stdout *and* messages, including foreach's
# `.verbose` bookkeeping) to <output_dir>/logs/console.txt while run(),
# export_summaries() and export_tables() are working.
#
# Regression: with `logs: yes` foreach runs with `.verbose = TRUE` and dumps all
# of its per-task bookkeeping in one ~21 KB burst when the loop returns. If
# stdout is a terminal whose output has been stopped by a stray ^S -- which
# survives detaching a screen session -- that burst blocks in write() forever.
# A Wales run lost 17h29m to exactly this on 2026-07-28. Writing to a file
# removes the terminal from the path entirely.
#
# The two methods are exercised through a minimal stub that binds the *real*
# implementations out of the Simulation generator, so no input data, synthpop
# or YAML is needed.

library(tinytest)
library(IMPACTncdEngland)

gen <- IMPACTncdEngland::Simulation
pm <- gen$private_methods

expect_true(is.function(pm$start_console_log),
            info = "start_console_log is a private method of Simulation")
expect_true(is.function(pm$stop_console_log),
            info = "stop_console_log is a private method of Simulation")

Stub <- R6::R6Class(
  "ConsoleLogStub",
  public = list(
    design = NULL,
    initialize = function(output_dir, logs) {
      self$design <- list(sim_prm = list(output_dir = output_dir, logs = logs))
    },
    # thin public shims so the test can drive the private methods
    start = function(phase = "phase") private$start_console_log(phase),
    stop = function() private$stop_console_log(),
    handle = function() private$console_log,
    # a body that writes and then fails, mirroring how run() uses on.exit()
    run_failing = function() {
      private$start_console_log("failing phase")
      on.exit(private$stop_console_log(), add = TRUE)
      cat("stdout inside the redirect\n")
      message("message inside the redirect")
      stop("boom")
    },
    # mirrors run()'s own safety cleanup, which pops sinks leaked by run_sim()
    # before stop_console_log() gets a chance to. Bounded to the depth seen on
    # entry, exactly as run() does it.
    run_with_safety_cleanup = function() {
      entry <- sink.number(type = "output")
      private$start_console_log("safety")
      on.exit(private$stop_console_log(), add = TRUE)
      cat("captured\n")
      if (sink.number(type = "message") > 2L) sink(type = "message")
      while (sink.number(type = "output") > entry) sink(type = "output")
      invisible(TRUE)
    }
  ),
  private = list(
    console_log = NULL,
    output_dir = pm$output_dir,
    start_console_log = pm$start_console_log,
    stop_console_log = pm$stop_console_log
  )
)

# Sink depths must be pristine before and after every block, otherwise a leak
# in one test silently swallows the output of the next.
clean_sinks <- function() {
  sink.number(type = "output") == 0L && sink.number(type = "message") == 2L
}
expect_true(clean_sinks(), info = "no sinks active at the start of the file")

# --- logs = FALSE: complete no-op ------------------------------------------
d <- file.path(tempdir(), paste0("clog_off_", as.integer(runif(1, 1, 1e6))))
s <- Stub$new(d, FALSE)
s$start("should not happen")
expect_null(s$handle(), info = "logs = FALSE leaves no console-log handle")
expect_false(file.exists(file.path(d, "logs", "console.txt")),
             info = "logs = FALSE writes no console.txt")
expect_true(clean_sinks(), info = "logs = FALSE touches no sinks")
expect_silent(s$stop())
expect_true(clean_sinks(), info = "stop() on a no-op redirect is harmless")

# --- logs = TRUE: captures stdout and messages ------------------------------
d <- file.path(tempdir(), paste0("clog_on_", as.integer(runif(1, 1, 1e6))))
s <- Stub$new(d, TRUE)
f <- file.path(d, "logs", "console.txt")

s$start("run(sc0)")
cat("STDOUT_MARKER\n")
message("MESSAGE_MARKER")
s$stop()

expect_true(clean_sinks(), info = "sinks fully unwound after stop()")
expect_true(file.exists(f), info = "console.txt is created under logs/")
txt <- readLines(f, warn = FALSE)
expect_true(any(grepl("STDOUT_MARKER", txt)), info = "stdout is captured")
expect_true(any(grepl("MESSAGE_MARKER", txt)), info = "messages are captured")
expect_true(any(grepl("===== run\\(sc0\\) at:", txt)),
            info = "the phase banner is written")
expect_null(s$handle(), info = "handle cleared after stop()")

# --- appends across phases rather than truncating ---------------------------
s$start("export_summaries()")
cat("SECOND_PHASE\n")
s$stop()
txt <- readLines(f, warn = FALSE)
expect_true(any(grepl("STDOUT_MARKER", txt)) && any(grepl("SECOND_PHASE", txt)),
            info = "a later phase appends, it does not truncate")
expect_equal(sum(grepl("^===== ", txt)), 2L,
             info = "one banner per phase")

# --- start() does not nest --------------------------------------------------
s$start("outer")
h1 <- s$handle()
depth_after_first <- sink.number(type = "output")
s$start("inner")
depth_after_second <- sink.number(type = "output")
expect_identical(h1$path, s$handle()$path,
                 info = "a nested start() keeps the outer handle")
expect_equal(depth_after_second, depth_after_first,
             info = "a nested start() pushes no extra sink")
s$stop()
expect_true(clean_sinks(), info = "a single stop() undoes a nested start()")

# --- an error still reaches the caller, and unwinds the sinks ---------------
# This is the safety property that makes the redirect acceptable: without it a
# fatal error would be written into the file and the terminal would go silent.
d <- file.path(tempdir(), paste0("clog_err_", as.integer(runif(1, 1, 1e6))))
s <- Stub$new(d, TRUE)
err <- tryCatch(s$run_failing(), error = function(e) conditionMessage(e))
expect_equal(err, "boom", info = "the error propagates to the caller")
expect_true(clean_sinks(), info = "sinks unwound even when the body throws")
txt <- readLines(file.path(d, "logs", "console.txt"), warn = FALSE)
expect_true(any(grepl("stdout inside the redirect", txt)),
            info = "output written before the error is still captured")

# --- tolerates run()'s safety cleanup having already popped the sinks -------
d <- file.path(tempdir(), paste0("clog_safety_", as.integer(runif(1, 1, 1e6))))
s <- Stub$new(d, TRUE)
s$run_with_safety_cleanup()
expect_true(clean_sinks(),
            info = "no over-popping after the safety cleanup already unwound")

# --- a CALLER's sink must survive ------------------------------------------
# run()'s safety cleanup used to unwind to depth 0, which destroys sinks the
# caller owns: run() inside capture.output()/knitr/expect_silent() would eat
# the rest of the caller's output. It must unwind only to its entry depth.
# Assertions are made *after* the outer sink is popped, or they would be
# captured by it.
outer_f <- file.path(tempdir(), paste0("clog_outer_", as.integer(runif(1, 1, 1e6))))
oc <- file(outer_f, open = "wt")
sink(oc, type = "output")
depth_with_caller_sink <- sink.number(type = "output")
Stub$new(file.path(tempdir(), "clog_inner"), TRUE)$run_with_safety_cleanup()
depth_after_run <- sink.number(type = "output")
sink(type = "output")
close(oc)

expect_equal(depth_after_run, depth_with_caller_sink,
             info = "the caller's sink survives run()'s safety cleanup")
expect_true(clean_sinks(), info = "back to a clean stack afterwards")

# --- stop() without start() is a no-op --------------------------------------
s <- Stub$new(file.path(tempdir(), "clog_never"), TRUE)
expect_silent(s$stop())
expect_true(clean_sinks(), info = "stop() without start() leaves sinks alone")

# --- an unwritable output_dir degrades to a warning, not an error -----------
s <- Stub$new(file.path("/proc/definitely_not_writable"), TRUE)
expect_warning(s$start("doomed"), pattern = "console log")
expect_null(s$handle(), info = "no handle when the log could not be opened")
expect_true(clean_sinks(), info = "a failed start() leaves no dangling sink")
