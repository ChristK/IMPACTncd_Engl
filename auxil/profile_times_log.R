## profile_times_log.R -------------------------------------------------------
## Parse an IMPACTncd `<output_dir>/logs/times.txt` and report where the wall
## clock actually went: per-phase durations, per-mc-iteration durations, the
## realised worker concurrency, parallel efficiency, and any dead gaps during
## which no worker was running.
##
## `times.txt` is written by Simulation$..$time_mark() (Simulation_class.R).
## The run phase logs a matched "Start/End mc iteration <n>" pair per synthpop
## iteration; the summary-export phase logs *starts only*, so export durations
## are inferred from a saturated-pool assumption (see `pool_export`).
##
## Usage:
##   Rscript auxil/profile_times_log.R <path/to/times.txt> [pool_export]
## e.g.
##   Rscript auxil/profile_times_log.R /mnt/storage_fast4/Wales/output/logs/times.txt 5

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
times_file <- if (length(args) >= 1L) args[[1L]] else
  stop("Supply the path to a times.txt")
pool_export <- if (length(args) >= 2L) as.integer(args[[2L]]) else NA_integer_

fmt_dur <- function(secs) {
  secs <- as.numeric(secs)
  h <- floor(secs / 3600)
  m <- floor((secs - 3600 * h) / 60)
  s <- secs - 3600 * h - 60 * m
  ifelse(h > 0, sprintf("%dh %02dm %02ds", h, m, round(s)),
    ifelse(m > 0, sprintf("%dm %02ds", m, round(s)), sprintf("%.1fs", s))
  )
}

# ---- parse ----------------------------------------------------------------
raw <- readLines(times_file, warn = FALSE)
raw <- raw[nzchar(trimws(raw))]

ev <- data.table(
  label = sub(" at: [0-9].*$", "", raw),
  ts = as.POSIXct(sub("^.* at: ", "", raw), tz = "", format = "%Y-%m-%d %H:%M:%OS")
)
ev <- ev[!is.na(ts)]
setorder(ev, ts)

ev[, kind := fifelse(grepl("^Start mc iteration [0-9]+$", label), "run_start",
  fifelse(grepl("^End mc iteration [0-9]+$", label), "run_end",
    fifelse(grepl("^Start mc iteration \\(summary export\\)", label), "exp_start",
      "phase"
    )
  )
)]
ev[kind != "phase", mc := as.integer(sub("^.*?([0-9]+)$", "\\1", label))]

# ---- top-level phase timeline ---------------------------------------------
ph <- ev[kind == "phase"]
ph[, dur_to_next := c(as.numeric(difftime(ts[-1], ts[-.N], units = "secs")), NA_real_)]

cat("\n=== PHASE TIMELINE ===\n")
print(ph[, .(phase = label, at = format(ts, "%Y-%m-%d %H:%M:%S"),
  until_next = fmt_dur(dur_to_next))], nrows = 100)

total <- as.numeric(difftime(max(ev$ts), min(ev$ts), units = "secs"))
cat("\nTotal wall clock (first -> last mark): ", fmt_dur(total), "\n", sep = "")

# ---- run phase: per-iteration durations ------------------------------------
runs <- merge(
  ev[kind == "run_start", .(mc, start = ts)],
  ev[kind == "run_end", .(mc, end = ts)],
  by = "mc", all = TRUE
)
runs[, dur := as.numeric(difftime(end, start, units = "secs"))]

if (nrow(runs)) {
  cat("\n=== RUN PHASE (per synthpop mc iteration) ===\n")
  cat("iterations completed : ", sum(!is.na(runs$dur)), " / ", nrow(runs), "\n", sep = "")
  cat("duration  min/med/max: ", fmt_dur(min(runs$dur, na.rm = TRUE)), " / ",
    fmt_dur(median(runs$dur, na.rm = TRUE)), " / ",
    fmt_dur(max(runs$dur, na.rm = TRUE)), "\n", sep = "")
  cat("total worker-time    : ", fmt_dur(sum(runs$dur, na.rm = TRUE)), "\n", sep = "")

  # realised concurrency as a step function over start/end events
  tl <- rbind(
    runs[!is.na(start), .(t = start, d = 1L)],
    runs[!is.na(end), .(t = end, d = -1L)]
  )
  setorder(tl, t, -d)
  tl[, n := cumsum(d)]
  tl[, span := c(as.numeric(difftime(t[-1], t[-.N], units = "secs")), 0)]

  busy <- tl[n > 0, sum(span)]
  span_run <- as.numeric(difftime(max(runs$end, na.rm = TRUE),
    min(runs$start, na.rm = TRUE), units = "secs"))
  peak <- max(tl$n)
  mean_conc <- tl[, sum(n * span)] / span_run

  cat("\nwall span of the loop: ", fmt_dur(span_run), "\n", sep = "")
  cat("peak concurrency     : ", peak, " workers\n", sep = "")
  cat("mean concurrency     : ", round(mean_conc, 2), " workers  (efficiency ",
    round(100 * mean_conc / peak, 1), "% of the pool)\n", sep = "")
  cat("idle worker-time lost: ", fmt_dur(peak * span_run - sum(runs$dur, na.rm = TRUE)),
    "\n", sep = "")

  # time spent with the pool NOT saturated, and fully empty
  cat("time with 0 workers running inside the loop: ",
    fmt_dur(tl[n == 0, sum(span)]), "\n", sep = "")

  # tail gap: last worker finished -> next phase mark
  last_end <- max(runs$end, na.rm = TRUE)
  nxt <- ph[ts > last_end][1L]
  if (!is.na(nxt$ts)) {
    cat("\n!! GAP after the last worker finished (", format(last_end, "%H:%M:%S"),
      ") until '", nxt$label, "' (", format(nxt$ts, "%H:%M:%S"), "): ",
      fmt_dur(as.numeric(difftime(nxt$ts, last_end, units = "secs"))),
      "  [", round(100 * as.numeric(difftime(nxt$ts, last_end, units = "secs")) / total, 1),
      "% of total wall clock]\n", sep = "")
  }

  # wave structure: with a saturated pool of size `peak`, iterations run in
  # waves; a wave costs max(), not mean(), of its members.
  setorder(runs, start)
  runs[, wave := ((seq_len(.N) - 1L) %/% peak) + 1L]
  wv <- runs[!is.na(dur), .(
    n = .N, wave_cost = max(dur), slack = sum(max(dur) - dur)
  ), keyby = wave]
  cat("\nstraggler cost (sum over waves of max-minus-each): ",
    fmt_dur(sum(wv$slack)), " of worker-time\n", sep = "")
}

# ---- summary export phase ---------------------------------------------------
exp_ev <- ev[kind == "exp_start"][order(ts)]
if (nrow(exp_ev)) {
  cat("\n=== SUMMARY EXPORT PHASE ===\n")
  end_mark <- ph[grepl("^End of exporting summaries", label), ts][1L]
  span_exp <- as.numeric(difftime(end_mark, min(exp_ev$ts), units = "secs"))
  cat("iterations           : ", nrow(exp_ev), "\n", sep = "")
  cat("wall span            : ", fmt_dur(span_exp), "\n", sep = "")

  # infer pool size from the initial burst of near-simultaneous starts, unless
  # the caller supplied it (clusternumber_export in the sim design yaml)
  gaps <- as.numeric(diff(exp_ev$ts), units = "secs")
  k <- if (!is.na(pool_export)) pool_export else sum(gaps[seq_len(min(10L, length(gaps)))] < 5) + 1L
  cat("assumed pool size    : ", k, " workers\n", sep = "")
  # with a saturated pool, task i ends when task i+k starts
  d <- as.numeric(difftime(exp_ev$ts[-seq_len(k)], exp_ev$ts[seq_len(nrow(exp_ev) - k)],
    units = "secs"))
  cat("inferred per-task    : median ", fmt_dur(median(d)), ", max ", fmt_dur(max(d)),
    "\n", sep = "")
  cat("implied worker-time  : ", fmt_dur(median(d) * nrow(exp_ev)), " over ", k,
    " workers = ", fmt_dur(median(d) * nrow(exp_ev) / k), " ideal wall\n", sep = "")
}

cat("\n")
