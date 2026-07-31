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
# `mc` is NOT unique across a multi-scenario run: every scenario logs its own
# "Start/End mc iteration 1..n", so a merge on `mc` alone goes cartesian (a
# 17-scenario x 40-chunk run gives 680 starts each joining 17 ends = 11,560
# rows, and data.table refuses). Segment the log into BLOCKS first -- one per
# parallel section, i.e. one per scenario -- and pair within a block.
#
# Block-awareness also matters for the statistics, not just the join: pooling
# every scenario would measure the "loop" from the first start to the last end
# across the whole run, reporting an 8-hour span at a few percent efficiency.
ev[, block := cumsum(kind == "phase" &
  grepl("^Start of (parallelisation|single-core run)", label))]

st <- ev[kind == "run_start", .(block, mc, start = ts)]
en <- ev[kind == "run_end", .(block, mc, end = ts)]
dup <- st[, .N, by = .(block, mc)][N > 1L]
if (nrow(dup)) {
  warning("mc repeats within a block (", nrow(dup),
          " cases); pairing by order of appearance instead")
  st[, k := seq_len(.N), by = .(block, mc)]
  en[, k := seq_len(.N), by = .(block, mc)]
  runs <- merge(st, en, by = c("block", "mc", "k"), all = TRUE)
} else {
  runs <- merge(st, en, by = c("block", "mc"), all = TRUE)
}
runs[, dur := as.numeric(difftime(end, start, units = "secs"))]

# one row of statistics per block
block_stats <- function(d) {
  d <- d[!is.na(start)]
  if (!nrow(d)) return(NULL)
  tl <- rbind(d[!is.na(start), .(t = start, dd = 1L)],
              d[!is.na(end), .(t = end, dd = -1L)])
  setorder(tl, t, -dd)
  tl[, n := cumsum(dd)]
  tl[, span := c(as.numeric(difftime(t[-1], t[-.N], units = "secs")), 0)]
  peak <- max(tl$n)
  span_run <- if (all(is.na(d$end))) NA_real_ else
    as.numeric(difftime(max(d$end, na.rm = TRUE), min(d$start), units = "secs"))
  setorder(d, start)
  d[, wave := ((seq_len(.N) - 1L) %/% max(1L, peak)) + 1L]
  slack <- d[!is.na(dur), .(s = sum(max(dur) - dur)), keyby = wave][, sum(s)]
  data.table(
    started = nrow(d), completed = sum(!is.na(d$dur)),
    began = min(d$start), span = span_run,
    worker_time = sum(d$dur, na.rm = TRUE),
    med = suppressWarnings(median(d$dur, na.rm = TRUE)),
    peak = peak,
    mean_conc = if (is.na(span_run)) NA_real_ else tl[, sum(n * span)] / span_run,
    idle = if (is.na(span_run)) NA_real_ else
      peak * span_run - sum(d$dur, na.rm = TRUE),
    empty = tl[n == 0, sum(span)], straggler = slack)
}

if (nrow(runs)) {
  bs <- rbindlist(lapply(split(runs, by = "block", keep.by = TRUE),
                         block_stats), idcol = "block")
  nblk <- nrow(bs)
  cat("\n=== RUN PHASE (", nblk,
      if (nblk != 1L) " parallel sections, one per scenario" else
        " parallel section", ") ===\n", sep = "")
  cat("iterations completed : ", sum(bs$completed), " / ", sum(bs$started), "\n", sep = "")

  if (nblk > 1L) {
    cat("\n-- per section --\n")
    print(bs[, .(block, began = format(began, "%m-%d %H:%M:%S"),
                 n = started, done = completed,
                 span = fmt_dur(span), median_iter = fmt_dur(med),
                 conc = round(mean_conc, 2),
                 straggler = fmt_dur(straggler))], nrows = 100)
    slow <- bs[which.max(span)]
    if (nrow(slow) && is.finite(slow$span) &&
          slow$span > 1.5 * median(bs$span, na.rm = TRUE))
      cat("\n!! section ", slow$block, " took ", fmt_dur(slow$span),
          " against a median of ", fmt_dur(median(bs$span, na.rm = TRUE)),
          " -- check for contention from other jobs\n", sep = "")
  }

  cat("\n-- overall --\n")
  cat("duration  min/med/max: ", fmt_dur(min(runs$dur, na.rm = TRUE)), " / ",
    fmt_dur(median(runs$dur, na.rm = TRUE)), " / ",
    fmt_dur(max(runs$dur, na.rm = TRUE)), "\n", sep = "")
  cat("total worker-time    : ", fmt_dur(sum(runs$dur, na.rm = TRUE)), "\n", sep = "")
  cat("summed section span  : ", fmt_dur(sum(bs$span, na.rm = TRUE)), "\n", sep = "")
  cat("peak concurrency     : ", max(bs$peak), " workers\n", sep = "")
  cat("mean concurrency     : ",
    round(sum(bs$mean_conc * bs$span, na.rm = TRUE) /
            sum(bs$span, na.rm = TRUE), 2),
    " workers  (efficiency ",
    round(100 * sum(runs$dur, na.rm = TRUE) /
            (max(bs$peak) * sum(bs$span, na.rm = TRUE)), 1),
    "% of the pool)\n", sep = "")
  cat("idle worker-time lost: ", fmt_dur(sum(bs$idle, na.rm = TRUE)), "\n", sep = "")
  cat("time with 0 workers running inside the sections: ",
    fmt_dur(sum(bs$empty, na.rm = TRUE)), "\n", sep = "")
  cat("straggler cost (sum over waves of max-minus-each): ",
    fmt_dur(sum(bs$straggler, na.rm = TRUE)), " of worker-time\n", sep = "")

  # tail gap per section: last worker finished -> next phase mark. This is the
  # one CLAUDE.md warns about (a Wales run lost 17h29m to a stopped tty here).
  gaps <- rbindlist(lapply(seq_len(nblk), function(b) {
    d <- runs[block == bs$block[b] & !is.na(end)]
    if (!nrow(d)) return(NULL)
    le <- max(d$end)
    nx <- ph[ts > le][1L]
    if (is.na(nx$ts)) return(NULL)
    data.table(block = bs$block[b], last_end = le, next_label = nx$label,
               next_ts = nx$ts,
               gap = as.numeric(difftime(nx$ts, le, units = "secs")))
  }))
  if (nrow(gaps)) {
    big <- gaps[gap > 60]
    if (nrow(big)) {
      cat("\n!! GAPS after the last worker finished (should be ~0s):\n")
      print(big[order(-gap), .(block, last_end = format(last_end, "%m-%d %H:%M:%S"),
                               until = next_label, gap = fmt_dur(gap))])
      cat("   total dead time: ", fmt_dur(sum(big$gap)), "  [",
          round(100 * sum(big$gap) / total, 1), "% of wall clock]\n", sep = "")
    } else {
      cat("\ntail gaps after each section: all under 60s (max ",
          fmt_dur(max(gaps$gap)), ") -- no dead time\n", sep = "")
    }
  }
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
  # With a saturated pool, task i ends when task i+k starts. That inference needs
  # MORE tasks than workers -- if the pool is at least as large as the task list
  # every task runs in a single wave, nothing is ever queued behind anything, and
  # there are no start-gaps to read a duration from. Say so rather than printing
  # the NA that median(numeric(0)) produces.
  if (nrow(exp_ev) > k) {
    d <- as.numeric(difftime(exp_ev$ts[-seq_len(k)],
      exp_ev$ts[seq_len(nrow(exp_ev) - k)], units = "secs"))
    cat("inferred per-task    : median ", fmt_dur(median(d)), ", max ",
      fmt_dur(max(d)), "\n", sep = "")
    cat("implied worker-time  : ", fmt_dur(median(d) * nrow(exp_ev)), " over ", k,
      " workers = ", fmt_dur(median(d) * nrow(exp_ev) / k), " ideal wall\n", sep = "")
  } else {
    cat("inferred per-task    : not inferable -- ", nrow(exp_ev), " task",
      if (nrow(exp_ev) != 1L) "s" else "", " on a pool of ", k,
      " runs in ONE wave, so no task waits behind another.\n", sep = "")
    cat("                       the wall span above IS the cost of that wave.\n")
  }
}

cat("\n")
