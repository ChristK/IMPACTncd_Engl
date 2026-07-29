## scan_times_gaps.R ---------------------------------------------------------
## Companion to auxil/profile_times_log.R.
##
## Scans many `<output_dir>/logs/times.txt` files and, for every
## "Start of parallelisation" ... "End of parallelisation" block, reports the
## TAIL GAP: the wall time between the last worker writing "End mc iteration N"
## and the parent writing "End of parallelisation".
##
## The parent does nothing between those two marks except return from
## `foreach(...) %dopar% {...}` (Simulation_class.R ~L505-531), so a tail gap is
## pure overhead in mclapply result collection / foreach accumulation. Comparing
## the gap across runs of different sizes tells us whether it scales with the
## number of tasks (systematic) or is a one-off environmental stall.
##
## Usage:
##   Rscript auxil/scan_times_gaps.R <times.txt> [<times.txt> ...]
##   Rscript auxil/scan_times_gaps.R $(find /data -name times.txt -path '*logs*')

library(data.table)

files <- commandArgs(trailingOnly = TRUE)
if (!length(files)) stop("Supply one or more times.txt paths")

fmt_dur <- function(secs) {
  secs <- as.numeric(secs)
  h <- floor(secs / 3600)
  m <- floor((secs - 3600 * h) / 60)
  s <- round(secs - 3600 * h - 60 * m)
  fifelse(is.na(secs), NA_character_,
    fifelse(h > 0, sprintf("%dh%02dm", h, m),
      fifelse(m > 0, sprintf("%dm%02ds", m, s), sprintf("%ds", s))))
}

scan_one <- function(f) {
  raw <- tryCatch(readLines(f, warn = FALSE), error = function(e) character())
  raw <- raw[nzchar(trimws(raw))]
  if (!length(raw)) return(NULL)

  d <- data.table(
    label = sub(" at: [0-9].*$", "", raw),
    ts = as.POSIXct(sub("^.* at: ", "", raw), tz = "",
      format = "%Y-%m-%d %H:%M:%OS")
  )
  d <- d[!is.na(ts)]
  if (!nrow(d)) return(NULL)

  # segment into parallelisation blocks; a file accumulates one block per
  # run() call (i.e. per scenario), appended over time
  d[, blk := cumsum(label == "Start of parallelisation")]
  d <- d[blk > 0L]
  if (!nrow(d)) return(NULL)

  rbindlist(lapply(split(d, by = "blk", keep.by = TRUE), function(b) {
    starts <- b[grepl("^Start mc iteration [0-9]+$", label)]
    ends <- b[grepl("^End mc iteration [0-9]+$", label)]
    fin <- b[label == "End of parallelisation", ts][1L]
    if (!nrow(ends) || is.na(fin)) return(NULL)

    # realised pool size = peak concurrency of the start/end step function
    tl <- rbind(starts[, .(t = ts, dlt = 1L)], ends[, .(t = ts, dlt = -1L)])
    setorder(tl, t, -dlt)
    pool <- max(cumsum(tl$dlt))

    last_end <- max(ends$ts)
    loop_span <- as.numeric(difftime(last_end, min(starts$ts), units = "secs"))
    gap <- as.numeric(difftime(fin, last_end, units = "secs"))
    data.table(
      file = f, run = b$blk[1L],
      n_iter = nrow(ends), pool = pool,
      loop_span_s = loop_span, gap_s = gap,
      gap_pct = 100 * gap / (loop_span + gap),
      finished = format(fin, "%Y-%m-%d %H:%M")
    )
  }), fill = TRUE)
}

res <- rbindlist(lapply(files, scan_one), fill = TRUE)
if (!nrow(res)) stop("No parsable parallelisation blocks found")

setorder(res, -gap_s)
res[, `:=`(
  loop_span = fmt_dur(loop_span_s),
  gap = fmt_dur(gap_s),
  gap_pct = round(gap_pct, 1),
  file = sub("/logs/times.txt$", "", file)
)]

cat("\nTail gap = last 'End mc iteration' -> 'End of parallelisation'\n")
cat("(the parent only returns from foreach in between; anything large is overhead)\n\n")
print(res[, .(file, run, n_iter, pool, loop_span, gap, gap_pct, finished)],
  nrows = 200, class = FALSE)

cat("\n--- gap vs number of tasks ---\n")
print(res[, .(n_iter, pool, gap_s = round(gap_s, 1),
  gap_s_per_task = round(gap_s / n_iter, 3))][order(n_iter)],
  nrows = 200, class = FALSE)
cat("\n")
