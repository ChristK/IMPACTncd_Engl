#!/bin/bash

#===============================================================================
# MEMORY PROFILER FOR R SCRIPTS
#===============================================================================
# Author: IMPACTncd Japan Team
# Purpose: Monitor memory usage (RSS and VSZ) of R scripts during execution
#
# This script runs an R script while continuously monitoring its memory usage
# and generates both a log file and a visualization plot of memory consumption
# over time.
#
# USAGE:
#   ./mem_profile.sh [OPTIONS] <R_script> <output_log_file>
#
# OPTIONS:
#   -p, --parallel  Enable parallel mode: track memory across ALL child processes
#                   spawned by the R script (e.g., from parallel::mclapply,
#                   foreach, future, etc.). Without this flag, only the main
#                   R process is monitored.
#   -h, --help      Show this help message and exit
#
# ARGUMENTS:
#   R_script        - Path to the R script to profile
#   output_log_file - Path where memory log will be saved (*.txt)
#
# OUTPUTS:
#   1. Memory log file (*.txt) - Timestamped memory usage data
#   2. Memory plot (*.png) - Visual representation of memory usage over time
#
# MEMORY METRICS:
#   RSS (Resident Set Size) - Physical memory currently used by the process
#   VSZ (Virtual Memory Size) - Total virtual memory used by the process
#   MEM% - Percentage of system memory used by the process
#   (In parallel mode, RSS and VSZ are summed across all child processes)
#
# EXAMPLES:
#   # Profile a single-threaded R script:
#   ./mem_profile.sh simulate.R memory_log.txt
#
#   # Profile a parallel R script (tracks all forked workers):
#   ./mem_profile.sh --parallel simulate.R memory_log.txt
#   ./mem_profile.sh -p calibrate.R /tmp/calib_memory.txt
#
# NOTES:
#   - Sampling interval: 0.1 seconds (100ms)
#   - Memory values are reported in GB for readability
#   - Script can be interrupted with Ctrl+C safely
#   - Requires R and basic Unix utilities (ps, awk, date)
#
# PARALLEL MODE DETAILS:
#   When -p/--parallel is enabled, the profiler tracks:
#   - The main R process
#   - All direct child processes (first-level workers)
#   - All descendant processes (nested parallelism)
#
#   This is essential for accurate memory profiling when using:
#   - parallel::mclapply(), parallel::parLapply()
#   - foreach with doParallel/doFuture backends
#   - future::plan(multisession/multicore)
#   - Any other fork-based or spawn-based parallelism
#
#   Note: Child process detection uses /proc filesystem on Linux. The script
#   will fall back to single-process tracking if child detection fails.
#===============================================================================

# Function to display help message
show_help() {
  echo "MEMORY PROFILER FOR R SCRIPTS"
  echo ""
  echo "USAGE: $0 [OPTIONS] <R_script> <output_log_file>"
  echo ""
  echo "OPTIONS:"
  echo "  -p, --parallel  Track memory across all child processes (for parallel R code)"
  echo "  -h, --help      Show this help message and exit"
  echo ""
  echo "ARGUMENTS:"
  echo "  R_script        Path to the R script to profile"
  echo "  output_log_file Path for memory log output (*.txt)"
  echo ""
  echo "EXAMPLES:"
  echo "  $0 simulate.R memory_usage.txt              # Single-process mode"
  echo "  $0 -p simulate.R memory_usage.txt           # Parallel mode"
  echo "  $0 --parallel calibrate.R /tmp/calib.txt    # Parallel mode (long flag)"
  echo ""
  echo "OUTPUT FILES:"
  echo "  - <output_log_file>     : Timestamped memory usage data"
  echo "  - <output_log_file>.png : Memory usage visualization"
  echo ""
  echo "PARALLEL MODE:"
  echo "  When -p/--parallel is enabled, memory is summed across the main R process"
  echo "  and ALL descendant processes. Use this when profiling code that uses:"
  echo "    - parallel::mclapply(), parallel::parLapply()"
  echo "    - foreach with doParallel/doFuture backends"
  echo "    - future::plan(multisession/multicore)"
  echo "    - Any fork-based or spawn-based parallelism"
}

# Parse command line options
PARALLEL_MODE=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    -p|--parallel)
      PARALLEL_MODE=true
      shift
      ;;
    -h|--help)
      show_help
      exit 0
      ;;
    -*)
      echo "ERROR: Unknown option: $1"
      echo "Use '$0 --help' for usage information"
      exit 1
      ;;
    *)
      break  # Stop parsing options, remaining args are positional
      ;;
  esac
done

# Check for required positional arguments
if [ "$#" -ne 2 ]; then
  echo "ERROR: Incorrect number of arguments"
  echo ""
  show_help
  exit 1
fi

# Parse positional arguments
SCRIPT="$1"
MEMLOG="$2"
PLOTFILE="${MEMLOG%.txt}.png"

# Validate input files and paths
if [ ! -f "$SCRIPT" ]; then
  echo "ERROR: R script '$SCRIPT' not found"
  exit 1
fi

if [ ! -d "$(dirname "$MEMLOG")" ]; then
  echo "ERROR: Directory '$(dirname "$MEMLOG")' does not exist"
  exit 1
fi

echo "==============================================================================="
echo "MEMORY PROFILER STARTING"
echo "==============================================================================="
echo "R Script:    $SCRIPT"
echo "Memory Log:  $MEMLOG"
echo "Plot File:   $PLOTFILE"
if [ "$PARALLEL_MODE" = true ]; then
  echo "Mode:        PARALLEL (tracking main process + all children)"
else
  echo "Mode:        SINGLE-PROCESS (main R process only)"
fi
echo "Start Time:  $(date)"
echo "==============================================================================="

# Record start time for duration calculation
START_TIME=$(date +%s.%N)

# Start R script in background and capture its PID
echo "Starting R script..."
Rscript "$SCRIPT" &
R_PID=$!

# Set up signal handler for graceful interruption
trap "kill $R_PID 2>/dev/null; echo ''; echo 'INTERRUPTED: Cleaning up...'; exit 1" INT TERM

# Initialize peak memory trackers
max_rss_kb=0  # Peak RSS in KB
max_vsz_kb=0  # Peak VSZ in KB

# Create memory log file with header
echo "# Memory Profile Log for: $SCRIPT" > "$MEMLOG"
echo "# Generated on: $(date)" >> "$MEMLOG"
if [ "$PARALLEL_MODE" = true ]; then
  echo "# Mode: PARALLEL (main process + all child processes)" >> "$MEMLOG"
else
  echo "# Mode: SINGLE-PROCESS (main R process only)" >> "$MEMLOG"
fi
echo "# Format: Timestamp RSS_GB VSZ_GB MEM_PCT NUM_PROCS" >> "$MEMLOG"
echo "# RSS_GB: Resident Set Size in Gigabytes" >> "$MEMLOG"
echo "# VSZ_GB: Virtual Memory Size in Gigabytes" >> "$MEMLOG"
echo "# MEM_PCT: Percentage of system memory used" >> "$MEMLOG"
echo "# NUM_PROCS: Number of processes being tracked" >> "$MEMLOG"

echo "Monitoring memory usage (sampling every 0.1 seconds)..."
sample_count=0

# Function to get all descendant PIDs of a process (recursive)
# Uses /proc filesystem for reliable child detection on Linux
get_all_descendants() {
  local parent_pid=$1
  local all_pids="$parent_pid"
  local children

  # Get direct children of this process
  children=$(ps --ppid "$parent_pid" -o pid= 2>/dev/null | tr -d ' ')

  for child in $children; do
    # Recursively get descendants of each child
    all_pids="$all_pids $(get_all_descendants $child)"
  done

  echo "$all_pids"
}

# Function to get aggregated memory for a list of PIDs
get_aggregated_memory() {
  local pids="$1"
  local total_rss=0
  local total_vsz=0
  local total_mem=0
  local proc_count=0

  for pid in $pids; do
    if ps -p "$pid" > /dev/null 2>&1; then
      mem_line=$(ps -o rss=,vsz=,%mem= -p "$pid" 2>/dev/null)
      if [ -n "$mem_line" ]; then
        rss=$(echo "$mem_line" | awk '{print $1}')
        vsz=$(echo "$mem_line" | awk '{print $2}')
        mem=$(echo "$mem_line" | awk '{print $3}')
        total_rss=$((total_rss + rss))
        total_vsz=$((total_vsz + vsz))
        total_mem=$(awk "BEGIN {printf \"%.1f\", $total_mem + $mem}")
        proc_count=$((proc_count + 1))
      fi
    fi
  done

  echo "$total_rss $total_vsz $total_mem $proc_count"
}

# Main monitoring loop - continues while R process is running
while ps -p $R_PID > /dev/null; do
  # Get current timestamp
  ts=$(date +"%H:%M:%S.%3N")

  if [ "$PARALLEL_MODE" = true ]; then
    # PARALLEL MODE: Get memory for main process + all descendants
    all_pids=$(get_all_descendants $R_PID)
    mem_result=$(get_aggregated_memory "$all_pids")

    rss_kb=$(echo "$mem_result" | awk '{print $1}')
    vsz_kb=$(echo "$mem_result" | awk '{print $2}')
    mem_pct=$(echo "$mem_result" | awk '{print $3}')
    num_procs=$(echo "$mem_result" | awk '{print $4}')
  else
    # SINGLE-PROCESS MODE: Get memory for main R process only
    mem_info=$(ps -o rss=,vsz=,%mem= -p $R_PID 2>/dev/null)

    rss_kb=$(echo "$mem_info" | awk '{print $1}')
    vsz_kb=$(echo "$mem_info" | awk '{print $2}')
    mem_pct=$(echo "$mem_info" | awk '{print $3}')
    num_procs=1
  fi

  # Handle empty values (process may have just exited)
  [ -z "$rss_kb" ] && rss_kb=0
  [ -z "$vsz_kb" ] && vsz_kb=0
  [ -z "$mem_pct" ] && mem_pct=0
  [ -z "$num_procs" ] && num_procs=0

  # Update peak memory trackers
  [ "$rss_kb" -gt "$max_rss_kb" ] 2>/dev/null && max_rss_kb=$rss_kb
  [ "$vsz_kb" -gt "$max_vsz_kb" ] 2>/dev/null && max_vsz_kb=$vsz_kb

  # Convert KB to GB for better readability
  rss_gb=$(awk "BEGIN {printf \"%.3f\", $rss_kb/1024/1024}")
  vsz_gb=$(awk "BEGIN {printf \"%.3f\", $vsz_kb/1024/1024}")

  # Write to log file (including process count)
  echo "$ts $rss_gb $vsz_gb $mem_pct $num_procs" >> "$MEMLOG"

  # Progress indicator (every 50 samples = 5 seconds)
  sample_count=$((sample_count + 1))
  if [ $((sample_count % 50)) -eq 0 ]; then
    if [ "$PARALLEL_MODE" = true ]; then
      echo "  Sample $sample_count: RSS=${rss_gb}GB, VSZ=${vsz_gb}GB, MEM%=${mem_pct}%, Procs=${num_procs}"
    else
      echo "  Sample $sample_count: RSS=${rss_gb}GB, VSZ=${vsz_gb}GB, MEM%=${mem_pct}%"
    fi
  fi

  # Wait before next sample
  sleep 0.1
done

# Wait for R script to complete and calculate execution time
wait $R_PID
R_EXIT_CODE=$?
END_TIME=$(date +%s.%N)
DURATION=$(awk "BEGIN {printf \"%.2f\", $END_TIME - $START_TIME}")

# Convert peak memory values to GB
MAX_RSS_GB=$(awk "BEGIN {printf \"%.3f\", $max_rss_kb/1024/1024}")
MAX_VSZ_GB=$(awk "BEGIN {printf \"%.3f\", $max_vsz_kb/1024/1024}")

# Print execution summary
echo ""
echo "==============================================================================="
echo "MEMORY PROFILING COMPLETED"
echo "==============================================================================="
echo "R Script:      $SCRIPT"
echo "Exit Code:     $R_EXIT_CODE"
echo "Duration:      ${DURATION}s"
echo "Total Samples: $sample_count"
if [ "$PARALLEL_MODE" = true ]; then
  echo "Mode:          PARALLEL (main process + all children)"
else
  echo "Mode:          SINGLE-PROCESS"
fi
echo ""
echo "PEAK MEMORY USAGE:"
echo "  Peak RSS:    ${MAX_RSS_GB} GB  (Physical memory)"
echo "  Peak VSZ:    ${MAX_VSZ_GB} GB  (Virtual memory)"
echo ""
echo "OUTPUT FILES:"
echo "  Memory Log:  $MEMLOG"
echo "  Plot File:   $PLOTFILE"
echo "==============================================================================="

# Generate memory usage visualization plot
echo ""
echo "Generating memory usage plot..."

# Create plot using R
Rscript - <<EOF
# Read memory log data
log <- read.table("$MEMLOG", header = FALSE, comment.char = "#")
colnames(log) <- c("Time", "RSS_GB", "VSZ_GB", "MEM_PCT", "NUM_PROCS")
log\$Index <- seq_len(nrow(log))

# Determine if parallel mode was used (more than 1 process at any point)
parallel_mode <- $( [ "$PARALLEL_MODE" = true ] && echo "TRUE" || echo "FALSE" )
max_procs <- max(log\$NUM_PROCS, na.rm = TRUE)

# Create PNG plot - taller if showing process count
plot_height <- ifelse(parallel_mode && max_procs > 1, 750, 600)
png("$PLOTFILE", width = 1000, height = plot_height, res = 100)

if (parallel_mode && max_procs > 1) {
  # Two-panel plot: memory on top, process count on bottom
  layout(matrix(c(1, 2), nrow = 2), heights = c(3, 1))
  par(mar = c(2, 4, 4, 2) + 0.1, cex.main = 1.2, cex.lab = 1.1)
} else {
  # Single panel plot
  par(mar = c(5, 4, 4, 2) + 0.1, cex.main = 1.2, cex.lab = 1.1)
}

# Create the main memory plot
plot(log\$Index, log\$RSS_GB, type = "l", col = "blue", lwd = 2,
     xlab = ifelse(parallel_mode && max_procs > 1, "", "Time (samples at 0.1 sec intervals)"),
     ylab = "Memory Usage (GB)",
     main = paste("Memory Profile:", basename("$SCRIPT"),
                  ifelse(parallel_mode, "(Parallel Mode)", "(Single Process)")),
     ylim = range(c(log\$RSS_GB, log\$VSZ_GB)))

# Add VSZ line
lines(log\$Index, log\$VSZ_GB, col = "red", lwd = 2)

# Add grid for better readability
grid(nx = NULL, ny = NULL, col = "gray90", lty = "dotted")

# Add legend with additional information
legend_text <- c(paste("RSS (Physical) - Peak:", sprintf("%.3f GB", max(log\$RSS_GB))),
                 paste("VSZ (Virtual) - Peak:", sprintf("%.3f GB", max(log\$VSZ_GB))))
if (parallel_mode) {
  legend_text <- c(legend_text, paste("Max Processes:", max_procs))
}
legend("topright",
       legend = legend_text,
       col = c("blue", "red", NA),
       lwd = c(2, 2, NA),
       cex = 0.9)

# Add process count panel if parallel mode with multiple processes
if (parallel_mode && max_procs > 1) {
  par(mar = c(5, 4, 1, 2) + 0.1)
  plot(log\$Index, log\$NUM_PROCS, type = "l", col = "darkgreen", lwd = 2,
       xlab = "Time (samples at 0.1 sec intervals)",
       ylab = "# Processes",
       ylim = c(0, max_procs + 1))
  grid(nx = NULL, ny = NULL, col = "gray90", lty = "dotted")
  abline(h = 1, col = "gray50", lty = "dashed")
}

# Add summary text at the bottom
mtext(paste("Duration:", "$DURATION", "seconds | Samples:", nrow(log)),
      side = 1, line = 4, cex = 0.8, col = "gray50", outer = FALSE)

dev.off()

# Print plot statistics
cat("Plot generated successfully:\\n")
cat("  File:", "$PLOTFILE", "\\n")
cat("  Dimensions: 1000x", plot_height, " pixels\\n", sep = "")
cat("  Data points:", nrow(log), "\\n")
if (parallel_mode) {
  cat("  Max concurrent processes:", max_procs, "\\n")
}
EOF

# Optional: Display plot file information
if [ -f "$PLOTFILE" ]; then
  echo ""
  echo "PLOT GENERATION SUCCESSFUL"
  echo "Plot saved to: $PLOTFILE"
  echo "File size: $(du -h "$PLOTFILE" | cut -f1)"
else
  echo ""
  echo "WARNING: Plot generation may have failed"
  echo "Check R installation and plotting capabilities"
fi

echo ""
echo "==============================================================================="
echo "MEMORY PROFILING SESSION COMPLETE"
echo "==============================================================================="

# Exit with the same code as the R script
exit $R_EXIT_CODE

# Optional plot viewing (commented out as it doesn't work over SSH)
# Uncomment these lines if running locally with GUI support:
# if command -v xdg-open &> /dev/null; then
#   echo "Opening plot with default viewer..."
#   xdg-open "$PLOTFILE" &> /dev/null &
# elif command -v open &> /dev/null; then  # for macOS
#   echo "Opening plot with default viewer..."
#   open "$PLOTFILE" &> /dev/null &
# fi