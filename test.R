source("./global.R")
IMPACTncd <- Simulation$new("testing/sim_design_testing.yaml")
IMPACTncd$
  del_logs()$
  del_outputs()$
  # del_synthpops()$
  # del_parfs()$
run(1:2, multicore = TRUE, "sc0")

# library(data.table)
# tbl <- read_parquet_dt(design$exposures$smok_dur_ex$file_path)
# tbl[, smok_status := as.integer(as.character(smok_status))]
# arrow::write_dataset(tbl, design$exposures$smok_dur_ex$file_path, partitioning = "year")
# kc <- sort(setdiff(
#   names(tbl),
#   c("mu", "sigma", "nu", "tau", "maxq", "minq", design$exposures$veg$thresholds)
# ))
# kc <- kc[order(match(kc, "year"))]
# setcolorder(tbl, kc)
# setkeyv(tbl, kc)
# pop <- tbl[, .(year, age, ethnicity, qimd, sex, sha, rankstat_smok_dur_curr = runif(.N))]
# design$exposures$smok_dur_curr$generate(pop, design, 1:10)

# set.seed(42)
# row <- 54000
# pop <- tbl[1:7e5, .(year, age, ethnicity, qimd, sex, sha, rank_veg = runif(.N))]
# while (row <= 7e5) {
#   if ((row %% 1e3) == 0) {
#     print(row)
#   }
#   t0 <- proc.time()[["elapsed"]]
#   out <- fqDEL(
#     tbl[row, minq] + pop[row, rank_veg] * (tbl[row, maxq] - tbl[row, minq]),
#     tbl[row, mu],
#     tbl[row, sigma],
#     tbl[row, nu]
#   )
#   t1 <- proc.time()[["elapsed"]] - t0
#   if (out > 100) {
#     print(paste0(row, " : ", out, " : ", t1))
#     stop("inf produced")
#   }
#   if (t1 > 0.5) {
#     print(paste0(row, " : ", t1))
#   }
#   row <- row + 1
# }

# row <- 56897
#  system.time(fqDEL(
#     0.99999999,
#     1.987479,
#     0.0001356873,
#     0.8305514
#   ))


# summary(tbl[630e2:640e2])
# tbl[630e2]
# fqDEL(1e-9, tbl[630e2, mu], tbl[630e2, sigma], tbl[630e2, nu])

# design$exposures$smok_cig_ex$add_qmin_qmax(
#   read_parquet_dt(design$exposures$smok_cig_ex$file_path),
#   force = FALSE,
#   write_to_disk = TRUE,
# )
# tt <- read_parquet_dt(design$exposures$smok_cig_ex$file_path)
# names(tt)
# system.time(design$exposures$smok_cig_ex$add_qmin_qmax(
#   tt,
#   force = FALSE,
#   write_to_disk = FALSE,
# ))


