source("./global.R")
e <- read_fst("./inputs/mortality/mrtl_clb.fst", as.data.table = TRUE)
e[, mrtl_clbr := 1] # Cancel previous calibrarion
write_fst(e, "./inputs/mortality/mrtl_clb.fst")

IMPACTncd <- Simulation$new("./inputs/mortality/sim_design_mrtl_clbr.yaml")

scenario_fn <- function(sp) NULL

IMPACTncd$
  # del_logs()$
  # del_outputs()$
  run(1:200, multicore = TRUE, "sc0")$
  export_summaries(multicore = TRUE, type = "allcause_mrtl")



lifetable_all <- fread("./inputs/mortality/lt.csv")
lifetable_all[, `:=` (year = year - 2000L)]
# mx to qx from https://rdrr.io/cran/MortalityLaws/src/R/LifeTable.R
# for (j in grep("mx_", names(lifetable_all), value = TRUE)) {
#   nam <- gsub("^mx_", "qx_", j)
#   lifetable_all[agegrp != "<1", (nam) :=  1 - exp(-5 * get(j))]
#   lifetable_all[agegrp == "<1", (nam) :=  1 - exp(-1 * get(j))]
# }


# get mortality (qx) results from a previous run of the model
tt <- fread("/mnt/storage_fast/output/hf_real_mrtl_clbr/summaries/mrtl_scaled_up.csv.gz"
)[scenario == "sc0"][, `:=` (dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]
e <- tt[, .(mrtl_rate = sum(all_cause_mrtl)/sum(popsize)), keyby = .(year, agegrp, sex, dimd)]
# outstrata <- c("mc", "year", "sex", "agegrp", "dimd", "scenario")
# e <- tt[, lapply(.SD, sum), .SDcols = patterns("_mrtl$|^popsize$"), keyby = eval(outstrata)
# ][, lapply(.SD, function(x) x/popsize), keyby = eval(outstrata)]
# e[, popsize := NULL]
# e <- e[, as.list(CKutils::fquantile(all_cause_mrtl, 0.5)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(e, c(setdiff(outstrata, "mc"), "mrtl_rate"))
# setkeyv(e, setdiff(outstrata, "mc"))

e[lifetable_all, on = c("year", "agegrp", "sex", "dimd"), mrtl_clbr := mx_total/mrtl_rate]
e[, summary(mrtl_clbr)]
e[is.infinite(mrtl_clbr), ]
e[is.infinite(mrtl_clbr), mrtl_clbr := 1]

tt <- data.table(age = 30:99, agegrp = CKutils::agegrp_name(30, 99, match_input = TRUE, match_input_max_age = 99))
e <- tt[e, on = .NATURAL, allow.cartesian = TRUE]
e[, c("mrtl_rate", "agegrp") := NULL]
e[, sex := factor(sex, c("men", "women"))]

# calibration to account for changes in prevalence because of mortality corrections in previous years
e[, mrtl_clbr := 0.95 * mrtl_clbr]
e[dimd == "10 least deprived", mrtl_clbr := (0.993 ^ (year - 13)) * mrtl_clbr]
e[dimd == "10 least deprived", mrtl_clbr := 0.997 * mrtl_clbr]
e[dimd == "10 least deprived" & between(age, 90, 99), mrtl_clbr := (0.996 ^ (year - 13)) * mrtl_clbr]

e[dimd == "9", mrtl_clbr := (0.998 ^ (year - 13)) * mrtl_clbr]

e[dimd == "8", mrtl_clbr := (0.9935 ^ (year - 13)) * mrtl_clbr]
e[dimd == "8" & between(age, 40, 69), mrtl_clbr := (1.002 ^ (year - 13)) * mrtl_clbr]

e[dimd == "7", mrtl_clbr := (0.9965 ^ (year - 13)) * mrtl_clbr]
e[dimd == "7" & between(age, 50, 59), mrtl_clbr := (1.0002 ^ (year - 13)) * mrtl_clbr]

e[dimd == "6", mrtl_clbr := (0.998 ^ (year - 13)) * mrtl_clbr]
e[dimd == "6" & between(age, 90, 99), mrtl_clbr := (1.0004 ^ (year - 13)) * mrtl_clbr]

e[dimd == "5", mrtl_clbr := (0.9945 ^ (year - 13)) * mrtl_clbr]

e[dimd == "4", mrtl_clbr := (0.9985 ^ (year - 13)) * mrtl_clbr]
e[dimd == "4" & between(age, 60, 69), mrtl_clbr := (0.9999 ^ (year - 13)) * mrtl_clbr]

e[dimd == "3", mrtl_clbr := (0.9965 ^ (year - 13)) * mrtl_clbr]
e[dimd == "3" & between(age, 30, 39), mrtl_clbr := (0.998 ^ (year - 13)) * mrtl_clbr]

e[dimd == "2", mrtl_clbr := (0.997 ^ (year - 13)) * mrtl_clbr]
e[dimd == "2" & between(age, 30, 39), mrtl_clbr := (0.998 ^ (year - 13)) * mrtl_clbr]

e[dimd == "1 most deprived", mrtl_clbr := (0.999 ^ (year - 13)) * mrtl_clbr]
e[dimd == "1 most deprived" & between(age, 30, 39), mrtl_clbr := (0.993 ^ (year - 13)) * mrtl_clbr]
e[dimd == "1 most deprived" & between(age, 50, 69), mrtl_clbr := 1.05 * mrtl_clbr]
e[dimd == "1 most deprived" & between(age, 80, 89), mrtl_clbr := 0.94 * mrtl_clbr]
e[dimd == "1 most deprived" & between(age, 90, 99), mrtl_clbr := 0.96 * mrtl_clbr]


setkeyv(e, c("year", "age", "dimd", "sex"))
IMPACTncdEngl::is_valid_lookup_tbl(e, c("year", "age", "sex", "dimd"))
write_fst(e, "./inputs/mortality/mrtl_clb.fst")



IMPACTncd <- Simulation$new("./inputs/sim_design.yaml")

scenario_fn <- function(sp) NULL

IMPACTncd$
  del_logs()$
  del_outputs()$
  run(1:200, multicore = TRUE, "sc0")$
  export_summaries(multicore = TRUE)

