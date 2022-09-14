onsmrtl <- read_fst("./inputs/mortality/mrtl_by_dimd.fst", as.data.table = TRUE)
onsmrtl[, `:=` (year = year - 2000L, mx = deaths/popsize)]
# mx to qx from https://rdrr.io/cran/MortalityLaws/src/R/LifeTable.R
onsmrtl[agegroup != "90+", qx :=  1 - exp(-5 * mx)]
onsmrtl[agegroup == "90+", qx :=  1 - exp(-10 * mx)]
setnames(onsmrtl, "agegroup", "agegrp")

# get mortality (qx) results from a previous run of the model
tt <- fread("/mnt/storage_fast/output/hf_real/summaries/mrtl_scaled_up.csv.gz"
)[scenario == "sc0" & year < 17,][, `:=` (sex = factor(sex, c("men", "women")),
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]
tt[agegrp %in% c("90-94", "95-99"), agegrp := "90+"]
tt <- tt[, .(qx_mdl = sum(all_cause_mrtl)/sum(popsize)),
         keyby = .(agegrp, year, sex, dimd)]

absorb_dt(tt, onsmrtl)

tt[, mrtl_clbr := mx/qx_mdl]

e <- data.table(age = 30:100, agegrp = agegrp_name(30, 90, match_input = TRUE, match_input_max_age = 100))
e <- e[tt, on = .NATURAL, allow.cartesian = TRUE]
e <- e[year > 13, .(mrtl_clbr = mean(mrtl_clbr)), keyby = .(age, sex, dimd)]
is_valid_lookup_tbl(e, c("age", "sex", "dimd"))
write_fst(e[, .(age, sex, dimd, mrtl_clbr)], "./inputs/mortality/mrtl_clb.fst")
