source("./global.R")
design <- Design$new("./inputs/sim_design.yaml")
# RR ----
# Create a named list of Exposure objects for the files in ./inputs/RR
fl <- list.files(path = "./inputs/RR", pattern = ".csvy$", full.names = TRUE)
# RR <- lapply(fl, Exposure$new, design)
# names(RR) <- sapply(RR, function(x) x$get_name())
# lapply(RR, function(x) {
#     x$gen_stochastic_effect(design, overwrite = FALSE, smooth = FALSE)
# })
RR <- future_lapply(fl, Exposure$new, design,future.seed = 950480304L)
names(RR) <- sapply(RR, function(x) x$get_name())
invisible(future_lapply(RR, function(x) {
    x$gen_stochastic_effect(design, overwrite = FALSE, smooth = FALSE)
},
future.seed = 627524136L))
# NOTE smooth cannot be exported to Design for now, because the first time
# this parameter changes we need logic to overwrite unsmoothed files
rm(fl)
#
# Generate diseases ----
diseases <- lapply(design$sim_prm$diseases, function(x) {
    print(x$name)
    x[["design_"]] <- design
    x[["RR"]] <- RR
    do.call(Disease$new, x)
})
names(diseases) <- sapply(design$sim_prm$diseases, `[[`, "name")

mk_scenario_init2 <- function(scenario_name, diseases_, sp, design_) {
    if (nzchar(scenario_name)) { # TODO get suffix from design
        scenario_suffix_for_pop <- paste0("_", scenario_name)
    } else {
        scenario_suffix_for_pop <- scenario_name
    }
    list(
        "exposures"          = design_$sim_prm$exposures,
        "scenarios"          = design_$sim_prm$scenarios, # to be generated programmatically
        "scenario"           = scenario_name,
        "kismet"             = design_$sim_prm$kismet, # If TRUE random numbers are the same for each scenario.
        "init_year"          = design_$sim_prm$init_year,
        "pids"               = "pid",
        "years"              = "year",
        "ages"               = "age",
        "sexs"               = "sex",
        "dimds"              = "dimd",
        "ageL"               = design_$sim_prm$ageL,
        "all_cause_mrtl"     = paste0("all_cause_mrtl", scenario_suffix_for_pop),
        "cms_score"          = paste0("cms_score", scenario_suffix_for_pop),
        "cms_count"          = paste0("cms_count", scenario_suffix_for_pop),
        "strata_for_outputs" = c("pid", "year", "age", "sex", "dimd"),
        "diseases"           = lapply(diseases_, function(x)
            x$to_cpp(sp, design_, scenario_suffix_for_pop))
    )
}

# sim <- SynthPop$new(0L, design)
# sim$write_synthpop(1:500)
# sim$delete_synthpop(NULL)
# ll <- sim$gen_synthpop_demog(design)
sp <- SynthPop$new(1L, design)

# diseases$ckd$gen_parf(sp, design)

nn <- "t2dm"
self <- diseases[[nn]]$.__enclos_env__$self
private <- diseases[[nn]]$.__enclos_env__$private
design_ <- design
diseases_ <- diseases
check = design_$sim_prm$logs
keep_intermediate_file = TRUE
mc_iter = mc_ = 1

# self    <- diseases[[nn]]$get_rr()[[1]]$.__enclos_env__$self
# private <-  diseases[[nn]]$get_rr()[[1]]$.__enclos_env__$private


# lapply(diseases, function(x) x$harmonise_epi_tables(sp))
lapply(diseases, function(x) {
    print(x)
    x$gen_parf_files(design)
})
lapply(diseases, function(x) {
    print(x)
    x$gen_parf(sp, design, diseases)
})
lapply(diseases, function(x) {
    print(x)
    x$set_init_prvl(sp, design)
})
lapply(diseases, function(x) {
    print(x)
    x$set_rr(sp, design)
})
lapply(diseases, function(x) {
    print(x)
    x$set_incd_prb(sp, design)
})
lapply(diseases, function(x) {
    print(x)
    x$set_dgns_prb(sp, design)
})
lapply(diseases, function(x) {
    print(x)
    x$set_mrtl_prb(sp, design)
})

# diseases$t2dm$harmonise_epi_tables(sp)
# diseases$t2dm$gen_parf(sp, design)
# diseases$t2dm$set_init_prvl(sp, design)
# diseases$t2dm$set_rr(sp, design)
# diseases$t2dm$set_incd_prb(sp, design)
# diseases$t2dm$set_dgns_prb(sp, design)
# diseases$t2dm$set_mrtl_prb(sp, design)


# diseases$af$harmonise_epi_tables(sp)
# diseases$af$gen_parf(sp, design)
# diseases$af$set_init_prvl(sp, design)
# diseases$af$set_rr(sp, design)
# diseases$af$set_incd_prb(sp, design)
# diseases$af$set_dgns_prb(sp, design)
# diseases$af$set_mrtl_prb(sp, design)

# diseases$asthma$harmonise_epi_tables(sp)
# diseases$asthma$gen_parf(sp, design)
# diseases$asthma$set_init_prvl(sp, design)
# diseases$asthma$set_rr(sp, design)
# diseases$asthma$set_incd_prb(sp, design)
# diseases$asthma$set_dgns_prb(sp, design)
# diseases$asthma$set_mrtl_prb(sp, design)

# self <- diseases$asthma$.__enclos_env__$self
# private <- diseases$asthma$.__enclos_env__$private

# diseases$copd$harmonise_epi_tables(sp)
# diseases$copd$gen_parf(sp, design)
# diseases$copd$set_init_prvl(sp, design)
# diseases$copd$set_rr(sp, design)
# diseases$copd$set_incd_prb(sp, design)
# diseases$copd$set_dgns_prb(sp, design)
# diseases$copd$set_mrtl_prb(sp, design)

# diseases$htn$harmonise_epi_tables(sp)
# diseases$htn$gen_parf(sp, design)
# diseases$htn$set_init_prvl(sp, design)
# diseases$htn$set_rr(sp, design)
# diseases$htn$set_incd_prb(sp, design)
# diseases$htn$set_dgns_prb(sp, design)
# diseases$htn$set_mrtl_prb(sp, design)

# diseases$obesity$harmonise_epi_tables(sp)
# diseases$obesity$gen_parf(sp, design)
# diseases$obesity$set_init_prvl(sp, design)
# diseases$obesity$set_rr(sp, design)
# diseases$obesity$set_incd_prb(sp, design)
# diseases$obesity$set_dgns_prb(sp, design)
# diseases$obesity$set_mrtl_prb(sp, design)

# diseases$dementia$harmonise_epi_tables(sp)
# diseases$dementia$gen_parf(sp, design)
# diseases$dementia$set_init_prvl(sp, design)
# diseases$dementia$set_rr(sp, design)
# diseases$dementia$set_incd_prb(sp, design)
# diseases$dementia$set_dgns_prb(sp, design)
# diseases$dementia$set_mrtl_prb(sp, design)

# diseases$chd$harmonise_epi_tables(sp)
# diseases$chd$gen_parf(sp, design)
# diseases$chd$set_init_prvl(sp, design)
# diseases$chd$set_rr(sp, design)
# diseases$chd$set_incd_prb(sp, design)
# diseases$chd$set_dgns_prb(sp, design)
# diseases$chd$set_mrtl_prb(sp, design)

# diseases$stroke$harmonise_epi_tables(sp)
# diseases$stroke$gen_parf(sp, design)
# diseases$stroke$set_init_prvl(sp, design)
# diseases$stroke$set_rr(sp, design)
# diseases$stroke$set_incd_prb(sp, design)
# diseases$stroke$set_dgns_prb(sp, design)
# diseases$stroke$set_mrtl_prb(sp, design)

# diseases$ckd$harmonise_epi_tables(sp) # bmi_rr already present in the data.
# diseases$ckd$gen_parf(sp, design)
# diseases$ckd$set_init_prvl(sp, design)
# diseases$ckd$set_rr(sp, design)
# diseases$ckd$set_incd_prb(sp, design)
# diseases$ckd$set_dgns_prb(sp, design)
# diseases$ckd$set_mrtl_prb(sp, design)

# diseases$lung_ca$harmonise_epi_tables(sp)
# diseases$lung_ca$gen_parf(sp, design)
# diseases$lung_ca$set_init_prvl(sp, design)
# diseases$lung_ca$set_rr(sp, design)
# diseases$lung_ca$set_incd_prb(sp, design)
# diseases$lung_ca$set_dgns_prb(sp, design)
# diseases$lung_ca$set_mrtl_prb(sp, design)

# diseases$colorect_ca$harmonise_epi_tables(sp)
# diseases$colorect_ca$gen_parf(sp, design)
# diseases$colorect_ca$set_init_prvl(sp, design)
# diseases$colorect_ca$set_rr(sp, design)
# diseases$colorect_ca$set_incd_prb(sp, design)
# diseases$colorect_ca$set_dgns_prb(sp, design)
# diseases$colorect_ca$set_mrtl_prb(sp, design)

# diseases$prostate_ca$harmonise_epi_tables(sp)
# diseases$prostate_ca$gen_parf(sp, design)
# diseases$prostate_ca$set_init_prvl(sp, design)
# diseases$prostate_ca$set_rr(sp, design)
# diseases$prostate_ca$set_incd_prb(sp, design)
# diseases$prostate_ca$set_dgns_prb(sp, design)
# diseases$prostate_ca$set_mrtl_prb(sp, design)

# diseases$breast_ca$harmonise_epi_tables(sp)
# diseases$breast_ca$gen_parf(sp, design)
# diseases$breast_ca$set_init_prvl(sp, design)
# diseases$breast_ca$set_rr(sp, design)
# diseases$breast_ca$set_incd_prb(sp, design)
# diseases$breast_ca$set_dgns_prb(sp, design)
# diseases$breast_ca$set_mrtl_prb(sp, design)

# diseases$andep$harmonise_epi_tables(sp)
# diseases$andep$gen_parf(sp, design)
# diseases$andep$set_init_prvl(sp, design)
# diseases$andep$set_rr(sp, design)
# diseases$andep$set_incd_prb(sp, design)
# diseases$andep$set_dgns_prb(sp, design)
# diseases$andep$set_mrtl_prb(sp, design)

# diseases$other_ca$harmonise_epi_tables(sp)
# diseases$other_ca$gen_parf(sp, design)
# diseases$other_ca$set_init_prvl(sp, design)
# diseases$other_ca$set_rr(sp, design)
# diseases$other_ca$set_incd_prb(sp, design)
# diseases$other_ca$set_dgns_prb(sp, design)
# diseases$other_ca$set_mrtl_prb(sp, design)

# diseases$hf$harmonise_epi_tables(sp)
# diseases$hf$gen_parf(sp, design)
# diseases$hf$set_init_prvl(sp, design)
# diseases$hf$set_rr(sp, design)
# diseases$hf$set_incd_prb(sp, design)
# diseases$hf$set_dgns_prb(sp, design)
# diseases$hf$set_mrtl_prb(sp, design)

# diseases$nonmodelled$harmonise_epi_tables(sp)
# diseases$nonmodelled$gen_parf(sp, design)
# diseases$nonmodelled$set_init_prvl(sp, design)
# diseases$nonmodelled$set_rr(sp, design)
# diseases$nonmodelled$set_incd_prb(sp, design)
# diseases$nonmodelled$set_dgns_prb(sp, design)
# diseases$nonmodelled$set_mrtl_prb(sp, design)



# lapply(diseases, function(x) {
#     print(x$name)
#     x$gen_parf(sp, design)$
#     set_init_prvl(sp, design)$
#     set_rr(sp, design)$
#     set_incd_prb(sp, design)$
#     set_dgns_prb(sp, design)$
#     set_mrtl_prb(sp, design)
# })

transpose(sp$pop[, lapply(.SD, anyNA)], keep.names = "rn")[(V1)]



# qsave(sp, "./simulation/tmp.qs", nthreads = 4)
# sp <- qread("./simulation/tmp.qs"); setDT(sp$pop)
l <- mk_scenario_init2("", diseases, sp, design)
simcpp(sp$pop, l, sp$mc)

lapply(diseases, function(x) {
    print(x)
    x$calibrate_incd_prb(sp, design)
})

simcpp(sp$pop, l, sp$mc)


sp$pop[year >= 13, sum(asthma_prvl > 0), keyby = year][, plot(year, V1)]
sp$pop[year >= 13, sum(asthma_prvl == 1), keyby = year][, plot(year, V1)]

sp$pop[year >= 13, sum(andep_prvl > 0), keyby = year][, plot(year, V1)]
sp$pop[year >= 13, sum(andep_prvl == 1), keyby = year][, plot(year, V1)]


sp$pop[year >= 13, sum(pain_prvl > 0), keyby = year][, plot(year, V1)]
sp$pop[year >= 13, sum(pain_prvl == 1), keyby = year][, plot(year, V1)]

sp$pop[year >= 13, sum(constipation_prvl > 0), keyby = year][, plot(year, V1)]
sp$pop[year >= 13, sum(constipation_prvl == 1), keyby = year][, plot(year, V1)]


sp$pop[year >= 13, sum(chd_prvl > 0), keyby = year][, plot(year, V1)]
sp$pop[year >= 13, sum(lung_ca_prvl > 0), keyby = year][, plot(year, V1)]

sp$pop[asthma_prvl > 0 & year == 13, mean(asthma_prvl)] # 4
sp$pop[andep_prvl > 0 & year == 13, mean(andep_prvl)] # 3
sp$pop[pain_prvl > 0 & year == 13, mean(pain_prvl)] # 7
sp$pop[constipation_prvl > 0 & year == 13, mean(constipation_prvl)]  #5


sp$pop[asthma_dgns != asthma_prvl, .N]
sp$pop[asthma_prvl > 0, mean(asthma_prvl)]
sp$pop[year >= 13, sum(asthma_prvl > 0), keyby = age][, plot(age, V1)]
sp$pop[year >= 13, sum(asthma_prvl > 0), keyby = year][, plot(year, V1)]
sp$pop[asthma_prvl > 0, hist(asthma_prvl)]
sp$pop[asthma_prvl > 0 & year == 13, table(asthma_prvl)]



sp$pop[asthma_prvl > 0, mean(asthma_prvl), keyby = year][, plot(year, V1)]
sp$pop[asthma_prvl > 0, max(asthma_prvl), keyby = year][, plot(year, V1)]

id <- sp$pop[asthma_prvl > 0, unique(pid)]
sp$pop[pid %in% id, sum(asthma_prvl > 0), by = pid][, table(V1)]
View(sp$pop[pid %in% id, .(as.character(pid), year, asthma_prvl)])

id <- sp$pop[constipation_prvl == 1 & year == 14L, unique(pid)]
sp$pop[pid %in% id, sum(constipation_prvl > 0), by = pid][, table(V1)]
View(sp$pop[pid %in% id, .(as.character(pid), year, constipation_prvl)])


disnam <- paste0(names(diseases), "_prvl")
disnam <- disnam[disnam != "nonmodelled_prvl"]
cmswt <- c(
    "htn_prvl"      = 0.08,
    "andep_prvl"    = 0.5,
    "pain_prvl"     = 0.92,
    "helo_prvl"     = 0.09, # hearing loss
    "ibs_prvl"      = 0.21,
    "asthma_prvl"   = 0.19,
    "dm_prvl"       = 0.75, # t2dm + t1dm
    "chd_prvl"      = 0.49,
    "ckd45_prvl"    = 0.53, # TODO change to ckd45
    "af_prvl"       = 1.34,
    "constip_prvl"  = 1.12,
    "stroke_prvl"   = 0.80,
    "copd_prvl"     = 1.46,
    "ctdra_prvl"    = 0.43, # connective tissue disorders + rheumatoid arthritis
    "cancer_prvl"   = 1.53,
    "alcpr_prvl"    = 0.65,
    "hf_prvl"       = 1.18,
    "dementia_prvl" = 2.50,
    "psychos_prvl"  = 0.64,
    "epilepsy_prvl" = 0.92
)

disnam <- c(grep("_ca_prvl$|^ctd_prvl$|^ra_prvl$|^t1dm_prvl$|^t2dm_prvl$",
                 disnam, value = TRUE, invert = TRUE))
setdiff(names(cmswt), disnam) # ideally should be empty. If not there are diseases in CMS that we do not model
for (i in setdiff(disnam, names(cmswt))) cmswt[[i]] <- 0 # Fill cmswt with 0 for conditions that are not in CMS

hlpfn <- function(disprvl, dt, cmswt = cmswt) clamp(dt[[disprvl]]) * cmswt[[disprvl]]
sp$pop[, cms_scoreR := Reduce(`+`, lapply(disnam, hlpfn, sp$pop, cmswt))]
sp$pop[year >= 13 & age >= 30, table(cms_score == cms_scoreR)]
tt <- sp$pop[year >= 13 & cms_score != cms_scoreR, unique(pid)]
View(sp$pop[pid %in% tt, .(as.character(pid), year, age, all_cause_mrtl,
                           cms_scoreR, cms_score, cms_count,
                       asthma_prvl, asthma_dgns,
                       ibs_prvl, ibs_dgns,
                       helo_prvl, helo_dgns)])
View(sp$pop[pid %in% tt, .SD, .SDcols = patterns("^pid$|^cms_|_mrtl$|_prvl$|_dgns$")])











sp$update_pop_weights()
sp$pop[, mc := sp$mc_aggr]

# export xps
dt <- copy(sp$pop)
mc_ <- sp$mc_aggr
export_xps <- function(mc_,
                       dt,
                       write_to_disk = TRUE,
                       filenam = "val_xps_output.csv") {
    to_agegrp(dt, 20L, 99L, "age", "agegrp20", min_age = 30, to_factor = TRUE) # TODO link max age to design

    dt[, smok_never_curr_xps := fifelse(smok_status_curr_xps == "1", 1L, 0L)]
    dt[, smok_active_curr_xps := fifelse(smok_status_curr_xps == "4", 1L, 0L)]

    xps <- grep("_curr_xps$", names(dt), value = TRUE)
    xps <- xps[-which(xps %in% c("smok_status_curr_xps", "met_curr_xps",
                                 "bpmed_curr_xps", "t2dm_prvl_curr_xps",
                                 "af_prvl_curr_xps"))]
    out_xps <- groupingsets(
        dt[all_cause_mrtl >= 0L & year >= 13, ], # TODO link to design
        j = lapply(.SD, weighted.mean, wt),
        by = c("year", "sex", "agegrp20", "qimd", "ethnicity", "sha"),
        .SDcols = xps,
        sets = list(
            c("year", "sex", "agegrp20", "qimd"),
            c("year", "sex"),
            c("year", "agegrp20"),
            c("year", "qimd"),
            c("year", "ethnicity"),
            c("year", "sha")
        )
    )[, `:=` (year = year + 2000L, mc = mc_)]
    for (j in seq_len(ncol(out_xps)))
        set(out_xps, which(is.na(out_xps[[j]])), j, "All")
    dt[, c(
        "agegrp20",
        "smok_never_curr_xps",
        "smok_active_curr_xps"
    ) := NULL]

    setkey(out_xps, year)

    fwrite_safe(out_xps, output_dir(filenam))

    invisible(out_xps)
}






nam <- c("mc", "pid", "year", "sex", "dimd", "ethnicity", "sha", grep("_prvl$|_mrtl$", names(sp$pop), value = TRUE))
fwrite_safe(sp$pop[all_cause_mrtl >= 0L, ..nam],
            file.path(design$sim_prm$output_dir, "lifecourse", paste0(sp$mc_aggr, "_lifecourse.csv")))

parf <- fread("/mnt/storage_fast/output/hf_real/parf/parf.csv")

fl <- list.files("/mnt/storage_fast/output/hf_real/lifecourse/",
                 "_lifecourse.csv$", full.names = TRUE)

out <- rbindlist(lapply(fl, fread))

sp$pop[!is.na(all_cause_mrtl), median(bmi_curr_xps), keyby = year][, plot(year, V1)]
sp$pop[!is.na(all_cause_mrtl), median(sbp_curr_xps), keyby = year][, plot(year, V1)]
sp$pop[!is.na(all_cause_mrtl), weighted.mean(smok_status_curr_xps == "4", wt, na.rm), keyby = year][, plot(year, V1)]
sp$pop[!is.na(all_cause_mrtl), mean(smok_cig_curr_xps), keyby = year][, plot(year, V1)]


sp$pop[!is.na(all_cause_mrtl), sum(chd_prvl == 1)/.N, keyby = year][, plot(year, V1)]
sp$pop[!is.na(all_cause_mrtl), sum(stroke_prvl == 1)/.N, keyby = year][, plot(year, V1)]
sp$pop[!is.na(all_cause_mrtl), sum(af_prvl == 1)/.N, keyby = year][, plot(year, V1)]
sp$pop[!is.na(all_cause_mrtl), sum(t2dm_prvl == 1)/.N, keyby = year][, plot(year, V1)]
sp$pop[!is.na(all_cause_mrtl), sum(obesity_prvl == 1)/.N, keyby = year][, plot(year, V1)]
sp$pop[!is.na(all_cause_mrtl), sum(htn_prvl == 1)/.N, keyby = year][, plot(year, V1)]
sp$pop[!is.na(all_cause_mrtl), sum(copd_prvl == 1)/.N, keyby = year][, plot(year, V1)]
sp$pop[!is.na(all_cause_mrtl), sum(lung_ca_prvl == 1)/.N, keyby = year][, plot(year, V1)]
sp$pop[!is.na(all_cause_mrtl), sum(breast_ca_prvl == 1)/.N, keyby = year][, plot(year, V1)]
sp$pop[!is.na(all_cause_mrtl), sum(colorect_ca_prvl == 1)/.N, keyby = year][, plot(year, V1)]
sp$pop[!is.na(all_cause_mrtl), sum(prostate_ca_prvl == 1)/.N, keyby = year][, plot(year, V1)]
sp$pop[!is.na(all_cause_mrtl), sum(hf_prvl == 1)/.N, keyby = year][, plot(year, V1)]
sp$pop[!is.na(all_cause_mrtl), sum(andep_prvl == 1)/.N, keyby = year][, plot(year, V1)]
sp$pop[!is.na(all_cause_mrtl), sum(other_ca_prvl == 1)/.N, keyby = year][, plot(year, V1)]

sp$pop[, sum(all_cause_mrtl > 0, na.rm = T)/.N, keyby = year][, plot(year, V1)]
sp$pop[, sum(all_cause_mrtl == 12, na.rm = T)/.N, keyby = year][, plot(year, V1)]

sp$pop[, table(all_cause_mrtl, useNA = "a")]
sp$pop[year == 13, table(chd_prvl)]

sp$pop[!is.na(all_cause_mrtl), sum(chd_prvl > 0)/.N, keyby = year][, plot(year, V1)]
sp$pop[!is.na(all_cause_mrtl), sum(stroke_prvl > 0)/.N, keyby = year][, plot(year, V1)]
sp$pop[!is.na(all_cause_mrtl), sum(af_prvl > 0)/.N, keyby = year][, plot(year, V1)]
sp$pop[!is.na(all_cause_mrtl), sum(t2dm_prvl > 0)/.N, keyby = year][, plot(year, V1)]
sp$pop[!is.na(all_cause_mrtl), sum(obesity_prvl > 0)/.N, keyby = year][, plot(year, V1)]
sp$pop[!is.na(all_cause_mrtl), sum(htn_prvl > 0)/.N, keyby = year][, plot(year, V1)]
sp$pop[!is.na(all_cause_mrtl), sum(copd_prvl > 0)/.N, keyby = year][, plot(year, V1)]
sp$pop[!is.na(all_cause_mrtl), sum(lung_ca_prvl > 0)/.N, keyby = year][, plot(year, V1)]
sp$pop[!is.na(all_cause_mrtl), sum(hf_prvl > 0)/.N, keyby = year][, plot(year, V1)]
sp$pop[!is.na(all_cause_mrtl), sum(andep_prvl > 0)/.N, keyby = year][, plot(year, V1)]


sp$pop[between(age, 60, 64), sum(af_prvl > 0)/.N, keyby = year][, plot(year, V1)]

sp$pop[, sum(sbp_curr_xps > 140) / .N, keyby = year]


fwrite_safe(sp$pop[1:10], "/mnt/storage_fast/output/hf_real/lifecourse/test.csv")
fwrite_safe(sp$pop[11:20], "/mnt/storage_fast/output/hf_real/lifecourse/test.csv")

tt <- list.files("/mnt/storage_fast/output/hf_real/lifecourse/", full.names = T)
lapply(tt, function(x) {
    print(x)
    fread(x)
})
f <- fread(tt[[3]])
p <- read_fst("./inputs/disease_burden/af_prvl.fst", as.data.table = T)

absorb_dt(f, p)
f[year == 13, sum(mu)]
f[year == 13, sum(af_prvl > 0)]
