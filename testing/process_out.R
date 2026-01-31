library(data.table)
library(yaml)
library(CKutils)
library(ggplot2)
library(ggthemes)
library(scales)
library(arrow)
library(IMPACTncdEngland)


prbl <- c(0.5, 0.025, 0.975, 0.1, 0.9)
baseline_year_for_change_outputs <- 2019L
theme_set(new = theme_economist())
theme_update(axis.text.x = element_text(size = 9), plot.title = element_text(hjust = 0.5))

# what <- "qalys"
# population <- "ons"
# strata <- "year"
# baseline_year <- baseline_year_for_change_outputs

# Prvl not standardised----
design <- Design$new("./testing/sim_design_testing.yaml")
sSummariesSubDirPath <- file.path(design$sim_prm$output_dir, "summaries")
sTablesSubDirPath <- file.path(design$sim_prm$output_dir, "tables")
output_dir <- design$sim_prm$output_dir

tbl_smmrs <- function(
    what = c(
            "prvl", "prvl_change", "incd", "incd_change",
            "ftlt", "ftlt_change", "mrtl", "mrtl_change",
            "dis_mrtl", "dis_mrtl_change",
            "cms_score", "cms_score_change", "cms_score_age",
            "cms_score_age_change", "cms_count", "cms_count_change",
            "pop", "qalys", "costs", "cypp", "cpp", "dpp",
            "net_qalys", "net_costs"
    ),
    population = c("ons", "esp"),
    strata,
    output_dir = output_dir,
    prbl = c(0.5, 0.025, 0.975, 0.1, 0.9),
    baseline_year = 2019L, # only used for prvl_change etc.
    comparator_scenario = "sc0",
    comparison_starting_year = baseline_year,
    two_agegrps = FALSE # if TRUE, agegrp is 30-64 and 65-99
    ) {
        strata <- lapply(strata, function(x) {
                c("mc", "scenario", x)
        })

        # construct file path to read from summaries
        str0 <- c(
                "prvl" = "prvl",
                "prvl_change" = "prvl",
                "incd" = "incd",
                "incd_change" = "incd",
                "ftlt" = "/dis_mrtl",
                "ftlt_change" = "/dis_mrtl",
                "mrtl" = "mrtl",
                "mrtl_change" = "mrtl",
                "dis_mrtl" = "dis_mrtl",
                "dis_mrtl_change" = "dis_mrtl",
                "cms_score" = "cms_score",
                "cms_score_change" = "cms_score",
                "cms_score_age" = "cms_score_by_age",
                "cms_score_age_change" = "cms_score_by_age",
                "cms_count" = "cms_count",
                "cms_count_change" = "cms_count",
                "qalys" = "qalys",
                "net_qalys" = "qalys",
                "costs" = "costs",
                "net_costs" = "costs",
                "cypp" = "prvl",
                "cpp" = "incd",
                "dpp" = "mrtl",
                "pop" = "prvl"
        )

        str1 <- c("ons" = "_scaled_up", "esp" = "_esp")

        # other useful strings
        str2 <- c(
                "prvl" = "_prvl$|^popsize$",
                "prvl_change" = "_prvl$|^popsize$",
                "incd" = "_incd$|^popsize$",
                "incd_change" = "_incd$|^popsize$",
                "ftlt" = "_deaths$|_prvl$",
                "ftlt_change" = "_deaths$|_prvl$",
                "mrtl" = "_mrtl$|^popsize$",
                "mrtl_change" = "_mrtl$|^popsize$",
                "dis_mrtl" = "^nonmodelled_deaths$|^chd_deaths$|^stroke_deaths$|^popsize$",
                "dis_mrtl_change" = "^nonmodelled_deaths$|^chd_deaths$|^stroke_deaths$|^popsize$",
                "cms_score" = "cms_score",
                "cms_score_change" = "cms_score",
                "cms_score_age" = "cms_score",
                "cms_score_age_change" = "cms_score",
                "cms_count" = "cms_count",
                "cms_count_change" = "cms_count",
                "qalys" = "^EQ5D5L$|^HUI3$",
                "net_qalys" = "^EQ5D5L$|^HUI3$",
                "costs" = "_costs$",
                "net_costs" = "_costs$",
                "cypp" = "_prvl$",
                "cpp" = "_incd$", 
                "dpp" = "_mrtl$",
                "pop" = "^popsize$"
        ) # used in grep
        str3 <- c(
                "prvl" = "prvl_rate_",
                "prvl_change" = "prct_change_",
                "incd" = "incd_rate_",
                "incd_change" = "prct_change_",
                "ftlt" = "ftlt_rate_",
                "ftlt_change" = "ftlt_rate_",
                "mrtl" = "mrtl_rate_",
                "mrtl_change" = "mrtl_change_",
                "dis_mrtl" = "disease_mrtl_rate_",
                "dis_mrtl_change" = "disease_mrtl_change_",                
                "cms_score" = "mean_cms_score_",
                "cms_score_change" = "mean_cms_score_",
                "cms_score_age" = "mean_cms_score_",
                "cms_score_age_change" = "mean_cms_score_",
                "cms_count" = "mean_cms_count_",
                "cms_count_change" = "mean_cms_count_",
                "qalys" = "qalys_",
                "net_qalys" = "net_qalys_",
                "costs" = "costs_",
                "net_costs" = "net_costs_",
                "cypp" = "cypp_",
                "cpp" = "cpp_",
                "dpp" = "dpp_",
                "pop" = "pop_size_"
        ) # used to col name output
        str4 <- c(
                "prvl" = "prevalence by ",
                "prvl_change" = "prevalence change by ",
                "incd" = "incidence by ",
                "incd_change" = "incidence change by ",
                "ftlt" = "case fatality by ",
                "ftlt_change" = "case fatality change by ",
                "mrtl" = "all-cause mortality by ",
                "mrtl_change" = "all-cause mortality change by ",
                "dis_mrtl" = "disease-specific mortality by ",
                "dis_mrtl_change" = "disease-specific mortality change by ",
                "cms_score" = "mean CMS score by ",
                "cms_score_change" = "mean CMS score change by ",
                "cms_score_age" = "mean CMS score by ",
                "cms_score_age_change" = "mean CMS score change by ",
                "cms_count" = "mean CMS count by ",
                "cms_count_change" = "mean CMS count change by ",
                "qalys" = "QALYs by ",
                "net_qalys" = "net QALYs by ",
                "costs" = "costs by ",
                "net_costs" = "net costs by ",
                "cypp" = "case-years prevented or postponed by ",
                "cpp" = "cases prevented or postponed by ",
                "dpp" = "deaths prevented or postponed by ",
                "pop" = "pop size by "
        ) # used in output filenames

        fpth <- file.path(output_dir, "summaries", paste0(str0[[what]], str1[[population]]))
        if (!file.exists(fpth)) {
                message(fpth, " doesn't exist")
                return(NULL)
        }

        tt <- as.data.table(open_dataset(fpth)) # numerator data

        # Check if comparison metrics are requested but only baseline scenario exists
        comparison_metrics <- c("cypp", "cpp", "dpp", "net_qalys", "net_costs")
        if (what %in% comparison_metrics) {
                available_scenarios <- unique(tt$scenario)
                non_comparator_scenarios <- setdiff(available_scenarios, comparator_scenario)
                if (length(non_comparator_scenarios) == 0) {
                        message("Skipping ", what, " because no intervention scenarios exist (only '",
                                comparator_scenario, "' found). Comparison metrics require at least one non-baseline scenario.")
                        return(NULL)
                }
        }

		if (two_agegrps) {
                sTablesSubDirPath <- file.path(design$sim_prm$output_dir, "tables2agegrps/")
                if (!dir.exists(sTablesSubDirPath)) dir.create(sTablesSubDirPath)
                if ("agegrp" %in% names(tt)) {
                        tt[agegrp %in% c("30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64"), agegrp := "30-64"]
                        tt[agegrp %in% c("65-69", "70-74", "75-79", "80-84", "85-89", "90-94", "95-99"), agegrp := "65-99"]
                }
        }

        # For ftlt I need prvl for the denominator
        if (grepl("^ftlt", what)) {
                fpth <- file.path(output_dir, "summaries", paste0(str0[["prvl"]], str1[[population]]))
                if (!file.exists(fpth)) stop(fpth, " doesn't exist")

                t1 <- as.data.table(open_dataset(fpth)) # denominator data  
                setnames(t1, "popsize", "nonmodelled_prvl")
				if (two_agegrps && "agegrp" %in% names(t1)) {
                        t1[agegrp %in% c("30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64"), agegrp := "30-64"]
                        t1[agegrp %in% c("65-69", "70-74", "75-79", "80-84", "85-89", "90-94", "95-99"), agegrp := "65-99"]
                }
                absorb_dt(tt, t1)
                tt <- tt[nonmodelled_prvl > 0] # This is the denom. Cannot be 0 and it is meaningless anyway
        }

        lapply(strata, function(x) {
                if (grepl("^cms_", what)) {
                        d <- tt[, .("value" = weighted.mean(get(str2[[what]]), popsize)),
                                keyby = eval(x)
                        ]

                        if (grepl("_change$", what)) { # when calculating change
                                d19 <- d[year == baseline_year][, year := NULL]
                                d[d19, on = c(setdiff(x, "year")), value := value / i.value]
                        }
                        d <- d[, as.list(fquantile(value, prbl)), keyby = eval(setdiff(x, "mc"))]
                        setnames(d, c(setdiff(x, "mc"), percent(prbl, prefix = str3[[what]])))
                        setkeyv(d, setdiff(x, "mc"))
                        setcolorder(d, setdiff(x, "mc"))
                } else if (grepl("^qalys$", what)) {
                        d <- tt[, .("EQ5D5L" = sum(EQ5D5L),
                                    "HUI3" = sum(HUI3)
                                   ),
                                keyby = eval(x)
                        ]
                        d <- melt(d, id.vars = x, variable.name = "scale", value.name = "QALYs")

                        # if (grepl("_change$", what)) { # when calculating change
                        #         d19 <- d[year == baseline_year][, year := NULL]
                        #         d[d19, on = c(setdiff(x, "year"), "scale"), QALYs := QALYs / i.QALYs]
                        # }
                        setkeyv(d, c(x[x != "year"], "scale", "year"))
                        d[, cumulative := cumsum(QALYs), keyby = c(setdiff(x, "year"), "scale")]
                        d <- melt(d, id.vars = c(x, "scale"), variable.name = "type")
                        d[, type := fifelse(type == "cumulative", "QALYs_cuml", "QALYs")]

                        setkey(d, "type", "scale")
                        d <-
                                d[, fquantile_byid(value, prbl, id = as.character(type), rounding = FALSE),
                                        keyby = eval(setdiff(c(x, "scale"), "mc"))
                                ]
                        setnames(d, c(setdiff(c(x, "scale"), "mc"), "type", percent(prbl, prefix = str3[[what]])))
                        setkeyv(d, c("type", setdiff(c(x, "scale"), "mc")))
                        setcolorder(d, setdiff(c(x, "scale"), "mc"))


                } else if (grepl("^net_qalys$", what)) {
                        d <- tt[, .("EQ5D5L" = sum(EQ5D5L),
                                    "HUI3" = sum(HUI3)
                                   ),
                                keyby = eval(x)
                        ]
                        d <- melt(d, id.vars = x, variable.name = "scale", value.name = "QALYs")
                        d_sc0 <- d[scenario == comparator_scenario & year >= comparison_starting_year][, scenario := NULL]
                        d <- d[scenario != comparator_scenario & year >= comparison_starting_year][d_sc0, on = c(setdiff(x, "scenario"), "scale"), net_QALYs := QALYs - i.QALYs] # positive numbers for prevention
                        d[, QALYs := NULL]
                        setkeyv(d, c(x[x != "year"], "scale", "year"))
                        d[, cumulative := cumsum(net_QALYs), keyby = c(setdiff(x, "year"), "scale")]
                        d <- melt(d, id.vars = c(x, "scale"), variable.name = "type")
                        d[type == "cumulative", type := "net_QALYs_cuml"]
                        setkey(d, "type", "scale")
                        d <-
                          d[, fquantile_byid(value, prbl, id = as.character(type), rounding = FALSE),
                          keyby = eval(setdiff(c(x, "scale"), "mc"))]
                        x <- c(x, "scale")
                        setnames(d, c(setdiff(x, "mc"), "type", percent(prbl, prefix = str3[[what]])))
                        setkeyv(d, c("type", setdiff(x, "mc")))
                        setcolorder(d, setdiff(x, "mc"))                        
                } else if (grepl("^costs", what)) {
                        d <- tt[, lapply(.SD, sum), .SDcols = patterns("_costs$"), keyby = eval(x)]
                        d <- melt(d, id.vars = x, variable.name = "costs_type", value.name = "costs")

                        # if (grepl("_change$", what)) { # when calculating change
                        #         d19 <- d[year == baseline_year][, year := NULL]
                        #         d[d19, on = c(setdiff(x, "year"), "type"), costs := costs / i.costs]
                        # }
                        d[, cumulative := cumsum(costs), keyby = c(setdiff(x, "year"), "costs_type")]
                        d <- melt(d, id.vars = c(x, "costs_type"), variable.name = "type")
                        d[type == "cumulative", type := "costs_cuml"]
                        setkey(d, "type", "costs_type")
                        d <-
                                d[, fquantile_byid(value, prbl, id = as.character(type), rounding = FALSE),
                                        keyby = eval(setdiff(c(x, "costs_type"), "mc"))
                                ]
                        setnames(d, c(
                                setdiff(c(x, "costs_type"), "mc"),
                                "type",
                                percent(prbl, prefix = str3[[what]])
                        ))
                        setkeyv(d, c("type", setdiff(c(x, "costs_type"), "mc")))
                        setcolorder(d, setdiff(c(x, "costs_type"), "mc"))

                } else if (grepl("^net_costs", what)) {
                        d <- tt[, lapply(.SD, sum), .SDcols = patterns("_costs$"), keyby = eval(x)]
                        d <- melt(d, id.vars = x, variable.name = "costs_type", value.name = "value")
                        d_sc0 <- d[scenario == comparator_scenario & year >= comparison_starting_year][, scenario := NULL]
                        d <- d[scenario != comparator_scenario & year >= comparison_starting_year][d_sc0, on = c(setdiff(x, "scenario"), "costs_type"), net_costs := value - i.value] # negative numbers for prevention
                        d[, value := NULL]
                        setkeyv(d, c(x[x != "year"], "costs_type", "year"))
                        d[, cumulative := cumsum(net_costs), keyby = c(setdiff(x, "year"), "costs_type")]
                        d <- melt(d, id.vars = c(x, "costs_type"), variable.name = "type")
                        d[type == "cumulative", type := "net_costs_cuml"]
                        setkey(d, "type", "costs_type")
                        d <-
                                d[, fquantile_byid(value, prbl, id = as.character(type), rounding = FALSE),
                                        keyby = eval(setdiff(c(x, "costs_type"), "mc"))
                                ]
                        x <- c(x, "costs_type")
                        setnames(d, c(
                                setdiff(x, "mc"),
                                "type",
                                percent(prbl, prefix = str3[[what]])
                        ))
                        setkeyv(d, c("type", setdiff(x, "mc")))
                        setcolorder(d, setdiff(x, "mc"))

                } else { # if not cms or qalys or costs...
                        d <- tt[, lapply(.SD, sum),
                                .SDcols = patterns(str2[[what]]),
                                keyby = x
                        ]
                        # convert int cols to numeric (avoids warning with melt())
                        is_int <- sapply(d[, .SD, .SDcols = -x], is.integer)
                        is_int <- names(is_int[is_int])
                        d[, (is_int) := lapply(.SD, as.numeric), .SDcols = is_int]
                        
                        if (grepl("^ftlt", what)) {
                                nm <- names(d)
                                nm <- grep("_deaths$", nm, value = TRUE)
                                nm <- gsub("_deaths$", "", nm)
                                nm <- setdiff(nm, "alive")
                                for (i in nm) {
                                        set(
                                                d, NULL, paste0(i, "_ftlt"),
                                                d[[paste0(i, "_deaths")]] / d[[paste0(i, "_prvl")]]
                                        )
                                }

                                nm <- names(d)
                                nm <- grep("_deaths$|_prvl$", nm, value = TRUE)
                                d[, (nm) := NULL]
                                setnafill(d, "const", 0, cols = grep("_ftlt$", names(d), value = TRUE))
                        } else if (!what %in% c("pop", "cypp", "cpp", "dpp")) { # avoid calculating rates for pop, cypp, cpp, dpp
                                d <- d[, lapply(.SD, function(y) {
                                        y / popsize
                                }), keyby = x]
                        }
                        

                        d <- melt(d, id.vars = x)

                        if (grepl("_change$", what)) { # when calculating change
                                d19 <- d[year == baseline_year][, year := NULL]
                                d[d19, on = c(setdiff(x, "year"), "variable"), value := value / i.value]
                        }

                        if (grepl("^cypp$|^cpp$|^dpp$", what)) {
                                d_sc0 <- d[scenario == comparator_scenario & year >= comparison_starting_year][, scenario := NULL]
                                d <- d[scenario != comparator_scenario & year >= comparison_starting_year][d_sc0, on = c(setdiff(x, "scenario"), "variable"), value := i.value - value] # positive numbers for prevention
                                d[, variable := gsub(paste0("_", str0[[what]]), "", variable)]
                                setkeyv(d, c(x[x != "year"], "variable", "year"))
                                d[, cumulative := cumsum(value), keyby = c(setdiff(x, "year"), "variable")]
                                d <- melt(d, id.vars = c(x, "variable"), variable.name = "type")
                                d[, type := fifelse(type == "cumulative", paste0(what, "_cuml"), what)]
                                setkey(d, "type", "variable")
                                d <-
                                        d[, fquantile_byid(value, prbl, id = as.character(variable), rounding = (what %in% c("pop", "cypp", "cpp", "dpp"))),
                                                keyby = eval(setdiff(c(x, "type"), "mc"))
                                        ]
                                x <- c(x, "type")
                                setnames(d, c(setdiff(x, "mc"), "disease", percent(prbl, prefix = str3[[what]])))
                        } else {
                                setkey(d, "variable")
                                d <-
                                        d[, fquantile_byid(value, prbl, id = as.character(variable), rounding = what == "pop"),
                                                keyby = eval(setdiff(x, "mc"))
                                        ]
                                setnames(d, c(setdiff(x, "mc"), "disease", percent(prbl, prefix = str3[[what]])))
                        }

                        if (what == "pop") {
                                d[, disease := NULL]
                        } else {
                             if ("popsize" %in% d$disease) d <- d[disease != "popsize"]
                        }
                  setkeyv(d, setdiff(x, "mc"))
                  setcolorder(d, setdiff(x, "mc"))
                }
                str5 <- c(
                        "ons" = " (not standardised).csv",
                        "esp" = paste0(" (", paste(setdiff(c("mc", "scenario", "year", "age", "sex"), x),
                                collapse = "-"
                        ), " standardised).csv")
                )
        
                str6 <- paste0(
                        str4[[what]],
                        paste(setdiff(x, c("mc", "scenario", "type", "scale", "costs_type")), collapse = "-"),
                        str5[[population]]
                ) # used for output file name/path
                fwrite(d, file.path(
                        output_dir, ifelse(two_agegrps, "tables2agegrps", "tables"), str6
                ))
        })
}


outperm <- expand.grid(
        what = c(
                "prvl", "prvl_change", "incd", "incd_change",
                "ftlt", "ftlt_change", "mrtl", "mrtl_change",
                "dis_mrtl", "dis_mrtl_change", "qalys", "costs",
                # "cms_score", "cms_score_change", "cms_score_age",
                # "cms_score_age_change", "cms_count", "cms_count_change",
                "cypp", "cpp", "dpp", "net_qalys", "net_costs", "pop"
        ),
        population = c("ons", "esp")
)

for (i in seq_len(nrow(outperm))) {
        what <- as.character(outperm$what[[i]])
        population <- as.character(outperm$population[[i]])
        if (population == "ons") {
                strata <- list(
                        "year",
                        c("year", "sex"),
                        c("year", "agegrp"),
                        c("year", "agegrp", "sex")
                )
        } else if (population == "esp") {
                strata <- list(
                        "year",
                        c("year", "sex")
                )
        } else {
                stop()
        }

        if (grepl("_age", what)) {
                strata <- lapply(strata, function(st) {
                        st[st == "agegrp"] <- "age"
                        st
                })
        }

        if ((what == "pop" && population == "esp") || (grepl("_age", what) && population == "esp")) next()

        print(paste0(what, "-", population))
        tbl_smmrs(what = what, population = population, strata = strata, output_dir = output_dir,
                prbl = prbl, baseline_year = baseline_year_for_change_outputs,
                comparator_scenario = "sc0",
                comparison_starting_year = baseline_year_for_change_outputs
        )
}

tbl_smmrs(what = "pop", population = "ons", list(
        "year",
        c("year", "sex"),
        c("year", "agegrp"),
        c("year", "agegrp", "sex")
), output_dir, prbl = prbl, baseline_year = baseline_year_for_change_outputs,
   comparator_scenario = "sc0",
   comparison_starting_year = baseline_year_for_change_outputs,
   two_agegrps = FALSE)

# 2 agegroups ----
outperm <- expand.grid(
        what = c(
                "prvl", "prvl_change", "incd", "incd_change",
                "ftlt", "ftlt_change", "mrtl", "mrtl_change",
                "dis_mrtl", "dis_mrtl_change", "qalys", "costs",
                # "cms_score", "cms_score_change", "cms_score_age",
                # "cms_score_age_change", "cms_count", "cms_count_change",
                "cypp", "cpp", "dpp", "net_qalys", "net_costs", "pop"
        ),
        population = "ons")
for (i in seq_len(nrow(outperm))) {
        what <- as.character(outperm$what[[i]])
        population <- as.character(outperm$population[[i]])
        if (population == "ons") {
                strata <- list(
                        c("year", "agegrp"),
                        c("year", "agegrp", "sex")
                )
        } else if (population == "esp") {
                strata <- list(
                        "year",
                        c("year", "sex")
                )
        } else {
                stop()
        }

        if (grepl("_age", what)) {
                strata <- lapply(strata, function(st) {
                        st[st == "agegrp"] <- "age"
                        st
                })
        }


		if ((what == "pop" && population == "esp") || (grepl("_age", what) && population == "esp")) next()


        print(paste0(what, "-", population))
        tbl_smmrs(what = what, population = population, strata = strata, output_dir = output_dir,
                prbl = prbl, baseline_year = baseline_year_for_change_outputs,
                comparator_scenario = "sc0",
                comparison_starting_year = baseline_year_for_change_outputs,
                two_agegrps = TRUE
        )
}

tbl_smmrs(what = "pop", population = "ons", list(
        "year",
        c("year", "sex"),
        c("year", "agegrp"),
        c("year", "agegrp", "sex")
), output_dir, prbl = prbl, baseline_year = baseline_year_for_change_outputs,
   comparator_scenario = "sc0",
   comparison_starting_year = baseline_year_for_change_outputs,
   two_agegrps = TRUE)


# All-cause mortality by disease not standardised ----
tt <- as.data.table(open_dataset(file.path(sSummariesSubDirPath, "all_cause_mrtl_by_dis_scaled_up")))

outstrata <- c("mc", "year", "scenario")
d <- tt[, lapply(.SD, sum),
        .SDcols = patterns("^deaths_|^cases_"),
        keyby = eval(outstrata)
]
d <- melt(d, id.vars = outstrata)
cases <- d[grep("^cases_", variable)][, variable := gsub("^cases_", "", variable)]
d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
d[cases, on = c(outstrata, "variable"), value := value / i.value]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, file.path(sTablesSubDirPath, "all-cause mortality given disease-year (not standardised).csv"))

outstrata <- c("mc", "year", "sex", "scenario")
d <- tt[, lapply(.SD, sum),
        .SDcols = patterns("^deaths_|^cases_"),
        keyby = eval(outstrata)
]
d <- melt(d, id.vars = outstrata)
cases <- d[grep("^cases_", variable)][, variable := gsub("^cases_", "", variable)]
d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
d[cases, on = c(outstrata, "variable"), value := value / i.value]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, file.path(sTablesSubDirPath, "all-cause mortality given disease-year-sex (not standardised).csv"))

outstrata <- c("mc", "year", "agegrp", "sex", "scenario")
d <- tt[, lapply(.SD, sum),
        .SDcols = patterns("^deaths_|^cases_"),
        keyby = eval(outstrata)
]
d <- melt(d, id.vars = outstrata)
cases <- d[grep("^cases_", variable)][, variable := gsub("^cases_", "", variable)]
d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
d[cases, on = c(outstrata, "variable"), value := value / i.value]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, file.path(sTablesSubDirPath, "all-cause mortality given disease-year-agegroup-sex (not standardised).csv"))

outstrata <- c("mc", "year", "agegrp", "sex", "scenario")
d <- tt[, lapply(.SD, sum),
        .SDcols = patterns("^deaths_|^cases_"),
        keyby = eval(outstrata)
]
d <- melt(d, id.vars = outstrata)
cases <- d[grep("^cases_", variable)][, variable := gsub("^cases_", "", variable)]
d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
d[cases, on = c(outstrata, "variable"), value := value / i.value]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, file.path(sTablesSubDirPath, "all-cause mortality given disease-year-agegroup-sex (not standardised).csv"))

# All-cause mortality by disease not standardised pop denominator----
tt <- as.data.table(open_dataset(file.path(sSummariesSubDirPath, "all_cause_mrtl_by_dis_scaled_up")))
pp <- as.data.table(open_dataset(file.path(sSummariesSubDirPath, "prvl_scaled_up")))


outstrata <- c("mc", "year", "scenario")
cases <- pp[, lapply(.SD, sum), .SDcols = patterns("^popsize$"), keyby = eval(outstrata)]
d <- tt[, lapply(.SD, sum),
        .SDcols = patterns("^deaths_|^cases_"),
        keyby = eval(outstrata)
]
d <- melt(d, id.vars = outstrata)
d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
d[cases, on = outstrata, value := value / popsize]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, file.path(sTablesSubDirPath, "all-cause mortality given disease-year popdenom (not standardised).csv"))

outstrata <- c("mc", "year", "sex", "scenario")
cases <- pp[, lapply(.SD, sum), .SDcols = patterns("^popsize$"), keyby = eval(outstrata)]
d <- tt[, lapply(.SD, sum),
        .SDcols = patterns("^deaths_|^cases_"),
        keyby = eval(outstrata)
]
d <- melt(d, id.vars = outstrata)
d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
d[cases, on = outstrata, value := value / popsize]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, file.path(sTablesSubDirPath, "all-cause mortality given disease-year-sex popdenom (not standardised).csv"))

outstrata <- c("mc", "year", "agegrp", "sex", "scenario")
cases <- pp[, lapply(.SD, sum), .SDcols = patterns("^popsize$"), keyby = eval(outstrata)]
d <- tt[, lapply(.SD, sum),
        .SDcols = patterns("^deaths_|^cases_"),
        keyby = eval(outstrata)
]
d <- melt(d, id.vars = outstrata)
d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
d[cases, on = outstrata, value := value / popsize]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, file.path(sTablesSubDirPath, "all-cause mortality given disease-year-agegroup-sex popdenom (not standardised).csv"))

outstrata <- c("mc", "year", "agegrp", "sex", "scenario")
cases <- pp[, lapply(.SD, sum), .SDcols = patterns("^popsize$"), keyby = eval(outstrata)]
d <- tt[, lapply(.SD, sum),
        .SDcols = patterns("^deaths_|^cases_"),
        keyby = eval(outstrata)
]
d <- melt(d, id.vars = outstrata)
d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
d[cases, on = outstrata, value := value / popsize]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, file.path(sTablesSubDirPath, "all-cause mortality given disease-year-agegroup-sex popdenom (not standardised).csv"))

rm(pp)

# All-cause mortality by disease standardised----
tt <- as.data.table(open_dataset(file.path(sSummariesSubDirPath, "all_cause_mrtl_by_dis_esp")))
outstrata <- c("mc", "year", "scenario")
d <- tt[, lapply(.SD, sum),
        .SDcols = patterns("^deaths_|^cases_"),
        keyby = eval(outstrata)
]
d <- melt(d, id.vars = outstrata)
cases <- d[grep("^cases_", variable)][, variable := gsub("^cases_", "", variable)]
d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
d[cases, on = c(outstrata, "variable"), value := value / i.value]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, file.path(sTablesSubDirPath, "all-cause mortality given disease-year (age-sex standardised).csv"))

outstrata <- c("mc", "year", "sex", "scenario")
d <- tt[, lapply(.SD, sum),
        .SDcols = patterns("^deaths_|^cases_"),
        keyby = eval(outstrata)
]
d <- melt(d, id.vars = outstrata)
cases <- d[grep("^cases_", variable)][, variable := gsub("^cases_", "", variable)]
d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
d[cases, on = c(outstrata, "variable"), value := value / i.value]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, file.path(sTablesSubDirPath, "all-cause mortality given disease-year-sex (age standardised).csv"))

outstrata <- c("mc", "year", "scenario")
d <- tt[, lapply(.SD, sum),
        .SDcols = patterns("^deaths_|^cases_"),
        keyby = eval(outstrata)
]
d <- melt(d, id.vars = outstrata)
cases <- d[grep("^cases_", variable)][, variable := gsub("^cases_", "", variable)]
d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
d[cases, on = c(outstrata, "variable"), value := value / i.value]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, file.path(sTablesSubDirPath, "all-cause mortality given disease-year (age-sex standardised).csv"))

outstrata <- c("mc", "year", "sex", "scenario")
d <- tt[, lapply(.SD, sum),
        .SDcols = patterns("^deaths_|^cases_"),
        keyby = eval(outstrata)
]
d <- melt(d, id.vars = outstrata)
cases <- d[grep("^cases_", variable)][, variable := gsub("^cases_", "", variable)]
d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
d[cases, on = c(outstrata, "variable"), value := value / i.value]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, file.path(sTablesSubDirPath, "all-cause mortality given disease-year-sex (age standardised).csv"))
rm(cases)


# Disease characteristics non standardised ----
tt <- as.data.table(open_dataset(file.path(sSummariesSubDirPath, "dis_characteristics_scaled_up")))

tt[, `:=`(mean_cms_count_cms1st_cont = as.numeric(mean_cms_count_cms1st_cont))]
d1 <- tt[, .SD, .SDcols = patterns("mc|scenario|year|sex|^cases_")]
d1 <- melt(d1, id.vars = c("mc", "year", "scenario", "sex"))
d1 <- unique(d1, by = c("mc", "year", "scenario", "sex", "variable"))
d1[, `:=`(disease = gsub("^cases_", "", variable), variable = NULL)]
tt <- tt[, .SD, .SDcols = patterns("mc|scenario|year|sex|^mean_duration_|^mean_age_incd_|^mean_age_1st_onset_|^mean_age_prvl_|^mean_cms_score_|^mean_cms_count_")]
tt[, mean_cms_count_cmsmm1 := as.double(mean_cms_count_cmsmm1)]
tt <- melt(tt, id.vars = c("mc", "year", "scenario", "sex"))
tt[, disease := gsub("^mean_duration_|^mean_age_incd_|^mean_age_1st_onset_|^mean_age_prvl_|^mean_cms_score_|^mean_cms_count_", "", variable)]
tt[d1, on = c("mc", "year", "scenario", "sex", "disease"), cases := i.value]
rm(d1) # NOTE mean_age_incd contains NAs

outstrata <- c("mc", "year", "scenario")
d <- tt[, weighted.mean(value, cases, na.rm = TRUE), keyby = c(outstrata, "variable")] # na.rm = TRUE for mean_age_incd
setkey(d, "variable")
d <- d[, fquantile_byid(V1, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "variable", percent(prbl, prefix = "value_")))
d[, disease := gsub("^mean_duration_|^mean_age_incd_|^mean_age_1st_onset_|^mean_age_prvl_|^mean_cms_score_|^mean_cms_count_", "", variable)]
d[grep("^mean_duration_", variable), type := "mean_duration"]
d[grep("^mean_age_incd_", variable), type := "mean_age_incd"]
d[grep("^mean_age_1st_onset_", variable), type := "mean_age_1st_onset"]
d[grep("^mean_age_prvl_", variable), type := "mean_age_prvl"]
d[grep("^mean_cms_score_", variable), type := "mean_cms_score"]
d[grep("^mean_cms_count_", variable), type := "mean_cms_count"]
d[, variable := NULL]
setkeyv(d, c(setdiff(outstrata, "mc"), "disease", "type"))
setcolorder(d)
fwrite(d, file.path(sTablesSubDirPath, "disease characteristics by year (not standardised).csv"))

outstrata <- c("mc", "year", "sex", "scenario")
d <- tt[, weighted.mean(value, cases, na.rm = TRUE), keyby = c(outstrata, "variable")]
setkey(d, "variable")
d <- d[, fquantile_byid(V1, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "variable", percent(prbl, prefix = "value_")))
d[, disease := gsub("^mean_duration_|^mean_age_incd_|^mean_age_1st_onset_|^mean_age_prvl_|^mean_cms_score_|^mean_cms_count_", "", variable)]
d[grep("^mean_duration_", variable), type := "mean_duration"]
d[grep("^mean_age_incd_", variable), type := "mean_age_incd"]
d[grep("^mean_age_1st_onset_", variable), type := "mean_age_1st_onset"]
d[grep("^mean_age_prvl_", variable), type := "mean_age_prvl"]
d[grep("^mean_cms_score_", variable), type := "mean_cms_score"]
d[grep("^mean_cms_count_", variable), type := "mean_cms_count"]
d[, variable := NULL]
setkeyv(d, c(setdiff(outstrata, "mc"), "disease", "type"))
setcolorder(d)
fwrite(d, file.path(sTablesSubDirPath, "disease characteristics by year-sex (not standardised).csv"))

outstrata <- c("mc", "year", "scenario")
d <- tt[, weighted.mean(value, cases, na.rm = TRUE), keyby = c(outstrata, "variable")]
setkey(d, "variable")
d <- d[, fquantile_byid(V1, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "variable", percent(prbl, prefix = "value_")))
d[, disease := gsub("^mean_duration_|^mean_age_incd_|^mean_age_1st_onset_|^mean_age_prvl_|^mean_cms_score_|^mean_cms_count_", "", variable)]
d[grep("^mean_duration_", variable), type := "mean_duration"]
d[grep("^mean_age_incd_", variable), type := "mean_age_incd"]
d[grep("^mean_age_1st_onset_", variable), type := "mean_age_1st_onset"]
d[grep("^mean_age_prvl_", variable), type := "mean_age_prvl"]
d[grep("^mean_cms_score_", variable), type := "mean_cms_score"]
d[grep("^mean_cms_count_", variable), type := "mean_cms_count"]
d[, variable := NULL]
setkeyv(d, c(setdiff(outstrata, "mc"), "disease", "type"))
setcolorder(d)
fwrite(d, file.path(sTablesSubDirPath, "disease characteristics by year (not standardised).csv"))

outstrata <- c("mc", "year", "sex", "scenario")
d <- tt[, weighted.mean(value, cases, na.rm = TRUE), keyby = c(outstrata, "variable")]
setkey(d, "variable")
d <- d[, fquantile_byid(V1, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "variable", percent(prbl, prefix = "value_")))
d[, disease := gsub("^mean_duration_|^mean_age_incd_|^mean_age_1st_onset_|^mean_age_prvl_|^mean_cms_score_|^mean_cms_count_", "", variable)]
d[grep("^mean_duration_", variable), type := "mean_duration"]
d[grep("^mean_age_incd_", variable), type := "mean_age_incd"]
d[grep("^mean_age_1st_onset_", variable), type := "mean_age_1st_onset"]
d[grep("^mean_age_prvl_", variable), type := "mean_age_prvl"]
d[grep("^mean_cms_score_", variable), type := "mean_cms_score"]
d[grep("^mean_cms_count_", variable), type := "mean_cms_count"]
d[, variable := NULL]
setkeyv(d, c(setdiff(outstrata, "mc"), "disease", "type"))
setcolorder(d)
fwrite(d, file.path(sTablesSubDirPath, "disease characteristics by year-sex (not standardised).csv"))

rm(d, tt)

# XPS ----
xps_tab <- as.data.table(open_dataset(file.path(design$sim_prm$output_dir, "xps", "xps20")))

xps_names <- grep("_curr_xps$", names(xps_tab), value = TRUE)

outstrata <- c("mc", "year", "agegrp20", "sex", "scenario")
d <- xps_tab[sex != "All" & agegrp20 != "All"] # This should depend on outstrata
d <- d[, lapply(.SD, mean), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "exposure", percent(prbl, prefix = "xps_mean_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, file.path(sTablesSubDirPath, "exposures by year-agegroup-sex (not standardised).csv"))

outstrata <- c("mc", "year", "scenario")
d <- xps_tab[sex == "All" & agegrp20 == "All"] # This should depend on outstrata
d <- d[, lapply(.SD, mean), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "exposure", percent(prbl, prefix = "xps_mean_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, file.path(sTablesSubDirPath, "exposures by year (not standardised).csv"))

outstrata <- c("mc", "year", "agegrp20", "scenario")
d <- xps_tab[sex == "All" & agegrp20 != "All"] # This should depend on outstrata
d <- d[, lapply(.SD, mean), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "exposure", percent(prbl, prefix = "xps_mean_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, file.path(sTablesSubDirPath, "exposures by year-agegroup (not standardised).csv"))

# Soshiro added programs from here 2023 12 22 due to missing the exposure file by "year-sex" (not standardised).vs 
outstrata <- c("mc", "year", "sex", "scenario")
d <- xps_tab[sex != "All"] # This should depend on outstrata
d <- d[, lapply(.SD, mean), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "exposure", percent(prbl, prefix = "xps_mean_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, file.path(sTablesSubDirPath, "exposures by year-sex (not standardised).csv"))


# XPS standardised ----
xps_tab <- as.data.table(open_dataset(file.path(design$sim_prm$output_dir, "xps", "xps5")))
xps_names <- grep("_curr_xps$", names(xps_tab), value = TRUE)

outstrata <- c("mc", "year", "sex", "scenario")
d <- xps_tab[sex != "All"] # This should depend on outstrata
d <- d[, lapply(.SD, mean), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "exposure", percent(prbl, prefix = "xps_mean_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, file.path(sTablesSubDirPath, "exposures by year-sex (age standardised).csv"))

outstrata <- c("mc", "year", "scenario")
d <- xps_tab[sex == "All"] # This should depend on outstrata
d <- d[, lapply(.SD, mean), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "exposure", percent(prbl, prefix = "xps_mean_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, file.path(sTablesSubDirPath, "exposures by year (age-sex standardised).csv"))

