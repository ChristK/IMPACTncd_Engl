## IMPACTncd_Engl is an implementation of the IMPACTncd framework, developed by Chris
## Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been
## funded by NIHR  HTA Project: 16/165/01 - IMPACTncd_Engl: Health Outcomes
## Research Simulation Environment.  The views expressed are those of the
## authors and not necessarily those of the NHS, the NIHR or the Department of
## Health.
##
## Copyright (C) 2018-2020 University of Liverpool, Chris Kypridemos
##
## IMPACTncd_Engl is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version. This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details. You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/> or write
## to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
## Boston, MA 02110-1301 USA.

cat("Validating IMPACTncd_Engl model...\n\n")
# setwd("~/My_Models/IMPACTncd_Engl/")
XPS <- TRUE
delete_previous_results <- FALSE # any change to IMPACTncdEnglmisc deletes output

# necessary for Rscript
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.0")

output_dir <-
  function(x = character(0))
    paste0("./validation/plot_trends/", x)

if (delete_previous_results)
  file.remove(list.files(output_dir(), full.names = TRUE, recursive = TRUE))

if (!require(CKutils)) {
  if (!require(remotes))
    install.packages("remotes")
  roxygen2::roxygenise("../CKutils/") # TODO remove before deployment
  remotes::install_github("ChristK/CKutils")
  library(CKutils)
}

# check if IMPACTncdEnglmisc has changed, reinstall package and set switch to_delete
# <- NULL to pass down below and delete existing synthpops
to_delete <- -1L
if (file.exists("./Rpackage/.version_snapshot.rds")) {
  snapshot <-
    readRDS("./Rpackage/.version_snapshot.rds")
  tt <- changedFiles(snapshot)
  if (sum(sapply(tt[1:3], length)) > 0) {
    roxygen2::roxygenise("./Rpackage/IMPACTncd_Engl_model_pkg/")
    # TODO remove before deployment
    remotes::install_local("./Rpackage/IMPACTncd_Engl_model_pkg/", force = TRUE)
    to_delete <- NULL
    file.remove(list.files(output_dir(), full.names = TRUE, recursive = TRUE))

    saveRDS(fileSnapshot(list.dirs("./Rpackage/IMPACTncd_Engl_model_pkg/"),
                         md5sum = TRUE),
            "./Rpackage/.version_snapshot.rds")
  }

} else {
  # if snapshot doesn't exists.
  saveRDS(fileSnapshot(list.dirs("./Rpackage/IMPACTncd_Engl_model_pkg/"),
                       md5sum = TRUE),
          "./Rpackage/.version_snapshot.rds")
}

library(IMPACTncdEnglmisc)

dependencies(
  c(
    # "gamlss", # only necessary when fitting the models
    # "gamlss.tr", # only necessary when fitting the models
    # "mc2d", # only necessary for generating fixed_mc
    "doParallel",
    "doRNG",
    "foreach",
    # "mc2d", # for rpert()
    "gamlss.dist",
    # For distr in prevalence.R
    "dqrng",
    "qs",
    "fst",
    "wrswoR",
    "ggplot2",
    "cowplot",
    "viridis",
    "dichromat",
    "future",
    "data.table"
  )
)
theme_set(theme_cowplot())
options(datatable.verbose = FALSE)
options(datatable.showProgress = FALSE)

design <-
  Design$new("./validation/sim_design_for_trends_validation.yaml")
parameters <- qread("./validation/parameters.qs")
parameters <- fromGUI_prune(parameters) # TODO delete for production
design$update_fromGUI(parameters)

options(future.fork.enable = TRUE) # enable fork in Rstudio TODO remove for production
plan(multiprocess, workers = design$sim_prm$clusternumber)


pr <-
  c(0.5, 0.025, 0.975, 0.1, 0.9, 0.2, 0.8) # percentiles for summaries

# Calculate quantiles TODO correct handling of rates
summarise_distr <-
  function(dt,
           probabilities,
           groupings = FALSE,
           rounding = FALSE) {
    stopifnot(is.data.table(dt))
    past_keys <- key(dt)
    setkeyv(dt, "variable")
    if (groupings) {
      # TODO dynamic (not hardcoded)
      # Note: xps_out already in groupings
      out <- rbind(dt[, fifelse(grepl("_rate$", variable),
                                weighted.mean(value, pop_size),
                                sum(value)),
                      by = c("scenario", "year", "variable", "mc")][, fquantile_byid(V1, probabilities,
                       as.character(variable), rounding),
                        by = c("scenario", "year")],
                   dt[, fifelse(grepl("_rate$", variable),
                                weighted.mean(value, pop_size),
                                sum(value)),
                      by =  c("scenario", "year", "sex", "variable", "mc")][, fquantile_byid(V1, probabilities,
                       as.character(variable), rounding),
                       by =  c("scenario", "year", "sex")],
                   dt[, fifelse(grepl("_rate$", variable),
                                weighted.mean(value, pop_size),
                                sum(value)),
                      by =  c("scenario", "year", "sex", "agegroup", "variable", "mc")][, fquantile_byid(V1, probabilities,
                        as.character(variable), rounding),
                        by =  c("scenario", "year", "sex", "agegroup")],
                   fill = TRUE)
    } else {
      out <-
        # xps_out already in groupings. I don't use weighted mean cause xps_out
        # have no pop_size
        dt[, mean(value),
           by = c("year",
                  "agegrp20", "sex", "qimd",
                  "mc", "variable")
           ][, fquantile_byid(V1, probabilities,
                                  as.character(variable),
                                  rounding),
               by = c("year", "agegrp20", "sex", "qimd")]
    }
    setkeyv(dt, past_keys)
    setnames(out, "V1", "variable")
    setnames(out, paste0("V", seq_along(pr) + 1L), scales::percent(pr))
    for (j in seq_len(ncol(out)))
      set(out, which(is.na(out[[j]])), j, "All")
    setkey(out, year)
    invisible(out)
  }


std_graph <-
  function(var,
           dt,
           subdir,
           y_axis_nam = var,
           filename = var,
           mltp_factor = 1) {
    dir.create(output_dir(subdir), FALSE, TRUE)
    gg <- ggplot(
      dt[variable == var],
      aes(
        x = year,
        y =    `50.0%` * mltp_factor,
        ymin = `2.5%` * mltp_factor,
        ymax = `97.5%` * mltp_factor,
        col = type,
        fill = type
      )
    ) +
      geom_point(size = 1,
                 alpha = 5 / 5,
                 show.legend = FALSE) +
      geom_line(size = 1, alpha = 5 / 5) +
      geom_ribbon(alpha = 1 / 5,
                  linetype = 0,
                  show.legend = FALSE) +
      scale_x_continuous(name = "Year", guide = guide_axis(n.dodge = 2)) +
      scale_y_continuous(
        name = fifelse(
          mltp_factor == 1,
          y_axis_nam,
          paste0(y_axis_nam, " per ", mltp_factor)
        ),
        expand = c(0, 0),
        limits = c(0, NA)
      ) +
      ggtitle(filename) +
      scale_color_viridis_d() +
      scale_fill_viridis_d() +
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1
        )
      )

    nam <- ""
    if (uniqueN(dt$sex) > 1L &&
        uniqueN(dt$agegrp20) == 1L &&
        unique(dt$agegrp20) == "All" &&
        uniqueN(dt$qimd) == 1L &&
        unique(dt$qimd) == "All") {
      gg <- gg + facet_grid(sex ~ ., scales = "free")
      nam <- "_sex"
    } else if (uniqueN(dt$sex) == 1L &&
               unique(dt$sex) == "All" &&
               uniqueN(dt$agegrp20) > 1L &&
               # unique(dt$agegrp20) == "All" &&
               uniqueN(dt$qimd) == 1L &&
               unique(dt$qimd) == "All") {
      gg <- gg + facet_grid(agegrp20 ~ ., scales = "free")
      nam <- "_agegroup"
    } else if (uniqueN(dt$sex) == 1L &&
               unique(dt$sex) == "All" &&
               uniqueN(dt$agegrp20) == 1L &&
               unique(dt$agegrp20) == "All" &&
               uniqueN(dt$qimd) > 1L # &&
               # unique(dt$qimd) == "All"
               ) {
               gg <- gg + facet_grid(qimd ~ ., scales = "free")
               nam <- "_qimd"
  } else if (uniqueN(dt$sex) > 1L &&
             uniqueN(dt$agegrp20) > 1L &&
             uniqueN(dt$qimd) > 1L) {
    gg <- gg + facet_grid(agegrp20 ~ sex + qimd, scales = "free")
    nam <- "_agegroup_sex_qimd"
  }


cowplot::ggsave2(
  filename = paste0(filename, nam, ".png"),
  gg,
  height = 9,
  width = 16,
  units = "cm",
  dpi = 300,
  scale = 2.5,
  path = output_dir(subdir)
)
}

# ******************************************************************
# IMPACT NCD IMPACTncd_Engl
# ******************************************************************

filenames <- c("val_incd_output.csv", "val_mrtl_output.csv")
if (XPS)
  filenames <- c(filenames, "val_xps_output_post.csv")

# Get trends ----
if (all(file.exists(output_dir(filenames)))) {
  if (XPS)  out_xps  <- fread(output_dir(filenames[[3]]))
  out_incd <- fread(output_dir(filenames[[1]]))
  out_mrtl <- fread(output_dir(filenames[[2]]))

} else {
  # delete old
  #   SynthPop$
  #   new(0, design)$
  # delete_synthpop(to_delete) #$
  # write_synthpop(1:design$sim_prm$iteration_n)$
  # count_synthpop()



  if (Sys.info()["sysname"] == "Windows") {
    cl <-
      makeCluster(design$sim_prm$clusternumber) # used for clustering. Windows compatible
    registerDoParallel(cl)
  } else {
    registerDoParallel(design$sim_prm$clusternumber)  # used for forking. Only Linux/OSX compatible
  }

  out <- foreach(
    mc_iter =  1:(design$sim_prm$iteration_n_final *
        design$sim_prm$n_synthpop_aggregation),
    # .combine = rbind,
    .inorder = FALSE,
    .options.multicore = list(preschedule = FALSE),
    .verbose = TRUE,
    .packages = c("data.table", "IMPACTncdEnglmisc")
    # .export = ls(),
    # .noexport = c("time_mark")
  ) %dopar%
    {

      if (is.null(to_delete)) {
        SynthPop$ # delete existing file but not primer
        new(0, design)$
        delete_synthpop(mc_iter)
      }

      POP <- SynthPop$new(mc_iter, design)

      nam <- grep("^rn_", names(POP$pop), value = TRUE)
      POP$pop[, (nam) := NULL]
      POP$pop <- POP[!(dead), ]

      if (XPS) export_xps(mc_iter, POP$pop, TRUE, "val_xps_output_post.csv")
      export_incd(mc_iter, POP$pop, TRUE)
      export_mrtl(mc_iter, POP$pop, TRUE)

      export_all_incd(mc_iter, POP$pop, TRUE)
      export_all_prvl(mc_iter, POP$pop, TRUE)
      rm(POP)
      NULL
    }

  if (exists("cl")) stopCluster(cl)

  if (XPS)
    out_xps  <- fread(output_dir(filenames[[3]]))
  out_incd <- fread(output_dir(filenames[[1]]))
  out_mrtl <- fread(output_dir(filenames[[2]]))

  if (design$sim_prm$n_synthpop_aggregation > 1L) {
    out_incd[, mc := as.integer(ceiling(mc / design$sim_prm$n_synthpop_aggregation))]
    out_incd <-
      out_incd[, lapply(.SD, sum), keyby = .(year, sex, agegrp20, qimd, mc)]
    fwrite(out_incd, output_dir(filenames[[1]]))
    out_mrtl[, mc := as.integer(ceiling(mc / design$sim_prm$n_synthpop_aggregation))]
    out_mrtl <-
      out_mrtl[, lapply(.SD, sum), keyby = .(year, sex, agegrp20, qimd, mc)]
    fwrite(out_mrtl, output_dir(filenames[[2]]))
  }
}



# XPS ----
if (XPS) {
  out_xps[, fruit_curr_xps := fruit_curr_xps / 80]
  out_xps[, veg_curr_xps := veg_curr_xps / 80]
  out_xps <-
    melt(out_xps, c("year",
                    "agegrp20", "sex", "qimd",
                    "mc"))
  out_xps <- summarise_distr(out_xps, pr, FALSE, FALSE)
  out_xps[, variable := gsub("_curr_xps$", "", variable)]
  out_xps[, type := "IMPACTncd_Engl"]

  if (file.exists("../IMPACTncd_Engl/preparatory_work/HSE_ts.fst")) {
    HSE_ts <-
      read_fst("../IMPACTncd_Engl/preparatory_work/HSE_ts.fst", as.data.table = TRUE)
  } else {
    source("../IMPACTncd_Engl/preparatory_work/preprocess_HSE.R", local = TRUE)
  }
  HSE <-
    HSE_ts[!is.na(age) &
             !is.na(sex) & !is.na(qimd)  & !is.na(ethnicity) &
             !is.na(sha) &
             between(age, design$sim_prm$ageL, design$sim_prm$ageH)]
  HSE[, alcohol := NULL]
  HSE[, year := year + 2000L]
  HSE[, smok_never := fifelse(smok_status == "1", 1, 0)]
  HSE[, smok_active := fifelse(smok_status == "4", 1, 0)]
  HSE[, agegrp20 := NULL]
  to_agegrp(HSE, 20L, 89L, "age", "agegrp20", to_factor = TRUE)
  setnames(
    HSE,
    c("frtpor", "vegpor", "dm", "bp_med", "totalwu"),
    c("fruit", "veg", "t2dm_prvl", "bpmed", "alcohol")
  )

  HSE[, smok_dur := 0]
  HSE[smok_status == "2", smok_dur := smok_dur_ex]
  HSE[smok_status == "3", smok_dur := smok_dur_ex]
  HSE[smok_status == "4", smok_dur := clamp(age - smok_init_age, 0L, 100L)]
  HSE[, smok_cig := 0]
  HSE[smok_status == "2", smok_cig := 1L]
  HSE[smok_status == "3", smok_cig := clamp(smok_cig_ex, 1L, 100L)]
  HSE[smok_status == "4", smok_cig := clamp(smok_cig_curr, 1L, 100L)]
  setnafill(HSE, "c", 0L, cols = "smok_quit_yrs")

  HSE[, ets := as.numeric(as.character(ets))]
  HSE[, t2dm_prvl := as.numeric(as.character(t2dm_prvl))]
  HSE[, bpmed := as.numeric(as.character(bpmed))]
  HSE[year == 2014, `:=`(fruit = NA, veg = NA)]
  xps <- unique(out_xps$variable)
  xps <- xps[-which(xps %in% c("af_dgn", "af_prvl", "ckd_prvl"))]
  xps_wt_int <- c(
    "active_days",
    "fruit",
    "veg",
    "smok_quit_yrs",
    "smok_dur",
    "smok_cig",
    "ets",
    "alcohol",
    "smok_never",
    "smok_active"
  )
  xps_wt_nurse <- c("bmi", "sbp", "bpmed", "statin_px")
  xps_wt_blood <- c("tchol", "t2dm_prvl")

  # HSE[, lapply(.SD, class), .SDcols = xps]
  out_xps2 <- groupingsets(
    HSE,
    j = lapply(.SD, weighted.mean, wt_int, na.rm = TRUE),
    by = c("year", "sex", "agegrp20", "qimd"),
    .SDcols = xps_wt_int,
    sets = list(
      c("year", "sex", "agegrp20", "qimd"),
      c("year", "sex"),
      c("year", "agegrp20"),
      c("year", "qimd")
    )
  )
  out_xps3 <- groupingsets(
    HSE,
    j = lapply(.SD, weighted.mean, wt_nurse, na.rm = TRUE),
    by = c("year", "sex", "agegrp20", "qimd"),
    .SDcols = xps_wt_nurse,
    sets = list(
      c("year", "sex", "agegrp20", "qimd"),
      c("year", "sex"),
      c("year", "agegrp20"),
      c("year", "qimd")
    )
  )
  out_xps4 <- groupingsets(
    HSE,
    j = lapply(.SD, weighted.mean, wt_blood, na.rm = TRUE),
    by = c("year", "sex", "agegrp20", "qimd"),
    .SDcols = xps_wt_blood,
    sets = list(
      c("year", "sex", "agegrp20", "qimd"),
      c("year", "sex"),
      c("year", "agegrp20"),
      c("year", "qimd")
    )
  )
  absorb_dt(out_xps2, out_xps3)
  absorb_dt(out_xps2, out_xps4)
  for (j in c("agegrp20", "sex", "qimd"))
    set(out_xps2, which(is.na(out_xps2[[j]])), j, "All")
  out_xps2 <-
    melt(out_xps2, c("year",
                     "agegrp20", "sex", "qimd"))
  out_xps2[, type := "HSE"]
  setnames(out_xps2, "value", "50.0%")
  xpsval <- rbind(out_xps, out_xps2, fill = TRUE)
  xpsval[, `:=` (
    sex = factor(sex),
    qimd = factor(
      qimd,
      levels = c("1 most deprived", "2", "3", "4", "5 least deprived", "All"),
      labels = c("1 most\ndeprived", "2", "3", "4", "5 least\ndeprived", "All")
    ),
    agegrp20 = factor(agegrp20)
  )]
  xpsval[, `2.5%`  := `50.0%`] # TODO better fix for UCI, LCI
  xpsval[, `97.5%` := `50.0%`] # TODO better fix for UCI, LCI

  for (nam in xps) {
    future({
      std_graph(nam,
                xpsval[sex != "All" &
                         agegrp20 == "All" & qimd == "All"],
                "validation_xps")
      std_graph(nam,
                xpsval[sex == "All" &
                         agegrp20 != "All" & qimd == "All"],
                "validation_xps")
      std_graph(nam,
                xpsval[sex == "All" &
                         agegrp20 == "All" & qimd != "All"],
                "validation_xps")
      std_graph(nam,
                xpsval[sex != "All" &
                         agegrp20 != "All" & qimd != "All"],
                "validation_xps")
    })
  }

  # std_graph("smok_never",
  #   xpsval[sex != "All" & agegrp20 == "All" & qimd == "All"],
  #   "validation_xps",
  #   y_axis_nam = "Never smokers proportion",
  #   filename = "Never smokers")
  # std_graph("smok_never",
  #   xpsval[sex == "All" & agegrp20 != "All" & qimd == "All"],
  #   "validation_xps",
  #   y_axis_nam = "Never smokers proportion",
  #   filename = "Never smokers")
  # std_graph("smok_never",
  #   xpsval[sex == "All" & agegrp20 == "All" & qimd != "All"],
  #   "validation_xps",
  #   y_axis_nam = "Never smokers proportion",
  #   filename = "Never smokers")
  # std_graph("smok_never",
  #   xpsval[sex != "All" & agegrp20 != "All" & qimd != "All"],
  #   "validation_xps",
  #   y_axis_nam = "Never smokers proportion",
  #   filename = "Never smokers")
}

# Mortality ---------------------------------------------------------------
# out_mrtl <- fread(output_dir(filenames[[2]]))

nam <- grep("^deaths_", names(out_mrtl), value = TRUE)
namMX <- gsub("^deaths_", "Mx_", nam)
nam3 <- gsub("^deaths_", "", nam)

for (i in nam)
  set(out_mrtl, NULL, namMX[which(nam == i)], out_mrtl[[i]] / out_mrtl$pops)
out_mrtl <-
  melt(out_mrtl, c("year",
                   "agegrp20", "sex", "qimd",
                   "mc"))
out_mrtl <- summarise_distr(out_mrtl, pr, FALSE, FALSE)
out_mrtl[, type := "IMPACTncd_Engl"]
out_mrtl <- out_mrtl[year > 2012L, ]


for (i in seq_along(nam)) {
  print(nam3[i])
  tt <- get_ons_mrtl(nam3[i], "rate", 20L)
  tt[, `:=` (pops = NULL,
             variable = namMX[i],
             type = "ONS")]
  setnames(tt, namMX[i], "50.0%")
  tt <- rbind(out_mrtl[variable == namMX[i]], tt, fill = TRUE)
  tt[, `:=` (
    sex = factor(sex),
    qimd = factor(
      qimd,
      levels = c("1 most deprived", "2", "3", "4", "5 least deprived", "All"),
      labels = c("1 most\ndeprived", "2", "3", "4", "5 least\ndeprived", "All")
    ),
    agegrp20 = factor(agegrp20)
  )]
  future({
    std_graph(namMX[i],
              tt[sex != "All" & agegrp20 != "All" & qimd != "All"],
              "validation_mrtl", mltp_factor = 1e5)
  })
  future({
    std_graph(namMX[i],
              tt[sex != "All" & agegrp20 == "All" & qimd == "All"],
              "validation_mrtl", mltp_factor = 1e5)
  })
  future({
    std_graph(namMX[i],
              tt[sex == "All" & agegrp20 != "All" & qimd == "All"],
              "validation_mrtl", mltp_factor = 1e5)
  })
  future({
    std_graph(namMX[i],
              tt[sex == "All" & agegrp20 == "All" & qimd != "All"],
              "validation_mrtl", mltp_factor = 1e5)
  })
}

# Incidence ---------------------------------------
# out_incd <- fread(output_dir(filenames[[1]]))

nam <- grep("^cases_", names(out_incd), value = TRUE)
namrate <- gsub("^cases_", "rate_", nam)
nam3 <- gsub("^cases_", "", nam)

for (i in nam)
  set(out_incd, NULL, namrate[which(nam == i)], out_incd[[i]] / out_incd$pops)
out_incd <-
  melt(out_incd, c("year",
                   "agegrp20", "sex", "qimd",
                   "mc"))
out_incd <- summarise_distr(out_incd, pr, FALSE, FALSE)
out_incd[, type := "IMPACTncd_Engl"]
out_incd <- out_incd[year > 2012,] # rate for 2003 is > 0

for (i in seq_along(nam)) {
  print(i)
  print(nam3[i])
  tt <- get_ons_incd(nam3[i], "rate", 20L)
  tt[, `:=` (pops = NULL,
             variable = namrate[i],
             type = "ONS")]
  setnames(tt, namrate[i], "50.0%")
  tt <- rbind(out_incd[variable == namrate[i]], tt, fill = TRUE)
  tt[, `:=` (
    sex = factor(sex),
    qimd = factor(
      qimd,
      levels = c("1 most deprived", "2", "3", "4", "5 least deprived", "All"),
      labels = c("1 most\ndeprived", "2", "3", "4", "5 least\ndeprived", "All")
    ),
    agegrp20 = factor(agegrp20)
  )]
  future({
    std_graph(namrate[i],
              tt[sex != "All" & agegrp20 != "All" & qimd != "All"],
              "validation_incd", mltp_factor = 1e5)
  })
  future({
    std_graph(namrate[i],
              tt[sex != "All" & agegrp20 == "All" & qimd == "All"],
              "validation_incd", mltp_factor = 1e5)
  })
  future({
    std_graph(namrate[i],
              tt[sex == "All" & agegrp20 != "All" & qimd == "All"],
              "validation_incd", mltp_factor = 1e5)
  })
  future({
    std_graph(namrate[i],
              tt[sex == "All" & agegrp20 == "All" & qimd != "All"],
              "validation_incd", mltp_factor = 1e5)
  })
}

# Incidence (ALL) ------------------------------------
out_incd <- fread(output_dir("val_all_incd_output.csv"))

nam <- grep("^cases_", names(out_incd), value = TRUE)
namrate <- gsub("^cases_", "rate_", nam)
nam3 <- gsub("^cases_", "", nam)

for (i in nam)
  set(out_incd, NULL, namrate[which(nam == i)], out_incd[[i]] / out_incd$pops)
out_incd <-
  melt(out_incd, c("year",
                   "agegrp20", "sex", "qimd",
                   "mc"))
out_incd <- summarise_distr(out_incd, pr, FALSE, FALSE)
out_incd[, type := "IMPACTncd_Engl"]



for (i in seq_along(nam)) {
  print(nam3[i])
  tt <- out_incd[variable == namrate[i]]
  tt[, `:=` (
    sex = factor(sex),
    qimd = factor(
      qimd,
      levels = c("1 most deprived", "2", "3", "4", "5 least deprived", "All"),
      labels = c("1 most\ndeprived", "2", "3", "4", "5 least\ndeprived", "All")
    ),
    agegrp20 = factor(agegrp20)
  )]
  future({
    std_graph(namrate[i],
              tt[sex != "All" & agegrp20 != "All" & qimd != "All"],
              "validation_incd_all", mltp_factor = 1e5)
  })
  future({
    std_graph(namrate[i],
              tt[sex != "All" & agegrp20 == "All" & qimd == "All"],
              "validation_incd_all", mltp_factor = 1e5)
  })
  future({
    std_graph(namrate[i],
              tt[sex == "All" & agegrp20 != "All" & qimd == "All"],
              "validation_incd_all", mltp_factor = 1e5)
  })
  future({
    std_graph(namrate[i],
              tt[sex == "All" & agegrp20 == "All" & qimd != "All"],
              "validation_incd_all", mltp_factor = 1e5)
  })
}

# Prevalence (ALL) ------------------------------
out_prvl <- fread(output_dir("val_all_prvl_output.csv"))

nam <- grep("^cases_", names(out_prvl), value = TRUE)
namrate <- gsub("^cases_", "rate_", nam)
nam3 <- gsub("^cases_", "", nam)

for (i in nam)
  set(out_prvl, NULL, namrate[which(nam == i)], out_prvl[[i]] / out_prvl$pops)
out_prvl <-
  melt(out_prvl, c("year",
                   "agegrp20", "sex", "qimd",
                   "mc"))
out_prvl <- summarise_distr(out_prvl, pr, FALSE, FALSE)
out_prvl[, type := "IMPACTncd_Engl"]



for (i in seq_along(nam)) {
  print(nam3[i])
  tt <- out_prvl[variable == namrate[i]]
  tt[, `:=` (
    sex = factor(sex),
    qimd = factor(
      qimd,
      levels = c("1 most deprived", "2", "3", "4", "5 least deprived", "All"),
      labels = c("1 most\ndeprived", "2", "3", "4", "5 least\ndeprived", "All")
    ),
    agegrp20 = factor(agegrp20)
  )]
  future({
    std_graph(namrate[i],
              tt[sex != "All" & agegrp20 != "All" & qimd != "All"],
              "validation_prvl_all", mltp_factor = 1e5)
  })
  future({
    std_graph(namrate[i],
              tt[sex != "All" & agegrp20 == "All" & qimd == "All"],
              "validation_prvl_all", mltp_factor = 1e5)
  })
  future({
    std_graph(namrate[i],
              tt[sex == "All" & agegrp20 != "All" & qimd == "All"],
              "validation_prvl_all", mltp_factor = 1e5)
  })
  future({
    std_graph(namrate[i],
              tt[sex == "All" & agegrp20 == "All" & qimd != "All"],
              "validation_prvl_all", mltp_factor = 1e5)
  })
}

