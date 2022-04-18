input_path <-
  function(x = character(0))
    paste0("/mnt/", Sys.info()[["user"]],
           "/UoL/CPRD2021/alhead_cprd/HFdemandmodeldata/processed_data/", x)
output_path <-
  function(x = character(0))
    paste0("/mnt/", Sys.info()[["user"]],
           "/UoL/CPRD2021/epi_models/", x)

harmonise <- function(dt) {
  setnames(dt, c("gender", "imd10", "region"), c("sex", "dimd", "sha"))
  dt[, `:=` (ethnicity = factor(ethnicity,
                                levels = c("white", "indian", "pakistani", "bangladeshi",
                                           "oth_asian", "bl_carib", "bl_afric", "chinese", "other"),
                                labels = c("white", "indian", "pakistani", "bangladeshi",
                                           "other asian", "black caribbean", "black african",
                                           "chinese", "other")),
             sha = factor(sha,
                          levels = c("North East", "North West", "Yorkshire And The Humber",
                                     "East Midlands", "West Midlands", "East of England",
                                     "London", "South East Coast", "South Central", "South West"),
                          labels = c("North East", "North West", "Yorkshire and the Humber",
                                     "East Midlands", "West Midlands", "East of England",
                                     "London", "South East Coast", "South Central", "South West")),
             sex = factor(sex, levels = c("M", "F"), labels = c("men", "women"))
  )
  ]
}

validate_plots <- function(dt, y, mod_max, suffix, disnm, strata) {
  y_pred <- predictAll(mod_max, dt, data = dt)$mu

  dt <- cbind(y, y_pred, dt)

  dir.create(output_path(paste0("validation/", disnm)), recursive = TRUE, showWarnings = FALSE)

  for(i in seq_along(strata)) {
    png(output_path(paste0(
      "validation/", disnm, "/", strata[[i]], suffix, ".png"
    )))

    dt[, .(incd = sum(V1) / sum(V1 + V2)), keyby = eval(strata[[i]])][, matplot(
      get(strata[[i]]),
      incd,
      pch = "o",
      xlab = strata[[i]],
      ylab = gsub("^_", "", suffix),
      main = strata[[i]],
      ylim = c(min(incd) * 0.8, max(incd) * 1.2)
    )]
    dt[, .(incd = sum(y_pred * (V1 + V2)) / sum(V1 + V2)), keyby = eval(strata[[i]])][, matlines(get(strata[[i]]), incd, col = "red")] # red is the pred
    legend(
      x = "bottomright",
      legend = c("CPRD", "GAMLSS"),
      lty = c(1, 1),
      col = c("black", "red"),
      lwd = 2
    )
    dev.off()
  }
}

mod_sel <- function(m, ncpus = 10L) stepGAIC(m, direction = "backward", parallel = "multicore", ncpus = ncpus)
