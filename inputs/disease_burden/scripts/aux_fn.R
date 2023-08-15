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



# For making RRs with CIs from a logit model. Base on the package logisticRR
RRfn <- function (fit, # a glm binomial("logit") model object
                  basecov = 0, fixcov = NULL,
                  data # a data.table
)
{
  setDF(data)
  #tmp <- strsplit(as.character(formula(fit))[3], "[+]")
  #if more than 500 characters, as.character inserts a new line - need to remove this
  tmp <- strsplit(gsub('[\n]', '', as.character(formula(fit))[3] ), "[+]")
  varnames <- gsub(" ", "", tmp[[1]])
  if (class(data[, names(data) == varnames[1]]) == "factor")
    return("Please use nominalRR")
  p <- length(varnames) - 1
  if (p == 0) {
    newfixcov = NULL
  }
  else if (p > 0) {
    newfixcov <- t(as.matrix(rep(0, p)))
    subdat = as.data.frame(data[, which(names(data) %in%
                                          varnames[-1])])
    # subdat = data[, varnames[-1], with = FALSE]
    tmp <- which(apply(subdat, 2, class) != "numeric")
    for (q in 1:p) {
      if (class(subdat[, q]) == "factor") {
        newfixcov[q] <- levels(as.factor(subdat[, q]))[1]
      }
      else {
        newfixcov[q] <- min(subdat[, q])
      }
    }
    newfixcov <- as.data.frame(newfixcov)
    names(newfixcov) = names(data)[which(names(data) %in%
                                           varnames[-1])]
  }
  if (sum(names(fixcov) %in% names(newfixcov)) > 0) {
    tmpind <- which(names(newfixcov) %in% names(fixcov))
    for (j in 1:length(tmpind)) {
      newfixcov[tmpind[j]] = eval(parse(text = paste0("fixcov$",
                                                      names(newfixcov[tmpind])[j])))
    }
  }
  fixcov = newfixcov
  expose.cov <- data.frame(basecov + 1)
  names(expose.cov) <- varnames[1]
  unexpose.cov <- data.frame(basecov)
  names(unexpose.cov) <- varnames[1]
  if (length(fixcov) > 0 & length(names(fixcov)) > 0 & length(fixcov) ==
      length(varnames) - 1) {
    expose.cov <- cbind(expose.cov, fixcov)
    unexpose.cov <- cbind(unexpose.cov, fixcov)
  }
  else if (length(names(fixcov)) == 0 & length(fixcov) > 0) {
    expose.cov <- cbind(expose.cov, fixcov)
    names(expose.cov)[2:length(expose.cov)] = varnames[2:length(varnames)]
    unexpose.cov <- cbind(unexpose.cov, fixcov)
    names(unexpose.cov)[2:length(unexpose.cov)] = varnames[2:length(varnames)]
  }
  else if (p > 0) {
    return("Invalid data frame for confounders")
  }
  for (i in 1:ncol(expose.cov)) {
    if (class(data[, names(data) == names(expose.cov)[i]]) !=
        "factor") {
      expose.cov[, i] <- as.numeric(expose.cov[, i])
      unexpose.cov[, i] <- as.numeric(unexpose.cov[, i])
    }
  }
  betas <- coefficients(fit)
  exposed <- exp(-predict(fit, expose.cov, type = "link"))
  unexposed <- exp(-predict(fit, unexpose.cov, type = "link"))
  RR <- (1 + unexposed)/(1 + exposed)
  n.par <- length(betas)
  B.vec <- rep(0, n.par)
  B.vec[1] <- (-exposed + unexposed)/(1 + exposed)^2
  B.vec[2] = (-(basecov + 1) * exposed * (1 + unexposed) +
                basecov * unexposed * (1 + exposed))/(1 + exposed)^2
  if (n.par > 2) {
    for (j in 3:n.par) {
      if (names(coefficients(fit))[j] %in% names(fixcov)) {
        tmp <- which(names(fixcov) %in% names(coefficients(fit))[j])
        B.vec[j] <- as.numeric(as.character(fixcov[tmp])) *
          (unexposed - exposed)/(1 + exposed)^2
      }
      else if (sum(startsWith(names(coefficients(fit))[j],
                              names(fixcov))) > 0) {
        tmp <- which(startsWith(names(coefficients(fit))[j],
                                names(fixcov)))
        if (gsub(names(fixcov)[tmp], "", names(coefficients(fit))[j]) ==
            as.character(fixcov[, tmp])) {
          B.vec[j] <- 1 * (unexposed - exposed)/(1 +
                                                   exposed)^2
        }
        else {
          B.vec[j] <- 0 * (unexposed - exposed)/(1 +
                                                   exposed)^2
        }
      }
    }
  }
  cov.mat <- summary(fit)$cov.unscaled
  deltavar <- 0
  for (i in 1:n.par) {
    for (j in 1:n.par) {
      deltavar <- deltavar + cov.mat[i, j] * B.vec[i] *
        B.vec[j]
    }
  }
  return(list(RR = RR, UCI_RR = RR + qnorm(0.975) * sqrt(deltavar), LCI_RR = RR - qnorm(0.975) * sqrt(deltavar), delta.var = deltavar, fix.cov = fixcov))
}








#For writing the .csvy files for disease dependencies
write_xps_tmplte_file <-
  function(dt, j, dpnds, disnm, type, file_path ) {
    exposure <- dpnds[j]
    write_xps_prm_file = function(dt, metadata, file_path) {
      y <- paste0("---\n", yaml::as.yaml(metadata), "---\n")
      con <- textConnection(y)
      on.exit(close(con))
      m <- readLines(con)
      y <- paste0("", m[-length(m)], collapse = "\n")
      y <- c(y, "\n")
      cat(y, file = file_path)
      fwrite(x = dt, file = file_path, append = TRUE, col.names = TRUE)
    }

    file_path <- normalizePath(file_path, mustWork = FALSE)
    metadata <- list(
      "xps_name"     = paste0(exposure, type), 
      "outcome"       = disnm,
      "lag"          = 1L,
      "distribution" =  "lognormal",
      "source"        = "CPRD",
      "notes"         = paste0("other covariate conditions: ", paste(dpnds[!dpnds %in% exposure], collapse = ","))
    )

    effect <- CJ(
      agegroup = factor(CKutils::agegrp_name(min_age = 0, max_age = 99, grp_width = 5)),
      sex = c("men", "women"),
      dimd = factor(c("1 most deprived", "2", "3", "4", "5", "6", "7", "8", "9", "10 least deprived"),
                    levels = c("1 most deprived", "2", "3", "4", "5", "6", "7", "8", "9", "10 least deprived"))
    )
    
    tmp <- dt[dpnd_on == exposure,
                   .(sex, agegroup, dimd,  rr = round(rr, digits = 2), ci_rr =round(ci_rr, digits = 2))]
    effect <- merge(effect, tmp, by = c("agegroup","sex", "dimd"), all = T)
    setkey(effect, sex, dimd, agegroup)
    
    effect[is.na(rr), `:=` (rr = 1, ci_rr =1)] #putting in RR = 1 for the missing age groups 
    
    write_xps_prm_file(effect, metadata, file_path)
  }
# 
# write_xps_tmplte_file <-
#   function(exposure, disnm, other, type, file_path ) {
#     write_xps_prm_file = function(dt, metadata, file_path) {
#       y <- paste0("---\n", yaml::as.yaml(metadata), "---\n")
#       con <- textConnection(y)
#       on.exit(close(con))
#       m <- readLines(con)
#       y <- paste0("", m[-length(m)], collapse = "\n")
#       y <- c(y, "\n")
#       cat(y, file = file_path)
#       fwrite(x = dt, file = file_path, append = TRUE, col.names = TRUE)
#     }
#     
#     file_path <- normalizePath(file_path, mustWork = FALSE)
#     metadata <- list(
#       "xps_name"     = paste0(exposure, type), 
#       "outcome"       = disnm,
#       "lag"          = 1L,
#       "distribution" =  "lognormal",
#       "source"        = "CPRD",
#       "notes"         = paste0("other covariate conditions: ", paste(other[!other %in% exposure], collapse = ","))
#     )
#     
#     effect <- CJ(
#       agegroup = factor(CKutils::agegrp_name(min_age = 0, max_age = 99, grp_width = 5)),
#       sex = c("men", "women"),
#       dimd = factor(c("1 most deprived", "2", "3", "4", "5", "6", "7", "8", "9", "10 least deprived"),
#                     levels = c("1 most deprived", "2", "3", "4", "5", "6", "7", "8", "9", "10 least deprived"))
#     )
#     
#     tmp <- results[dpnd_on == exposure,
#                    .(sex, agegroup, dimd,  rr = round(rr, digits = 2), ci_rr =round(ci_rr, digits = 2))]
#     effect <- merge(effect, tmp, by = c("agegroup","sex", "dimd"), all = T)
#     setkey(effect, sex, dimd, agegroup)
#     
#     effect[is.na(rr), `:=` (rr = 1, ci_rr =1)] #putting in RR = 1 for the missing age groups 
#     
#     write_xps_prm_file(effect, metadata, file_path)
#   }
