library(data.table)
library(fst)
library(dfoptim)
library(CKutils)
library(doParallel)
library(parallelly)
library(foreach)
validation <- TRUE
# used for forking. Only Linux/OSX compatible
# registerDoParallel(15L)

cl <-
    makeClusterPSOCK(
        15L, # cores
        dryrun = FALSE,
        quiet = FALSE,
        rscript_startup = quote(local({
            library(data.table)
            library(fst)
            library(dfoptim)
            library(CKutils)
        })),
        rscript_args = c(
            "--no-init-file",
            "--no-site-file",
            "--no-environ"
        ),
        setup_strategy = "parallel"
    ) # used for clustering. Windows compatible
registerDoParallel(cl)
on.exit(if (exists("cl")) stopCluster(cl))

poproj <- read_fst("inputs/pop_projections/lad17_proj.fst", as.data.table = TRUE)
bydimd <- read_fst("inputs/pop_estimates_lsoa/lsoa_to_locality_indx.fst", as.data.table = TRUE)
poplsoa <- read_fst("inputs/pop_estimates_lsoa/LSOA_mid_year_population_estimates.fst", as.data.table = TRUE)
poplsoa <- melt(poplsoa, 1:3, variable.name = "age", value.name = "pops")
poplsoa <- poplsoa[, age := as.integer(as.character(age))]
poplsoa[bydimd, on = "LSOA11CD", `:=`(dimd = i.dimd, LAD17CD = i.LAD17CD)]
poplsoa <- poplsoa[, .(pops = sum(pops)), by = .(year, age, sex, LAD17CD, dimd)]
# poplsoa[, sum(pops), by = .(year, age, sex, LAD17CD)][V1 == 0, ] # what to do with 0s?
poplsoa[, prop := pops / sum(pops), by = .(year, age, sex, LAD17CD)]
poplsoa[, sum(prop), by = .(year, age, sex, LAD17CD)][, table(V1, useNA = "ifany")] # Check if prop sums to 1
poplsoa <- dcast(poplsoa, year + age + sex + LAD17CD ~ dimd, value.var = "prop") # distribution of dimd in lad17 by year, age, sex
setnafill(poplsoa, fill = 0, cols = c("1 most deprived", as.character(2:9), "10 least deprived"))
absorb_dt(poproj, poplsoa)
rm(bydimd, poplsoa)
poprojOriginal <- copy(poproj)
if (validation) poproj[year > min(year), c("1 most deprived", as.character(2:9), "10 least deprived") := NA]



f <- function(x = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1), ds = p) { # function to be minimised
    # x is a vector of growth rates by dimd
    err <- ds$popsYnext - (ds$pops * ds$`1 most deprived` * x[1] +
        ds$pops * ds$`2` * x[2] +
        ds$pops * ds$`3` * x[3] +
        ds$pops * ds$`4` * x[4] +
        ds$pops * ds$`5` * x[5] +
        ds$pops * ds$`6` * x[6] +
        ds$pops * ds$`7` * x[7] +
        ds$pops * ds$`8` * x[8] +
        ds$pops * ds$`9` * x[9] +
        ds$pops * ds$`10 least deprived` * x[10])
    sum(err^2)
}

ff <- function(x = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1), ds = p, growthRates = r$par, x_init) {
    # x is the dimd distribution of the next year that we try to estimate
    # x_init is the dimd distribution of the current year for starting values
    err <- ds$popsYnext - (ds$pops * growthRates[1] * x[1] +
        ds$pops * growthRates[2] * x[2] +
        ds$pops * growthRates[3] * x[3] +
        ds$pops * growthRates[4] * x[4] +
        ds$pops * growthRates[5] * x[5] +
        ds$pops * growthRates[6] * x[6] +
        ds$pops * growthRates[7] * x[7] +
        ds$pops * growthRates[8] * x[8] +
        ds$pops * growthRates[9] * x[9] +
        ds$pops * growthRates[10] * x[10])
    sum(err^2) +
    sum((x - x_init)^2) * 1e5 + # trick to 'force' x to be close to init values
        (sum(x) - 1)^2 * 1e7 # trick to 'force' sum(x) = 1
}

for (year_ in sort(unique(poproj$year))) {
    print(paste0("Year: ", year_, " ", Sys.time()))
    if (year_ == max(poproj$year)) next # skip the last year because we don't have popsYnext
    for (age_ in sort(unique(poproj$age))) {
        print(paste0("Age: ", age_, " ", Sys.time()))
        for (sex_ in sort(levels(poproj$sex))) {
            print(paste0("Sex: ", sex_))

            # find growth factor by dimd for all LADs (r$par)
            p <- poproj[year == year_ & age == age_ & sex == sex_] # p is the year, age slice of poproj
            if (age_ == 100L) {
                pp <- poproj[year == year_ + 1L & age == age_ & sex == sex_] # pp is the year+1, age slice of poplsoa
            } else {
                pp <- poproj[year == year_ + 1L & age == age_ + 1L & sex == sex_] # pp is the year+1, age+1 slice of poplsoa
            }
            p[pp, on = "LAD17CD", popsYnext := i.pops]
            r <- nmkb(rep(1, 10), f,
                ds = p,
                lower = 0, upper = Inf, control = list(maxfeval = 1e6, restarts.max = 9, trace = FALSE)
            )
            if (r$convergence != 0) stop(paste0(r$convergence, "\n", r$message))

            # apply the growth factor to dimd distribution of current year to get dimd distribution of next year, if missing
            if (pp[, all(sapply(.SD, is.na)), .SDcols = 7:16]) { # if distribution of dimd is missing
            # system.time({
                out <- foreach(
                    i = sort(unique(p$LAD17CD)),
                    .combine = rbind,
                    .inorder = FALSE,
                    .options.multicore = list(preschedule = FALSE),
                    .verbose = FALSE,
                    .packages = c("CKutils", "dfoptim", "fst", "data.table"),
                    .export = NULL, # ls(envir = globalenv()),
                    .noexport = c("poprojOriginal", "poproj", "poplsoa", "bydimd")
                ) %dopar% {
                    cn <- transpose(p[LAD17CD == i, 7:16])$V1
                    cn[cn == 0] <- 1e-10 # replace 0 values with a small number to enable convergence
                    cn <- cn / sum(cn) # normalise
                    cn <- rep(0.1, 10) # TESTING!!!
                    rr <- nmkb(cn, ff,
                        ds = p[LAD17CD == i], growthRates = r$par, x_init = transpose(p[LAD17CD == i, 7:16])$V1,
                        lower = 0, upper = 1, control = list(maxfeval = 1e6, restarts.max = 9, trace = FALSE)
                    )
                    # print(paste0("Sum 1: ", sum(rr$par)))
                    # print(paste0("Diff: ", paste(round(rr$par - transpose(p[LAD17CD == i, 7:16])$V1, 2), collapse = ", ")))
                    # print(paste0("Pop : ", round(p[LAD17CD == i]$popsYnext - (p[LAD17CD == i]$pops * r$par[1] * rr$par[1] +
                    #                      p[LAD17CD == i]$pops * r$par[2] * rr$par[2] +
                    #                      p[LAD17CD == i]$pops * r$par[3] * rr$par[3] +
                    #                      p[LAD17CD == i]$pops * r$par[4] * rr$par[4] +
                    #                      p[LAD17CD == i]$pops * r$par[5] * rr$par[5] +
                    #                      p[LAD17CD == i]$pops * r$par[6] * rr$par[6] +
                    #                      p[LAD17CD == i]$pops * r$par[7] * rr$par[7] +
                    #                      p[LAD17CD == i]$pops * r$par[8] * rr$par[8] +
                    #                      p[LAD17CD == i]$pops * r$par[9] * rr$par[9] +
                    #                      p[LAD17CD == i]$pops * r$par[10] *rr$par[10]))))
                    if (rr$convergence != 0) stop(paste0(rr$convergence, "\n", rr$message))
                    if (abs(sum(rr$par) - 1) > 1e-2) print(paste0("sum of dimd is ", sum(rr$par)))
                    if (sum(rr$par) != 1) rr$par <- rr$par / sum(rr$par)
                    data.table(year = year_ + 1L, age = age_ + 1L, sex = sex_, LAD17CD = as.character(i), "1 most deprived" = rr$par[1], "2" = rr$par[2], "3" = rr$par[3], "4" = rr$par[4], "5" = rr$par[5], "6" = rr$par[6], "7" = rr$par[7], "8" = rr$par[8], "9" = rr$par[9], "10 least deprived" = rr$par[10])
                } # end foreach
                # }) # system.time
                poproj[out, on = c("year", "age", "sex", "LAD17CD"), `:=`(
                    `1 most deprived` = `i.1 most deprived`,
                    `2` = `i.2`,
                    `3` = `i.3`,
                    `4` = `i.4`,
                    `5` = `i.5`,
                    `6` = `i.6`,
                    `7` = `i.7`,
                    `8` = `i.8`,
                    `9` = `i.9`,
                    `10 least deprived` = `i.10 least deprived`
                )]

                #   Special case for age 0 wher age0 dimd distribution of next year depends on dimd distribution of current year
                if (age_ == 0L) {
                    # TODO create function to avoid code duplication
                    pp <- poproj[year == year_ + 1L & age == age_ & sex == sex_] # pp is the year+1, age slice of poplsoa
                    p[pp, on = "LAD17CD", popsYnext := i.pops]
                    r <- nmkb(rep(1, 10), f,
                        ds = p,
                        lower = 0, upper = Inf, control = list(maxfeval = 1e6, restarts.max = 9, trace = FALSE)
                    )
                    if (r$convergence != 0) stop(paste0(r$convergence, "\n", r$message))

                    out <- foreach(
                        i = sort(unique(p$LAD17CD)),
                        .combine = rbind,
                        .inorder = FALSE,
                        .options.multicore = list(preschedule = FALSE),
                        .verbose = FALSE,
                        .packages = c("CKutils", "dfoptim", "fst", "data.table"),
                        .export = NULL, # ls(envir = globalenv()),
                        .noexport = c("poprojOriginal", "poproj", "poplsoa", "bydimd")
                    ) %dopar% {
                        cn <- transpose(p[LAD17CD == i, 7:16])$V1
                        cn[cn == 0] <- 1e-10 # replace 0 values with a small number to enable convergence
                        cn <- cn / sum(cn) # normalise
                        cn <- rep(0.1, 10) # TESTING!!!
                        rr <- nmkb(cn, ff,
                            ds = p[LAD17CD == i], growthRates = r$par,  x_init = transpose(p[LAD17CD == i, 7:16])$V1,
                            lower = 0, upper = 1, control = list(maxfeval = 1e6, restarts.max = 9, trace = FALSE)
                        )
                        if (rr$convergence != 0) stop(paste0(rr$convergence, "\n", rr$message))
                        if (abs(sum(rr$par) - 1) > 1e-2) print(paste0("sum of dimd is ", sum(rr$par)))
                        if (sum(rr$par) != 1) rr$par <- rr$par / sum(rr$par)
                        data.table(year = year_ + 1L, age = age_, sex = sex_, LAD17CD = as.character(i), "1 most deprived" = rr$par[1], "2" = rr$par[2], "3" = rr$par[3], "4" = rr$par[4], "5" = rr$par[5], "6" = rr$par[6], "7" = rr$par[7], "8" = rr$par[8], "9" = rr$par[9], "10 least deprived" = rr$par[10])
                    }
                    poproj[out, on = c("year", "age", "sex", "LAD17CD"), `:=`(
                        `1 most deprived` = `i.1 most deprived`,
                        `2` = `i.2`,
                        `3` = `i.3`,
                        `4` = `i.4`,
                        `5` = `i.5`,
                        `6` = `i.6`,
                        `7` = `i.7`,
                        `8` = `i.8`,
                        `9` = `i.9`,
                        `10 least deprived` = `i.10 least deprived`
                    )]
                } # end if age_ == 0
            } # end if distribution of dimd is missing
        } # end sex loop
    } # end age loop
    if (validation) {
        write_fst(poproj, "./inputs/pop_projections/lad17_proj_byDIMD_validation.fst", 100) # if validation
    } else {
        write_fst(poproj, "./inputs/pop_projections/lad17_proj_byDIMD.fst", 100) # if not validation
    }
} # end year loop

if (exists("cl")) stopCluster(cl)

# check validation
if (FALSE) {
library(hexbin)
tt <- read_fst("inputs/pop_projections/lad17_proj_byDIMD_validation.fst", as.data.table = TRUE)
setkey(tt, year, age, sex, LAD17CD)
setkey(poprojOriginal, year, age, sex, LAD17CD)

tt[year == 20, ]

year_ <- 3L
all.equal(tt[year == year_, .SD, .SDcols = 1:6],poprojOriginal[year == year_, .SD, .SDcols = 1:6])
summary(tt[year == year_, .SD, .SDcols = 7:16] - poprojOriginal[year == year_, .SD, .SDcols = 7:16])
summary(tt[year == year_, 7] - poprojOriginal[year == year_, 7])
for (dimd_ in c("1 most deprived", as.character(2:9), "10 least deprived")) {
    plot(
        hexbin(poprojOriginal[age >= 0 & year == year_, pops * get(dimd_)],
            tt[age >= 0 & year == year_, pops * get(dimd_)] - poprojOriginal[age >= 0 & year == year_, pops * get(dimd_)],
            xbins = 200
        ),
        ylab = "bias", xlab = "Observed pop size",
        main = paste0("Year: ", year_, ", DIMD: ", dimd_)
    )
}

dimd_ <- "1 most deprived"
for (year_ in 3:year_) {
   
    plot(
        hexbin(poprojOriginal[age >= 0 & year == year_, pops * get(dimd_)],
            tt[age >= 0 & year == year_, pops * get(dimd_)] - poprojOriginal[age >= 0 & year == year_, pops * get(dimd_)],
            xbins = 200
        ),
        ylab = "bias", xlab = "Observed pop size",
        main = paste0("Year: ", year_, ", DIMD: ", dimd_)
    )
    # abline(h = 0)
}





View(cbind(tt[year == year_, 2:7], poprojOriginal[year == year_, 7],
 tt[year == year_, 7] - poprojOriginal[year == year_, 7], 
 tt[year == year_, 7] / poprojOriginal[year == year_, 7])
 )

 View(cbind(
     poprojOriginal[year == 2 & age == 30 & sex == "men", 1:7],
      poprojOriginal[year == 3 & age == 31 & sex == "men", 7],
     poprojOriginal[year == 2 & age == 30 & sex == "men", 7] -
      poprojOriginal[year == 3 & age == 31 & sex == "men", 7]
 ))

year_ <- 2L
age_ <- 60L
 summary(poprojOriginal[year == year_ & age == age_, 7:16] -
     poprojOriginal[year == year_ + 1L & age == age_ + 1, 7:16])

}