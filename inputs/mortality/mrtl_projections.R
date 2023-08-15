library(data.table)
library(demography)
library(fst)

# dt <- fread("P:/My Datasets/ONS mortality/all-cause by dimd 2001-2019/mortality_dimd_2001-19.csv")
dt <- fread("./inputs/mortality/mortality_dimd_2001-19.csv")
# strata <- c("year", "sex", "age", "dimd")

hor <- 70L # maximum simulation horizon
lifetable_all <- data.table(NULL)

for (l in unique(dt$dimd)) {
  for (k in unique(dt$sex)) {
    print(paste0(l, " ", k))

    # Decompose mortality
    x1 <- dcast(dt[sex == k & dimd == l, ],
                age ~ year, value.var = "mx")
    x1[, age := NULL]

    x2 <- dcast(dt[sex == k & dimd == l, ],
                age ~ year, value.var = "pop_size")
    x2[, age := NULL]

    nam <- paste0("England_", k, "_", l, "_allcause")
    rate <- as.matrix(x1)
    pop <- as.matrix(x2)

    xx <- demogdata(
      rate,
      pop,
      c(0:89, 95), # list of the mean of our age groups
      sort(unique(dt$year)),
      "mortality",
      paste0("England"),
      paste0("England_", k, "_", l, "_allcause"),
      lambda = 0
    )

    xx <- smooth.demogdata(xx, age.grid = 0:99, obs.var = "empirical")

    mort.fit <-
      fdm(extract.years(xx, 2001:2015),
          method = "rapca", # M crashes least deprived men
          weight = FALSE, # weight is the most important argument
          # beta = 0.9,
          max.age = 99)
    mortf   <- forecast(mort.fit, h = 4, max.d = 1, level = 99)
    rmse <- sqrt(mean((rate[1:90, 16:19] - mortf$rate[[1]][1:90, ])^2))

    mort.fit <-
      fdm(xx,
          method = "rapca", # M crashes least deprived men
          weight = FALSE, # weight is the most important argument
          # beta = 0.9,
          max.age = 99)

    for (i in seq(0.01, 0.99, 0.01)) {
      mort.fit2 <-
        fdm(extract.years(xx, 2001:2015),
            method = "rapca", # M crashes least deprived men
            weight = TRUE, # weight is the most important argument
            beta = i,
            max.age = 99)
      mortf   <- forecast(mort.fit2, h = 4, max.d = 1, level = 99)
      rmse2 <- sqrt(mean((rate[1:90, 16:19] - mortf$rate[[1]][1:90, ])^2))


      if (rmse2 < rmse) {
        print(paste0("SUCCESS! ", "weight ", i, " RMSE: ", rmse2))
        rmse <- rmse2
        keep <- i
        mort.fit <- fdm(xx,
                        method = "rapca", # M crashes least deprived men
                        weight = TRUE, # weight is the most important argument
                        beta = i,
                        max.age = 99)
      } else print(paste0("weight ", i, " RMSE: ", rmse2))

    }

    # if (exists("keep")) {
    # for (i in seq(keep - 0.09, keep + 0.09, 0.01)) {
    #   mort.fit2 <-
    #     fdm(extract.years(xx, 2001:2015),
    #         method = "rapca", # M crashes least deprived men
    #         weight = TRUE, # weight is the most important argument
    #         beta = i,
    #         max.age = 99)
    #   mortf   <- forecast(mort.fit2, h = 4, max.d = 1, level = 99)
    #   rmse2 <- sqrt(mean((rate[1:90, 16:19] - mortf$rate[[1]][1:90, ])^2))
    #
    #   if (rmse2 < rmse) {
    #     print(paste0("weight ", i, " RMSE: ", rmse2))
    #     rmse <- rmse2
    #     mort.fit <- fdm(xx,
    #                     method = "rapca", # M crashes least deprived men
    #                     weight = TRUE, # weight is the most important argument
    #                     beta = i,
    #                     max.age = 99)
    #   } else print(i)
    # }
    # }

    mortf   <- forecast(mort.fit, h = hor, max.d = 1, level = 99)
    mortf60 <- forecast(mort.fit, h = hor, max.d = 1, level = 60)
    mortf70 <- forecast(mort.fit, h = hor, max.d = 1, level = 70)
    mortf80 <- forecast(mort.fit, h = hor, max.d = 1, level = 80)
    mortf90 <- forecast(mort.fit, h = hor, max.d = 1, level = 90)

    # produce lui & uui
    mortf.1 <- mortf.99 <- mortf
    mortf.40 <- mortf.60 <- mortf60
    mortf.30 <- mortf.70 <- mortf70
    mortf.20 <- mortf.80 <- mortf80
    mortf.10 <- mortf.90 <- mortf90

    mortf.1$rate  <- mortf$rate$lower
    mortf.99$rate <- mortf$rate$upper
    mortf.40$rate <- mortf60$rate$lower
    mortf.60$rate <- mortf60$rate$upper
    mortf.30$rate <- mortf70$rate$lower
    mortf.70$rate <- mortf70$rate$upper
    mortf.20$rate <- mortf80$rate$lower
    mortf.80$rate <- mortf80$rate$upper
    mortf.10$rate <- mortf90$rate$lower
    mortf.90$rate <- mortf90$rate$upper

    strata <- c("sex", "dimd", "type")
    output <-
      as.data.table(mortf$rate[[1]],
                    keep.rownames = TRUE)[, `:=`(type = names(mortf$rate[1]))]
    output[, (strata) :=
             tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
    output.1 <-
      as.data.table(mortf.1$rate,
                    keep.rownames = TRUE)[, `:=`(type = names(mortf$rate[1]))]
    output.1[, (strata) :=
               tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
    output.99 <-
      as.data.table(mortf.99$rate,
                    keep.rownames = TRUE)[, `:=`(type = names(mortf$rate[1]))]
    output.99[, (strata) :=
                tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
    output.40 <-
      as.data.table(mortf.40$rate,
                    keep.rownames = TRUE)[, `:=`(type = names(mortf$rate[1]))]
    output.40[, (strata) :=
                tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
    output.60 <-
      as.data.table(mortf.60$rate,
                    keep.rownames = TRUE)[, `:=`(type = names(mortf$rate[1]))]
    output.60[, (strata) :=
                tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
    output.30 <-
      as.data.table(mortf.30$rate,
                    keep.rownames = TRUE)[, `:=`(type = names(mortf$rate[1]))]
    output.30[, (strata) :=
                tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
    output.70 <-
      as.data.table(mortf.70$rate,
                    keep.rownames = TRUE)[, `:=`(type = names(mortf$rate[1]))]
    output.70[, (strata) :=
                tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
    output.20 <-
      as.data.table(mortf.20$rate,
                    keep.rownames = TRUE)[, `:=`(type = names(mortf$rate[1]))]
    output.20[, (strata) :=
                tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
    output.80 <-
      as.data.table(mortf.80$rate,
                    keep.rownames = TRUE)[, `:=`(type = names(mortf$rate[1]))]
    output.80[, (strata) :=
                tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
    output.10 <-
      as.data.table(mortf.10$rate,
                    keep.rownames = TRUE)[, `:=`(type = names(mortf$rate[1]))]
    output.10[, (strata) :=
                tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
    output.90 <-
      as.data.table(mortf.90$rate,
                    keep.rownames = TRUE)[, `:=`(type = names(mortf$rate[1]))]
    output.90[, (strata) :=
                tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]


    strata <- c("rn", "sex", "dimd", "type")
    output <- melt(output,
                   id.vars = strata,
                   value.name = "mx_total")
    output.1 <-
      melt(output.1,
           id.vars = strata,
           value.name = "mx_total_1")
    output.99 <-
      melt(output.99,
           id.vars = strata,
           value.name = "mx_total_99")
    output.10 <-
      melt(output.10,
           id.vars = strata,
           value.name = "mx_total_10")
    output.20 <-
      melt(output.20,
           id.vars = strata,
           value.name = "mx_total_20")
    output.30 <-
      melt(output.30,
           id.vars = strata,
           value.name = "mx_total_30")
    output.40 <-
      melt(output.40,
           id.vars = strata,
           value.name = "mx_total_40")
    output.60 <-
      melt(output.60,
           id.vars = strata,
           value.name = "mx_total_60")
    output.70 <-
      melt(output.70,
           id.vars = strata,
           value.name = "mx_total_70")
    output.80 <-
      melt(output.80,
           id.vars = strata,
           value.name = "mx_total_80")
    output.90 <-
      melt(output.90,
           id.vars = strata,
           value.name = "mx_total_90")

    strata <- c("rn", "sex", "dimd", "type", "variable")
    output[output.1, mx_total_1 := i.mx_total_1, on = strata]
    output[output.99, mx_total_99 := i.mx_total_99, on = strata]
    output[output.10, mx_total_10 := i.mx_total_10, on = strata]
    output[output.20, mx_total_20 := i.mx_total_20, on = strata]
    output[output.30, mx_total_30 := i.mx_total_30, on = strata]
    output[output.40, mx_total_40 := i.mx_total_40, on = strata]
    output[output.60, mx_total_60 := i.mx_total_60, on = strata]
    output[output.70, mx_total_70 := i.mx_total_70, on = strata]
    output[output.80, mx_total_80 := i.mx_total_80, on = strata]
    output[output.90, mx_total_90 := i.mx_total_90, on = strata]

    original <-
      as.data.table(xx$rate[[1]], keep.rownames = TRUE)[, `:=`(type = names(xx$rate))]
    original[, c("sex", "dimd", "type") := tstrsplit(type, "_", fixed = TRUE, keep = 2:4)]

    original <-
      melt(original,
           id.vars = c("rn", "sex", "dimd", "type"),
           value.name = "mx_total")
    original[, paste0("mx_total_",
                      c(1, 99, 10, 20, 30, 40, 60, 70, 80, 90)) :=
               mx_total]

    out <-
      rbind(original, output, use.names = TRUE, fill = TRUE)
    out[, `:=` (age = as.integer(rn),
                year = as.integer(as.character(variable)),
                rn = NULL,
                variable = NULL
                )]

    lifetable_all <- rbind(lifetable_all, out)

  }
}

fwrite(lifetable_all, "./inputs/mortality/lt_by_age.csv")

setkey(lifetable_all, year)
# Convert from mid-year (1st-2nd July) to 1st of April (3 months)
foo <- function(x) {
  3/4 * x + 1/4 * shift(x)
}
lifetable_all <- lifetable_all[, c(list("year" = year), lapply(.SD, foo)), by = .(age, sex, dimd), .SDcols = patterns("^mx_")]

# by 5-year age groups
setkey(lifetable_all, age)
CKutils::to_agegrp(lifetable_all, max_age = 99, min_age = 0)
lifetable_all <- lifetable_all[year > 2001, lapply(.SD, function(x) last(1 - cumprod(1 - x))), keyby = .(year, agegrp, sex, dimd), .SDcols = patterns("^mx_")]
# The formula above estimates qx by age group.

# Convert qx to mx
# mx to qx from https://rdrr.io/cran/MortalityLaws/src/R/LifeTable.R
for (j in grep("mx_", names(lifetable_all), value = TRUE)) {
  lifetable_all[agegrp != "<1", (j) := (-log(1 - get(j))/5)]
  lifetable_all[agegrp == "<1", (j) := (-log(1 - get(j))/1)]
}


fwrite(lifetable_all, "./inputs/mortality/lt.csv")
lifetable_all <- fread("./inputs/mortality/lt.csv")


#test graph
library(ggplot2)

ggplot(lifetable_all[year < 2045, ],
  aes(
    y = mx_total,
    x = year,
    col = dimd
  )) +
  geom_line() +
  facet_grid(agegrp~sex , scales = "free")


lifetable_all[, quan := runif(1)]
# lifetable_all[year >= 2030, quan := runif(1)]
 ggplot(lifetable_all[year < 2045, ],
   aes(
     y = mx_total,
     ymin = mx_total_1,
     ymax = mx_total_99,
     x = year
   )) +
   geom_pointrange() +
   facet_grid(agegrp~sex + dimd, scales = "free")





# validation
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

gg <- ggplot(
  lifetable_all,
  aes(
    x = year,
    y = mx_total,
    ymin = mx_total_10,
    ymax = mx_total_90,
    col = factor(age),
    fill = factor(age)
  )
) +
  geom_point(size = 1,
             alpha = 5 / 5,
             show.legend = F) +
  geom_line(linewidth = 1, alpha = 5 / 5) +
  geom_ribbon(alpha = 1 / 5,
              linetype = 0,
              show.legend = F) +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "All-cause mortality") +
  facet_grid(dimd ~ .)
ggsave2(
  filename = "All_cause_mortality_projections_mx.png",
  gg, height = 9, width = 16, units = "in",
  path = "C:/Users/Zoe/OneDrive - The University of Liverpool/ZOE PROJECT 2/Datas/Results/Population_statistics"
)

