library(data.table)
library(demography)
library(fst)

# dt <- fread("P:/My Datasets/ONS mortality/all-cause by dimd 2001-2019/mortality_dimd_2001-19.csv")
dt <- fread("./inputs/mortality/mortality_dimd_2001-19.csv")
strata <- c("year", "sex", "age", "dimd")

hor <- 70L # maximum simulation horizon
rate <- vector("list", 0)
pop <- vector("list", 0)

for (l in unique(dt$dimd)) {
  for (k in unique(dt$sex)) {
    # Decompose mortality
    x1 <- dcast(dt[sex == k & dimd == l, ],
                age ~ year, value.var = "mx")
    x1[, age := NULL]

    x2 <- dcast(dt[sex == k & dimd == l, ],
                age ~ year, value.var = "pop_size")
    x2[, age := NULL]

    nam <- paste0("England_", k, "_", l, "_allcause")
    rate[[nam]] <- as.matrix(x1)
    pop[[nam]] <- as.matrix(x2)
  }
}

# demog data doesn't work on lists of matrices
xx <- demogdata(
  rate[[1]],
  pop[[1]],
  c(0:89, 95), # list of the mean of our age groups
  sort(unique(dt$year)),
  "mortality",
  paste0("England"),
  names(rate[1]),
  lambda = 0
)

# work around of above limitation
xx$rate <- rate
xx$pop  <- pop
xx$name <- names(rate)

xx <- smooth.demogdata(xx, age.grid = 0:99, obs.var = "empirical")

mort.fit <- # or coherentfdm
  coherentfdm(xx,
              3,
              6,
              method = "M",
              # level = TRUE,
              # weight = TRUE, # weight is the most important argument
              # beta = 0.9,
              max.age = 99)

# mort.fit2 <- # or coherentfdm
#   coherentfdm(xx,
#               3,
#               6,
#               method = "M",
#               level = TRUE,
#               weight = TRUE, # weight is the most important argument
#               beta = 0.99,
#               max.age = 99)
#
# summary(mort.fit$product)
# summary(mort.fit2$product)
#
# summary(mort.fit$ratio$`England_men_1 most deprived_allcause`)
# summary(mort.fit2$ratio$`England_men_1 most deprived_allcause`)
#
# summary(mort.fit$ratio$`England_men_10 least deprived_allcause`)
# summary(mort.fit2$ratio$`England_men_10 least deprived_allcause`)

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

output    <- vector("list", length(mortf) - 2)
output.1  <- vector("list", length(mortf) - 2)
output.99 <- vector("list", length(mortf) - 2)
output.40 <- vector("list", length(mortf) - 2)
output.60 <- vector("list", length(mortf) - 2)
output.30 <- vector("list", length(mortf) - 2)
output.70 <- vector("list", length(mortf) - 2)
output.20 <- vector("list", length(mortf) - 2)
output.80 <- vector("list", length(mortf) - 2)
output.10 <- vector("list", length(mortf) - 2)
output.90 <- vector("list", length(mortf) - 2)

strata <- c("sex", "dimd", "type")
for (ii in 1:(length(mortf) - 2)) {
  mortf.1[[ii]]$rate[[1]]  <- mortf[[ii]]$rate$lower
  mortf.99[[ii]]$rate[[1]] <- mortf[[ii]]$rate$upper
  mortf.40[[ii]]$rate[[1]] <- mortf60[[ii]]$rate$lower
  mortf.60[[ii]]$rate[[1]] <- mortf60[[ii]]$rate$upper
  mortf.30[[ii]]$rate[[1]] <- mortf70[[ii]]$rate$lower
  mortf.70[[ii]]$rate[[1]] <- mortf70[[ii]]$rate$upper
  mortf.20[[ii]]$rate[[1]] <- mortf80[[ii]]$rate$lower
  mortf.80[[ii]]$rate[[1]] <- mortf80[[ii]]$rate$upper
  mortf.10[[ii]]$rate[[1]] <- mortf90[[ii]]$rate$lower
  mortf.90[[ii]]$rate[[1]] <- mortf90[[ii]]$rate$upper

  output[[ii]] <-
    as.data.table(mortf[[ii]]$rate[[1]],
                  keep.rownames = T)[, `:=`(type = names(mortf[ii]))]
  output[[ii]][, (strata) :=
                 tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
  output.1[[ii]] <-
    as.data.table(mortf.1[[ii]]$rate[[1]],
                  keep.rownames = T)[, `:=`(type = names(mortf.1[ii]))]
  output.1[[ii]][, (strata) :=
                   tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
  output.99[[ii]] <-
    as.data.table(mortf.99[[ii]]$rate[[1]],
                  keep.rownames = T)[, `:=`(type = names(mortf.99[ii]))]
  output.99[[ii]][, (strata) :=
                    tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
  output.40[[ii]] <-
    as.data.table(mortf.40[[ii]]$rate[[1]],
                  keep.rownames = T)[, `:=`(type = names(mortf.40[ii]))]
  output.40[[ii]][, (strata) :=
                    tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
  output.60[[ii]] <-
    as.data.table(mortf.60[[ii]]$rate[[1]],
                  keep.rownames = T)[, `:=`(type = names(mortf.60[ii]))]
  output.60[[ii]][, (strata) :=
                    tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
  output.30[[ii]] <-
    as.data.table(mortf.30[[ii]]$rate[[1]],
                  keep.rownames = T)[, `:=`(type = names(mortf.30[ii]))]
  output.30[[ii]][, (strata) :=
                    tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
  output.70[[ii]] <-
    as.data.table(mortf.70[[ii]]$rate[[1]],
                  keep.rownames = T)[, `:=`(type = names(mortf.70[ii]))]
  output.70[[ii]][, (strata) :=
                    tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
  output.20[[ii]] <-
    as.data.table(mortf.20[[ii]]$rate[[1]],
                  keep.rownames = T)[, `:=`(type = names(mortf.20[ii]))]
  output.20[[ii]][, (strata) :=
                    tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
  output.80[[ii]] <-
    as.data.table(mortf.80[[ii]]$rate[[1]],
                  keep.rownames = T)[, `:=`(type = names(mortf.80[ii]))]
  output.80[[ii]][, (strata) :=
                    tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
  output.10[[ii]] <-
    as.data.table(mortf.10[[ii]]$rate[[1]],
                  keep.rownames = T)[, `:=`(type = names(mortf.10[ii]))]
  output.10[[ii]][, (strata) :=
                    tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
  output.90[[ii]] <-
    as.data.table(mortf.90[[ii]]$rate[[1]],
                  keep.rownames = T)[, `:=`(type = names(mortf.90[ii]))]
  output.90[[ii]][, (strata) :=
                    tstrsplit(type, "_", fixed = TRUE, keep = 2:(length(strata) +1L))]
}

output    <- rbindlist(output)
output.1  <- rbindlist(output.1)
output.99 <- rbindlist(output.99)
output.40 <- rbindlist(output.40)
output.60 <- rbindlist(output.60)
output.30 <- rbindlist(output.30)
output.70 <- rbindlist(output.70)
output.20 <- rbindlist(output.20)
output.80 <- rbindlist(output.80)
output.10 <- rbindlist(output.10)
output.90 <- rbindlist(output.90)

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

test <- copy(xx$rate)
original <- vector("list", length(test))

for (ii in 1:(length(test))) {
  original[[ii]] <-
    as.data.table(test[[ii]], keep.rownames = T)[, `:=`(type = names(test[ii]))]
  original[[ii]][, c("sex", "dimd", "type") := tstrsplit(type, "_", fixed = TRUE, keep = 2:4)]
}
original <- rbindlist(original)
original <-
  melt(original,
       id.vars = c("rn", "sex", "dimd", "type"),
       value.name = "mx_total")
original[, paste0("mx_total_",
                  c(1, 99, 10, 20, 30, 40, 60, 70, 80, 90)) :=
           mx_total]

lifetable_all <-
  rbind(original, output, use.names = TRUE, fill = TRUE)
lifetable_all[, age := as.integer(rn)]
lifetable_all[, year := as.integer(as.character(variable))]


#test coherence graph
library(ggplot2)

ggplot(lifetable_all[
  (age %% 10) == 5 & between(age, 0, 99) & year > 1000, ],
  aes(
    y = mx_total,
    x = year,
    col = dimd
  )) +
  geom_line() +
  facet_grid(age~sex , scales = "free")


lifetable_all[, quan := runif(1)]
# lifetable_all[year >= 2030, quan := runif(1)]
 ggplot(lifetable_all[
   (age %% 10) == 5 & between(age, 20, 99) & year > 1000, ],
   aes(
     y = mx_total,
     ymin = mx_total_1,
     ymax = mx_total_99,
     x = year
   )) +
   geom_pointrange() +
   facet_grid(age~sex + dimd, scales = "free")





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
  geom_line(size = 1, alpha = 5 / 5) +
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



data(apc)
cases
model1 <- bamp(t(rate[[1]]), t(pop[[1]]), age="rw2", period="rw2", cohort="rw2",
               periods_per_agegroup = 1)
