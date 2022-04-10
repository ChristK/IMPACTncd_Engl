## IMPACTncdEngl is an implementation of the IMPACTncd framework, developed by Chris
## Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been
## funded by NIHR  HTA Project: 16/165/01 - IMPACTncdEngl: Health Outcomes
## Research Simulation Environment.  The views expressed are those of the
## authors and not necessarily those of the NHS, the NIHR or the Department of
## Health.
##
## Copyright (C) 2018-2020 University of Liverpool, Chris Kypridemos
##
## IMPACTncdEngl is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version. This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details. You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/> or write
## to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
## Boston, MA 02110-1301 USA.

mortality_projection <-
  function(dt,
           cause_nam,
           cause_nam2,
           hor,
           pca_method = c("M", "rapca", "classical"),
           b = 0.5,
           weighting = TRUE,
           min_age = 30L,
           max_age = 90L,
           power = 0.4
  ) {
    # b closer to 0 give more weight to most recent observations
    # if more than one age sex combinationas input, then forecast is coherent
    # Otherwise, not
    dt <- dt[!agegrp %in% agegrp_name(0, min_age - 1L)]

    rate <- vector("list", 0)
    pop <- vector("list", 0)

    for (k in unique(dt$sex)) {
      for (l in unique(dt$qimd)) {
        # Decompose mortality
        x1 <- dcast(dt[sex == k & qimd == l, ],
                    agegrp ~ year, value.var = cause_nam)[, agegrp := NULL]
        x2 <- dcast(dt[sex == k & qimd == l, ],
                    agegrp ~ year, value.var = "pops")[, agegrp := NULL]
        nam <-
          gsub(" ", "_", paste0("England__", k, "__", l, "__", cause_nam2))
        rate[[nam]] <- as.matrix(x1)
        pop[[nam]] <- as.matrix(x2)
      }
    }

    # demog data doesn't work on lists of matrices
    xx <- demogdata(
      rate[[1]],
      pop[[1]],
      c(seq(min_age, 85, 5) + 2L, 95),
      # c(0, 1, seq(min_age, 85, 5) + 2L, 95),
      # 0:100,
      sort(unique(dt$year)),
      "mortality",
      "England",
      names(rate)[1],
      lambda = 0
    )

    if (length(rate) > 1L) {
      # work around of above limitation
      xx$rate <- rate
      xx$pop  <- pop
      xx$name <- names(rate)
    }

    xx <-
      # without smoothing forecast doesn't fit
      smooth.demogdata(xx, age.grid = min_age:max_age, power = power)
    # print(plot(xx))

    if (length(rate) > 1L) {
      mort.fit <-
        coherentfdm(
          xx, 6L, 6L,
          method = pca_method,
          weight = weighting,
          # weight is the most important arguement. Geometrically decaying weights is applied to the decentralized data
          beta = b,
          # the speed of geometric decay (default 0.1)
          level = TRUE,
          # include an additional (intercept) term that depends on t but not on x.
          max.age = max_age
        )
    } else {
      mort.fit <-
        fdm(
          xx,
          method = pca_method,
          weight = weighting,
          # weight is the most important arguement. Geometrically decaying weights is applied to the decentralized data
          beta = b,
          # the speed of geometric decay (default 0.1)
          level = TRUE,
          # include an additional (intercept) term that depends on t but not on x.
          max.age = max_age
        )
    }
    mortf <- forecast(mort.fit,
                      h = hor,
                      jumpchoice = "actual",
                      level = 99)
    mortf60 <- forecast(mort.fit,
                        h = hor,
                        jumpchoice = "actual",
                        level = 60)
    mortf70 <- forecast(mort.fit,
                        h = hor,
                        jumpchoice = "actual",
                        level = 70)
    mortf80 <- forecast(mort.fit,
                        h = hor,
                        jumpchoice = "actual",
                        level = 80)
    mortf90 <- forecast(mort.fit,
                        h = hor,
                        jumpchoice = "actual",
                        level = 90)

    # produce lui & uui
    mortf.1 <- mortf.99 <- mortf
    mortf.40 <- mortf.60 <- mortf60
    mortf.30 <- mortf.70 <- mortf70
    mortf.20 <- mortf.80 <- mortf80
    mortf.10 <- mortf.90 <- mortf90


    if (length(rate) > 1L) {
      output <- vector("list", length(mortf) - 2)
      output.1 <- vector("list", length(mortf) - 2)
      output.99 <- vector("list", length(mortf) - 2)
      output.40 <- vector("list", length(mortf) - 2)
      output.60 <- vector("list", length(mortf) - 2)
      output.30 <- vector("list", length(mortf) - 2)
      output.70 <- vector("list", length(mortf) - 2)
      output.20 <- vector("list", length(mortf) - 2)
      output.80 <- vector("list", length(mortf) - 2)
      output.10 <- vector("list", length(mortf) - 2)
      output.90 <- vector("list", length(mortf) - 2)

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
                        keep.rownames = TRUE)[, `:=`(type = names(mortf[ii]))]
        output[[ii]][, c("sex", "qimd", "type") :=
                       tstrsplit(type, "__", fixed = TRUE, keep = c(2, 3, 4))]
        output.1[[ii]] <-
          as.data.table(mortf.1[[ii]]$rate[[1]],
                        keep.rownames = TRUE)[, `:=`(type = names(mortf.1[ii]))]
        output.1[[ii]][, c("sex", "qimd", "type") :=
                         tstrsplit(type, "__", fixed = TRUE, keep = c(2, 3, 4))]
        output.99[[ii]] <-
          as.data.table(mortf.99[[ii]]$rate[[1]],
                        keep.rownames = TRUE)[, `:=`(type = names(mortf.99[ii]))]
        output.99[[ii]][, c("sex", "qimd", "type") :=
                          tstrsplit(type, "__", fixed = TRUE, keep = c(2, 3, 4))]
        output.40[[ii]] <-
          as.data.table(mortf.40[[ii]]$rate[[1]],
                        keep.rownames = TRUE)[, `:=`(type = names(mortf.40[ii]))]
        output.40[[ii]][, c("sex", "qimd", "type") :=
                          tstrsplit(type, "__", fixed = TRUE, keep = c(2, 3, 4))]
        output.60[[ii]] <-
          as.data.table(mortf.60[[ii]]$rate[[1]],
                        keep.rownames = TRUE)[, `:=`(type = names(mortf.60[ii]))]
        output.60[[ii]][, c("sex", "qimd", "type") :=
                          tstrsplit(type, "__", fixed = TRUE, keep = c(2, 3, 4))]
        output.30[[ii]] <-
          as.data.table(mortf.30[[ii]]$rate[[1]],
                        keep.rownames = TRUE)[, `:=`(type = names(mortf.30[ii]))]
        output.30[[ii]][, c("sex", "qimd", "type") :=
                          tstrsplit(type, "__", fixed = TRUE, keep = c(2, 3, 4))]
        output.70[[ii]] <-
          as.data.table(mortf.70[[ii]]$rate[[1]],
                        keep.rownames = TRUE)[, `:=`(type = names(mortf.70[ii]))]
        output.70[[ii]][, c("sex", "qimd", "type") :=
                          tstrsplit(type, "__", fixed = TRUE, keep = c(2, 3, 4))]
        output.20[[ii]] <-
          as.data.table(mortf.20[[ii]]$rate[[1]],
                        keep.rownames = TRUE)[, `:=`(type = names(mortf.20[ii]))]
        output.20[[ii]][, c("sex", "qimd", "type") :=
                          tstrsplit(type, "__", fixed = TRUE, keep = c(2, 3, 4))]
        output.80[[ii]] <-
          as.data.table(mortf.80[[ii]]$rate[[1]],
                        keep.rownames = TRUE)[, `:=`(type = names(mortf.80[ii]))]
        output.80[[ii]][, c("sex", "qimd", "type") :=
                          tstrsplit(type, "__", fixed = TRUE, keep = c(2, 3, 4))]
        output.10[[ii]] <-
          as.data.table(mortf.10[[ii]]$rate[[1]],
                        keep.rownames = TRUE)[, `:=`(type = names(mortf.10[ii]))]
        output.10[[ii]][, c("sex", "qimd", "type") :=
                          tstrsplit(type, "__", fixed = TRUE, keep = c(2, 3, 4))]
        output.90[[ii]] <-
          as.data.table(mortf.90[[ii]]$rate[[1]],
                        keep.rownames = TRUE)[, `:=`(type = names(mortf.90[ii]))]
        output.90[[ii]][, c("sex", "qimd", "type") :=
                          tstrsplit(type, "__", fixed = TRUE, keep = c(2, 3, 4))]
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
    } else {
      mortf.1$rate[[1]]  <- mortf$rate$lower
      mortf.99$rate[[1]] <- mortf$rate$upper
      mortf.40$rate[[1]] <- mortf60$rate$lower
      mortf.60$rate[[1]] <- mortf60$rate$upper
      mortf.30$rate[[1]] <- mortf70$rate$lower
      mortf.70$rate[[1]] <- mortf70$rate$upper
      mortf.20$rate[[1]] <- mortf80$rate$lower
      mortf.80$rate[[1]] <- mortf80$rate$upper
      mortf.10$rate[[1]] <- mortf90$rate$lower
      mortf.90$rate[[1]] <- mortf90$rate$upper

      output <-
        as.data.table(mortf$rate[[1]],
                      keep.rownames = TRUE)[, `:=`(type = names(rate))]
      output[, c("sex", "qimd", "type") :=
               tstrsplit(type, "__", fixed = TRUE, keep = c(2, 3, 4))]
      output.1 <-
        as.data.table(mortf.1$rate[[1]],
                      keep.rownames = TRUE)[, `:=`(type = names(rate))]
      output.1[, c("sex", "qimd", "type") :=
                 tstrsplit(type, "__", fixed = TRUE, keep = c(2, 3, 4))]
      output.99 <-
        as.data.table(mortf.99$rate[[1]],
                      keep.rownames = TRUE)[, `:=`(type = names(rate))]
      output.99[, c("sex", "qimd", "type") :=
                  tstrsplit(type, "__", fixed = TRUE, keep = c(2, 3, 4))]
      output.40 <-
        as.data.table(mortf.40$rate[[1]],
                      keep.rownames = TRUE)[, `:=`(type = names(rate))]
      output.40[, c("sex", "qimd", "type") :=
                  tstrsplit(type, "__", fixed = TRUE, keep = c(2, 3, 4))]
      output.60 <-
        as.data.table(mortf.60$rate[[1]],
                      keep.rownames = TRUE)[, `:=`(type = names(rate))]
      output.60[, c("sex", "qimd", "type") :=
                  tstrsplit(type, "__", fixed = TRUE, keep = c(2, 3, 4))]
      output.30 <-
        as.data.table(mortf.30$rate[[1]],
                      keep.rownames = TRUE)[, `:=`(type = names(rate))]
      output.30[, c("sex", "qimd", "type") :=
                  tstrsplit(type, "__", fixed = TRUE, keep = c(2, 3, 4))]
      output.70 <-
        as.data.table(mortf.70$rate[[1]],
                      keep.rownames = TRUE)[, `:=`(type = names(rate))]
      output.70[, c("sex", "qimd", "type") :=
                  tstrsplit(type, "__", fixed = TRUE, keep = c(2, 3, 4))]
      output.20 <-
        as.data.table(mortf.20$rate[[1]],
                      keep.rownames = TRUE)[, `:=`(type = names(rate))]
      output.20[, c("sex", "qimd", "type") :=
                  tstrsplit(type, "__", fixed = TRUE, keep = c(2, 3, 4))]
      output.80 <-
        as.data.table(mortf.80$rate[[1]],
                      keep.rownames = TRUE)[, `:=`(type = names(rate))]
      output.80[, c("sex", "qimd", "type") :=
                  tstrsplit(type, "__", fixed = TRUE, keep = c(2, 3, 4))]
      output.10 <-
        as.data.table(mortf.10$rate[[1]],
                      keep.rownames = TRUE)[, `:=`(type = names(rate))]
      output.10[, c("sex", "qimd", "type") :=
                  tstrsplit(type, "__", fixed = TRUE, keep = c(2, 3, 4))]
      output.90 <-
        as.data.table(mortf.90$rate[[1]],
                      keep.rownames = TRUE)[, `:=`(type = names(rate))]
      output.90[, c("sex", "qimd", "type") :=
                  tstrsplit(type, "__", fixed = TRUE, keep = c(2, 3, 4))]
    }



    output <- melt(
      output,
      id.vars = c("rn", "sex", "qimd", "type"),
      value.name = "mx_total",
      variable.name = "year"
    )
    output.1 <-
      melt(
        output.1,
        id.vars = c("rn", "sex", "qimd", "type"),
        value.name = "mx_total_1",
        variable.name = "year"
      )
    output.99 <-
      melt(
        output.99,
        id.vars = c("rn", "sex", "qimd", "type"),
        value.name = "mx_total_99",
        variable.name = "year"
      )
    output.10 <-
      melt(
        output.10,
        id.vars = c("rn", "sex", "qimd", "type"),
        value.name = "mx_total_10",
        variable.name = "year"
      )
    output.20 <-
      melt(
        output.20,
        id.vars = c("rn", "sex", "qimd", "type"),
        value.name = "mx_total_20",
        variable.name = "year"
      )
    output.30 <-
      melt(
        output.30,
        id.vars = c("rn", "sex", "qimd", "type"),
        value.name = "mx_total_30",
        variable.name = "year"
      )
    output.40 <-
      melt(
        output.40,
        id.vars = c("rn", "sex", "qimd", "type"),
        value.name = "mx_total_40",
        variable.name = "year"
      )
    output.60 <-
      melt(
        output.60,
        id.vars = c("rn", "sex", "qimd", "type"),
        value.name = "mx_total_60",
        variable.name = "year"
      )
    output.70 <-
      melt(
        output.70,
        id.vars = c("rn", "sex", "qimd", "type"),
        value.name = "mx_total_70",
        variable.name = "year"
      )
    output.80 <-
      melt(
        output.80,
        id.vars = c("rn", "sex", "qimd", "type"),
        value.name = "mx_total_80",
        variable.name = "year"
      )
    output.90 <-
      melt(
        output.90,
        id.vars = c("rn", "sex", "qimd", "type"),
        value.name = "mx_total_90",
        variable.name = "year"
      )
    absorb_dt(output, output.1)
    absorb_dt(output, output.10)
    absorb_dt(output, output.20)
    absorb_dt(output, output.30)
    absorb_dt(output, output.40)
    absorb_dt(output, output.60)
    absorb_dt(output, output.70)
    absorb_dt(output, output.80)
    absorb_dt(output, output.90)
    absorb_dt(output, output.99)

    test <- copy(xx$rate)
    original <- vector("list", length(test))

    for (ii in 1:(length(test))) {
      original[[ii]] <-
        as.data.table(test[[ii]], keep.rownames = T)[, `:=`(type = names(test[ii]))]
      original[[ii]][, c("sex", "qimd", "type") := tstrsplit(type, "__", fixed = TRUE, keep =
                                                               c(2, 3, 4))]
    }
    original <- rbindlist(original)
    original <-
      melt(
        original,
        id.vars = c("rn", "sex", "qimd", "type"),
        value.name = "mx_total",
        variable.name = "year"
      )
    original[, paste0("mx_total_",
                      c(1, 99, 10, 20, 30, 40, 60, 70, 80, 90)) :=
               mx_total]

    lifetable_all <-
      rbind(original, output)
    lifetable_all[, `:=` (age = as.integer(rn),
                          rn = NULL,
                          year = as.integer(as.character(year)))]

    lifetable_all[, lapply(.SD, clamp, 0, 1, TRUE), .SDcols = patterns("^mx_")]
    return(lifetable_all)
  }

plot_projection <- function(dt, ci = TRUE, min_age = 30L, max_age = 90L, max_observed_year = 2016L) {
  if (ci) {
    p <- ggplot(dt[(age %% 10) == 0 & between(age, min_age, max_age),],
                aes(
                  y = mx_total,
                  ymin = mx_total_20,
                  ymax = mx_total_80,
                  x = year,
                  col = year <= max_observed_year,
                  fil = year <= max_observed_year
                )) +
      geom_line() + geom_ribbon(alpha = 1/5) +
      # geom_pointrange() +
      facet_grid(age ~ sex + qimd, scales = "free") + # , scales = "free"
      theme(legend.position = "none")
  } else {
    p <- ggplot(dt[(age %% 10) == 0 & between(age, min_age, max_age),],
                aes(y = mx_total,
                    x = year,
                    col = year <= max_observed_year)) +
      geom_line() +
      facet_grid(age ~ ., scales = "free") + # , scales = "free"
      theme(legend.position = "none")
  }
  p <- p + ggtitle(paste0(unique(tt$type), "_", unique(tt$sex), "_", gsub(" ", "_", unique(tt$qimd))))

}
