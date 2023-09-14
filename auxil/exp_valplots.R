
source("./global.R")

library(ggthemes)
theme_set(new = theme_economist_white(gray_bg = FALSE))
theme_update(legend.position = "bottom",
             legend.title= element_blank(),
             axis.title = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10), size = 10),
             axis.text = element_text(size = 10),
             strip.text = element_text(size = 10)
)
#  axis.text = element_text(size = 12),
#              axis.title = element_text(size = 12,
#                                        margin = margin(t = 10, r = 10, b = 10, l = 10)),
#              legend.position = "bottom",
#              legend.text=element_text(size=12),
#              #plot.title = element_text(hjust = 0.5),
# )

plot_synthpop_val <-
  function(dt, # data
           x, # Name of risk factor not in inverted commas
           grp, # strata
           wt, # weights
           title,
           x_label,
           breaks,
           xangle = 0,
           standardised_to_grp = c(FALSE, TRUE),
           print_to_screen = c(FALSE, TRUE)) {
    tmp <- paste0(grp, collapse = "_") #this is collapsing the group for saving the plot name
    bks <- if(breaks == "short"){
      c(0,0.5,1)}else{
        c(0,0.25,0.5,0.75,1)}

    if (standardised_to_grp) {
      dt[, weight := get(wt) / sum(get(wt)), by = c("type", grp)]
    } else {
      dt[, weight := get(wt) / sum(get(wt)), by = c("type")]
    }

    x <- enquo(x)
    p <- ggplot(dt, aes(
      !!x,
      colour = type,
      weight = weight,
      linetype = type
    )) +
      geom_density() +
      facet_grid(grp) +
      xlab(x_label) +
      ylab("Probability") +
#      ggtitle(title) +
      scale_color_brewer(palette = "Set1") +
      scale_y_continuous(breaks = bks) +
      theme(axis.text.x = element_text(angle = xangle))

    if (print_to_screen)
      print(p)
    suppressWarnings(
      cowplot::ggsave2(
#        paste0(gsub(" ", "_", title), "_density.png"),
        paste0(title, "_", tmp, "_density.png"),
        p,
        width = 16,
        height = 9,
        units = "cm",
        scale = 2,
        dpi = 300,
        path = "./validation/synthpop_models"
      )
    )

    p <- ggplot(dt, aes(
      !!x,
      colour = type,
      weight = weight,
      linetype = type
    )) +
      stat_ecdf() +
      facet_grid(grp) +
      xlab(x_label) +
      ylab("Probability") +
#      ggtitle(title) +
      scale_color_brewer(palette = "Set1") +
      scale_y_continuous(breaks = bks) +
      theme(axis.text.x = element_text(angle = xangle))

        if (print_to_screen)
      print(p)


    suppressWarnings(
      cowplot::ggsave2(
#        paste0(gsub(" ", "_", title), "_cdf.png"),
        paste0(title, "_", tmp, "_cdf.png"),
        p,
        width = 16,
        height = 9,
        units = "cm",
        scale = 2,
        dpi = 300,
        path = "./validation/synthpop_models"
      )
    )
  }





#Alcohol
xps_mdl <- qread("./secure_data/lifecourse_models/alcohol_model.qs")
dt <- xps_mdl$data
distr_nam <- xps_mdl$family[[1]]

xlab_nam <- expression(bold(Alcohol ~ (g / d)))
dt[, age := age * 17.5 + 51.3] # Need to tell you where this is coming from
alcohol_model_tbl <- read_fst("./inputs/exposure_distributions/alcohol_table.fst", as.data.table =  TRUE)
zz <-
  validate_gamlss_tbl(dt, alcohol_model_tbl, 50, "totalwu", paste0("q", distr_nam)
  )[between(totalwu, quantile(totalwu, 0.01), quantile(totalwu, 0.99))]
zz[, type := ifelse(type == "Modelled", "IMPACTncd", "HSE")]
zz[, year := year + 2000]
zz[, qimd := factor(qimd,
                    levels = c("1 most deprived", "2", "3", "4", "5 least deprived"),
                    labels = c("1\nmost deprived", "2", "3", "4", "5\nleast deprived"))]



plot_synthpop_val(zz, totalwu, c("year","agegrp10"), "wt_int", "Alcohol", "long", xlab_nam, xangle = 0, FALSE, TRUE)
plot_synthpop_val(zz, totalwu, c("qimd","agegrp10"), "wt_int", "Alcohol", "long", xlab_nam, xangle = 0, FALSE, TRUE)






#Active days - this is a diff distribution. Don't know how to use validate_gamlss_tbl with this
xps_mdl <- qread("./secure_data/lifecourse_models/active_days_model.qs")
pa_model_tbl <- read_fst("./inputs/exposure_distributions/active_days_table.fst", as.data.table =  TRUE)


zz <- clone_dt(xps_mdl$data, 10)
zz[, active_days := NULL]
zz[, age := round(age*17+49.6)]
zz[, rank_pa := runif(.N)]
nam <- intersect(names(zz), names(pa_model_tbl))
zz[pa_model_tbl, active_days := (rank_pa > pa0) + (rank_pa > pa1) + (rank_pa > pa2) +
     (rank_pa > pa3) + (rank_pa > pa4) + (rank_pa > pa5) + (rank_pa > pa6),
   on = nam]
zz[, `:=` (
  type = "Modelled",
  active_days = factor(
    active_days,
    levels = 0:7,
    labels = 0:7,
    ordered = TRUE
  ),
  .id = NULL
)]
zz[, rank_pa := NULL]
zz <- rbind(zz, xps_mdl$data[, type := "Observed"])
zz[, active_days := as.integer(as.character(active_days))]
zz[, type := ifelse(type == "Modelled", "IMPACTncd", "HSE")]
zz[, year := year + 2000]
zz[, qimd := factor(qimd,
                    levels = c("1 most deprived", "2", "3", "4", "5 least deprived"),
                    labels = c("1\nmost deprived", "2", "3", "4", "5\nleast deprived"))]

future({
 # dir.create("./validation/synthpop_models", FALSE)
  zz[, weight := wt_int/sum(wt_int), by = type]
  png(
    "./validation/synthpop_models/Active_days_rel_dist.png",
    3840,
    2160,
    pointsize = 48

  )
  reldist_diagnostics(zz[type == "Observed", active_days],
                      zz[type == "Modelled", active_days],
                      zz[type == "Observed", weight],
                      zz[type == "Modelled", weight],
                      main = xlab_nam,
                      discrete = TRUE)
  dev.off()
})


plot_synthpop_val(zz, active_days, c("year","agegrp10"), "wt_int", "Active_days",  xlab_nam,  "long", xangle = 0, FALSE, TRUE)
plot_synthpop_val(zz, active_days, c("qimd","agegrp10"), "wt_int", "Active_days",  xlab_nam, "long", xangle = 0, FALSE, TRUE)




#bmi
xps_mdl <- qread("./secure_data/lifecourse_models/bmi_model.qs")
dt <- xps_mdl$data
distr_nam <- xps_mdl$family[[1]]

xlab_nam <- expression(bold(BMI ~ (kg / m)))
dt[, age := round(age*16.8+50.8)]

bmi_tbl <- read_fst("./inputs/exposure_distributions/bmi_table.fst", as.data.table =  TRUE)
zz <-
  validate_gamlss_tbl(dt, bmi_tbl, 50, "bmi", paste0("q", distr_nam)
  )[between(bmi, quantile(bmi, 0.01), quantile(bmi, 0.99))]
zz[, type := ifelse(type == "Modelled", "IMPACTncd", "HSE")]
zz[, year := year + 2000]
zz[, qimd := factor(qimd,
                    levels = c("1 most deprived", "2", "3", "4", "5 least deprived"),
                    labels = c("1\nmost deprived", "2", "3", "4", "5\nleast deprived"))]


plot_synthpop_val(zz, bmi, c("year","agegrp10"), "wt_nurse", "BMI", xlab_nam, "short", xangle = 90, FALSE, TRUE)
plot_synthpop_val(zz, bmi, c("qimd","agegrp10"), "wt_nurse", "BMI", xlab_nam, "long", xangle = 90,FALSE, TRUE)


#ets #what is the unit
xps_mdl <- qread("./secure_data/lifecourse_models/ets_model.qs")
dt <- xps_mdl$data
distr_nam <- xps_mdl$family[[1]]

xlab_nam <- expression(bold("Environmental tobacco smoke"))
dt[, age := round(age*17.4+50.8)]

ets_tbl <- read_fst("./inputs/exposure_distributions/ets_table.fst", as.data.table =  TRUE)
zz <-
  validate_gamlss_tbl(dt, ets_tbl, 50, "ets", paste0("q", distr_nam)
  )[between(ets, quantile(ets, 0.01), quantile(ets, 0.99))]
zz[, type := ifelse(type == "Modelled", "IMPACTncd", "HSE")]
zz[, year := year + 2000]
zz[, qimd := factor(qimd,
                    levels = c("1 most deprived", "2", "3", "4", "5 least deprived"),
                    labels = c("1\nmost deprived", "2", "3", "4", "5\nleast deprived"))]


plot_synthpop_val(zz, ets, c("year","agegrp10"), "wt_int", "ets", xlab_nam, "short", xangle = 90, FALSE, TRUE)
plot_synthpop_val(zz, ets, c("qimd","agegrp10"), "wt_int", "ets", xlab_nam, "long", xangle = 90, FALSE, TRUE)



#fruit
xps_mdl <- qread("./secure_data/lifecourse_models/frtpor_model.qs")
dt <- xps_mdl$data
distr_nam <- xps_mdl$family[[1]]

xlab_nam <- expression(bold("Fruit intake (portions/day)"))
dt[, age := round(age*17.4+50.6)]

fruit_tbl <- read_fst("./inputs/exposure_distributions/frtpor_table.fst", as.data.table =  TRUE)
zz <-
  validate_gamlss_tbl(dt, fruit_tbl, 50, "frtpor", paste0("q", distr_nam)
  )[between(frtpor, quantile(frtpor, 0.01), quantile(frtpor, 0.99))]
zz[, type := ifelse(type == "Modelled", "IMPACTncd", "HSE")]
zz[, year := year + 2000]
zz[, qimd := factor(qimd,
                    levels = c("1 most deprived", "2", "3", "4", "5 least deprived"),
                    labels = c("1\nmost deprived", "2", "3", "4", "5\nleast deprived"))]


plot_synthpop_val(zz, frtpor, c("year","agegrp10"), "wt_int", "frtpor", xlab_nam, "short", xangle = 0, FALSE, TRUE)
plot_synthpop_val(zz, frtpor, c("qimd","agegrp10"), "wt_int", "frtpor", xlab_nam, "long", xangle = 0, FALSE, TRUE)



#veg
xps_mdl <- qread("./secure_data/lifecourse_models/vegpor_model.qs")
dt <- xps_mdl$data
distr_nam <- xps_mdl$family[[1]]

xlab_nam <- expression(bold("Vegetable intake (portions/day)"))
dt[, age := round(age*17.4+50.6)]

veg_tbl <- read_fst("./inputs/exposure_distributions/vegpor_table.fst", as.data.table =  TRUE)
zz <-
  validate_gamlss_tbl(dt, veg_tbl, 50, "vegpor", paste0("q", distr_nam)
  )[between(vegpor, quantile(vegpor, 0.01), quantile(vegpor, 0.99))]
zz[, type := ifelse(type == "Modelled", "IMPACTncd", "HSE")]
zz[, year := year + 2000]
zz[, qimd := factor(qimd,
                    levels = c("1 most deprived", "2", "3", "4", "5 least deprived"),
                    labels = c("1\nmost deprived", "2", "3", "4", "5\nleast deprived"))]


plot_synthpop_val(zz, vegpor, c("year","agegrp10"), "wt_int", "vegpor", xlab_nam, "short", xangle = 0, FALSE, TRUE)
plot_synthpop_val(zz, vegpor, c("qimd","agegrp10"), "wt_int", "vegpor", xlab_nam, "long", xangle = 0, FALSE, TRUE)


#tchol
xps_mdl <- qread("./secure_data/lifecourse_models/tchol_model.qs")
dt <- xps_mdl$data
distr_nam <- xps_mdl$family[[1]]

xlab_nam <- expression(bold("Total cholesterol ~ (mmol/l)"))
dt[, age := round(age*16.6+52)]

tchol_tbl <- read_fst("./inputs/exposure_distributions/tchol_table.fst", as.data.table =  TRUE)
zz <-
  validate_gamlss_tbl(dt, tchol_tbl, 50, "tchol", paste0("q", distr_nam)
  )[between(tchol, quantile(tchol, 0.01), quantile(tchol, 0.99))]
zz[, type := ifelse(type == "Modelled", "IMPACTncd", "HSE")]
zz[, year := year + 2000]
zz[, qimd := factor(qimd,
                    levels = c("1 most deprived", "2", "3", "4", "5 least deprived"),
                    labels = c("1\nmost deprived", "2", "3", "4", "5\nleast deprived"))]


plot_synthpop_val(zz, tchol, c("year","agegrp10"), "wt_blood", "tchol", xlab_nam, "short", xangle = 0, FALSE, TRUE)
plot_synthpop_val(zz, tchol, c("qimd","agegrp10"), "wt_blood", "tchol", xlab_nam, "long", xangle = 0, FALSE, TRUE)



#sbp
xps_mdl <- qread("./secure_data/lifecourse_models/sbp_model.qs")
dt <- xps_mdl$data
distr_nam <- xps_mdl$family[[1]]

xlab_nam <- expression(bold("SBP ~ (mmHg)"))
dt[, age := round(age*17.1+52.1)]

sbp_tbl <- read_fst("./inputs/exposure_distributions/sbp_table.fst", as.data.table =  TRUE)
zz <-
  validate_gamlss_tbl(dt, sbp_tbl, 50, "sbp", paste0("q", distr_nam)
  )[between(sbp, quantile(sbp, 0.01), quantile(sbp, 0.99))]
zz[, type := ifelse(type == "Modelled", "IMPACTncd", "HSE")]
zz[, year := year + 2000]
zz[, qimd := factor(qimd,
                    levels = c("1 most deprived", "2", "3", "4", "5 least deprived"),
                    labels = c("1\nmost deprived", "2", "3", "4", "5\nleast deprived"))]


plot_synthpop_val(zz, sbp, c("year","agegrp10"), "wt_nurse", "sbp", xlab_nam, "short", xangle = 90, FALSE, TRUE)
plot_synthpop_val(zz, sbp, c("qimd","agegrp10"), "wt_nurse", "sbp", xlab_nam, "long", xangle = 90, FALSE, TRUE)



#smok_cess
xps_mdl <- qread("./secure_data/lifecourse_models/smok_cess_model.qs")
dt <- xps_mdl$data
distr_nam <- xps_mdl$family[[1]]

xlab_nam <- expression(bold(Smoking~Cessation))
dt[, age := round(age*15.8+43.8)]

smok_cess_tbl <- read_fst("./inputs/exposure_distributions/smok_cess_table.fst", as.data.table =  TRUE)
zz <-
  validate_gamlss_tbl(dt[age >= 20L], smok_cess_tbl, 50, "event", paste0("q", distr_nam)
  )[between(event, quantile(event, 0.01), quantile(event, 0.99))]
zz[, type := ifelse(type == "Modelled", "IMPACTncd", "HSE")]
zz[, year := year + 2000]
zz[, qimd := factor(qimd,
                    levels = c("1 most deprived", "2", "3", "4", "5 least deprived"),
                    labels = c("1\nmost deprived", "2", "3", "4", "5\nleast deprived"))]


plot_synthpop_val(zz, event, c("year","agegrp10"), "wt_int", "smok_cess", xlab_nam, "short", xangle = 90, FALSE, TRUE)
plot_synthpop_val(zz, event, c("qimd","agegrp10"), "wt_int", "smok_cess", xlab_nam, "long", xangle = 90, FALSE, TRUE)





#smok_incid - Is this the right one? these look wrong....
xps_mdl <- qread("./secure_data/lifecourse_models/smok_incid_model.qs")
dt <- xps_mdl$data
distr_nam <- xps_mdl$family[[1]]

xlab_nam <- expression(bold("Smoking Incidence"))
dt[, age := round(age*17.4 + 51.7)]

smok_incid_tbl <- read_fst("./inputs/exposure_distributions/smok_incid_table.fst", as.data.table =  TRUE)
zz <-
  validate_gamlss_tbl(dt[age >= 20], smok_incid_tbl, 50, "event", paste0("q", distr_nam)
  )[between(event, quantile(event, 0.01), quantile(event, 0.99))]
zz[, type := ifelse(type == "Modelled", "IMPACTncd", "HSE")]
zz[, year := year + 2000]
zz[, qimd := factor(qimd,
                    levels = c("1 most deprived", "2", "3", "4", "5 least deprived"),
                    labels = c("1\nmost deprived", "2", "3", "4", "5\nleast deprived"))]


plot_synthpop_val(zz, event, c("year","agegrp10"), "wt_int", "smok_incid", xlab_nam, "short", xangle = 90, FALSE, TRUE)
plot_synthpop_val(zz, event, c("qimd","agegrp10"), "wt_int", "smok_incid", xlab_nam, "long", xangle = 90, FALSE, TRUE)


#smok_status - why doesn't it join...
xps_mdl <- qread("./secure_data/lifecourse_models/smok_status_model.qs")
dt <- xps_mdl$data
distr_nam <- xps_mdl$family[[1]]

xlab_nam <- expression(bold(Smoking~Status))
dt[, age := round(age*17.4 + 50.7)]

smok_status_tbl <- read_fst("./inputs/exposure_distributions/smok_status_table.fst", as.data.table =  TRUE)
zz <-
  validate_gamlss_tbl(dt, smok_status_tbl, 50, "smok_status", paste0("q", distr_nam)
  )[between(smok_status, quantile(smok_status, 0.01), quantile(smok_status, 0.99))]
zz[, type := ifelse(type == "Modelled", "IMPACTncd", "HSE")]
zz[, year := year + 2000]
zz[, qimd := factor(qimd,
                    levels = c("1 most deprived", "2", "3", "4", "5 least deprived"),
                    labels = c("1\nmost deprived", "2", "3", "4", "5\nleast deprived"))]


plot_synthpop_val(zz, smok_status, c("year","agegrp10"), "wt_int", "smok_status", xlab_nam, "short", xangle = 0, FALSE, TRUE)
plot_synthpop_val(zz, smok_status, c("qimd","agegrp10"), "wt_int", "smok_status", xlab_nam, "long", xangle = 0, FALSE, TRUE)


#statin_px - the validate gives a weird error
xps_mdl <- qread("./secure_data/lifecourse_models/statin_px_model.qs")
dt <- xps_mdl$data
distr_nam <- xps_mdl$family[[1]]

xlab_nam <- expression(bold("Statins prescriptions"))
dt[, age := round(age* 16.6 + 52.4)]

statins_tbl <- read_fst("./inputs/exposure_distributions/statin_px_table.fst", as.data.table =  TRUE)
dt[, tchol:= round(clamp(tchol, 2, 12), 0)]
zz <-
  validate_gamlss_tbl(dt[age >= 20L], statins_tbl, 50, "statin_px", paste0("q", distr_nam)
  )[between(statin_px, quantile(statin_px, 0.01), quantile(statin_px, 0.99))]
zz[, type := ifelse(type == "Modelled", "IMPACTncd", "HSE")]
zz[, year := year + 2000]
zz[, qimd := factor(qimd,
                    levels = c("1 most deprived", "2", "3", "4", "5 least deprived"),
                    labels = c("1\nmost deprived", "2", "3", "4", "5\nleast deprived"))]


plot_synthpop_val(zz, statin_px, c("year","agegrp10"), "wt_blood", "statin_px", xlab_nam, "long", xangle = 90, FALSE, TRUE)
plot_synthpop_val(zz, statin_px, c("qimd","agegrp10"), "wt_blood", "statin_px", xlab_nam, "long", xangle = 90, FALSE, TRUE)



#smok_dur_ex
xps_mdl <- qread("./secure_data/lifecourse_models/smok_dur_ex_model.qs")
dt <- xps_mdl$data
distr_nam <- xps_mdl$family[[1]]

xlab_nam <- expression(bold("Years~of~smoking~(exsmokers)"))
dt[, age := round(age*17+56.2)]

smok_dur_ex_tbl <- read_fst("./inputs/exposure_distributions/smok_dur_ex_table.fst", as.data.table =  TRUE)
zz <-
  validate_gamlss_tbl(dt, smok_dur_ex_tbl, 50, "smok_dur_ex", paste0("q", distr_nam)
  )[between(smok_dur_ex, quantile(smok_dur_ex, 0.01), quantile(smok_dur_ex, 0.99))]
zz[, type := ifelse(type == "Modelled", "IMPACTncd", "HSE")]
zz[, year := year + 2000]
zz[, qimd := factor(qimd,
                    levels = c("1 most deprived", "2", "3", "4", "5 least deprived"),
                    labels = c("1\nmost deprived", "2", "3", "4", "5\nleast deprived"))]


plot_synthpop_val(zz, smok_dur_ex, c("year","agegrp10"), "wt_int", "smok_dur_ex", xlab_nam, "short", xangle = 0, FALSE, TRUE)
plot_synthpop_val(zz, smok_dur_ex, c("qimd","agegrp10"), "wt_int", "smok_dur_ex", xlab_nam, "long", xangle = 0, FALSE, TRUE)



#smok_dur_curr
xps_mdl <- qread("./secure_data/lifecourse_models/smok_dur_curr_model.qs")
dt <- xps_mdl$data
distr_nam <- xps_mdl$family[[1]]

xlab_nam <- expression(bold(Years~of~smoking~(current~smokers)))
dt[, age := round(age*15.3 + 45)]

smok_dur_curr_tbl <- read_fst("./inputs/exposure_distributions/smok_dur_curr_table.fst", as.data.table =  TRUE)
zz <-
  validate_gamlss_tbl(dt, smok_dur_curr_tbl, 50, "smok_dur_curr", paste0("q", distr_nam)
  )[between(smok_dur_curr, quantile(smok_dur_curr, 0.01), quantile(smok_dur_curr, 0.99))]
zz[, type := ifelse(type == "Modelled", "IMPACTncd", "HSE")]
zz[, year := year + 2000]
zz[, qimd := factor(qimd,
                    levels = c("1 most deprived", "2", "3", "4", "5 least deprived"),
                    labels = c("1\nmost deprived", "2", "3", "4", "5\nleast deprived"))]


plot_synthpop_val(zz, smok_dur_curr, c("year","agegrp10"), "wt_int", "smok_dur_curr", xlab_nam, "short", xangle = 0, FALSE, TRUE)
plot_synthpop_val(zz, smok_dur_curr, c("qimd","agegrp10"), "wt_int", "smok_dur_curr", xlab_nam, "long", xangle = 0, FALSE, TRUE)




#smok_cig_ex
xps_mdl <- qread("./secure_data/lifecourse_models/smok_cig_ex_model.qs")
dt <- xps_mdl$data
distr_nam <- xps_mdl$family[[1]]

xlab_nam <- expression(bold("Cigarettes per day (ex-smokers)"))
dt[, age := round(age*15.3 + 45)]

smok_cig_ex_tbl <- read_fst("./inputs/exposure_distributions/smok_cig_curr_table.fst", as.data.table =  TRUE)
zz <-
  validate_gamlss_tbl(dt[age >=20], smok_cig_ex_tbl, 50, "smok_cig_ex", paste0("q", distr_nam)
  )[between(smok_cig_ex, quantile(smok_cig_ex, 0.01), quantile(smok_cig_ex, 0.99))]
zz[, type := ifelse(type == "Modelled", "IMPACTncd", "HSE")]
zz[, year := year + 2000]
zz[, qimd := factor(qimd,
                    levels = c("1 most deprived", "2", "3", "4", "5 least deprived"),
                    labels = c("1\nmost deprived", "2", "3", "4", "5\nleast deprived"))]


plot_synthpop_val(zz, smok_cig_ex, c("year","agegrp10"), "wt_int", "smok_cig_ex", xlab_nam, "short", xangle = 0, FALSE, TRUE)
plot_synthpop_val(zz, smok_cig_ex, c("qimd","agegrp10"), "wt_int", "smok_cig_ex", xlab_nam, "long", xangle = 0, FALSE, TRUE)





#smok_cig_curr
xps_mdl <- qread("./secure_data/lifecourse_models/smok_cig_curr_model.qs")
dt <- xps_mdl$data
distr_nam <- xps_mdl$family[[1]]

xlab_nam <- expression(bold("Cigarettes per day (current~smokers)"))
dt[, age := round(age*15.3 + 45)]

smok_cig_curr_tbl <- read_fst("./inputs/exposure_distributions/smok_cig_curr_table.fst", as.data.table =  TRUE)
zz <-
  validate_gamlss_tbl(dt, smok_cig_curr_tbl, 50, "smok_cig_curr", paste0("q", distr_nam)
  )[between(smok_cig_curr, quantile(smok_cig_curr, 0.01), quantile(smok_cig_curr, 0.99))]
zz[, type := ifelse(type == "Modelled", "IMPACTncd", "HSE")]
zz[, year := year + 2000]
zz[, qimd := factor(qimd,
                    levels = c("1 most deprived", "2", "3", "4", "5 least deprived"),
                    labels = c("1\nmost deprived", "2", "3", "4", "5\nleast deprived"))]


plot_synthpop_val(zz, smok_cig_curr, c("year","agegrp10"), "wt_int", "smok_cig_curr", xlab_nam, "short", xangle = 0, FALSE, TRUE)
plot_synthpop_val(zz, smok_cig_curr, c("qimd","agegrp10"), "wt_int", "smok_cig_curr", xlab_nam, "long", xangle = 0, FALSE, TRUE)



#smok_quit_yrs
xps_mdl <- qread("./secure_data/lifecourse_models/smok_quit_yrs_model.qs")
dt <- xps_mdl$data
distr_nam <- xps_mdl$family[[1]]

xlab_nam <- expression(bold("Years since cessation"))
dt[, age := round(age*17 + 56.2)]

smok_quit_yrs_tbl <- read_fst("./inputs/exposure_distributions/smok_quit_yrs_table.fst", as.data.table =  TRUE)
zz <-
  validate_gamlss_tbl(dt, smok_quit_yrs_tbl, 50, "smok_quit_yrs", paste0("q", distr_nam)
  )[between(smok_quit_yrs, quantile(smok_quit_yrs, 0.01), quantile(smok_quit_yrs, 0.99))]
zz[, type := ifelse(type == "Modelled", "IMPACTncd", "HSE")]
zz[, year := year + 2000]
zz[, qimd := factor(qimd,
                    levels = c("1 most deprived", "2", "3", "4", "5 least deprived"),
                    labels = c("1\nmost deprived", "2", "3", "4", "5\nleast deprived"))]


plot_synthpop_val(zz, smok_quit_yrs, c("year","agegrp10"), "wt_int", "smok_quit_yrs", xlab_nam, "short", xangle = 0, FALSE, TRUE)
plot_synthpop_val(zz, smok_quit_yrs, c("qimd","agegrp10"), "wt_int", "smok_quit_yrs", xlab_nam, "long", xangle = 0, FALSE, TRUE)
