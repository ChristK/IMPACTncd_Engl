source("./global.R")

IMPACTncd <- Simulation$new("./inputs/sim_design.yaml")

# g <- IMPACTncd$get_causal_structure(print_plot = TRUE)
# g <- IMPACTncd$get_causal_structure(processed = FALSE, print_plot = TRUE, focus = "chd")

#plot(igraph::make_ego_graph(g, order = 1, c("pain"), "in")[[1]])

g <- IMPACTncd$get_causal_structure(processed = FALSE, print_plot = TRUE)
c(g)[[9]][[3]]$name  # this gives all the names of the vertices - these are what we want to rename
# all <- c(g)[[9]][[3]]$name
# outcome <- c("breast_ca", "chd", "colorectal_ca", "nonmodelled", "stroke" , "af" , "constipation" ,
#          "t2dm"  ,        "alcpr"    ,     "psychosis",     "andep"   ,      "asthma"  ,
#        "ckd", "dementia" , "other_ca"    ,  "pain" ,         "hf" ,
#        "copd", "ctd", "epilepsy", "lung_ca" ,      "helo" ,
#        "htn"   ,        "ibs"       ,       "prostate_ca"  , "ra"     , "t1dm" )
# exposure <- all[!all %in% outcome ]
#
# g <- set_vertex_attr(g, name = "type", index = outcome, value = "outcome")
# g <- set_vertex_attr(g, name = "type", index = exposure, value = "exposure")
#
# all_edge <- E(g)
# exposure_edge <- E(g)[.inc(V(g)[name %in% exposure])]
# other_edge <- all_edge[!all_edge %in% exposure_edge]
# g <- set_edge_attr(g, name = "rel_type", index = other_edge, value = "outcome_edge")
# g <- set_edge_attr(g, name = "rel_type", index = exposure_edge, value = "exposure_edge")

names <- c("Active days", "Breast\ncancer", "CHD", "Colorectal\ncancer", "Nonmodelled\nmortality", "Stroke",
              "Atrial\nfibrillation", "Constipation", "Alcohol\nintake", "Diabetes\nT2", "Alcohol\nproblems",
              "Psychosis", "Anxiety &\nDepression", "Asthma", "BMI", "CKD", "Dementia", "Obesity", "Other\ncancers",
              "Pain", "Heart failure", "COPD", "Connective Tissue\nDisorders", "Epilepsy", "Second-hand\nsmoke",
              "Lung\ncancer", "Fruit\nintake", "Hearing\nloss", "Hypertension", "IBS", "Metabolic equivalent\ntime",
           "Prostate\ncancer", "Rheumatoid\nArthritis", "SBP", "Smoking", "Statins", "Diabetes\nT1", "Total\ncholesterol",
              "Vegetable\nintake")
V(g)$name <- names

library(RColorBrewer)
pal <- brewer.pal(length(unique(V(g)$type)), "Set1")
pal[3] <- "black"


#Making them more translucent for the edges
pal_e <- alpha(pal[1:2], 0.5)

all <- c(g)[[9]][[3]]$name
outcome <- c("Breast\ncancer", "CHD", "Colorectal\ncancer", "Nonmodelled\nmortality",
             "Stroke" , "Atrial\nfibrillation" , "Constipation" ,
             "Diabetes\nT2"  ,        "Alcohol\nproblems"    ,
             "Psychosis",     "Anxiety &\nDepression"   ,      "Asthma"  ,
             "CKD", "Dementia" , "Other\ncancers"    ,  "Pain" ,         "Heart failure" ,
             "COPD", "Connective Tissue\nDisorders", "Epilepsy", "Lung\ncancer" ,      "Hearing\nloss" ,
             "Hypertension"   ,        "IBS"       ,       "Prostate\ncancer"  , "Rheumatoid\nArthritis"
             , "Diabetes\nT1" )
exposure <- all[!all %in% outcome ]

g <- set_vertex_attr(g, name = "type", index = outcome, value = "outcome")
g <- set_vertex_attr(g, name = "type", index = exposure, value = "exposure")
V(g)[V(g)$type == "outcome"]$color <- pal[2]
V(g)[V(g)$type == "exposure"]$color <- pal[1]

all_edge <- E(g)
exposure_edge <- E(g)[.inc(V(g)[name %in% exposure])]
other_edge <- all_edge[!all_edge %in% exposure_edge]
g <- set_edge_attr(g, name = "rel_type", index = other_edge, value = "outcome_edge")
g <- set_edge_attr(g, name = "rel_type", index = exposure_edge, value = "exposure_edge")
E(g)[E(g)$rel_type == "outcome_edge"]$color <- pal_e[2]
E(g)[E(g)$rel_type == "exposure_edge"]$color <- pal_e[1]

#I manually changed the igraph function to adjust the self-arrow size
# trace("plot.igraph",edit=TRUE)
# #changed the cp to this:
# cp <- matrix(c(x0, y0, x0 + 0.2, y0 + 0.1, x0 + 0.2,
#                y0 - 0.1, x0, y0), ncol = 2, byrow = TRUE)


mod_plot <- function(x, focus = FALSE, circlesize = 35, fontsize = 1.5){
  if (missing(focus)) {
    graph <- x
  }
  else{
    if (length(focus) > 1L)
      stop("focus need to be scalar string.")
    if (!focus %in% names)
      stop("focus need to be an assigned name.")
    graph <- make_ego_graph(x,
                            order = 1, nodes = focus, mode = "in")[[1]]
    graph <- set_vertex_attr(graph, name = "type", index = focus, value = "Z")
    V(graph)[V(graph)$type == "Z"]$color <- pal[3]

  }
  plot.igraph(graph,
       vertex.shape = "none",
       #vertex.color = "NA",
       #vertex.frame.color = alpha(c(graph)[[9]][[3]]$color, 0.3),
       vertex.size = circlesize,
       vertex.label.font = 2,
       vertex.label.cex = fontsize,
       vertex.label.color = c(graph)[[9]][[3]]$color,
       edge.arrow.width = 1,
       edge.arrow.size = 0.5,
       edge.arrow.color = alpha(c(graph)[[9]][[3]]$color, 0.3),
       edge.width = 0.5,
       edge.lty = 3,
       #edge.color = pal[as.numeric(as.factor(edge_attr(graph, "rel_type")))],
       layout = layout_components,
       loop.angle = 0.3)
  }


# disnames <- c("Breast cancer", "CHD", "Colorectal cancer", "Nonmodelled mortality", "Stroke",
#            "Atrial fibrillation", "Constipation", "Diabetes T2", "Alcohol problems",
#            "Psychosis", "Anxiety & Depression", "Asthma", "CKD", "Dementia",  "Other cancers",
#            "Pain", "Heart failure", "COPD", "Connective Tissue Disorders", "Epilepsy",
#            "Lung cancer", "Hearing loss", "Hypertension", "IBS",  "Prostate cancer",
#            "Rheumatoid Arthritis", "Diabetes T1")


out_pth_plot <-
  function(x = character(0))
    paste0("/mnt/", Sys.info()[["user"]],
           "/UoL/CPRD2021/epi_models/structure_plots/", x)

# png(file = out_pth_plot(paste0("overall.png")),
#     width = 800, height = 800)
# mod_plot(g)
# dev.off()


# this might be a bit better?
png(file = out_pth_plot(paste0("overall.png")),
    width = 1800, height = 1600, res = 400, pointsize = 3)
mod_plot(g, circlesize = 20, fontsize = 1.2)
dev.off()

# SVG files for if want to edit
svg(file = out_pth_plot(paste0("overall.svg")),
    pointsize = 5)
print(mod_plot(g))
dev.off()


# for (i in outcome){ #not all of these diseases have relationships.
#   #There should be a way to discard conditions that don't, but I can't work out how
#   png(file = out_pth_plot(paste0(i, ".png")),
#       width = 1600, height = 1600)
#   mod_plot(g, i)
#   dev.off()
# }

for (i in outcome){ #not all of these diseases have relationships.
#This might be a bit better?
png(file = out_pth_plot(paste0(i, ".png")),
    width = 1600, height = 1600, res = 600, pointsize = 3)
mod_plot(g, i)
 dev.off()
}


#Some need different settings
i <- "Psychosis"
  png(file = out_pth_plot(paste0(i, ".png")),
      width = 1800, height = 1600, res = 600, pointsize = 3)
  mod_plot(g, i)
  dev.off()

i <- "Constipation"
png(file = out_pth_plot(paste0(i, ".png")),
    width = 1800, height = 1600, res = 600, pointsize = 3)
    mod_plot(g, i,  circlesize = 30, fontsize = 1.2)
  dev.off()

  i <- "Pain"
  png(file = out_pth_plot(paste0(i, ".png")),
      width = 1600, height = 1600, res = 600, pointsize = 3)
  mod_plot(g, i,  circlesize = 30, fontsize = 1.2)
  dev.off()

#svg files for if want to edit
  for (i in outcome){ #not all of these diseases have relationships.
    #This might be a bit better?
    svg(file = out_pth_plot(paste0(i, ".svg")), pointsize = 5)
    print(mod_plot(g, i))
    dev.off()
  }



scenario_fn <- function(sp) NULL


IMPACTncd$
  del_logs()$
  del_outputs()$
  run(1:2, multicore = TRUE, "sc0")

# IMPACTncd$export_summaries(multicore = TRUE)
# source("./auxil/CPRD_sim_validation_plots.R")

scenario_fn <- function(sp) {

  # The major risk factors included in our model in three cases have “healthy”
  # ranges according to government guidelines, body mass index (BMI), systolic
  # blood pressure and cholesterol: for example the NHS advises an adults BMI in
  # the range 18.5-24.9. In this setting, the 20% improvement of risk factors is
  # relative to the distance from the middle of the range (21.7): someone with a
  # BMI of 31.7 is 10kg/m2 above the middle of the range, so a 20% improvement
  # is a 2kg/m2 reduction to 29.7. The NHS advised ranges for blood pressure and
  # cholesterol are 90/60mmHg - 120/80mmHg, for total cholesterol the
  # recommendation is below 5mmol/L, but above 1mmol/L for “good cholesterol”.
  # We model improvements in the other risk factors as follows: fruit and
  # vegetable intake is 20% higher than the base case for each person; alcohol
  # intake is 20% lower; 20% of the population become more active by one active
  # day per week; smoking prevalence drops by 20% (relative); smoking
  # consumption for those that still smoke drops by 20%; passive smoking
  # prevalence drops by 20% (relative).
  sc_year <- 23L # The year the change starts
  change <- 0.2 # positive means health improvement. Do not set to 0

  bmi_target <- mean(c(18.5, 24.9))
  sp$pop[year >= sc_year & bmi_curr_xps > bmi_target, bmi_curr_xps := bmi_curr_xps + (bmi_target - bmi_curr_xps) * change]

  sbp_target <- mean(c(90, 120))
  sp$pop[year >= sc_year & sbp_curr_xps > sbp_target, sbp_curr_xps := sbp_curr_xps + (sbp_target - sbp_curr_xps) * change]

  tchol_target <- mean(c(3, 5)) # Lower bound arbitrary but high enough to ensure HDL > 1mmol/L possible
  sp$pop[year >= sc_year & tchol_curr_xps > tchol_target, tchol_curr_xps := tchol_curr_xps + (tchol_target - tchol_curr_xps) * change]



  sp$pop[year >= sc_year, alcohol_curr_xps := as.integer(round(alcohol_curr_xps * (1 - change)))]
  sp$pop[year >= sc_year, fruit_curr_xps := as.integer(round(fruit_curr_xps * (1 + change)))]
  sp$pop[year >= sc_year, veg_curr_xps := as.integer(round(veg_curr_xps * (1 + change)))]

  # (20%) of the population increases physical activity by 1 day and met by 20%
  tt <- sp$pop[year >= sc_year, .(unique(pid))]
  tt <- tt[as.logical(rbinom(.N, 1, abs(change))), V1]
  sp$pop[pid %in% tt & year >= sc_year, `:=` (
    active_days_curr_xps = clamp(x = active_days_curr_xps + as.integer(sign(change)), a = 0L, b = 7L, inplace = FALSE),
    met_curr_xps = as.integer(round(met_curr_xps * (1 + change)))
    )]

  # smoking
  if (change > 0) { # Note positive means better health

    sp$pop[year >= sc_year & smok_status_curr_xps == "4",
       hc_eff := rbinom(.N, 1L, change)]

    sp$pop[, (c("smok_status_curr_xps", "smok_quit_yrs_curr_xps", "smok_dur_curr_xps", "smok_cig_curr_xps")) :=
         simsmok_policy_impact_decr(
           smok_status_curr_xps,
           smok_quit_yrs_curr_xps,
           smok_dur_curr_xps,
           smok_cig_curr_xps,
           pid_mrk,
           hc_eff
         )]
  } else {
    # calculate policy effect with those quit smoking recently be more
    # likely to relapse
    tt <- sp$pop[year >= sc_year, .("ex"   = sum(smok_status_curr_xps == "3"),
                                "curr" = sum(smok_status_curr_xps == "4")), keyby = year]
    tt[, impacted := round(curr * (1 - change))] # Note change here is -ve so 1 - change is an increase

    # Make change to add up every year (for Vincy's SCC abstract)
    # tt[, impacted := round(curr * scenario_parms$sc_str_smk_change *
    #     (year - min(year) + 1L))]

    sp$pop[tt, `:=`(impacted = i.impacted,
                ex = i.ex), on = "year"]
    sp$pop[year >= sc_year & smok_status_curr_xps == "3",
       rid := 1:.N, by = year]
    sp$pop[, hc_eff := 0L]
    tt <- sp$pop[year >= sc_year & smok_status_curr_xps == "3",
             .(rid = sample_int_expj(first(ex), first(impacted),
                                     (smok_quit_yrs_curr_xps + 1L) ^
                                       -1)),
             keyby = year]
    sp$pop[tt, hc_eff := 1L, on = .(year, rid)]
    sp$pop[, c("impacted", "ex", "rid") := NULL]

    sp$pop[, (c("smok_status_curr_xps", "smok_quit_yrs_curr_xps", "smok_dur_curr_xps")) :=
         simsmok_policy_impact_incr(
           smok_status_curr_xps,
           smok_quit_yrs_curr_xps,
           smok_dur_curr_xps,
           pid_mrk,
           hc_eff
         )]
  }
  sp$pop[, hc_eff := NULL]

  sp$pop[year >= sc_year, smok_cig_curr_xps := as.integer(round(smok_cig_curr_xps * (1 - change)))]

  sp$pop[, smok_packyrs_curr_xps := as.integer(round(smok_cig_curr_xps * smok_dur_curr_xps / 20))]


  # ets need to be after smoking
  if (change > 0) {
    sp$pop[year >= sc_year & ets_curr_xps == 1L & smok_status_curr_xps != "4" & rbinom(.N, 1, abs(change)),
           ets_curr_xps := 0L]
  } else {
    sp$pop[year >= sc_year & ets_curr_xps == 0L & smok_status_curr_xps != "4" & rbinom(.N, 1, abs(change)),
           ets_curr_xps := 1L]
  }

}


IMPACTncd$
  run(1:2, multicore = TRUE, "sc1")$
  export_summaries(multicore = TRUE)

# IMPACTncd$export_summaries(multicore = TRUE)
source("./auxil/process_out_for_HF.R")
source("./auxil/CPRD_sim_validation_plots_CK.R")
# Bus error (core dumped)
