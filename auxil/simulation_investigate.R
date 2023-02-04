source("./global.R")

IMPACTncd <- Simulation$new("./inputs/sim_design.yaml")

# g <- IMPACTncd$get_causal_structure(print_plot = TRUE)

scenario_fn <- function(sp) NULL


IMPACTncd$
  del_logs()$
  del_outputs()$
  run(1:10, multicore = TRUE, "sc0")$
  export_summaries(multicore = TRUE, type = c("prvl", "incd", "mrtl", "dis_mrtl")) # "cms"

source("./auxil/process_out_investigate.R")
source("./auxil/CPRD_sim_validation_plots_CK.R")



# library(data.table)

# tt <- fread("/mnt/storage_fast/output/hf_real/tables/prevalence change by year (age-sex-dimd standardised).csv")
# plot(tt[disease == "alcpr_prvl", year], tt[disease == "alcpr_prvl", `prvl_rate_50.0%`], title("Alcohol Problems"))
# plot(tt[disease == "pain_prvl", year], tt[disease == "pain_prvl", `prvl_rate_50.0%`], title("Pain"))
# plot(tt[disease == "constipation_prvl", year], tt[disease == "constipation_prvl", `prvl_rate_50.0%`], title("Constipation"))
# plot(tt[disease == "asthma_prvl", year], tt[disease == "asthma_prvl", `prvl_rate_50.0%`], title("asthma"))
# plot(tt[disease == "andep_prvl", year], tt[disease == "andep_prvl", `prvl_rate_50.0%`], title("Anxiety/Depression"))
# plot(tt[disease == "cmsmm1.5_prvl", year], tt[disease == "cmsmm1.5_prvl", `prvl_rate_50.0%`], title("CMM"))



# tt <- fread("/mnt/storage_fast/output/hf_real/tables/mean CMS score change by year (age-sex-dimd standardised).csv")
# plot(tt$year, tt$`mean_cms_score_50.0%`)

# tt <- fread("/mnt/storage_fast/output/hf_real/tables/prevalence by year (age-sex-dimd standardised).csv")
# plot(tt[disease == "alcpr_prvl", year], tt[disease == "alcpr_prvl", `prvl_rate_50.0%`], title("Alcohol Problems"))
# plot(tt[disease == "pain_prvl", year], tt[disease == "pain_prvl", `prvl_rate_50.0%`], title("Pain"))
# plot(tt[disease == "cmsmm1.5_prvl", year], tt[disease == "cmsmm1.5_prvl", `prvl_rate_50.0%`], title("CMM"))

# tt <- fread("/mnt/storage_fast/output/hf_real/tables/prevalence by year (age-sex-dimd standardised).csv")
# plot(tt[disease == "alcpr_prvl", year], tt[disease == "alcpr_prvl", `prvl_rate_50.0%`], title("Alcohol Problems"))
# plot(tt[disease == "pain_prvl", year], tt[disease == "pain_prvl", `prvl_rate_50.0%`], title("Pain"))
# plot(tt[disease == "constipation_prvl", year], tt[disease == "constipation_prvl", `prvl_rate_50.0%`], title("Constipation"))
# plot(tt[disease == "asthma_prvl", year], tt[disease == "asthma_prvl", `prvl_rate_50.0%`], title("asthma"))
# plot(tt[disease == "andep_prvl", year], tt[disease == "andep_prvl", `prvl_rate_50.0%`], title("Anxiety/Depression"))
# plot(tt[disease == "cmsmm1.5_prvl", year], tt[disease == "cmsmm1.5_prvl", `prvl_rate_50.0%`], title("CMM"))
