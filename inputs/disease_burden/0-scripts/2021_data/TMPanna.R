write_xps_tmplte_file <-
  function(exposure, exposurenm, disnm, type, file_path = output_path("/template.csvy")) {
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
      "xps_name"     = paste0(exposurenm, type), 
      "outcome"       = disnm,
      "lag"          = 1L,
      "distribution" =  "lognormal",
      "source"        = "CPRD",
      "notes"         = paste0("other covariate conditions: ", paste(dpnd2[!dpnd2 %in% exposurenm], collapse = ","))
    )
    
    effect <- CJ(
      agegroup = factor(CKutils::agegrp_name(min_age = 0, max_age = 99, grp_width = 5)),
      sex = c("men", "women"),
      dimd = factor(c("1 most deprived", "2", "3", "4", "5", "6", "7", "8", "9", "10 least deprived"),
                    levels = c("1 most deprived", "2", "3", "4", "5", "6", "7", "8", "9", "10 least deprived"))
    )
    
    tmp <- results[dpnd_on == exposure,
                      .(sex, agegroup, dimd,  rr = round(rr, digits = 2), ci_rr =round(ci_rr, digits = 2))]
    effect <- merge(effect, tmp, by = c("agegroup","sex", "dimd"), all = T)
    setkey(effect, sex, dimd, agegroup)
    
    effect[is.na(rr), `:=` (rr = 1, ci_rr =1)] #putting in RR = 1 for the missing age groups 
    
    write_xps_prm_file(effect, metadata, file_path)
  }

results <-   qread(output_path(paste0("dependencies/",disnm, "_dpndRR.qs")), nthreads = 10)
results <- results[, .(sex, agegroup, dimd, rr = round(rr, digits = 2), ci_rr = round(uci_rr, digits = 2), dpnd_on )]

type <- "_prev"


for(j in 1:length(dpnd2)){
  write_xps_tmplte_file(dpnd[j],dpnd2[j], disnm2, type, output_path(paste0("dependencies/",dpnd2[j],"~",disnm2,".csvy")))
}
rm(j)
