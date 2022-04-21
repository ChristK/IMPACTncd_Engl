library(piggyback)
tag <- "v0.0.2"
repo <- "ChristK/IMPACTncd_Engl"

if (interactive()) {
  print("interactive")
  pb_download(repo = repo, tag = tag, show_progress = TRUE)
} else { # used with Rscript
  # i.e. Rscript /root/IMPACTncd_Engl/gh_deploy.R "/root/IMPACTncd_Engl/"
  args <- commandArgs(TRUE)
  pb_download(dest = args, repo = repo, tag = tag, show_progress = FALSE)
}


#
