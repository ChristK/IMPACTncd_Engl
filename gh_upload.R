library(piggyback)
tag <- "v0.0.2"
repo <- "ChristK/IMPACTncd_Engl"
# pb_new_release(repo, tag) # Only need to run the first time a github repo is created
fl <- list.files(getwd(), pattern = ".fst$|.xls$|.xlsx$.qs$", full.names = TRUE, recursive = TRUE)
fl <- grep("secure_data", fl, value = TRUE, invert = TRUE)

# pb_delete(file = fl,
#           repo = repo,
#           tag = tag)

pb_upload(fl, repo = repo, tag = tag)

pb_list(repo = repo, tag = tag)
