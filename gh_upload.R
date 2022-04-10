library(piggyback)
tag <- "v0.0.1"
repo <- "ChristK/IMPACTncd_Engl"
# pb_new_release(repo, tag)
fl <- list.files(getwd(), pattern = ".fst$", full.names = TRUE, recursive = TRUE)
fl <- grep("secure_data", fl, value = TRUE, invert = TRUE)
pb_track(fl)

pb_upload(pb_track(), repo = repo, tag = tag)

pb_list(repo = repo, tag = tag)
