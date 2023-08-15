if (!require(piggyback)) {
  install.packages("piggyback")
  library(piggyback)
}
if (!require(data.table)) {
  install.packages("data.table")
  library(data.table)
}

# Note consider making this a method for Simulation class

tag <- "v0.0.3"
repo <- "ChristK/IMPACTncd_Engl"
# pb_new_release(repo, tag) # Only need to run the first time a github repo is created
# pb_release_delete(repo, tag)


fl <- list.files(getwd(), pattern = ".fst$|.xls$|.xlsx$|.qs$", full.names = TRUE, recursive = TRUE)
fl <- grep("secure_data", fl, value = TRUE, invert = TRUE)
fl <- grep("tmp.qs$", fl, value = TRUE, invert = TRUE)

filindx <- data.table(orig_file = basename(fl),
                        sanit_file = gsub("[^[:alnum:]&&^\\.]", ".", basename(fl)), # replace all non alphanumerics by ... except .
                        abs_dir = dirname(fl),
                        rel_dir = gsub(getwd(), "", dirname(fl)), key = "orig_file")

if (any(duplicated(filindx$sanit_file))) stop("Duplicated filenames found")

fwrite(filindx, "./aux/filindx.csv")

# NOTE due to the renaming the following 3 lines are risky if the process get
# interrupted before the final rename to original names. Should be easy to
# recover id this happens though and I prefer it rather than copy the files
# before renaming.
file.rename(file.path(filindx$abs_dir, filindx$orig_file), file.path(filindx$abs_dir, filindx$sanit_file))
pb_upload(file.path(filindx$abs_dir, filindx$sanit_file), repo = repo, tag = tag)
file.rename(file.path(filindx$abs_dir, filindx$sanit_file), file.path(filindx$abs_dir, filindx$orig_file))


# pb_list(repo = repo, tag = tag)

# pb_delete(file = filindx["lsoa_to_locality_indx.fst", sanit_file][[2]],
#           repo = repo,
#           tag = tag)

