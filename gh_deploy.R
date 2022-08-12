if (!require(piggyback)) {
  install.packages("piggyback")
  library(piggyback)
}
if (!require(data.table)) {
  install.packages("data.table")
  library(data.table)
}


tag <- "v0.0.3"
repo <- "ChristK/IMPACTncd_Engl"

args <- ifelse(interactive(), getwd(), commandArgs(TRUE))
filindx <- fread(file.path(args, "/aux/filindx.csv"), key = "orig_file")

sapply(unique(file.path(args, filindx$rel_dir)), dir.create,  showWarnings = FALSE, recursive = TRUE)

# Note there is a bug in the function. Ideally dest = file.path(args,
# filindx$rel_dir, filindx$sanit_file) should be dest = file.path(args,
# filindx$rel_dir). The latter only works if scalar
pb_download(file = filindx$sanit_file,
            dest = file.path(args, filindx$rel_dir, filindx$sanit_file),
            repo = repo, tag = tag, overwrite = interactive(),
            use_timestamps = FALSE)
file.rename(file.path(args, filindx$rel_dir, filindx$sanit_file), file.path(args, filindx$rel_dir, filindx$orig_file))

