library(piggyback)
library(data.table)

tag <- "v0.0.3"
repo <- "ChristK/IMPACTncd_Engl"

args <- ifelse(interactive(), getwd(), commandArgs(TRUE))
filindx <- fread(file.path(args, "/aux/filindx.csv"), key = "orig_file")

pb_download(filindx$sanit_file, dest = file.path(args, filindx$rel_dir), repo = repo, tag = tag, overwrite = !interactive())
file.rename(file.path(args, filindx$rel_dir, filindx$sanit_file), file.path(args, filindx$rel_dir, filindx$orig_file))

