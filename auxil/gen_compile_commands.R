#!/usr/bin/env Rscript

# Regenerate a Clang-style compile_commands.json for the package C/C++ sources
# so VS Code cpptools uses the same flags as the R build: C++17, the LinkingTo
# include dirs, -DNDEBUG, -fpic and the platform hardening flags. This is more
# accurate and lower-maintenance than listing includePath entries by hand.
#
# It stays faithful by parsing a dry run of the real build (R CMD SHLIB -n)
# rather than assembling flags manually. Re-run after any change to LinkingTo,
# Makevars, the R version, or a linked package's headers:
#
#   Rscript auxil/gen_compile_commands.R
#
# Output: Rpackage/IMPACTncd_England_model_pkg/src/compile_commands.json

pkg_dir <- "Rpackage/IMPACTncd_England_model_pkg"
src_dir <- normalizePath(file.path(pkg_dir, "src"), mustWork = TRUE)
desc <- file.path(pkg_dir, "DESCRIPTION")

# LinkingTo gives the -I flags that R CMD INSTALL injects but SHLIB does not.
lt <- tryCatch(read.dcf(desc, "LinkingTo")[1, 1],
               error = function(e) NA_character_)
clink <- ""
if (!is.na(lt)) {
  pkgs <- sub("\\s*\\(.*\\)", "", trimws(strsplit(lt, ",")[[1]]))
  inc <- vapply(pkgs, function(p) system.file("include", package = p),
                character(1))
  names(inc) <- pkgs
  miss <- names(inc)[!nzchar(inc)]
  if (length(miss)) {
    warning("LinkingTo package(s) without an include dir: ",
            paste(miss, collapse = ", "), call. = FALSE)
  }
  inc <- inc[nzchar(inc)]
  clink <- paste(sprintf('-I"%s"', inc), collapse = " ")
}

# Dry-run the build on a throwaway copy. With no .o present, make must emit a
# compile line for every source (in src/ the up-to-date .o files hide them).
tmp <- tempfile("ccjson_")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
to_copy <- list.files(
  src_dir, full.names = TRUE,
  pattern = "\\.(c|cc|cpp|h|hpp|hh)$|^Makevars$"
)
invisible(file.copy(to_copy, tmp))
srcs <- list.files(tmp, pattern = "\\.(c|cc|cpp)$")
if (!length(srcs)) {
  stop("No C/C++ sources found in ", src_dir, call. = FALSE)
}

# Pass CLINK_CPPFLAGS via R's env so its spaces/quotes are not re-parsed.
old_clink <- Sys.getenv("CLINK_CPPFLAGS", unset = NA)
Sys.setenv(CLINK_CPPFLAGS = clink)
on.exit({
  if (is.na(old_clink)) {
    Sys.unsetenv("CLINK_CPPFLAGS")
  } else {
    Sys.setenv(CLINK_CPPFLAGS = old_clink)
  }
}, add = TRUE)

old_wd <- setwd(tmp)
on.exit(setwd(old_wd), add = TRUE)
lines <- system2("R", c("CMD", "SHLIB", "-n", shQuote(srcs)),
                 stdout = TRUE, stderr = TRUE)
setwd(old_wd)

# Keep only the per-file compile steps (" -c file -o file.o ").
comp <- grep("[[:space:]]-c[[:space:]]", lines, value = TRUE)
if (!length(comp)) {
  stop("Could not capture compile commands from `R CMD SHLIB -n`.\n",
       paste(lines, collapse = "\n"), call. = FALSE)
}

entries <- lapply(comp, function(cmd) {
  src <- sub(".*[[:space:]]-c[[:space:]]+([^[:space:]]+).*", "\\1", cmd)
  cmd <- sub("^\\s*(ccache|sccache)\\s+", "", cmd)
  file <- normalizePath(file.path(src_dir, src), mustWork = FALSE)
  list(file = file, command = trimws(cmd))
})

# Emit compile_commands.json (base R, with JSON string escaping).
esc <- function(s) {
  s <- gsub("\\\\", "\\\\\\\\", s)
  s <- gsub('"', '\\\\"', s)
  gsub("\t", "\\\\t", s)
}
blocks <- vapply(entries, function(e) {
  paste(
    "  {",
    sprintf('    "directory": "%s",', esc(src_dir)),
    sprintf('    "file": "%s",', esc(e$file)),
    sprintf('    "command": "%s"', esc(e$command)),
    "  }",
    sep = "\n"
  )
}, character(1))
out <- file.path(src_dir, "compile_commands.json")
writeLines(paste0("[\n", paste(blocks, collapse = ",\n"), "\n]"), out)

message(sprintf("Wrote %d compile entr%s to %s",
                length(blocks),
                if (length(blocks) == 1) "y" else "ies", out))
