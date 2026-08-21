get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)))
  }
  getwd()
}

source(file.path(get_script_dir(), "zenodo_utils.R"))

args <- commandArgs(trailingOnly = TRUE)

zenodo_preview_url <- paste0(
  "https://zenodo.org/records/22044327?preview=1&token=",
  "eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjBjOWNkZmU0LTQ0YWYtNGM1NC1iMmVjLWQ5ZjFkNWUyNTBmNSIs",
  "ImRhdGEiOnt9LCJyYW5kb20iOiI1YTk4MjI3NjQ5NjcxOTM3ZTc2NGVjODE5ZDA1NjFmOCJ9.",
  "jEVYwx-yi_9xN8WLCdwnBe7zvyHErlPygNNLPLk9NNYubLpb0I4sLPb6Dn6ZIplN8oY1ZPQy9COVF0m8BMXDzg"
)

usage <- paste(
  "Usage:",
  "  Rscript download_zenodo_inputs.R [destination-dir] [--overwrite]",
  "",
  "Examples:",
  "  Rscript download_zenodo_inputs.R",
  "  Rscript download_zenodo_inputs.R input_folder --overwrite",
  sep = "\n"
)

overwrite <- "--overwrite" %in% args
args <- args[args != "--overwrite"]

if (length(args) > 1) {
  stop(usage, call. = FALSE)
}

destination_dir <- if (length(args) == 1) args[[1]] else "input_folder"

download_zenodo_inputs(
  identifier = zenodo_preview_url,
  destination_dir = destination_dir,
  overwrite = overwrite
)
