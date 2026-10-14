#
# readme2yaml.R -- write project.yaml for the bicdelivery permissions
# database (init_impact_project_permissions.py)
#
# usage:
#   Rscript readme2yaml.R PROJNO ROOT RUNFOLDER GENOME
#
#   PROJNO      project number, without Proj_ (e.g. 18143)
#   ROOT        the project's delivery folder, the parent of RUNFOLDER
#   RUNFOLDER   r_NNN
#   GENOME      build tag of the run (GENOME= in 00.RUNSTATUS.txt)
#
# The PI and the investigator are the user names (the part before @,
# lowercased) of "Lab head / PI email" and "Your email" in the first
# README.txt found going up from the current directory to Proj_<PROJNO>.
# Called by deliverResults.sh. Adapted from adagio's bin/readme2yaml.R,
# which looks only in . and .., derives the root from README.txt and
# the genome from its Species.
#

cat("## START: readme2yaml.R\n")

suppressPackageStartupMessages(library(tidyverse))

PIPELINE <- "atacseq"

fatal <- function(msg) {
  cat(str_glue("\n    FATAL ERROR: {msg}\n\n"))
  quit(status = 1)
}

argv <- commandArgs(trailingOnly = TRUE)
if (len(argv) != 4) {
  cat("\nusage: Rscript readme2yaml.R PROJNO ROOT RUNFOLDER GENOME\n\n")
  quit(status = 1)
}
projno <- argv[1]
root <- argv[2]
runfolder <- argv[3]
genome <- argv[4]

if (!str_detect(runfolder, "^r_\\d+$")) {
  fatal(str_glue("RUNFOLDER [{runfolder}] is not r_NNN"))
}
if (!str_starts(root, "/")) {
  fatal(str_glue("ROOT [{root}] is not an absolute path"))
}

path_parts <- fs::path_split(getwd())[[1]]
proj_idx <- which(path_parts == str_glue("Proj_{projno}"))
if (len(proj_idx) == 0) {
  fatal(str_glue("no Proj_{projno} directory in {getwd()}"))
}
readme_file <- rev(max(proj_idx):len(path_parts)) |>
  map_chr(\(i) fs::path_join(path_parts[seq_len(i)])) |>
  fs::path("README.txt") |>
  keep(fs::file_exists) |>
  head(1)
if (len(readme_file) == 0) {
  fatal(str_glue("no README.txt from {getwd()} up to Proj_{projno}"))
}
cat(str_glue("Processing {readme_file} ...\n\n"))

readme <- tibble(line = read_lines(readme_file)) |>
  filter(str_detect(line, ": ")) |>
  separate_wider_delim(line, ": ", names = c("key", "value"), too_many = "merge") |>
  mutate(across(everything(), str_trim))

#' The user name of the email address in one README.txt field
#'
#' @param field README.txt key, e.g. "Your email"
#' @return the part before @, lowercased
email_user <- function(field) {
  value <- readme |>
    filter(key == field, value != "") |>
    pull(value)
  if (len(value) != 1) {
    fatal(str_glue("{readme_file} has {len(value)} values for [{field}], not 1"))
  }
  value |>
    str_remove("@.*$") |>
    str_to_lower()
}

project <- list(
  root = root,
  runfolder = runfolder,
  project = projno,
  pi = email_user("Lab head / PI email"),
  invest = email_user("Your email"),
  genome = genome,
  pipeline = PIPELINE
)

yaml::write_yaml(project, "project.yaml")
cat(yaml::as.yaml(project))

cat("## END: readme2yaml.R\n")
