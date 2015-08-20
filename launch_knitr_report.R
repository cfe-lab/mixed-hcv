# We can't pass arguments directly to an R script that we want to spin to a knitr HTML doc.
# But the spun R script will persist existing R environment variables.
library(knitr)
args <- commandArgs(TRUE)
if (length(args) < 2) stop(paste0("Bad args, usage:", 
                                  "Rscript  'launch_knitr_report.R' [R script to spin] [filepath of output HTML report] [args to pass to R script]"))

knitr_script <- args[1]
output_html <- args[2]

knit_script_prefix <- unlist(strsplit(knitr_script, "\\.[rR]$", perl=TRUE))[1]

opts_chunk$set(progress = TRUE, verbose = TRUE, width=1500, tidy = FALSE, error= TRUE, warning = FALSE, message=FALSE, echo=FALSE)

spin(knitr_script, knit=FALSE)
knit2html(paste0(knit_script_prefix, ".Rmd"), stylesheet="markdown_bigwidth.css")
file.copy(paste0(knit_script_prefix, ".html"), output_html, overwrite=TRUE)