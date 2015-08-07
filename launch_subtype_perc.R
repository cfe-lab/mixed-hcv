# We can't pass arguments directly to an R script that we want to spin to a knitr HTML doc.
# But the spun R script will persist existing R environment variables.
library(knitr)
args <- commandArgs(TRUE)
if (length(args) < 2) stop("Bad args, usage refdir cmpdir")

subtype_hits_csv <- args[1]
expected_mixture_csv <- args[2]
runname <- args[3]

opts_chunk$set(progress = TRUE, verbose = TRUE, width=1500, tidy = FALSE, error= TRUE, warning = FALSE, message=FALSE, echo=FALSE)

spin("subtype_perc.R", knit=FALSE)
knit2html("subtype_perc.Rmd", stylesheet="markdown_bigwidth.css")
file.copy("subtype_perc.html", paste0("/media/macdatafile/mixed-hcv/mixture_reports/", runname, ".mixture.html"), overwrite=TRUE)