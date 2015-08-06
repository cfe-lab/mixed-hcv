#!/usr/bin/Rscript --vanilla

## Collect some information on mapping quality scores in a MiSeq run.

## TRUE prints only the arguments we specify; FALSE specifies lots of
## other things like --slave, --no-restore, etc.
## Template:
## ./map_quality_analysis.R [map quality CSV] [output plot file] [sample name]
args <- commandArgs(TRUE)

## This has columns qname, flag, and mapq.  We want to make a histogram
## of mapq.
mapqs <- read.csv(args[1])

middle.quantiles <- quantile(mapqs$mapq, c(0.25, 0.5, 0.75))

pdf(args[2])
map.quals <- hist(mapqs$mapq, breaks=-0.5:(ceiling(max(mapqs$mapq))+0.5),
                  main=paste("Mapping quality for sample", args[3]),
                  ylab="frequency", xlab="score")

quantile.colours <- rainbow(3, alpha=0.5)
abline(v=middle.quantiles, col=quantile.colours)
legend("topright",
       fill=quantile.colours,
       legend=c(paste("25th percentile:", middle.quantiles[1]),
         paste("median:", middle.quantiles[2]),
         paste("75th percentile:", middle.quantiles[3])))
dev.off()
