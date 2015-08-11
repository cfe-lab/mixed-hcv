# Find any difference between subtype mixture calls when we use any hit against HCV subtype full genome or 
# when the alignments are filtered for regions that are predetermined for best genotyping
# from http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0122082#pone-0122082-t002
#


library(plyr)
library(knitr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)
print(args)


opts_chunk$set(progress = TRUE, verbose = TRUE, width=1500, tidy = FALSE, error= TRUE, warning = FALSE, message=FALSE, echo=FALSE)
#options(width=200)

# fullgen_hits_csv <- "/media/macdatafile/mixed-hcv/gb-ref+hg38_v2/150720_M01841_0148_000000000-AE93J.gb-ref+hg38_v2.csv"
# slice_hits_csv <- "/home/thuy/gitrepo/mixed-hcv/out/150720_M01841_0148_000000000-AE93J/150720_M01841_0148_000000000-AE93J.deli.csv"
# expected_mixture_csv <- "/media/macdatafile/mixed-hcv/expected_mixture/150720_M01841_0148_000000000-AE93J.expected_mixture.csv"
# runname <- "150720_M01841_0148_000000000-AE93J"


fullgen_hits_csv <- args[2]
slice_hits_csv <- args[3]
expected_mixture_csv <- args[4]
runname <- args[5]

#'
#' Compare HCV Full Genome Hits Vs Target Region Hits  for Run `r runname`
#' -------------------------------------------------------
#' 

expected_mixture <- read.table(expected_mixture_csv, sep=",", header=TRUE)
expected_mixture$gtype <- as.factor(as.character(expected_mixture$gtype))
colnames(expected_mixture)[grep("perc", colnames(expected_mixture))] <- paste0(colnames(expected_mixture)[grep("perc", colnames(expected_mixture))], ".exp")

fullgen_hits <- read.table(fullgen_hits_csv, sep=",", header=TRUE)
fullgen_hits <- subset(fullgen_hits, fullgen_hits$subtype != "" & grepl("hg38", fullgen_hits$subtype) == FALSE, select=-perc)

# find total reads that hit HCV
fullgen_total_hcv <- ddply(.data=fullgen_hits,
                   .variables=c("runname", "sample", "snum", "total"),
                   .fun=function(x) {
                     data.frame(
                       totalhcv.fullgen=sum(x$count[x$subtype != "" & grepl("hg38", x$subtype) == FALSE], na.rm=TRUE)
                     )
                   })


# runname,sample,snum,gene,start,end,rank,count,subtype,nucseq,aaseq
slice_hits <- read.table(slice_hits_csv, sep=",", header=TRUE, na.strings=c(""))

# Defensive programming - aggregate over the subtypes, 
# in case we encounter an earlier version of the CSV that split count entries 
# for the same sample-subtype-targetregion by full read sequence
slice_hits <- ddply(.data=slice_hits,
                        .variables=c("runname", "sample", "snum", "subtype", "gene", "start", "end"),
                        .fun=function(x) {
                          data.frame(
                            count=sum(x$count[x$subtype != "" & grepl("hg38", x$subtype) == FALSE], na.rm=TRUE)
                          )
                        })


# find total reads that hit HCV in the target regions
slice_total_hcv <- ddply(.data=slice_hits,
                   .variables=c("runname", "sample", "snum"),
                   .fun=function(x) {
                     data.frame(
                       totalhcv.slice=sum(x$count[x$subtype != "" & grepl("hg38", x$subtype) == FALSE], na.rm=TRUE)
                     )
                   })

# Aggregate total reads that hit subtype
slice_subtype_hits <- ddply(.data=slice_hits,
                         .variables=c("runname", "sample", "snum", "subtype"),
                         .fun=function(x) {
                           data.frame(
                             totalsubtype.slice=sum(x$count[x$subtype != "" & grepl("hg38", x$subtype) == FALSE], na.rm=TRUE)
                           )
                         })




# Merge in the full genome and slice data frames to compare subtype % according to each method
cmp_subtype <- merge(x=subset(fullgen_hits, select=-total),
                     y=slice_subtype_hits,
                     all=TRUE)
colnames(cmp_subtype)[grep("count", colnames(cmp_subtype))] <- "totalsubtype.fullgen"
cmp_subtype$totalsubtype.fullgen[is.na(cmp_subtype$totalsubtype.fullgen)] <- 0
cmp_subtype$totalsubtype.slice[is.na(cmp_subtype$totalsubtype.slice)] <- 0


# Merge in the total hcv hits for each method so that we can calculate % subtype for each method
cmp_subtype <- merge(x=cmp_subtype,
                     y=fullgen_total_hcv,
                     all=TRUE)

cmp_subtype <- merge(x=cmp_subtype,
                     y=slice_total_hcv,
                     all=TRUE)

cmp_subtype$perc.fullgen <- cmp_subtype$totalsubtype.fullgen * 100 / cmp_subtype$totalhcv.fullgen
cmp_subtype$perc.slice <- cmp_subtype$totalsubtype.slice * 100 / cmp_subtype$totalhcv.slice

# hack to get same column names to plot diff dataframes on same ggplot
expected_mixture_hack <-  expected_mixture
expected_mixture_hack$perc.fullgen <- expected_mixture$subtype_perc.exp 
expected_mixture_hack$perc.slice <- expected_mixture$subtype_perc.exp 


#' **Circles are actual subtype %.  Triangles are expected subtype %.**
#' **Black line is y=x**
#' 
#+ fig.width=15, fig.height=12
fig <- ggplot(cmp_subtype, aes(x=perc.fullgen, y=perc.slice)) + 
  geom_abline(slope=1) + 
  geom_point(shape=1, size=10, color="red") +   
  geom_point(data=expected_mixture_hack, aes(x=perc.fullgen, y=perc.slice), 
             shape=2, size=10, color="Blue") + 
  xlab("Subtype % [Genome HCV Hits]") + 
  ylab("Subtype % [Target Region HCV Hits]") + 
  ggtitle("Compare Subtype % Calculated Via Full Genome Hits Vs Target Region Hits") + 
  #scale_shape_discrete(guide=FALSE) + 
  #scale_color_discrete(guide=guide_legend(nrow=3) ) + 
  theme(strip.text = element_text(size=rel(1.2)), 
        axis.title=element_text(size=rel(2)), 
        axis.text=element_text(size=rel(1.2), angle=45, hjust=1),
        plot.title=element_text(size=rel(2)),
        legend.title=element_text(size=rel(2)),
        legend.text=element_text(size=rel(1.5)),
        legend.position="top") +   
  facet_wrap(~sample)
print(fig)

