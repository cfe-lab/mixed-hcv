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
# 
# fullgen_hits_csv <- "/media/macdatafile/mixed-hcv/gb-ref+hg38_v2/150731_M01841_0153_000000000-AE96E.csv"
# slice_hits_csv <- "/media/macdatafile/mixed-hcv/gb-ref+hg38_v2/deli/150731_M01841_0153_000000000-AE96E.csv"
# expected_mixture_csv <- "/media/macdatafile/mixed-hcv/expected_mixture/150731_M01841_0153_000000000-AE96E.expected_mixture.csv"
# runname <- "150731_M01841_0153_000000000-AE96E"


fullgen_hits_csv <- args[2]
slice_hits_csv <- args[3]
expected_mixture_csv <- args[4]
runname <- args[5]

#'
#' Run `r runname`  :  What is the Difference in HCV Subtype % When Full Genome Aligments Are Unfiltered vs Filtered for Target Regions? 
#' -------------------------------------------------------
#' 
#' Target Regions are those determined by a previous study to differentiate between genotypes:  
#' http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0122082#pone-0122082-t002
#' 
#' Only NS5a, NS5b, NS3 target regions are used as filters.
#' 
expected_mixture <- read.table(expected_mixture_csv, sep=",", header=TRUE, na.strings=c(""))
expected_mixture$gtype <- as.character(expected_mixture$gtype)
colnames(expected_mixture)[grep("perc", colnames(expected_mixture))] <- paste0(colnames(expected_mixture)[grep("perc", colnames(expected_mixture))], ".exp")
expected_mixture$resolved_subtype <-  ifelse(is.na(expected_mixture$subtype) | expected_mixture$subtype == "", 
                                                 as.character(expected_mixture$gtype), 
                                                 as.character(expected_mixture$subtype))
expected_mixture$resolved_subtype_perc.exp <- ifelse(is.na(expected_mixture$subtype) | expected_mixture$subtype == "", 
                                                expected_mixture$gtype_perc, 
                                                expected_mixture$subtype_perc.exp)

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

# Put in genotype info
cmp_subtype$gtype <- substr(cmp_subtype$subtype, 1, 1)

total_gtype <- ddply(.data=cmp_subtype,
                     .variables=c("runname", "sample", "snum", "gtype"),
                     .fun=function(x) {
                       data.frame(
                         totalgtype.fullgen = sum(x$totalsubtype.fullgen, na.rm=TRUE),
                         totalgtype.slice = sum(x$totalsubtype.slice, na.rm=TRUE)
                       )                       
                     })

cmp_subtype <- merge(x=cmp_subtype, y=total_gtype, all=TRUE)


# calc subtype, genotype percentages of HCV population
cmp_subtype$subtype_perc.fullgen <- cmp_subtype$totalsubtype.fullgen * 100 / cmp_subtype$totalhcv.fullgen
cmp_subtype$subtype_perc.slice <- cmp_subtype$totalsubtype.slice * 100 / cmp_subtype$totalhcv.slice
cmp_subtype$gtype_perc.fullgen <- cmp_subtype$totalgtype.fullgen * 100 / cmp_subtype$totalhcv.fullgen
cmp_subtype$gtype_perc.slice <- cmp_subtype$totalgtype.slice * 100 / cmp_subtype$totalhcv.slice


# Merge in expected to compare them
cmp_subtype_exp <- adply(.data=cmp_subtype,
                         .margins=1,
                         .fun=function(x) {
                           exp_sample <- subset(expected_mixture, expected_mixture$runname == x$runname &
                                                                    expected_mixture$sample == x$sample)
                           
                           is_exp_subtype <- FALSE
                           resolved_subtype <- "Unexpected"
                           if (x$subtype %in% exp_sample$resolved_subtype) {
                             is_exp_subtype <- TRUE
                             resolved_subtype <- as.character(x$subtype)
                           } else if (x$gtype %in% exp_sample$resolved_subtype) {
                             is_exp_subtype <- TRUE
                             resolved_subtype <- as.character(x$gtype)
                           } 
                           
                            # If the expected subtype is left empty, then any subtype is OK as long as it falls under the expected genotype
                           resolved_subtype = 
                           data.frame(
                             is_exp_gtype = x$gtype %in% exp_sample$gtype,                             
                             is_exp_subtype = is_exp_subtype,
                             resolved_subtype = resolved_subtype
                           )
                         })
cmp_subtype_exp$resolved_subtype <- as.factor(cmp_subtype_exp$resolved_subtype)

cmp_subtype_resolved_exp <- ddply(.data=cmp_subtype_exp,
                                  .variables=c("runname", "sample", "snum", "resolved_subtype", "totalhcv.fullgen", "totalhcv.slice"),
                                  .fun=function(x) {
                                    data.frame(
                                      total_resolved_subtype.fullgen = sum(x$totalsubtype.fullgen, na.rm=TRUE),
                                      total_resolved_subtype.slice = sum(x$totalsubtype.slice, na.rm=TRUE),
                                      resolved_subtype_perc.fullgen = 100 * sum(x$totalsubtype.fullgen, na.rm=TRUE) / x$totalhcv.fullgen[1],
                                      resolved_subtype_perc.slice = 100 * sum(x$totalsubtype.slice, na.rm=TRUE) / x$totalhcv.slice[1]
                                    )
                                  })

  
# hack to get same column names to plot diff dataframes on same ggplot
expected_mixture_hack <-  expected_mixture
expected_mixture_hack$resolved_subtype_perc.fullgen <- expected_mixture$resolved_subtype_perc.exp
expected_mixture_hack$resolved_subtype_perc.slice <- expected_mixture$resolved_subtype_perc.exp


#' **Circles are actual subtype %.  Triangles are expected subtype %.**
#' **Black line is y=x**
#' 
#+ fig.width=15, fig.height=12
colourCount <-  length(unique(cmp_subtype_resolved_exp$resolved_subtype))
getPalette <- colorRampPalette(brewer.pal(min(colourCount, 8), "Accent"))  # returns a function that takes number of colors as argument
fig <- ggplot(cmp_subtype_resolved_exp, 
              aes(x=resolved_subtype_perc.fullgen, y=resolved_subtype_perc.slice, color=resolved_subtype)) + 
  geom_abline(slope=1) + 
  geom_point(data=expected_mixture_hack, aes(x=resolved_subtype_perc.fullgen, y=resolved_subtype_perc.slice), 
             shape=17, size=8) + 
  geom_point(size=8, shape=1, pch=34) +   
  xlab("Subtype % [Full Genome HCV Hits]") + 
  ylab("Subtype % [Target Region HCV Hits]") + 
  ggtitle("Compare Subtype % Derived From Full Genome Hits Vs Target Region Hits") + 
  scale_color_manual(name="Expected Type", values = getPalette(colourCount)) +
  theme(strip.text = element_text(size=rel(1.2)), 
        axis.title=element_text(size=rel(2)), 
        axis.text=element_text(size=rel(1.2), angle=45, hjust=1),
        plot.title=element_text(size=rel(2)),
        legend.title=element_text(size=rel(2)),
        legend.text=element_text(size=rel(1.5)),
        legend.position="top") +   
  facet_wrap(~sample)
print(fig)

