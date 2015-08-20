# when the alignments are filtered for regions that are predetermined for best genotyping
# from http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0122082#pone-0122082-t002
#
# Usage:
# Rscript -e "library(knitr);  spin('subtype_slice_perc.R', knit=FALSE); knit2html('subtype_slice_perc.Rmd', stylesheet='markdown_bigwidth.css')" [subtype slice hits csv] [expected mix csv] [runname]


library(plyr)
library(knitr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)
print(args)

MIN_SUBTYPE_REPORT_PERC <- 1
MIN_GTYPE_REPORT_PERC <- 1

MAX_GTYPE_DIFF_PERC <- 30
MAX_SUBTYPE_DIFF_PERC <- 30

opts_chunk$set(progress = TRUE, verbose = TRUE, width=1500, tidy = FALSE, error= TRUE, warning = FALSE, message=FALSE, echo=FALSE)
#options(width=200)


# subtype_slice_hits_csv <- "/home/thuy/gitrepo/mixed-hcv/out/150720_M01841_0148_000000000-AE93J/150720_M01841_0148_000000000-AE93J.deli.csv"
# expected_mixture_csv <- "/media/macdatafile/mixed-hcv/expected_mixture/150720_M01841_0148_000000000-AE93J.expected_mixture.csv"
# runname <- "150720_M01841_0148_000000000-AE93J"



subtype_slice_hits_csv <- args[3]
expected_mixture_csv <- args[4]
runname <- args[5]



expected_mixture <- read.table(expected_mixture_csv, sep=",", header=TRUE, na.strings=c(""))
expected_mixture$gtype <- as.factor(as.character(expected_mixture$gtype))



# runname,sample,snum,gene,start,end,rank,count,subtype,nucseq,aaseq
subtype_slice_hits <- read.table(subtype_slice_hits_csv, sep=",", header=TRUE, na.strings=c(""))


total_subtypes <- ddply(.data=subtype_slice_hits,
                   .variables=c("runname", "sample", "snum", "subtype"),
                   .fun=function(x) {
                     data.frame(
                       subtype_count=sum(x$count[x$subtype != "" & grepl("hg38", x$subtype) == FALSE], na.rm=TRUE)
                     )
                   })


# find total reads that hit HCV in the target regions
total_hcv <- ddply(.data=subtype_slice_hits,
                   .variables=c("runname", "sample", "snum"),
                   .fun=function(x) {
                     data.frame(
                       sample_slice_count=sum(x$count[x$subtype != "" & grepl("hg38", x$subtype) == FALSE], na.rm=TRUE)
                     )
                   })

slice_hits <- merge(x=total_subtypes,
                     y=total_hcv,
                     by=c("runname", "sample", "snum"),
                     all=TRUE)


slice_hits$subtype_perc <- ifelse(slice_hits$subtype == "" | grepl("hg38", slice_hits$subtype) == TRUE,
                       NA,
                       slice_hits$subtype_count * 100/slice_hits$sample_slice_count)


# Collapse all subtypes with percentage < MIN_REPORT_PERC into "other" major subtype
# Collapse all human hits into "Human" category, all unaligned reads to "Unaligned" category, hcv hits to "HCV" category
slice_hits <- adply(.data=slice_hits,
              .margins=1,  # by row
              .fun=function(x) {
                                
                act_hcv_subtype_categ <- ""  # adply will convert NA to a numeric for some reason.  Workaround set to empty string.     
                
                if (x$subtype == "" | grepl("hg38", x$subtype) == TRUE) {
                  act_hcv_subtype_categ <- ""
                } else if (x$subtype_perc < MIN_SUBTYPE_REPORT_PERC) {
                  act_hcv_subtype_categ <- "Any Subtype < 1% in Sample"                  
                } else {
                  act_hcv_subtype_categ <- as.character(x$subtype)
                
                }
                                
                act_hcv_gtype <- ""
                if (x$subtype == "" | grepl("hg38", x$subtype) == TRUE) {                  
                  act_hcv_gtype <- ""
                } else if (x$subtype_perc < MIN_GTYPE_REPORT_PERC) {                  
                  act_hcv_gtype <- substr(x$subtype, 1, 1)
                } else {                  
                  act_hcv_gtype <- substr(x$subtype, 1, 1)
                }
                
                data.frame(act_hcv_gtype=act_hcv_gtype,                                                      
                           act_hcv_subtype_categ=act_hcv_subtype_categ
                )
              })

# Reorder levels, adply screws up level ordering
slice_hits$act_hcv_subtype_categ[slice_hits$act_hcv_subtype_categ == ""] <- NA
slice_hits$act_hcv_subtype_categ <- as.factor(slice_hits$act_hcv_subtype_categ)
slice_hits$act_hcv_subtype_categ <- factor(slice_hits$act_hcv_subtype_categ,
                                     levels=sort(levels(slice_hits$act_hcv_subtype_categ)))

slice_hits$act_hcv_gtype[slice_hits$act_hcv_gtype == ""] <- NA
slice_hits$act_hcv_gtype <- as.factor(slice_hits$act_hcv_gtype)
slice_hits$act_hcv_gtype <- factor(slice_hits$act_hcv_gtype,
                             levels=sort(levels(slice_hits$act_hcv_gtype)))


#' Target Region HCV Mixture Results For Run `r runname`
#' ===========================================================
#' 
#' 

#' Target Region HCV Population Genotype Percentages For Run `r runname`
#' -------------------------------------------------------
#'
#'
#' **Thresholds for determining if genotype mixture call is accurate:** 
#' 
#' **Max difference between expected and actual genotype % = `r MAX_GTYPE_DIFF_PERC`**
#' 
#' **Total % of unexpected genotypes that appear in sample must be less than the smallest expected genotype %**
#' 
#' **All expected genotypes must be detected at non-zero percentages.**
#' 
#' 
# Aggregate by major hcv genotype  (ignore nonhcv hits).  Genotypes < MIN_REPORT_PERC % are collapsed into "other" genotype
major_gtypes <- ddply(.data=slice_hits[!is.na(slice_hits$act_hcv_gtype),],
                      .variables=c("runname", "sample", "snum", "act_hcv_gtype"),
                      .fun=function(x) {
                        data.frame(
                          act_hcv_gtype_perc = sum(x$subtype_perc, na.rm=TRUE),
                          act_hcv_gtype_categ = ifelse (sum(x$subtype_perc, na.rm=TRUE) < MIN_GTYPE_REPORT_PERC, 
                                                        "Any Genotype < 1% in Sample", 
                                                        substr(x$subtype, 1, 1))
                        )
                      })

major_gtypes$act_hcv_gtype_categ[major_gtypes$act_hcv_gtype_categ == ""] <- NA
major_gtypes$act_hcv_gtype_categ <- as.factor(major_gtypes$act_hcv_gtype_categ)
major_gtypes$act_hcv_gtype_categ <- factor(major_gtypes$act_hcv_gtype_categ,
                                           levels=sort(levels(major_gtypes$act_hcv_gtype_categ)))


# Check if the actual genotype percentage is within +/- MAX_GTYPE_DIFF_PERC of the expected
check_gtype <- merge(x=unique(expected_mixture[, c("runname", "sample", "gtype", "gtype_perc")]),
                     y=subset(major_gtypes, !is.na(major_gtypes$act_hcv_gtype)),
                     by.x=c("runname", "sample", "gtype"),
                     by.y=c("runname", "sample", "act_hcv_gtype"),
                     all.x=TRUE, 
                     all.y=FALSE)

check_gtype$act_hcv_gtype_perc[is.na(check_gtype$act_hcv_gtype_perc)] <- 0
check_gtype$diff_gtype_perc <- abs(check_gtype$gtype_perc - check_gtype$act_hcv_gtype_perc)

sample_gtype_result <- ddply(.data=check_gtype,
                             .variables=c("runname", "sample"),
                             .fun=function(x) {
                               # each expected genotype must be within MAX_GTYPE_DIFF_PERC of actual genotype
                               # total unexpected genotype must be less than the smallest expected genotype
                               unexpected_gtype_perc <- 100 - sum(x$act_hcv_gtype_perc, na.rm=TRUE)
                               total_gtype_offrange <- sum(x$diff_gtype_perc > MAX_GTYPE_DIFF_PERC, na.rm=TRUE)
                               
                               total_nonexist_expected_gtype <- sum(!is.na(x$gtype_perc) & x$gtype_perc != "" & x$gtype_perc > 0 &
                                                                      (is.na(x$act_hcv_gtype_perc) | x$act_hcv_gtype_perc == 0))
                               
                               data.frame(
                                 unexpected_gtype_perc = unexpected_gtype_perc,
                                 total_gtype_offrange = total_gtype_offrange,
                                 total_nonexist_expected_gtype = total_nonexist_expected_gtype,
                                 is_gtype_ok=ifelse(total_gtype_offrange > 0 | 
                                                      unexpected_gtype_perc > min(x$gtype_perc) | 
                                                      total_nonexist_expected_gtype > 0,
                                                    "Wrong", 
                                                    "Right"))
                             })
sample_gtype_result$is_gtype_ok <- factor(sample_gtype_result$is_gtype_ok, levels=c("Right", "Wrong"))

# hacks to get ggplot working with 2 diff datasets.  Columns must match.
sample_gtype_result$act_hcv_gtype_categ <- NA
sample_gtype_result$act_hcv_gtype_perc <- NA

#+ fig.width=20, fig.height=10
colourCount <-  length(unique(major_gtypes$act_hcv_gtype_categ))
getPalette <- colorRampPalette(brewer.pal(min(8, colourCount) , "Accent"))  # returns a function that takes number of colors as argument
fig <- ggplot(major_gtypes,   
              aes(x=sample, weight=act_hcv_gtype_perc, fill=act_hcv_gtype_categ)) + 
  geom_bar(color="black") + 
  geom_text(data=sample_gtype_result, aes(x=sample, y=-10, label=is_gtype_ok, color=is_gtype_ok)) + 
  guides(color=FALSE) + 
  scale_color_manual(values=c("darkgreen", "Red")) + 
  scale_fill_manual(name="HCV Genotype", values = getPalette(colourCount)) +
  ggtitle("Sample HCV Population Genotype %") + 
  xlab("\nSample") + 
  ylab("% HCV Population") + 
  theme_bw(base_size=12) + 
  theme(strip.text.x = element_text(size=rel(1)), 
        axis.title=element_text(size=rel(2)), 
        axis.text=element_text(size=rel(1.2), angle=45, hjust=1),
        plot.title=element_text(size=rel(2)),
        legend.title=element_text(size=rel(2)),
        legend.text=element_text(size=rel(1.5)),
        legend.position="top")
print(fig)

major_gtypes_categ <- ddply(.data=major_gtypes,
                            .variables=c("runname", "sample", "snum", "act_hcv_gtype_categ"),
                            .fun=function(x) {
                              data.frame(act_hcv_gtype_categ_perc=sum(x$act_hcv_gtype_perc, na.rm=TRUE))
                            })
sort_major_gtypes_categ <- major_gtypes_categ[order(major_gtypes_categ$sample, -major_gtypes_categ$act_hcv_gtype_categ_perc),]
sort_major_gtypes_categ <- subset(sort_major_gtypes_categ, select=-c(runname))
#+ results="asis"
kable(sort_major_gtypes_categ, 
      format="html",
      caption="Sample Genotype Breakdown",
      row.names=FALSE,
      col.names=c("Sample", "Sample Num", "Genotype", "% HCV Population"))

#'
#' Target Region HCV Population Subtype Percentages for Run `r runname`
#' -------------------------------------------------------
#' 
#' **Thresholds for determining if subtype mixture call is accurate:** 
#' 
#' **Max difference between expected and actual genotype % = `r MAX_GTYPE_DIFF_PERC`**
#' 
#' **Total % of unexpected genotypes that appear in sample must be less than the smallest expected genotype %**
#' 
#' **All expected genotypes must be detected at non-zero percentages.**
#' 
#' **Max difference between expected and actual subtypes % = `r MAX_SUBTYPE_DIFF_PERC`**
#' 
#' **Total % of unexpected subtypes that appear in sample must be less than the smallest expected subtype or genotype %**
#' 
#' **All expected subtypes must be detected at non-zero percentages.**
#' 
#' 

# Aggregate by major hcv subtype  (ignore nonhcv hits).  Subtypes < 1 % are collapsed into "other" subtype.
major_subtypes <- ddply(.data=slice_hits[!is.na(slice_hits$act_hcv_subtype_categ),],
                        .variables=c("runname", "sample", "snum", "subtype", "act_hcv_subtype_categ"),
                        .fun=function(x) {
                          data.frame(
                            act_hcv_subtype_perc = sum(x$subtype_perc, na.rm=TRUE)
                          )
                        })



# Check if the actual subtype percentage is within +/- MAX_SUBTYPE_DIFF_PERC of the expected
check_subtype <- merge(x=expected_mixture,
                       y=subset(major_subtypes, !is.na(major_subtypes$subtype)),
                       by.x=c("runname", "sample", "subtype"),
                       by.y=c("runname", "sample", "subtype"),
                       all.x=TRUE, all.y=FALSE)

check_gtype_subtype <- merge(x=subset(check_subtype, select=-snum),
                             y=subset(check_gtype, select=-snum),
                             all=TRUE)


# If an expected subtype is NA, then that means that any subtype of parent genotype is OK.
check_gtype_subtype$act_hcv_subtype_perc[is.na(check_gtype_subtype$act_hcv_subtype_perc)] <- 0
check_gtype_subtype$diff_subtype_perc <- abs(check_gtype_subtype$subtype_perc - check_gtype_subtype$act_hcv_subtype_perc)


sample_subtype_result <- ddply(.data=check_gtype_subtype,
                               .variables=c("runname", "sample"),
                               .fun=function(x) {
                                 
                                 # if there is an expected subtype, then add the actual subtype %
                                 # If there is any subtype for a specific genotype is allowed, then add the actual gtype %
                                 any_subtype_ok_x <- unique(x[is.na(x$subtype) | x$subtype == "", ])
                                 specific_subtype_ok_x <- x[!is.na(x$subtype) & x$subtype != "", ]
                                 
                                 unexpected_subtype_gtype_perc <- 100 - 
                                   sum(specific_subtype_ok_x$act_hcv_subtype_perc, na.rm=TRUE)  - 
                                   sum(any_subtype_ok_x$act_hcv_gtype_perc, na.rm=TRUE)
                                 
                                 total_subtype_offrange <- sum(x$diff_subtype_perc > MAX_SUBTYPE_DIFF_PERC, na.rm=TRUE)
                                 
                                 total_nonexist_expected_subtype <- sum(!is.na(x$subtype_perc) & x$subtype_perc != "" & x$subtype_perc > 0 &
                                                                          (is.na(x$act_hcv_subtype_perc) | x$act_hcv_subtype_perc == 0))
                                 data.frame(
                                   
                                   unexpected_subtype_gtype_perc = unexpected_subtype_gtype_perc,
                                   total_subtype_offrange = total_subtype_offrange,
                                   total_nonexist_expected_subtype = total_nonexist_expected_subtype,
                                   
                                   is_subtype_ok=ifelse(total_subtype_offrange > 0 | 
                                                          unexpected_subtype_gtype_perc > min(x$gtype_perc, x$subtype_perc, na.rm=TRUE) |
                                                          total_nonexist_expected_subtype > 0,
                                                        "Wrong", 
                                                        "Right"))
                               })

sample_gtype_subtype_result <- merge(x=sample_subtype_result,
                                     y=sample_gtype_result,
                                     all=TRUE)

sample_gtype_subtype_result$is_ok <- as.factor(ifelse(sample_gtype_subtype_result$is_subtype_ok == "Wrong" | 
                                                        sample_gtype_subtype_result$is_gtype_ok == "Wrong",
                                                      "Wrong",
                                                      "Right"))
sample_gtype_subtype_result$is_ok <- factor(sample_gtype_subtype_result$is_ok, levels=c("Right", "Wrong"))  # ensure ordering for colors

# hack to get ggplot working witn 2 diff dataframes
sample_gtype_subtype_result$act_hcv_subtype_categ <- NA
sample_gtype_subtype_result$act_hcv_subtype_perc <- NA

#+ fig.width=20, fig.height=10
colourCount <-  length(unique(major_subtypes$act_hcv_subtype_categ))
getPalette <- colorRampPalette(brewer.pal(min(8, colourCount) , "Accent"))  # returns a function that takes number of colors as argument
fig <- ggplot(major_subtypes,
              aes(x=sample, weight=act_hcv_subtype_perc, fill=act_hcv_subtype_categ)) + 
  geom_bar(color="black") + 
  scale_fill_manual(name="HCV Subtype", values = getPalette(colourCount)) +
  geom_text(data=sample_gtype_subtype_result, aes(x=sample, y=-10, label=is_ok, color=is_ok)) + 
  guides(color=FALSE) + 
  scale_color_manual(values=c("darkgreen", "Red")) + 
  ggtitle("Sample HCV Population Subtype %") + 
  xlab("\nSample") + 
  ylab("% HCV Population") + 
  theme_bw(base_size=12) + 
  theme(strip.text.x = element_text(size=rel(1)), 
        axis.title=element_text(size=rel(2)), 
        axis.text=element_text(size=rel(1.2), angle=45, hjust=1),
        plot.title=element_text(size=rel(2)),
        legend.title=element_text(size=rel(2)),
        legend.text=element_text(size=rel(1.5)),
        legend.position="top")
print(fig)

major_subtypes_categ <- ddply(major_subtypes,
                              .variables=c("runname", "sample", "snum", "act_hcv_subtype_categ"),
                              .fun=function(x) {
                                data.frame(act_hcv_subtype_categ_perc=sum(x$act_hcv_subtype_perc, na.rm=TRUE))
                              })
sort_major_subtypes_categ <- major_subtypes_categ[order(major_subtypes_categ$sample, -major_subtypes_categ$act_hcv_subtype_categ_perc),]
sort_major_subtypes_categ <- subset(sort_major_subtypes_categ, select=-c(runname))
#+ results="asis"
kable(sort_major_subtypes_categ, 
      format="html",
      caption="Sample Subtype Breakdown",
      row.names=FALSE,
      col.names=c("Sample", "Sample Num", "Subtype", "% HCV Population"))


#' Coverage Per Target Region
#' ==================================
#' 
#' **Read only included as a hit if:**
#' 
#' **- both mates covered the entire the target region**
#' 
#' **- merged mates contain less than 20% N's after masking for low quality and discordant bases**
#' 
#' **- alignment was the best (primary) alignment for the read and didn't hit multiple references**
#'
#' 
# 1-based coordinates for labelling
#+ fig.width=20, fig.height=10
coord_cover <- ddply(.data=subtype_slice_hits,
                     .variables=c("runname", "sample", "snum", "gene", "start", "end"),
                     .fun=function(x) {
                       data.frame(cover=sum(x$count, na.rm=TRUE),
                                  coord=paste0(x$gene, ":", x$start+1, "-", x$end))
                     })

colourCount <-  length(unique(coord_cover$coord))
getPalette <- colorRampPalette(brewer.pal(min(8, colourCount) , "Accent"))  # returns a function that takes number of colors as argument
fig <- ggplot(coord_cover, aes(x=sample, y=cover, color=coord)) + 
  geom_point(size=10) + 
  scale_color_manual(name="Target Coordinates", values = getPalette(colourCount)) + 
  ggtitle("Coverage for Target Coordinates") + 
  xlab("Sample") + 
  ylab("Coverage") +   
  theme(strip.text.x = element_text(size=rel(1)), 
        axis.title=element_text(size=rel(2)), 
        axis.text=element_text(size=rel(1.2), angle=45, hjust=1),
        plot.title=element_text(size=rel(2)),
        legend.title=element_text(size=rel(2)),
        legend.text=element_text(size=rel(1.5)),
        legend.position="top")
print(fig)

