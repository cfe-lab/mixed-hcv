# Calculates the HCV subtype percentages for each sample and compares them to the expected mixture

# ENVIRONMENT VARIABLES must be set beforehand:
# - subtype_hits_csv
# - expected_mixture_csv

library(plyr)
library(knitr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

MIN_SUBTYPE_REPORT_PERC <- 1
MIN_GTYPE_REPORT_PERC <- 1

MAX_GTYPE_DIFF_PERC <- 10
MAX_SUBTYPE_DIFF_PERC <- 10

opts_chunk$set(progress = TRUE, verbose = TRUE, width=1500, tidy = FALSE, error= TRUE, warning = FALSE, message=FALSE, echo=FALSE)
options(width=200)


# subtype_hits_csv <- "/media/macdatafile/mixed-hcv/gb-ref+hg38_v2/150731_M01841_0153_000000000-AE96E.csv"
# expected_mixture_csv <- "/media/macdatafile/mixed-hcv/expected_mixture/150731_M01841_0153_000000000-AE96E.expected_mixture.csv"
# runname <- "150731_M01841_0153_000000000-AE96E"


subtype_hits <- read.table(subtype_hits_csv, sep=",", header=TRUE)
expected_mixture <- read.table(expected_mixture_csv, sep=",", header=TRUE)
expected_mixture$gtype <- as.factor(as.character(expected_mixture$gtype))





# find total reads that hit HCV
total_hcv <- ddply(.data=subtype_hits,
                    .variables=c("runname", "sample", "snum", "total"),
                    .fun=function(x) {
                      data.frame(
                        totalhcv=sum(x$count[x$subtype != "" & grepl("hg38", x$subtype) == FALSE], na.rm=TRUE)
                      )
                    })

hits <- merge(x=subtype_hits,
             y=total_hcv,
             all=TRUE)
hits$hcvperc <- ifelse(hits$subtype == "" | grepl("hg38", hits$subtype) == TRUE,
                        NA,
                        hits$count * 100/hits$totalhcv)


# Collapse all subtypes with percentage < MIN_REPORT_PERC into "other" major subtype
# Collapse all human hits into "Human" category, all unaligned reads to "Unaligned" category, hcv hits to "HCV" category
hits <- adply(.data=hits,
              .margins=1,  # by row
              .fun=function(x) {
                category <- NA
                if (x$subtype == "") {
                  category <- "Unaligned"
                } else if (grepl("hg38", x$subtype) == TRUE) {
                  category <- "Human"
                } else  {
                  category <- "HCV"
                }
                
                act_hcv_subtype_categ <- ""  # adply will convert NA to a numeric for some reason.  Workaround set to empty string.
                if (x$subtype == "" | grepl("hg38", x$subtype) == TRUE) {
                  act_hcv_subtype_categ <- ""
                } else if (x$hcvperc < MIN_SUBTYPE_REPORT_PERC) {
                  act_hcv_subtype_categ <- "Other"
                } else {
                  act_hcv_subtype_categ <- as.character(x$subtype)
                }
                
                act_hcv_gtype_categ <- ""
                if (x$subtype == "" | grepl("hg38", x$subtype) == TRUE) {
                  act_hcv_gtype_categ <- ""
                } else if (x$hcvperc < MIN_GTYPE_REPORT_PERC) {
                  act_hcv_gtype_categ <- "Other"
                } else {
                  act_hcv_gtype_categ <- substr(x$subtype, 1, 1)
                }
                
                data.frame(act_hcv_gtype_categ=act_hcv_gtype_categ,
                           act_hcv_subtype_categ=act_hcv_subtype_categ,
                           category=category
                           )
              })
hits$category <- as.factor(hits$category)

# Reorder levels, adply screws up level ordering
hits$act_hcv_subtype_categ[hits$act_hcv_subtype_categ == ""] <- NA
hits$act_hcv_subtype_categ <- as.factor(hits$act_hcv_subtype_categ)
hits$act_hcv_subtype_categ <- factor(hits$act_hcv_subtype_categ,
                                                  levels=sort(levels(hits$act_hcv_subtype_categ)))

hits$act_hcv_gtype_categ[hits$act_hcv_gtype_categ == ""] <- NA
hits$act_hcv_gtype_categ <- as.factor(hits$act_hcv_gtype_categ)
hits$act_hcv_gtype_categ <- factor(hits$act_hcv_gtype_categ,
                                                  levels=sort(levels(hits$act_hcv_gtype_categ)))
#' HCV Mixture Results For Run `r runname`
#' ===========================================================
#' 
#' 

#' HCV Population Genotype Percentages For Run `r runname`
#' -------------------------------------------------------
#'
#'
#' **Thresholds for determining if genotype mixture call is accurate:** 
#' 
#' **Max difference between expected and actual genotype % = `r MAX_GTYPE_DIFF_PERC`**
#' 
#' **Total % of unexpected genotypes that appear in sample must be less than the smallest expected genotype %**
#' 
#' 
# Aggregate by major hcv genotype  (ignore nonhcv hits).  Genotypes < MIN_REPORT_PERC % are collapsed into "other" genotype
major_gtypes <- ddply(.data=hits[!is.na(hits$act_hcv_gtype_categ),],
                              .variables=c("runname", "sample", "snum", "act_hcv_gtype_categ"),
                              .fun=function(x) {
                                data.frame(
                                  act_hcv_gtype_categ_perc = sum(x$hcvperc, na.rm=TRUE)
                                )
                              })

# Check if the actual genotype percentage is within +/- MAX_GTYPE_DIFF_PERC of the expected
check_gtype <- merge(x=unique(expected_mixture[, c("runname", "sample", "gtype", "gtype_perc")]),
                     y=subset(major_gtypes, major_gtypes$act_hcv_gtype_categ != "Other"),
                     by.x=c("runname", "sample", "gtype"),
                     by.y=c("runname", "sample", "act_hcv_gtype_categ"),
                     all.x=TRUE, 
                     all.y=FALSE)

check_gtype$act_hcv_gtype_categ_perc[is.na(check_gtype$act_hcv_gtype_categ_perc)] <- 0
check_gtype$diff_gtype_perc <- abs(check_gtype$gtype_perc - check_gtype$act_hcv_gtype_categ_perc)

sample_gtype_result <- ddply(.data=check_gtype,
                              .variables=c("runname", "sample"),
                              .fun=function(x) {
                                # each expected genotype must be within MAX_GTYPE_DIFF_PERC of actual genotype
                                # total unexpected genotype must be less than the smallest expected genotype
                                unexpected_gtype_perc <- 100 - sum(x$act_hcv_gtype_categ_perc, na.rm=TRUE)
                                total_gtype_offrange <- sum(x$diff_gtype_perc > MAX_GTYPE_DIFF_PERC, na.rm=TRUE)
                                
                                data.frame(
                                  unexpected_gtype_perc = unexpected_gtype_perc,
                                  total_gtype_offrange = total_gtype_offrange,
                                  is_gtype_ok=ifelse(total_gtype_offrange > 0 | unexpected_gtype_perc > min(x$gtype_perc),
                                                        "Wrong", 
                                                        "Right"))
                              })
sample_gtype_result$is_gtype_ok <- as.factor(sample_gtype_result$is_gtype_ok)

# hacks to get ggplot working with 2 diff datasets.  Columns must match.
sample_gtype_result$act_hcv_gtype_categ <- NA
sample_gtype_result$act_hcv_gtype_categ_perc <- NA

#+ fig.width=20, fig.height=10
colourCount <-  length(unique(major_gtypes$act_hcv_gtype_categ))
getPalette <- colorRampPalette(brewer.pal(8, "Accent"))  # returns a function that takes number of colors as argument
fig <- ggplot(major_gtypes,   
              aes(x=sample, weight=act_hcv_gtype_categ_perc, fill=act_hcv_gtype_categ)) + 
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

sort_major_gtypes <- major_gtypes[order(major_gtypes$sample, -major_gtypes$act_hcv_gtype_categ_perc),]
sort_major_gtypes <- subset(sort_major_gtypes, select=-c(runname))
#+ results="asis"
kable(sort_major_gtypes, 
      format="html",
      caption="Sample Genotype Breakdown",
      row.names=FALSE,
      col.names=c("Sample", "Sample Num", "Genotype", "% HCV Population"))

#'
#' HCV Population Subtype Percentages for Run `r runname`
#' -------------------------------------------------------
#' 
#' **Thresholds for determining if subtype mixture call is accurate:** 
#' 
#' **Max difference between expected and actual genotype % = `r MAX_GTYPE_DIFF_PERC`**
#' 
#' **Total % of unexpected genotypes that appear in sample must be less than the smallest expected genotype %**
#' 
#' 
#' **Max difference between expected and actual subtypes % = `r MAX_SUBTYPE_DIFF_PERC`**
#' 
#' **Total % of unexpected subtypes that appear in sample must be less than the smallest expected subtype or genotype %**
#' 

# Aggregate by major hcv subtype  (ignore nonhcv hits).  Subtypes < 1 % are collapsed into "other" subtype.
major_subtypes <- ddply(.data=hits[!is.na(hits$act_hcv_subtype_categ),],
                              .variables=c("runname", "sample", "snum", "act_hcv_subtype_categ"),
                              .fun=function(x) {
                                data.frame(
                                  act_hcv_subtype_categ_perc = sum(x$hcvperc, na.rm=TRUE)
                                )
                              })



# Check if the actual subtype percentage is within +/- MAX_SUBTYPE_DIFF_PERC of the expected
check_subtype <- merge(x=expected_mixture,
                     y=subset(major_subtypes, major_subtypes$act_hcv_subtype_categ != "Other"),
                     by.x=c("runname", "sample", "subtype"),
                     by.y=c("runname", "sample", "act_hcv_subtype_categ"),
                     all.x=TRUE, all.y=FALSE)

check_gtype_subtype <- merge(x=subset(check_subtype, select=-snum),
                        y=subset(check_gtype, select=-snum),
                        all=TRUE)


# If an expected subtype is NA, then that means that any subtype of parent genotype is OK.
check_gtype_subtype$act_hcv_subtype_categ[is.na(check_gtype_subtype$act_hcv_subtype_categ)] <- 0
check_gtype_subtype$diff_subtype_perc <- abs(check_gtype_subtype$subtype_perc - check_gtype_subtype$act_hcv_subtype_categ)


sample_subtype_result <- ddply(.data=check_gtype_subtype,
                       .variables=c("runname", "sample"),
                       .fun=function(x) {
                         
                         # if there is an expected subtype, then add the actual subtype %
                         # If there is any subtype for a specific genotype is allowed, then add the actual gtype %
                         any_subtype_ok_x <- unique(x[is.na(x$subtype) | x$subtype == "", ])
                         specific_subtype_ok_x <- x[!is.na(x$subtype) & x$subtype != "", ]
                         
                         unexpected_subtype_gtype_perc <- 100 - 
                           sum(specific_subtype_ok_x$act_hcv_subtype_categ_perc, na.rm=TRUE)  - 
                           sum(any_subtype_ok_x$act_hcv_gtype_categ_perc, na.rm=TRUE)
                        
                         total_subtype_offrange <- sum(x$diff_subtype_perc > MAX_SUBTYPE_DIFF_PERC, na.rm=TRUE)
                         
                         data.frame(
                           
                           unexpected_subtype_gtype_perc = unexpected_subtype_gtype_perc,
                           total_subtype_offrange = total_subtype_offrange,
                           
                           is_subtype_ok=ifelse(total_subtype_offrange > 0 | unexpected_subtype_gtype_perc > min(x$gtype_perc, x$subtype_perc, na.rm=TRUE),
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

# hack to get ggplot working witn 2 diff dataframes
sample_gtype_subtype_result$act_hcv_subtype_categ <- NA
sample_gtype_subtype_result$act_hcv_subtype_categ_perc <- NA

#+ fig.width=20, fig.height=10
colourCount <-  length(unique(major_subtypes$act_hcv_subtype_categ))
getPalette <- colorRampPalette(brewer.pal(8, "Accent"))  # returns a function that takes number of colors as argument
fig <- ggplot(major_subtypes,
            aes(x=sample, weight=act_hcv_subtype_categ_perc, fill=act_hcv_subtype_categ)) + 
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


sort_major_subtypes <- major_subtypes[order(major_subtypes$sample, -major_subtypes$act_hcv_subtype_categ_perc),]
sort_major_subtypes <- subset(sort_major_subtypes, select=-c(runname))
#+ results="asis"
kable(sort_major_subtypes, 
      format="html",
      caption="Sample Subtype Breakdown",
      row.names=FALSE,
      col.names=c("Sample", "Sample Num", "Subtype", "% HCV Population"))

#'
#' Species Breakdown For Run `r runname`
#' ---------------------
#'

#'
#+ fig.width=20, fig.height=10
fig <- ggplot(hits, aes(x=sample, weight=count, fill=category)) + 
  geom_bar(position="dodge", color="black") + 
  scale_fill_discrete(name="Species") + 
  ggtitle("Alignment Breakdown by Species") + 
  xlab("Sample") + 
  ylab("Total Hits") + 
  theme_bw(base_size=12) + 
  theme(strip.text.x = element_text(size=rel(1)), 
        axis.title=element_text(size=rel(2)), 
        axis.text=element_text(size=rel(1.2), angle=45, hjust=1),
        plot.title=element_text(size=rel(2)),
        legend.title=element_text(size=rel(2)),
        legend.text=element_text(size=rel(1.5)),
        legend.position="top")
print(fig)
  