library(plyr)
library(knitr)
library(ggplot2)
library(reshape2)

MAX_HCV_PERC_DIFF <- 1
MIN_REPORT_PERC <- 1

opts_chunk$set(progress = TRUE, verbose = TRUE, width=1500, tidy = FALSE, error= TRUE, warning = FALSE, message=FALSE, echo=FALSE)
options(width=200)


jul20_hcvhuman <- read.table("/media/macdatafile/mixed-hcv/gb-ref+hg38_v2/150720_M01841_0148_000000000-AE93J.gb-ref+hg38_v2.csv",
                        sep=",", header=TRUE)
# summary(jul20_hcvhuman)
# head(jul20_hcvhuman)
# dim(jul20_hcvhuman)

jul20_hcv <- read.table("/media/macdatafile/mixed-hcv/gb-ref/150720_M01841_0148_000000000-AE93J.gb-ref.csv",
                        sep=",", header=TRUE)
# summary(jul20_hcv)
# head(jul20_hcv)
# dim(jul20_hcv)

cmp20 <- merge(x=jul20_hcvhuman,
             y=jul20_hcv,
             by.x=c("runname", "sample", "snum", "subtype", "total"),
             by.y=c("runname", "sample", "snum", "subtype", "total"),
             all=TRUE,
             suffixes=c(".20HcvHum", ".20Hcv"))
cmp20$count.20Hcv[is.na(cmp20$count.20Hcv)] <- 0
cmp20$count.20HcvHum[is.na(cmp20$count.20HcvHum)] <- 0


popperc.20 <- ddply(.data=cmp20,
                          .variables=c("runname", "sample", "snum", "total"),
                          .fun=function(x) {
                            data.frame(
                              totalhcv.20HcvHum=sum(x$count.20HcvHum[x$subtype != "" & grepl("hg38", x$subtype) == FALSE], na.rm=TRUE),
                              totalhcv.20Hcv=sum(x$count.20Hcv[x$subtype != "" & grepl("hg38", x$subtype) == FALSE], na.rm=TRUE)
                              )
                          })

# summary(popperc.20)
# head(popperc.20)
# dim(popperc.20)

cmp20 <- merge(x=cmp20,
               y=popperc.20,
               all=TRUE)
cmp20$hcvperc.20HcvHum <- ifelse(is.na(cmp20$count.20HcvHum*100/cmp20$totalhcv.20HcvHum), 0, cmp20$count.20HcvHum*100/cmp20$totalhcv.20HcvHum)
cmp20$hcvperc.20Hcv <- ifelse(is.na(cmp20$count.20Hcv*100/cmp20$totalhcv.20Hcv), 0, cmp20$count.20Hcv*100/cmp20$totalhcv.20Hcv)
cmp20$hcvperc.20HcvHum[cmp20$subtype == "" | grepl("hg38", cmp20$subtype) == TRUE] <- NA
cmp20$hcvperc.20Hcv[cmp20$subtype == "" | grepl("hg38", cmp20$subtype) == TRUE] <- NA
cmp20$hcvperc_diff_20HcvHum_20Hcv <- cmp20$hcvperc.20HcvHum - cmp20$hcvperc.20Hcv
#summary(cmp20)
#head(cmp20)
#dim(cmp20)



write.table(cmp20, "/media/macdatafile/mixed-hcv/cmp_jul20_Hcv_HcvHuman.csv", sep=",", row.names=FALSE, na="")

#' July 20: Can We Recover Expected HCV Population Mixtures From Sequencing?  (Ref=HCV+Human)
#' ===================================================================
#' 
#' **Mostly Yes.  Most samples match expected predominant mixtures, except for 56585A-HCV which is expected 1b/2b but is only 2b**
#' 
nice_cmp20 <- cmp20[grepl("hg38", cmp20$subtype) == FALSE & cmp20$subtype != "",]
nice_cmp20$nice_subtype <- as.character(nice_cmp20$subtype)
nice_cmp20$nice_subtype[nice_cmp20$hcvperc.20HcvHum < MIN_REPORT_PERC] <- "other"
nice_cmp20$nice_subtype <- as.factor(nice_cmp20$nice_subtype)
nice_cmp20 <- ddply(.data=nice_cmp20,
                    .variables=c("runname", "sample", "snum", "total", "nice_subtype"),
                    .fun=function(x) {
                      data.frame(hcvperc.20HcvHum=sum(x$hcvperc.20HcvHum, na.rm=TRUE))
                    })


#+ fig.width=20, fig.height=8
fig <- ggplot(nice_cmp20, 
              aes(x=sample, y=hcvperc.20HcvHum, color=nice_subtype)) + 
  geom_point(size=10) + 
  scale_color_brewer(name="HCV Subtype", type="qual") + 
  ggtitle("July 20: HCV Population Subtype % Per Sample") + 
  xlab("\nSample") + 
  ylab("% HCV Population") + 
  theme(strip.text.x = element_text(size=rel(1)), 
        axis.title=element_text(size=rel(2)), 
        axis.text=element_text(size=rel(1.5)),
        plot.title=element_text(size=rel(2)),
        legend.title=element_text(size=rel(2)),
        legend.text=element_text(size=rel(1.5)))
print(fig)

#+ results="asis"
kable(subset(nice_cmp20, select=-c(runname, total))[order(nice_cmp20$sample, -nice_cmp20$hcvperc.20HcvHum),], 
      format="html",
      caption="July 20: HCV Population Subtype % Per Sample",
      row.names=FALSE,
      col.names=c("Sample", "Sample Num", "Subtype", "% HCV Population"))


#' July 20: Does Reported Subtype % Change Between Reference=HCV+Human Vs Reference=HCV?
#' =====================================================================================
#' 
#' **Reference Choice Makes a Big Difference in Subtype Percentage in HCV Population for Subtypes 1a, 1b, 2a**
#' 
#' **The ratio of 1a:2a decreases drastically for all the 5LOG samples.**
#' 
#' **Including Human in the reference reduces the reads aligned to HCV, which means that bowtie (alignment software) thought
#' certain Human Reads were actually HCV**



#+ fig.width=20, fig.height=8
# IQR in % coverage change due to change in reference by subtype
fig <- ggplot(cmp20[grepl("hg38", cmp20$subtype) == FALSE & cmp20$subtype != "",], 
              aes(x=subtype, y=hcvperc_diff_20HcvHum_20Hcv)) + 
  geom_boxplot(fill="red") + 
  ggtitle("July 20: Reference Makes Big Difference in Subtype % in HCV Population") + 
  xlab("HCV Subtype") + 
  ylab("Difference in % HCV Reads Aligned When \nRef=HCV+Human Vs Ref=HCV") + 
  theme_bw(base_size = 12) + 
  theme(panel.grid.major = element_line(size = .3, color = "lightgray")) + 
  theme(strip.text.x = element_text(size=rel(1)), 
        axis.title=element_text(size=rel(2)), 
        axis.text=element_text(size=rel(1.5), angle=45, hjust=1),
        plot.title=element_text(size=rel(2)),
        legend.title=element_text(size=rel(2)),
        legend.text=element_text(size=rel(1.5)))
print(fig)


# 
# bigchange_jul20 <- cmp20[abs(cmp20$hcvperc_diff_20HcvHum_20Hcv) > MAX_HCV_PERC_DIFF & 
#                            grepl("hg38", cmp20$subtype)==FALSE & cmp20$subtype != "" ,]
# 
# #+ results="asis"
# kable(subset(bigchange_jul20, select=-c(runname, perc.20HcvHum, perc.20Hcv)), 
#       format="html", 
#       caption=paste0("July 20 Sample Subtypes That Change > ", MAX_HCV_PERC_DIFF, "% of HCV Population, Depending on Reference"), 
#       row.names=FALSE, padding=2, 
#       col.names=c("Sample", "SampleNum", "Subtype", "Sample Reads", "Read Hits: Ref=HCV+Human", "Read Hits:Ref=HCV",
#                   "HCV Reads:Ref=HCV+Human", "HCV Reads:Ref=HCV", "% HCV Reads:Ref=HCV+Human", "% HCV Reads:Ref=HCV", 
#                   "% Diff HCV Reads"))


cmp20_melt <- reshape2:::melt(cmp20[grepl("hg38", cmp20$subtype) == FALSE & cmp20$subtype != "" &
                                      !(cmp20$hcvperc.20HcvHum < MIN_REPORT_PERC &
                                          cmp20$hcvperc.20Hcv < MIN_REPORT_PERC),], 
                              id.vars=c("runname", "sample", "snum", "total", "subtype"),
                              measure.vars=c("hcvperc.20HcvHum", "hcvperc.20Hcv"),
                              variable.name="ref",
                              value.name="hcvperc")

#+ fig.width=20, fig.height=10
fig <- ggplot(cmp20_melt, 
              aes(x=subtype, y=hcvperc, color=ref)) + 
  facet_wrap(~sample) + 
  geom_point(size=10) + 
  scale_color_discrete(name="Reference", labels=c("HCV+Human", "HCV")) + 
  ggtitle("July 20: Reference choice causes drastic change in 1a:2a ratio in 5LOG samples") + 
  xlab("\nPredominant HCV Subtypes") + 
  ylab("% HCV Reads Aligned") +
  theme_bw(base_size = 12) + 
  theme(panel.grid.major = element_line(size = .3, color = "lightgray")) + 
  theme(strip.text.x = element_text(size=rel(2)), 
        axis.title=element_text(size=rel(2)), 
        axis.text=element_text(size=rel(1.5), angle=45, hjust=1),
        plot.title=element_text(size=rel(2)),
        legend.title=element_text(size=rel(2)),
        legend.text=element_text(size=rel(1.5)))
print(fig)


cmp20_total_melt <- reshape2:::melt(unique(cmp20[grepl("hg38", cmp20$subtype) == FALSE & cmp20$subtype != "",
                                                 c("runname", "sample", "snum", "totalhcv.20HcvHum", "totalhcv.20Hcv")]), 
                              id.vars=c("runname", "sample", "snum"),
                              measure.vars=c("totalhcv.20HcvHum", "totalhcv.20Hcv"),
                              variable.name="ref",
                              value.name="totalhcv")

fig <- ggplot(cmp20_total_melt, 
              aes(x=sample, weight=totalhcv, fill=ref)) +   
  geom_bar(position="dodge", color="black") + 
  scale_fill_discrete(name="Reference", labels=c("HCV+Human", "HCV")) + 
  ggtitle("July 20: Including Human In Reference Reduces Reads Aligned to HCV") + 
  xlab("Sample") + 
  ylab("Total Reads Aligned to HCV") +
  theme_bw(base_size = 12) + 
  theme(panel.grid.major = element_line(size = .3, color = "lightgray")) + 
  theme(strip.text.x = element_text(size=rel(2)), 
        axis.title=element_text(size=rel(2)), 
        axis.text=element_text(size=rel(1.5)),
        plot.title=element_text(size=rel(2)),
        legend.title=element_text(size=rel(2)),
        legend.text=element_text(size=rel(1.5)))
print(fig)

#  JUL 24
#########################

jul24_hcvhuman <- read.table("/media/macdatafile/mixed-hcv/gb-ref+hg38_v2/150724_M01841_0150_000000000-AE8E0.csv",
                             sep=",", header=TRUE)
# summary(jul24_hcvhuman)
# head(jul24_hcvhuman)
# dim(jul24_hcvhuman)


jul24_hcv <- read.table("/media/macdatafile/mixed-hcv/gb-ref/150724_M01841_0150_000000000-AE8E0.csv",
                        sep=",", header=TRUE)
# summary(jul24_hcv)
# head(jul24_hcv)
# dim(jul24_hcv)

cmp24 <- merge(x=jul24_hcvhuman,
               y=jul24_hcv,
               by.x=c("runname", "sample", "snum", "subtype", "total"),
               by.y=c("runname", "sample", "snum", "subtype", "total"),
               all=TRUE,
               suffixes=c(".24HcvHum", ".24Hcv"))
cmp24$count.24Hcv[is.na(cmp24$count.24Hcv)] <- 0
cmp24$count.24HcvHum[is.na(cmp24$count.24HcvHum)] <- 0

popperc.24 <- ddply(.data=cmp24,
                          .variables=c("runname", "sample", "snum", "total"),
                          .fun=function(x) {
                            data.frame(
                              totalhcv.24HcvHum=sum(x$count.24HcvHum[x$subtype != "" & grepl("hg38", x$subtype) == FALSE], na.rm=TRUE),
                              totalhcv.24Hcv=sum(x$count.24Hcv[x$subtype != "" & grepl("hg38", x$subtype) == FALSE], na.rm=TRUE)
                            )
                          })
# summary(popperc.24)
# head(popperc.24)
# dim(popperc.24)

cmp24 <- merge(x=cmp24,
               y=popperc.24,
               all=TRUE)
cmp24$hcvperc.24HcvHum <- ifelse(is.na(cmp24$count.24HcvHum * 100 / cmp24$totalhcv.24HcvHum), 0, cmp24$count.24HcvHum * 100 / cmp24$totalhcv.24HcvHum)
cmp24$hcvperc.24Hcv <- ifelse(is.na(cmp24$count.24Hcv * 100/ cmp24$totalhcv.24Hcv), 0, cmp24$count.24Hcv * 100/ cmp24$totalhcv.24Hcv)
cmp24$hcvperc.24HcvHum[cmp24$subtype == "" | grepl("hg38", cmp24$subtype) == TRUE] <- NA
cmp24$hcvperc.24Hcv[cmp24$subtype == "" | grepl("hg38", cmp24$subtype) == TRUE] <- NA
cmp24$hcvperc_diff_24HcvHum_24Hcv <- cmp24$hcvperc.24HcvHum - cmp24$hcvperc.24Hcv
# summary(cmp24)
# head(cmp24)
# dim(cmp24)

write.table(cmp24, "/media/macdatafile/mixed-hcv/cmp_jul24_Hcv_HcvHuman.csv", sep=",", row.names=FALSE,  na="")


#' July 24: Can We Recover Expected HCV Population Mixtures Sample From Sequencing?  (Ref=HCV+Human)
#' ===================================================================
#' 
#' **Mostly Yes.  Most samples match expected predominant mixtures, except for 56585A-HCV which should have 1b/2b but instead only has 2b.**
#' 
#' 
nice_cmp24 <- cmp24[grepl("hg38", cmp24$subtype) == FALSE & cmp24$subtype != "",]
nice_cmp24$nice_subtype <- as.character(nice_cmp24$subtype)
nice_cmp24$nice_subtype[nice_cmp24$hcvperc.24HcvHum < MIN_REPORT_PERC] <- "other"
nice_cmp24$nice_subtype <- as.factor(nice_cmp24$nice_subtype)
nice_cmp24 <- ddply(.data=nice_cmp24,
                    .variables=c("runname", "sample", "snum", "total", "nice_subtype"),
                    .fun=function(x) {
                      data.frame(hcvperc.24HcvHum=sum(x$hcvperc.24HcvHum, na.rm=TRUE))
                    })


#+ fig.width=20, fig.height=8
fig <- ggplot(nice_cmp24, 
              aes(x=sample, y=hcvperc.24HcvHum, color=nice_subtype)) + 
  geom_point(size=10) + 
  scale_color_brewer(name="HCV Subtype", type="qual") + 
  ggtitle("July 24: HCV Population Subtype % Per Sample") + 
  xlab("\nSample") + 
  ylab("% HCV Population") + 
  theme(strip.text.x = element_text(size=rel(1)), 
        axis.title=element_text(size=rel(2)), 
        axis.text=element_text(size=rel(1.5)),
        plot.title=element_text(size=rel(2)),
        legend.title=element_text(size=rel(2)),
        legend.text=element_text(size=rel(1.5)))
print(fig)

#+ results="asis"
kable(subset(nice_cmp24, select=-c(runname, total))[order(nice_cmp24$sample, -nice_cmp24$hcvperc.24HcvHum),], 
      format="html",
      caption="July 24: HCV Population Subtype % Per Sample",
      row.names=FALSE,
      col.names=c("Sample", "Sample Num", "Subtype", "% HCV Population"))

#' July 24 Run:  Does Reported Subtype % Change Between Reference=HCV+Human Vs Reference=HCV?
#' =====================================================================================
#' 
#' 
#' **Yes!  Effects Are Similar to Those Found in Jul 20 Run.**
#' **Including Human Genome in Reference Makes a Big Difference for Sample Subtype HCV Population Percentage 1a, 1b, 2a.**
#' 
#' **The ratio of 1a:2a decreases drastically for all the 5LOG samples.**
#' 
#' **Including Human in the reference reduces the total reads aligned to HCV, which means that 
#' there are human reads that are similar enough to hcv to confuse bowtie**

#'


#' 
#+ fig.width=20, fig.height=8
# IQR in % coverage change due to change in reference by subtype
fig <- ggplot(cmp24[grepl("hg38", cmp24$subtype) == FALSE & cmp24$subtype != "",], 
              aes(x=subtype, y=hcvperc_diff_24HcvHum_24Hcv)) + 
  geom_boxplot(fill="red") + 
  ggtitle("July 24: Reference Makes Big Difference in Subtype % in HCV Population") + 
  xlab("HCV Subtype") + 
  ylab("Difference in % HCV Reads Aligned When \nRef=HCV+Human Vs Ref=HCV") + 
  theme_bw(base_size = 12) + 
  theme(panel.grid.major = element_line(size = .3, color = "lightgray")) + 
  theme(strip.text.x = element_text(size=rel(1)), 
        axis.title=element_text(size=rel(2)), 
        axis.text=element_text(size=rel(1.5), angle=45, hjust=1),
        plot.title=element_text(size=rel(2)),
        legend.title=element_text(size=rel(2)),
        legend.text=element_text(size=rel(1.5)))
print(fig)

# bigchange_jul24 <- cmp24[abs(cmp24$hcvperc_diff_24HcvHum_24Hcv) > MAX_HCV_PERC_DIFF & grepl("hg38", cmp24$subtype)==FALSE & cmp24$subtype != "" ,]
# 
# #+ results="asis"
# kable(subset(bigchange_jul24, select=-c(runname, perc.24HcvHum, perc.24Hcv)), 
#       format="html", 
#       caption=paste0("July 24 Sample Subtypes That Change > ", MAX_HCV_PERC_DIFF, "% of HCV Population, Depending on Reference"), 
#       row.names=FALSE, padding=2, 
#       col.names=c("Sample", "SampleNum", "Subtype", "Sample Reads", "Read Hits: Ref=HCV+Human", "Read Hits:Ref=HCV",
#                   "HCV Reads:Ref=HCV+Human", "HCV Reads:Ref=HCV", "% HCV Reads:Ref=HCV+Human", "% HCV Reads:Ref=HCV", 
#                   "% Diff HCV Reads"))


cmp24_melt <- reshape2:::melt(cmp24[grepl("hg38", cmp24$subtype) == FALSE & cmp24$subtype != "" &
                                      !(cmp24$hcvperc.24HcvHum < MIN_REPORT_PERC &
                                          cmp24$hcvperc.24Hcv < MIN_REPORT_PERC),], 
                              id.vars=c("runname", "sample", "snum", "total", "subtype"),
                              measure.vars=c("hcvperc.24HcvHum", "hcvperc.24Hcv"),
                              variable.name="ref",
                              value.name="hcvperc")


#' 
#+ fig.width=20, fig.height=10
fig <- ggplot(cmp24_melt, 
              aes(x=subtype, y=hcvperc, color=ref)) + 
  facet_wrap(~sample) + 
  geom_point(size=10) + 
  scale_color_discrete(name="Reference", labels=c("HCV+Human", "HCV")) + 
  ggtitle("July 20: Reference choice causes drastic change in 1a:2a ratio in 5LOG samples") + 
  xlab("\nPredominant HCV Subtypes") + 
  ylab("% HCV Reads Aligned") +
  theme_bw(base_size = 12) + 
  theme(panel.grid.major = element_line(size = .3, color = "lightgray")) + 
  theme(strip.text.x = element_text(size=rel(2)), 
        axis.title=element_text(size=rel(2)), 
        axis.text=element_text(size=rel(1.5), angle=45, hjust=1),
        plot.title=element_text(size=rel(2)),
        legend.title=element_text(size=rel(2)),
        legend.text=element_text(size=rel(1.5)))
print(fig)


cmp24_total_melt <- reshape2:::melt(unique(cmp24[grepl("hg38", cmp24$subtype) == FALSE & cmp24$subtype != "",
                                                 c("runname", "sample", "snum", "totalhcv.24HcvHum", "totalhcv.24Hcv")]), 
                                    id.vars=c("runname", "sample", "snum"),
                                    measure.vars=c("totalhcv.24HcvHum", "totalhcv.24Hcv"),
                                    variable.name="ref",
                                    value.name="totalhcv")

fig <- ggplot(cmp24_total_melt, 
              aes(x=sample, weight=totalhcv, fill=ref)) +   
  geom_bar(position="dodge", color="black") + 
  scale_fill_discrete(name="Reference", labels=c("HCV+Human", "HCV")) + 
  ggtitle("July 24: Including Human In Reference Reduces Reads Aligned to HCV") + 
  xlab("Sample") + 
  ylab("Total Reads Aligned to HCV") +
  theme_bw(base_size = 12) + 
  theme(panel.grid.major = element_line(size = .3, color = "lightgray")) + 
  theme(strip.text.x = element_text(size=rel(2)), 
        axis.title=element_text(size=rel(2)), 
        axis.text=element_text(size=rel(1.5)),
        plot.title=element_text(size=rel(2)),
        legend.title=element_text(size=rel(2)),
        legend.text=element_text(size=rel(1.5)))
print(fig)

#' Does Higher Cluster Density in July 20 Run Cause Different Results From July 24 Run?  (Ref=HCV+Human)
#' =====================================================================================
#' 
#' **No.  Higher cluster density in July 20 Run means each July 20 sample has more reads than the July 24 sample.**
#' 
#' **However, HCV population subtype % remain similar between Jul 20 and Jul 24 runs.**
#' 
cmp <- merge(x=subset(cmp20, select=-c(runname, total, perc.20Hcv, perc.20HcvHum, totalhcv.20HcvHum, totalhcv.20Hcv)),
             y=subset(cmp24, select=-c(runname, total, perc.24Hcv, perc.24HcvHum, totalhcv.24HcvHum, totalhcv.24Hcv)),
             by.x=c("sample", "snum", "subtype"),
             by.y=c("sample", "snum", "subtype"),
             all=TRUE,
             suffixes=c(".20", ".24"))

# merge in the totals, since alignments that exist in Jul 20 but not Jul 24 and vice versa will have empty total columns.
cmp <- merge(x=cmp, 
             y=subset(popperc.20, select=-c(runname, totalhcv.20Hcv)), 
             by=c("sample", "snum"), all=TRUE)

cmp <- merge(x=cmp, 
             y=subset(popperc.24, select=-c(runname, totalhcv.24Hcv)), 
             by=c("sample", "snum"), all=TRUE, suffixes=c(".20", ".24"))

# For all the alignments that exist in Jul 20 but not Jul 24, Fill in the counts, percentages with 0
cmp$count.20Hcv[is.na(cmp$count.20Hcv)] <- 0
cmp$count.20HcvHum[is.na(cmp$count.20HcvHum)] <- 0
cmp$count.24Hcv[is.na(cmp$count.24Hcv)] <- 0
cmp$count.24HcvHum[is.na(cmp$count.24HcvHum)] <- 0

cmp$hcvperc.24HcvHum[is.na(cmp$hcvperc.24HcvHum)] <- 0
cmp$hcvperc.24Hcv[is.na(cmp$hcvperc.24Hcv)] <- 0
cmp$hcvperc.20HcvHum[is.na(cmp$hcvperc.20HcvHum)] <- 0
cmp$hcvperc.20Hcv[is.na(cmp$hcvperc.20Hcv)] <- 0

cmp$hcvperc_diff_20HcvHum_20Hcv[is.na(cmp$hcvperc_diff_20HcvHum_20Hcv)] <- 0
cmp$hcvperc_diff_24HcvHum_24Hcv[is.na(cmp$hcvperc_diff_24HcvHum_24Hcv)] <- 0

cmp$hcvperc_diff_24HcvHuman_20HcvHuman <- cmp$hcvperc.24HcvHum - cmp$hcvperc.20HcvHum

# summary(cmp)
# head(cmp)
# dim(cmp)
#write.table(cmp, "/media/macdatafile/mixed-hcv/cmp_jul20_jul24_Hcv_HcvHuman.csv", sep=",", row.names=FALSE,  na="")


con <- file("/media/macdatafile/mixed-hcv/cmp_jul20_jul24_Hcv_HcvHuman.csv", open="wt")
writeLines("# sample: Miseq sample name", con)
writeLines("# snum: Miseq sample number", con)
writeLines("# subtype: HCV Subtype", con)
writeLines("# count.20HcvHum: Total Reads Aligned to Subtype in Jul 20 run where reference contained Human+HCV", con)
writeLines("# count.20Hcv: Total Reads Aligned to Subtype in Jul 20 run where reference contained HCV", con)
writeLines("# hcvperc.20HcvHum: HCV Population subtype % in Jul 20 run where reference contained Human+HCV.  ie 100*total reads hit subtype/total reads hit any HCV subtype", con)
writeLines("# hcvperc.20Hcv: HCV Population subtype % in Jul 20 run where reference contained HCV.  ie 100*total reads hit subtype/total reads hit any HCV subtype", con)
writeLines("# hcvperc_diff_20HcvHum_20Hcv: Difference in Jul 20 HCV Population Subtype % when reference contains human+hcv vs reference contains hcv", con)
writeLines("# count.24HcvHum: Total Reads Aligned to Subtype in Jul 24 run where reference contained Human+HCV", con)
writeLines("# count.24Hcv: Total Reads Aligned to Subtype in Jul 24 run where reference contained HCV", con)
writeLines("# hcvperc.24HcvHum: HCV Population subtype % in Jul 24 run where reference contained Human+HCV.  ie 100*total reads hit subtype/total reads hit any HCV subtype", con)
writeLines("# hcvperc.24Hcv: HCV Population subtype % in Jul 24 run where reference contained HCV.  ie 100*total reads hit subtype/total reads hit any HCV subtype", con)
writeLines("# hcvperc_diff_24HcvHum_24Hcv: Difference in Jul 24 run HCV Population Subtype % when reference contains human+hcv vs reference contains hcv", con)
writeLines("# total.20: Total Sample Reads in Jul 20 run", con)
writeLines("# total.24: Total Sample Reads in Jul 24 run", con)
writeLines("# totalhcv.20HcvHum: Total Reads Aligned to HCV in Jul 20 run against reference containing Human+HCV", con)
writeLines("# totalhcv.24HcvHum: Total Reads Aligned to HCV in Jul 24 run against reference containing Human+HCV", con)
writeLines("# hcvperc_diff_24HcvHuman_20HcvHuman: Difference in HCV Population Subtype % between Jul 24 run and Jul 20 run. Reference contains human+hcv", con)
write.table(cmp, con, sep=",", row.names=FALSE,  na="")
close(con)

# Find the effect on read count
cmp_total_melt <- reshape2:::melt(unique(cmp[grepl("hg38", cmp$subtype) == FALSE & cmp$subtype != "",
                                      c("sample", "snum", "total.20", "total.24")]), 
                              id.vars=c("sample", "snum"),
                              measure.vars=c("total.20", "total.24"),
                              variable.name="ref",
                              value.name="total")
#+ fig.width=20, fig.height=10
fig <- ggplot(cmp_total_melt, 
              aes(x=sample, weight=total, fill=ref)) + 
  geom_bar(position="dodge")+ 
  scale_fill_discrete(name="Run", labels=c("Jul 20", "Jul 24")) + 
  ggtitle("Each July 20 Sample Has More Reads Than July 24 Sample") + 
  xlab("Sample") + 
  ylab("Total Sample Reads") + 
  theme_bw(base_size = 12) + 
  theme(panel.grid.major = element_line(size = .3, color = "lightgray")) + 
  theme(strip.text.x = element_text(size=rel(1)), 
        axis.title=element_text(size=rel(2)), 
        axis.text=element_text(size=rel(1.5)),
        plot.title=element_text(size=rel(2)),
        legend.title=element_text(size=rel(2)),
        legend.text=element_text(size=rel(1.5)))
print(fig)

#' 
#' 
#+ fig.width=20, fig.height=8
# IQR in % coverage change due to change in reference by subtype
fig <- ggplot(cmp[grepl("hg38", cmp$subtype) == FALSE & cmp$subtype != "",], 
              aes(x=subtype, y=hcvperc_diff_24HcvHuman_20HcvHuman)) + 
  geom_boxplot(fill="red") + 
  ggtitle("Small Per-Sample Difference Between Jul 24 and Jul 20 Run: Ref=HCV+Human") + 
  xlab("HCV Subtype") + 
  ylab("Per Sample Difference in % HCV Reads Hits \nBetween Jul 24 and Jul 20") + 
  theme_bw(base_size = 12) + 
  theme(panel.grid.major = element_line(size = .3, color = "lightgray")) + 
  theme(strip.text.x = element_text(size=rel(1)), 
        axis.title=element_text(size=rel(2)), 
        axis.text=element_text(size=rel(1.5), angle=45, hjust=1),
        plot.title=element_text(size=rel(2)),
        legend.title=element_text(size=rel(2)),
        legend.text=element_text(size=rel(1.5)))
print(fig)


fig <- ggplot(cmp[grepl("hg38", cmp$subtype) == FALSE & cmp$subtype != "" &
                    abs(cmp$hcvperc.24HcvHum) >= MIN_REPORT_PERC &
                    abs(cmp$hcvperc.20HcvHum) >= MIN_REPORT_PERC,], 
              aes(x=subtype, y=hcvperc_diff_24HcvHuman_20HcvHuman, color=sample)) + 
  facet_wrap(~sample) + 
  guides(color=FALSE) + 
  geom_point(fill="red", size=10) + 
  ggtitle("Differences in Predominant Subtypes By Sample") + 
  xlab("Predominant HCV Subtypes") + 
  ylab("Difference in % HCV Reads Hits \nBetween Jul 24 and Jul 20") + 
  theme_bw(base_size = 12) + 
  theme(panel.grid.major = element_line(size = .3, color = "lightgray")) + 
  theme(strip.text.x = element_text(size=rel(2)), 
        axis.title=element_text(size=rel(2)), 
        axis.text=element_text(size=rel(1.5), angle=45, hjust=1),
        plot.title=element_text(size=rel(2)),
        legend.title=element_text(size=rel(2)),
        legend.text=element_text(size=rel(1.5)))
print(fig)
