library(plyr)
library(knitr)

MAX_PERC_DIFF <- 0.01

opts_chunk$set(progress = TRUE, verbose = TRUE, width=1500, tidy = FALSE, error= TRUE, warning = FALSE, message=FALSE, echo=FALSE)
options(width=200)

#' Compare July 20 Run:  Reference Contains HCV+Human Vs Reference Contains HCV Only
#' =====================================================================================
#' 
#' 
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
             by.x=c("runname", "sample", "snum", "subtype"),
             by.y=c("runname", "sample", "snum", "subtype"),
             all=TRUE,
             suffixes=c(".20HcvHum", ".20Hcv"))

cmp20$perc.20HcvHum[is.na(cmp20$perc.20HcvHum)] <- 0
cmp20$perc.20Hcv[is.na(cmp20$perc.20Hcv)] <- 0
cmp20$perc_diff_20HcvHum_20Hcv <- cmp20$perc.20HcvHum - cmp20$perc.20Hcv

# summary(cmp20)
# head(cmp20)
# dim(cmp20)

write.table(cmp20, "/media/macdatafile/mixed-hcv/cmp_jul20_Hcv_HcvHuman.csv", sep=",", row.names=FALSE, na="")

#' 
#' 
# summary(cmp20[abs(cmp20$perc_diff_20HcvHum_20Hcv) > MAX_PERC_DIFF & grepl("hg38", cmp20$subtype)==FALSE & cmp20$subtype != "" ,])
# dim(cmp20[abs(cmp20$perc_diff_20HcvHum_20Hcv) > MAX_PERC_DIFF & grepl("hg38", cmp20$subtype)==FALSE & cmp20$subtype != "" ,])

bigchange_jul20 <- cmp20[abs(cmp20$perc_diff_20HcvHum_20Hcv) > MAX_PERC_DIFF & grepl("hg38", cmp20$subtype)==FALSE & cmp20$subtype != "" ,]

#+ results="asis"
kable(subset(bigchange_jul20, select=-runname), 
      format="html", caption="July 20 Sample Subtypes That Change > 0.01% Depending on Reference", row.names=FALSE, padding=2)


#'
#' **Including Human Genome in Reference Makes a Small Difference for Sample Subtype Percentage 1a, 1b, 2a**
#' 
#' 
#+ fig.width=20, fig.height=8
# IQR in % coverage change due to change in reference by subtype
fig <- ggplot(cmp[grepl("hg38", cmp$subtype) == FALSE & cmp$subtype != "",], 
              aes(x=subtype, y=perc_diff)) + 
  geom_boxplot(fill="red") + 
  ggtitle("July 20: Reference Makes Small Difference in Subtype %") + 
  xlab("HCV Subtype") + 
  ylab("Difference in % Reads Aligned When \nRef=HCV+Human Vs Ref=HCV") + 
  theme_bw(base_size = 12) + 
  theme(panel.grid.major = element_line(size = .3, color = "lightgray")) + 
  theme(strip.text.x = element_text(size=rel(1)), 
        axis.title=element_text(size=rel(2)), 
        axis.text=element_text(size=rel(1.5), angle=45, hjust=1),
        plot.title=element_text(size=rel(2)),
        legend.title=element_text(size=rel(2)),
        legend.text=element_text(size=rel(1.5)))
print(fig)
                
# #' Compare July 24 Run:  Reference Contains HCV+Human Vs Reference Contains HCV Only
# #' =====================================================================================
# #' 
# #' 
# 
# 
# jul24_hcvhuman <- read.table("/media/macdatafile/mixed-hcv/gb-ref+hg38_v2/150724_M01841_0150_000000000-AE8E0.gb-ref+hg38_v2.csv",
#                              sep=",", header=TRUE)
# summary(jul24_hcvhuman)
# head(jul24_hcvhuman)
# dim(jul24_hcvhuman)
# 
# 
# jul24_hcv <- read.table("/media/macdatafile/mixed-hcv/gb-ref/150724_M01841_0150_000000000-AE8E0.gb-ref.csv",
#                         sep=",", header=TRUE)
# summary(jul24_hcv)
# head(jul24_hcv)
# dim(jul24_hcv)
# 
# cmp24 <- merge(x=jul24_hcvhuman,
#                y=jul24_hcv,
#                by.x=c("runname", "sample", "snum", "subtype"),
#                by.y=c("runname", "sample", "snum", "subtype"),
#                all=TRUE,
#                suffixes=c(".24HcvHum", ".24Hcv"))
# cmp24$perc.24HcvHum[is.na(cmp24$perc.24HcvHum)] <- 0
# cmp24$perc.24Hcv[is.na(cmp24$perc.24Hcv)] <- 0
# cmp24$perc_diff_24HcvHum_24Hcv <- cmp24$perc.24HcvHum - cmp24$perc.24Hcv
# 
# summary(cmp24)
# head(cmp24)
# dim(cmp24)
# 
# 
# summary(cmp24[abs(cmp24$perc_diff_24HcvHum_24Hcv) > MAX_PERC_DIFF & grepl("hg38", cmp24$subtype)==FALSE & cmp24$subtype != "" ,])
# dim(cmp24[abs(cmp24$perc_diff_24HcvHum_24Hcv) > MAX_PERC_DIFF & grepl("hg38", cmp24$subtype)==FALSE & cmp24$subtype != "" ,])
# 
# bigchange_jul24 <- cmp24[abs(cmp24$perc_diff_24HcvHum_24Hcv) > MAX_PERC_DIFF & grepl("hg38", cmp24$subtype)==FALSE & cmp24$subtype != "" ,]
# 
# #+ results="asis"
# kable(bigchange_jul24, format="html", caption="July 24 Sample Subtypes Affected By Reference", row.names=FALSE)
# 
# 
# #'
# 
# #' **July 24: Including Human Genome in Reference Makes a Small Difference for Sample Subtype Percentage 1a, 1b, 2a**
# #' 
# #' 
# #+ fig.width=15, fig.height=8
# # IQR in % coverage change due to change in reference by subtype
# fig <- ggplot(cmp[grepl("hg38", cmp24$subtype) == FALSE & cmp$subtype != "",], 
#               aes(x=subtype, y=perc_diff)) + 
#   geom_boxplot(fill="red") + 
#   ggtitle("July 24: Reference Makes Small Difference in Subtype %") + 
#   xlab("HCV Subtype") + 
#   ylab("Difference in % Reads Aligned When \nRef=HCV+Human Vs Ref=HCV") + 
#   theme_bw(base_size = 12) + 
#   theme(panel.grid.major = element_line(size = .3, color = "lightgray")) + 
#   theme(strip.text.x = element_text(size=rel(1)), 
#         axis.title=element_text(size=rel(2)), 
#         axis.text=element_text(size=rel(1.5), angle=45, hjust=1),
#         plot.title=element_text(size=rel(2)),
#         legend.title=element_text(size=rel(2)),
#         legend.text=element_text(size=rel(1.5)))
# print(fig)
# 
# 
# #' Compare July 24 Run - HCV+Human Reference Vs July 20 Run - HCV+Human Reference
# #' =====================================================================================
# #' 
# #' 
# cmp <- merge(x=subset(cmp20, select=-runname),
#              y=subset(cmp24, select=-runname),
#              by.x=c("sample", "snum", "subtype"),
#              by.y=c("sample", "snum", "subtype"),
#              all=TRUE)
# cmp$perc_diff_24HcvHuman_20HcvHuman
# 
# summary(cmp)
# head(cmp)
# dim(cmp)
# write.table(cmp, "/media/macdatafile/mixed-hcv/cmp_jul20_jul24_Hcv_HcvHuman.csv", sep=",", row.names=FALSE,  na="")
