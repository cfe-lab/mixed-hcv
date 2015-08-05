library(plyr)

jul20_hcvhuman <- read.table("/media/macdatafile/mixed-hcv/gb-ref+hg38_v2/150720_M01841_0148_000000000-AE93J.csv",
                        sep=",", header=TRUE)
summary(jul20_hcvhuman)
head(jul20_hcvhuman)
dim(jul20_hcvhuman)

sum_hcv_only <- ddply(.data=jul20_hcvhuman,
                      .variables=c("runname", "sample", "snum"),
                      .fun=function(x) {
                        data.frame(total_hcv=sum(x$count[grepl("hg38")]))
                      })

jul20_hcv <- read.table("/media/macdatafile/mixed-hcv/gb-ref/150720_M01841_0148_000000000-AE93J.csv",
                        sep=",", header=TRUE)
summary(jul20_hcv)
head(jul20_hcv)
dim(jul20_hcv)


cmp <- merge(x=jul20_hcvhuman,
             y=jul20_hcv,
             by.x=c("runname", "sample", "snum", "subtype"),
             by.y=c("runname", "sample", "snum", "subtype"),
             all=TRUE,
             suffixes=c(".20HcvHum", ".20Hcv"))
cmp$perc.20HcvHum[is.na(cmp$perc.20HcvHum)] <- 0
cmp$perc.20Hcv[is.na(cmp$perc.20Hcv)] <- 0
summary(cmp)
head(cmp)
dim(cmp)


cmp$perc_diff <- cmp$perc.20HcvHum - cmp$perc.20Hcv

summary(cmp)
summary(cmp[abs(cmp$perc_diff) > 1 & grepl("hg38", cmp$subtype)==FALSE & cmp$subtype != "" ,])
dim(cmp[abs(cmp$perc_diff) > 1 & grepl("hg38", cmp$subtype)==FALSE & cmp$subtype != "" ,])

cmp[abs(cmp$perc_diff) > 1 & grepl("hg38", cmp$subtype)==FALSE & cmp$subtype != "" ,]

