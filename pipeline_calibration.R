#!/usr/bin/Rscript --vanilla

## Compare the results of the mixed-hcv pipeline under different settings:
## - reference sequences: HCV only or HCV + human
## - regions considered: full genome or "deli-sliced"
## - mapping quality >= 0 and >= 10 (and maybe others)

## ##
## DATA ACQUISITION

## Directories that contain the data.  This is tailored to the
## directory structure we use to store it all.
ref.dirs <- list(HCVonly="/Volumes/MACDATAFILE/mixed-hcv/gb-ref",
                 HCVonly2="/Volumes/MACDATAFILE/mixed-hcv/gb-ref2",
                 HCVhuman="/Volumes/MACDATAFILE/mixed-hcv/gb-ref+hg38_v2")
mapq.cutoffs <- c(0, 10)
deli.subdir <- "deli"

## We will need to extract the run name from the filename for the deli-sliced
## data, as the runname in the files is wrong.
deli.runname.re <- "q(?:[[:digit:]]+)\\.(.+)\\.csv"

full.counts.table <- NULL
deli.counts.table <- NULL
for (ref.used in names(ref.dirs))
  {
    data.dir <- ref.dirs[[ref.used]]

    regions.to.consider <- "full"
    if (ref.used != "HCVonly")
      {
        regions.to.consider <- c(regions.to.consider, "deli")
      }
    
    for (regions.considered in regions.to.consider)
      {
        deli.prefix <- ""
        if (regions.considered == "deli")
          {
            data.dir <- paste(data.dir, "/deli", sep="")
          }
        
        for (mapq.cutoff in mapq.cutoffs)
          {
            ## Data using this mapping quality cutoff has "q[cutoff]." appended
            ## to the start.
            mapq.prefix <- paste("q", mapq.cutoff, ".", sep="")

            curr.run.csvs <- dir(data.dir,
                                 pattern=paste(deli.prefix, mapq.prefix,
                                   ".+\\.csv", sep=""),
                                 full.names=FALSE)

            for (curr.run.csv in curr.run.csvs)
              {
                cat("Reading ", curr.run.csv, " (ref=", ref.used,
                    ", regions=", regions.considered, ")\n", sep="")
                curr.run.table <- read.csv(paste(data.dir, "/", curr.run.csv, sep=""))

                curr.run.table$ref <- ref.used
                curr.run.table$mapq.cutoff <- mapq.cutoff                
                
                if (regions.considered == "deli")
                  {
                    curr.runname <- sub(deli.runname.re, "\\1", curr.run.csv)
                    curr.run.table$runname <- curr.runname
                    deli.counts.table <- rbind(deli.counts.table, curr.run.table)
                  }
                else
                  {
                    full.counts.table <- rbind(full.counts.table, curr.run.table)
                  }
              }
          }
      }
  }

## END DATA ACQUISITION
## ##


## full.counts.table looks like:
## runname,sample,snum,subtype,count,total,perc
full.counts.table <-
  full.counts.table[full.counts.table$subtype != "" &
                    grepl("hg38", full.counts.table$subtype) == FALSE,]

full.subtype.counts <- full.counts.table[,c("runname", "sample", "snum", "ref",
                                            "mapq.cutoff", "subtype", "count")]
full.subtype.counts <-
  full.subtype.counts[with(full.subtype.counts,
                           order(runname, sample, snum, ref, mapq.cutoff, subtype)),]
full.HCV.totals <-
  with(full.subtype.counts,
       aggregate(count,
                 by=list(runname=runname,
                   sample=sample,
                   snum=snum,
                   ref=ref,
                   mapq.cutoff=mapq.cutoff),
                 sum))
names(full.HCV.totals)[which(names(full.HCV.totals) == "x")] <- "total.HCV"
full.subtype.counts <- merge(full.subtype.counts, full.HCV.totals)

## The deli tables look like:
## runname,sample,snum,gene,start,end,rank,count,subtype,nucseq,aaseq
deli.subtype.counts <-
  with(deli.counts.table,
       aggregate(count,
                 by=list(runname=runname, sample=sample, snum=snum, gene=gene,
                   ref=ref, mapq.cutoff=mapq.cutoff, subtype=subtype),
                 sum))
## Rename the "x" column to "count".
names(deli.subtype.counts)[which(names(deli.subtype.counts) == "x")] <- "count"
deli.subtype.counts <-
  deli.subtype.counts[with(deli.subtype.counts,
                           order(runname, sample, snum, gene, ref,
                                 mapq.cutoff, subtype)),]

deli.HCV.totals <- with(deli.subtype.counts,
                        aggregate(count,
                                  by=list(runname=runname,
                                    sample=sample,
                                    snum=snum,
                                    gene=gene,
                                    ref=ref,
                                    mapq.cutoff=mapq.cutoff),
                                  sum))
names(deli.HCV.totals)[which(names(deli.HCV.totals) == "x")] <- "total.HCV"

deli.subtype.counts <- merge(deli.subtype.counts, deli.HCV.totals)

## Now stick them together.
full.subtype.counts$gene <- NA
full.subtype.counts$regions.considered <- "full"
deli.subtype.counts$regions.considered <- "deli"

subtype.counts <- rbind(full.subtype.counts, deli.subtype.counts)
subtype.counts <-
  subtype.counts[,c("runname", "sample", "snum", "ref", "mapq.cutoff",
                    "regions.considered", "gene", "subtype", "count", "total.HCV")]

write.csv(subtype.counts, file="HCV_subtype_data.csv", row.names=FALSE)

## ## We'd like to parse out the mixture types and proportions from the sample
## ## name, eg.
## ## HCVmixed9-RP1A1B-9010-HCV
## ## means 90% 1A and 10% 1B.
## mixture.re <- "HCVmixed.+-RP(?:SD)?([[:digit:]][[:alpha:]]?)([[:digit:]][[:alpha:]]?)-([[:digit:]]{2})([[:digit:]]{2})(?:.*)-HCV"

## counts.table$first.genotype <- sub(mixture.re, "\\1", counts.table$sample)
## counts.table$second.genotype <- sub(mixture.re, "\\2", counts.table$sample)
## counts.table$first.prop <- sub(mixture.re, "\\3", counts.table$sample)
## counts.table$second.prop <- sub(mixture.re, "\\4", counts.table$sample)
