#!/usr/bin/Rscript --vanilla

## Compare the results of the mixed-hcv pipeline under different settings:
## - reference sequences: HCV only or HCV + human
## - regions considered: full genome or "deli-sliced"
## - mapping quality >= 0 and >= 10 (and maybe others)

## ##
## DATA ACQUISITION

already.saved <- FALSE
save.filename <- "counts_data.RData"

## Directories that contain the data.  This is tailored to the
## directory structure we use to store it all.
ref.dirs <- list(HCVonly="/Volumes/MACDATAFILE/mixed-hcv/gb-ref",
                 HCVonly2="/Volumes/MACDATAFILE/mixed-hcv/gb-ref2",
                 HCVhuman="/Volumes/MACDATAFILE/mixed-hcv/gb-ref+hg38_v2")
mapq.cutoffs <- c(0, 10)
coverages <- c(1, 0.5)
deli.subdir <- "deli"

## We will need to extract the run name from the filename for the deli-sliced
## data, as the runname in the files is wrong.
deli.runname.re <- "q(?:[[:digit:]]+)\\.(?:.+minwidth.+\\.)?(15.+)\\.csv"

if (!already.saved)
  {
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
            if (regions.considered == "deli")
              {
                data.dir <- paste(data.dir, "/deli", sep="")
              }
            
            for (mapq.cutoff in mapq.cutoffs)
              {
                ## Data using this mapping quality cutoff has "q[cutoff]." appended
                ## to the start.
                mapq.prefix <- paste("q", mapq.cutoff, "\\.", sep="")

                for (coverage in coverages)
                  {
                    coverage.prefix <- ""
                    if (coverage == 0.5)
                      {
                        coverage.prefix <- "gb-ref\\+hg38_v2\\.minwidth0\\.5\\."
                      }
                    
                    curr.run.csvs <- dir(data.dir,
                                         pattern=paste(mapq.prefix,
                                           coverage.prefix, "15.+\\.csv", sep=""),
                                         full.names=FALSE)
                    
                    for (curr.run.csv in curr.run.csvs)
                      {
                        ## ## FIXME revisit this tomorrow!
                        ## if (grepl("150729", curr.run.csv))
                        ##   {
                        ##     next
                        ##   }
                        
                        cat("Reading ", curr.run.csv, " (ref=", ref.used,
                            ", regions=", regions.considered, ")\n", sep="")
                        curr.run.table <- read.csv(paste(data.dir, "/",
                                                         curr.run.csv, sep=""))
                        
                        curr.run.table$ref <- ref.used
                        curr.run.table$mapq.cutoff <- mapq.cutoff
                        curr.run.table$coverage <- coverage
                        
                        if (regions.considered == "deli")
                          {
                            curr.runname <- sub(deli.runname.re, "\\1", curr.run.csv)
                            curr.run.table$runname <- curr.runname
                            deli.counts.table <-
                              rbind(deli.counts.table, curr.run.table)
                          }
                        else
                          {
                            full.counts.table <-
                              rbind(full.counts.table, curr.run.table)
                          }
                      }
                  }
              }
          }
      }
    save(full.counts.table, deli.counts.table, file=save.filename)
  } else {
    load(save.filename)
  }
    
## END DATA ACQUISITION
## ##

## ##
## CREATING CSV OUTPUT

## full.counts.table looks like:
## runname,sample,snum,subtype,count,total,perc
full.counts.table <-
  full.counts.table[full.counts.table$subtype != "" &
                    grepl("hg38", full.counts.table$subtype) == FALSE,]

full.subtype.counts <- full.counts.table[,c("runname", "sample", "snum", "ref",
                                            "mapq.cutoff", "coverage", "subtype",
                                            "count")]
full.subtype.counts <-
  full.subtype.counts[with(full.subtype.counts,
                           order(runname, sample, snum, ref, mapq.cutoff, coverage,
                                 subtype)),]
full.HCV.totals <-
  with(full.subtype.counts,
       aggregate(count,
                 by=list(runname=runname,
                   sample=sample,
                   snum=snum,
                   ref=ref,
                   mapq.cutoff=mapq.cutoff,
                   coverage=coverage),
                 sum))
names(full.HCV.totals)[which(names(full.HCV.totals) == "x")] <- "total.HCV"
full.subtype.counts <- merge(full.subtype.counts, full.HCV.totals)

## The deli tables look like:
## runname,sample,snum,gene,start,end,rank,count,subtype,nucseq,aaseq
deli.subtype.counts <-
  with(deli.counts.table,
       aggregate(count,
                 by=list(runname=runname, sample=sample, snum=snum, gene=gene,
                   ref=ref, mapq.cutoff=mapq.cutoff, coverage=coverage,
                   subtype=subtype),
                 sum))
## Rename the "x" column to "count".

names(deli.subtype.counts)[which(names(deli.subtype.counts) == "x")] <- "count"
deli.subtype.counts <-
  deli.subtype.counts[with(deli.subtype.counts,
                           order(runname, sample, snum, gene, ref,
                                 mapq.cutoff, coverage, subtype)),]

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
  subtype.counts[,c("runname", "sample", "snum", "regions.considered",
                    "ref", "mapq.cutoff", "coverage",
                    "gene", "subtype", "count", "total.HCV")]

write.csv(subtype.counts, file="HCV_subtype_data.csv", row.names=FALSE)

## END CREATING CSV OUTPUT
## ##

## We wish to compare the total number of reads to the number of HCV reads.
total.reads.by.sample <- unique(full.counts.table[,c(1,2,3,6)])
total.reads.by.sample <- merge(total.reads.by.sample, full.HCV.totals)
total.reads.by.sample$proportion.HCV <- with(total.reads.by.sample, total.HCV/total)

## We parse out the mixture types and proportions from the sample name, eg.
## HCVmixed9-RP1A1B-9010-HCV
## means 90% 1A and 10% 1B.
mixture.re <- "HCVmixed.+-RP(?:SD)?([[:digit:]][[:alpha:]]?)([[:digit:]][[:alpha:]]?)-([[:digit:]]{2})([[:digit:]]{1,2})(?:.*)-HCV"

mixtures <- subtype.counts[grepl(mixture.re, subtype.counts$sample),]

mixtures$genotype <-
  sapply(mixtures$subtype,
         function (subtype)
         {
           result <- ""
           if (subtype == "1a")
             {
               result <- "1A"
             } else if (subtype == "1b") {
               result <- "1B"
             } else {
               result <- substring(subtype, 1, 1)
             }
           return(result)
         })

mixtures.by.genotype.full <-
  with(mixtures[mixtures$regions.considered == "full",],
       aggregate(count,
                 by=list(runname=runname, sample=sample, snum=snum,
                   ref=ref, mapq.cutoff=mapq.cutoff,
                   regions.considered=regions.considered,
                   coverage=coverage,
                   genotype=genotype),
                 sum))
mixtures.by.genotype.full$gene <- NA
names(mixtures.by.genotype.full)[which(names(mixtures.by.genotype.full) == "x")] <-
  "count"

mixtures.by.genotype.deli <-
  with(mixtures[mixtures$regions.considered == "deli",],
       aggregate(count,
                 by=list(runname=runname, sample=sample, snum=snum,
                   ref=ref, mapq.cutoff=mapq.cutoff,
                   regions.considered=regions.considered,
                   coverage=coverage,
                   gene=gene, genotype=genotype),
                 sum))
names(mixtures.by.genotype.deli)[which(names(mixtures.by.genotype.deli) == "x")] <-
  "count"

mixtures.by.genotype <- rbind(mixtures.by.genotype.full, mixtures.by.genotype.deli)
mixtures.by.genotype <-
  mixtures.by.genotype[,c("runname", "sample", "snum", "regions.considered", "ref",
                          "mapq.cutoff", "coverage",
                          "gene", "genotype", "count")]
mixtures.by.genotype <-
  mixtures.by.genotype[with(mixtures.by.genotype,
                            order(runname, sample, snum, regions.considered,
                                  ref, mapq.cutoff,
                                  coverage, gene, genotype)),]

mixed.samples.raw <- NULL
unique.settings <- unique(mixtures.by.genotype[,1:8])
genotypes <- sort(unique(mixtures.by.genotype$genotype))
for (i in 1:nrow(unique.settings))
  {
    curr.runname <- unique.settings$runname[i]
    curr.sample <- unique.settings$sample[i]
    curr.snum <- unique.settings$snum[i]
    curr.ref <- unique.settings$ref[i]
    curr.mapq.cutoff <- unique.settings$mapq.cutoff[i]
    curr.regions.considered <- unique.settings$regions.considered[i]
    curr.coverage <- unique.settings$coverage[i]
    curr.gene <- unique.settings$gene[i]

    curr.data <- mixtures[mixtures$runname == curr.runname &
                          mixtures$sample == curr.sample &
                          mixtures$snum == curr.snum &
                          mixtures$ref == curr.ref &
                          mixtures$mapq.cutoff == curr.mapq.cutoff &
                          mixtures$regions.considered == curr.regions.considered &
                          mixtures$coverage == curr.coverage &
                          ((mixtures$gene == curr.gene) |
                           (is.na(mixtures$gene) & is.na(curr.gene))),]

    curr.setting.data <- unique.settings[i,]
    curr.total <- 0
    for (genotype in genotypes)
      {
        curr.count <- sum(curr.data$count[genotype == curr.data$genotype])
        curr.total <- curr.total + curr.count
        curr.setting.data[1, genotype] <- curr.count
      }
    curr.setting.data$total <- curr.total
    mixed.samples.raw <- rbind(mixed.samples.raw, curr.setting.data)
  }

## For legibility let's drop runname and snum.
mixed.samples.raw$runname <- NULL
mixed.samples.raw$snum <- NULL

## Collapse the "deli" reads together.
collapse.deli.genes <-
  unique(mixed.samples.raw[mixed.samples.raw$regions.considered == "deli",
                           c("sample", "ref", "mapq.cutoff", "coverage")])

mixed.samples.deli <- NULL
for (i in 1:nrow(collapse.deli.genes))
  {
    new.row <- collapse.deli.genes[i,]

    curr.indices <- with(mixed.samples.raw,
                         sample == new.row$sample &
                         ref == new.row$ref &
                         mapq.cutoff == new.row$mapq.cutoff &
                         coverage == new.row$coverage &
                         regions.considered == "deli")
    curr.data <-
      mixed.samples.raw[curr.indices &
                        mixed.samples.raw$regions.considered == "deli",]

    for (genotype in genotypes)
      {
        new.row[, genotype] <- sum(curr.data[, genotype])
      }
    new.row$total <- sum(curr.data$total)
    mixed.samples.deli <- rbind(mixed.samples.deli, new.row)
  }
mixed.samples.deli$regions.considered <- "deli"

## This table drops the "gene" column.
mixed.samples <-
  rbind(mixed.samples.raw[mixed.samples.raw$regions.considered == "full",
                          c(1, 2, 3, 4, 5, 7:ncol(mixed.samples.raw))],
        mixed.samples.deli)
mixed.samples <-
  mixed.samples[with(mixed.samples,
                     order(sample, regions.considered, ref, mapq.cutoff,
                           coverage)),]
                      
mixed.samples.props <- mixed.samples[,1:5]
for (genotype in genotypes)
  {
    mixed.samples.props[, genotype] <-
      mixed.samples[, genotype]/mixed.samples$total
  }

## Now put together an analogous table with the *expected* values so we can
## make starplots of the expected against the actual.

expected.mixtures <- data.frame(sample=unique(mixed.samples$sample))
for (genotype in genotypes)
  {
    expected.mixtures[,genotype] <- 0
  }
for (i in 1:nrow(expected.mixtures))
  {
    curr.sample <- expected.mixtures$sample[i]
    curr.first.geno <- sub(mixture.re, "\\1", curr.sample)
    curr.second.geno <- sub(mixture.re, "\\2", curr.sample)
    curr.first.prop <- as.numeric(sub(mixture.re, "\\3", curr.sample))
    curr.second.prop <- as.numeric(sub(mixture.re, "\\4", curr.sample))
    
    expected.mixtures[i, curr.first.geno] <- curr.first.prop/100
    expected.mixtures[i, curr.second.geno] <- curr.second.prop/100
  }

## ## For each sample, we make a spider plot of:
## ## - the expected distribution
## ## - deli, HCVhuman, mapq.cutoff == 0, coverage == 0.5
## ## - deli, HCVhuman, mapq.cutoff == 0, coverage == 1
## ## - full, HCVhuman, mapq.cutoff == 0 and 10
## ## - full, HCVonly2, mapq.cutoff == 0 and 10

## cases.to.plot <-
##   data.frame(regions.considered=c(rep("deli", 2), rep("full", 4)),
##              ref=c(rep("HCVhuman", 4), rep("HCVonly2", 2)),
##              mapq.cutoff=c(0, 0, 0, 10, 0, 10),
##              coverage=c(0.5, rep(1, 5)))
## spiderplot.colours <- heat.colors(nrow(cases.to.plot)+1, alpha=0.25)
                            
## ## for (i in 1:nrow(expected.mixtures))
## for (i in 1)
##   {
##     curr.sample.name <- sample.names[i]
##     curr.sample.data <- 
##       mixed.samples.props[grepl(curr.sample.name, mixed.samples.props$sample),]

##     relevant.cases <- merge(cases.to.plot, curr.sample.data)

##     star.data <- rbind(expected.mixtures[i, 2:ncol(expected.mixtures)],
##                        relevant.cases[, genotypes])
##     ## stars(star.data, scale=FALSE, location=c(0,0),
##     ##       key.loc=c(0,0), main=sample.name, lty=2,
##     ##       col.lines=spiderplot.colours)
    
##     stars(star.data[1,], scale=FALSE, location=c(0,0),
##           key.loc=c(0,0), main=sample.name, lty=2,
##           col.lines=rgb(0.5, 0.5, 0.5, 0.5))
##     stars(star.data[1,], scale=FALSE, location=c(0,0),
##           lwd=2, add=TRUE)
##     stars(star.data[2:nrow(star.data),],
##           scale=FALSE, location=c(0,0), lty=1, add=TRUE,
##           col.lines=heat.colors(nrow(cases.to.plot), alpha=0.25))
##   }
