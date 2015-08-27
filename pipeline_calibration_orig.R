#!/usr/bin/Rscript --vanilla

## Compare the results of the mixed-hcv pipeline under different settings:
## - reference sequences: HCV only or HCV + human
## - regions considered: full genome or "deli-sliced"
## - mapping quality >= 0 and >= 10 (and maybe others)

## ##
## DATA ACQUISITION

already.saved <- TRUE
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
                                    mapq.cutoff=mapq.cutoff,
                                    coverage=coverage),
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
mixture.re <- "HCVmixed.+-RP(?:SD)?([[:digit:]][[:alpha:]]?)([[:digit:]][[:alpha:]]?)-([[:digit:]]{2})([[:digit:]]{1,2})(?:.*)"

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

## For each sample, we make a spider plot of:
## - the expected distribution
## - deli, HCVhuman, mapq.cutoff == 0, coverage == 0.5
## - deli, HCVhuman, mapq.cutoff == 0, coverage == 1
## - full, HCVhuman, mapq.cutoff == 0 and 10
## - full, HCVonly, mapq.cutoff == 0 and 10
## - full, HCVonly2, mapq.cutoff == 0 and 10

cases.to.plot <-
  data.frame(regions.considered=c(rep("deli", 2), rep("full", 6)),
             ref=c(rep("HCVhuman", 2), rep("HCVonly", 2), rep("HCVonly2", 2),
               rep("HCVhuman", 2)),
             mapq.cutoff=c(0, 0, 0, 10, 0, 10, 0, 10),
             coverage=c(0.5, rep(1, 7)))

radii.scale.factor <- 0.8

pdf("mixture_starplots.pdf")
for (i in 1:nrow(expected.mixtures))
## for (i in 3)
  {
    curr.sample.name <- expected.mixtures$sample[i]
    curr.sample.data <- 
      mixed.samples.props[grepl(curr.sample.name, mixed.samples.props$sample),]

    relevant.cases <- merge(cases.to.plot, curr.sample.data, sort=FALSE)

    expected.star.data <- expected.mixtures[i, 2:ncol(expected.mixtures)]
    star.data <- relevant.cases[, genotypes]
    combined.star.data <- rbind(expected.star.data, star.data)

    labels <- sapply(1:nrow(relevant.cases),
                     function (idx)
                     {
                       paste(relevant.cases$regions.considered[idx],
                             ", ",
                             relevant.cases$ref[idx],
                             "\nmapq=", relevant.cases$mapq.cutoff[idx],
                             ", coverage=", relevant.cases$coverage[idx],
                             sep="")
                     })
    
    star.locations <- stars(rbind(expected.star.data, star.data),
                            scale=FALSE, len=radii.scale.factor,
                            lwd=2, labels=c("expected", labels),
                            flip.labels=TRUE,
                            main=curr.sample.name,
                            radius=FALSE,
                            cex=0.5)
    
    for (i in 1:nrow(star.locations))
      {
        stars(combined.star.data[i,], scale=FALSE,
              len=radii.scale.factor,
              locations=c(star.locations$Var1[i], star.locations$Var2[i]),
              key.loc=c(star.locations$Var1[i], star.locations$Var2[i]),
              lty=2, add=TRUE)
      }
  }
dev.off()


## Let's get some assessments for how certain our assessments of
## the multinomial parameters are.
library(nnet)

## To solve for p_1, add up p_X/p_1 for all genotypes X to get (1-p1)/p1, then
## solve for p1:
## (1-p1)/p1 == X
## (1-p1) == p1*X
## p1(X + 1) == 1
## p1 == 1/(X+1)
probs <- function(coefs)
  {
    odds.ratios <- exp(coefs)
    not.p1.over.p1 <- sum(odds.ratios)
    p1 <- 1/(not.p1.over.p1 + 1)
    return(c(p1, odds.ratios*p1))
  }

## A function that calculates crude confidence intervals for a
## multinomial regression's fitted values.
mult.regr.bounds <- function(multinomial.regression)
  {
    std.err <- summary(multinomial.regression)$standard.errors
    estimates <- fitted(multinomial.regression)

    ## Each element of all.bounds is a list with two elements: the first
    ## is a two-vector of bounds coming from the intercept.
    all.bounds <-
      lapply(1:length(coef(multinomial.regression)),
             function (idx)
             {
               bound.coefs.1 <- coef(multinomial.regression)
               bound.coefs.2 <- bound.coefs.1
               bound.coefs.1[idx] <- bound.coefs.1[idx] + 2*std.err[idx]
               bound.coefs.2[idx] <- bound.coefs.2[idx] - 2*std.err[idx]
               
               bound.probs.1 <- probs(bound.coefs.1)
               bound.probs.2 <- probs(bound.coefs.2)
               return(list(bounds=range(bound.probs.1[idx+1],
                             bound.probs.2[idx+1]),
                           int.bounds=range(bound.probs.1[1],
                             bound.probs.2[1])))
             })
    intercept.bounds <-
      range(unlist(lapply(all.bounds, function (curr.bound) curr.bound[[2]])))
    param.bounds <- lapply(all.bounds, function(curr.bound) curr.bound[[1]])
    
    return(list(lower=c(intercept.bounds[1],
                  sapply(param.bounds, function (x) min(x))),
                upper=c(intercept.bounds[2],
                  sapply(param.bounds, function (x) max(x)))))
  }

analyze.trial <- function(curr.trial)
  {
    mult.regression <- multinom(as.matrix(curr.trial[, 6:14]) ~ 1)
    
    ## Translate the coefficients to probabilities.  The coefficients specify
    ## the log-odds of being in a specific genotype vs genotype 1 (non-AB).
    ## For example: log(odds of being in genotype 1A vs 1) == coef1A
    ## so p_{1A}/p_1 = exp(coef1A)
    
    ## We need some criteria for how certain the result would have to be
    ## before we consider it.
    std.err <- summary(mult.regression)$standard.errors
    estimates <- fitted(mult.regression)
    bounds <- mult.regr.bounds(mult.regression)

    details <-
      with(curr.trial,
           paste(regions.considered, ", ", ref, "\nmapq cutoff=",
                 mapq.cutoff, ", coverage=", coverage, ", N=", total, sep=""))
    
    x.coords <- barplot(estimates, ylim=c(0,1), col="lightgrey",
                        xlab="genotype", ylab="proportion",
                        main=details, cex.main=0.75)

    arrows(x0=x.coords, y0=bounds$lower,
           y1=bounds$upper, code=3,
           angle=90, length=0.05, col="red")
  }

pdf("mixture_confidences.pdf", width=8.5, height=11)
for (i in 1:nrow(expected.mixtures))
  {
    curr.sample.name <- expected.mixtures$sample[i]
    curr.sample.data <-
      mixed.samples[grepl(curr.sample.name, mixed.samples$sample),]

    par(mfrow=c(4,2), oma=c(0,0,3,0))
    for (idx in 1:nrow(cases.to.plot))
      {
        relevant.cases <- merge(cases.to.plot[idx,], curr.sample.data, sort=FALSE)
        if (nrow(relevant.cases) == 0)
          {
            plot.new()
          } else {
            analyze.trial(relevant.cases)
          }
      }
    mtext(curr.sample.name, outer=TRUE, cex=1.5)
  }
dev.off()

## ##
## Produce some plots of the estimates of each proportion, broken down
## by the mixture.
ingredient.re <- "HCVmixed[^-]+-([^-]+)-.+"
ingredient.strings <- sub(ingredient.re, "\\1", unique(mixed.samples$sample))

mixed.samples$first.geno <- sub(mixture.re, "\\1", mixed.samples$sample)
mixed.samples$second.geno <- sub(mixture.re, "\\2", mixed.samples$sample)
mixed.samples$first.prop <-
  as.numeric(sub(mixture.re, "\\3", mixed.samples$sample))/100
mixed.samples$second.prop <-
  as.numeric(sub(mixture.re, "\\4", mixed.samples$sample))/100
                  
## For each trial, perform a multinomial regression.
regressions <-
  lapply(1:nrow(mixed.samples),
         function (idx)
         {
           first.geno <- mixed.samples$first.geno[idx]
           second.geno <- mixed.samples$second.geno[idx]

           first.geno.count <- mixed.samples[idx, first.geno]
           second.geno.count <- mixed.samples[idx, second.geno]
           rest.count <- mixed.samples$total[idx] - first.geno.count -
             second.geno.count

           regr.matrix <-
             matrix(c(rest.count, first.geno.count, second.geno.count),
                    dimnames=list(c(""), c("rest", first.geno, second.geno)),
                    nrow=1)
           
           regr <- multinom(regr.matrix ~ 1)

           return(regr)
         })

## For each set of "ingredients", let's make plots of the estimated
## proportions of the first genotype, second genotype, and others.

## For each sample, we make a bar plot of:
## - deli, HCVhuman, mapq.cutoff == 0, coverage == 0.5
## - deli, HCVhuman, mapq.cutoff == 0, coverage == 1
## - full, HCVhuman, mapq.cutoff == 0 and 10
## - full, HCVonly2, mapq.cutoff == 0 and 10

cases.to.barplot <-
  data.frame(regions.considered=c(rep("deli", 2), rep("full", 4)),
             ref=c(rep("HCVhuman", 2), rep("HCVonly2", 2),
               rep("HCVhuman", 2)),
             mapq.cutoff=c(0, 0, 0, 10, 0, 10),
             coverage=c(0.5, rep(1, 5)))

pipeline.setting.descs <-
  sapply(1:nrow(cases.to.barplot),
         function (idx)
         {
           paste(cases.to.barplot$regions.considered[idx],
                 ", ", cases.to.barplot$ref[idx],
                 ", mapq >=", cases.to.barplot$mapq.cutoff[idx],
                 ", coverage=", cases.to.barplot$coverage[idx],
                 sep="")
         })

barplot.palette <- rainbow(nrow(cases.to.barplot), s=0.5, v=0.8)
ten.ninety.two.re <- "[^-]+-RP[^-]+-1090-2(?:-HCV)?"

pdf("mixture_barplots.pdf", width=11, height=8.5)
for (ingredient.string in unique(ingredient.strings))
  { 
    ## Note: in our dataset, this will never return a 0-length vector.
    relevant.indices <- which(grepl(ingredient.string, mixed.samples$sample))
    curr.samples <- mixed.samples[relevant.indices,]
    curr.samples$orig.idx <- relevant.indices

    genotypes <- c(first=curr.samples$first.geno[1],
                   second=curr.samples$second.geno[1],
                   rest="rest")

    expected <- list(first=c(0.1, 0.1, 0.5, 0.9))
    expected$second <- 1 - expected$first
    if (grepl("RPSD", ingredient.string))
      {
        expected$first <- c(0.75, 0.90, 0.95, 0.98)
        expected$second <- 1 - expected$first
      }
    expected$rest <- expected$first

    ## Let's build some bar charts of our data, with one set of bars
    ## per case we care about.  Thus we will need to build matrices
    ## whose columns are different expected proportions and whose rows
    ## represent different pipelines.
    
    estimates <- list(rest=NULL, first=NULL, second=NULL)
    lower <- list(rest=NULL, first=NULL, second=NULL)
    upper <- list(rest=NULL, first=NULL, second=NULL)
    for (idx in 1:nrow(cases.to.barplot))
      {
        ## This loops over rows of the matrix we're building:
        ## each loop deals with one pipeline case.
        
        curr.data <- merge(cases.to.barplot[idx,], curr.samples)
        ## We'll also order the 10/90 splits by replicate number.
        curr.data$replicate <- 1
        curr.data$replicate[grepl(ten.ninety.two.re, curr.data$sample)] <- 2
        curr.data <- curr.data[order(curr.data$first.prop,
                                     curr.data$replicate),]

        curr.estimates <- list(rest=rep(NA, 4), first=rep(NA, 4),
                               second=rep(NA, 4))
        curr.lower <- list(rest=rep(NA, 4), first=rep(NA, 4),
                           second=rep(NA, 4))
        curr.upper <- list(rest=rep(NA, 4), first=rep(NA, 4),
                           second=rep(NA, 4))
        for (curr.settings.idx in 1:nrow(curr.data))
          {
            curr.row <- curr.data[curr.settings.idx,]
            ## Now we fill in each entry in the row.  For
            ## a non-serial dilution case, the columns should be, from
            ## left to right, 10%, 10%, 50%, and 90%.
            curr.first.prop <- curr.row$first.prop
            
            if (grepl("RPSD", curr.row$sample))
              {
                entry.idx <- match(curr.row$first.prop, c(0.75, 0.9, 0.95, 0.98))
              } else {
                entry.idx <- NA
                if (curr.row$first.prop == 0.1)
                  {
                    if (curr.row$replicate == 1)
                      {
                        entry.idx <- 1
                      } else {
                        entry.idx <- 2
                      }
                  } else if (curr.row$first.prop == 0.5) {
                    entry.idx <- 3
                  } else {
                    entry.idx <- 4
                  }
              }
            
            curr.regr <- regressions[[curr.row$orig.idx]]
            regr.estimates <- fitted(curr.regr)
            regr.bounds <- mult.regr.bounds(curr.regr)

            for (j in 1:3)
              {
                curr.estimates[[j]][entry.idx] <- regr.estimates[j]
                curr.lower[[j]][entry.idx] <- regr.bounds$lower[j]
                curr.upper[[j]][entry.idx] <- regr.bounds$upper[j]
              }
          }

        ## Now add the newly constructed row to the matrices we're building.
        for (param in c("rest", "first", "second"))
          {
            estimates[[param]] <- rbind(estimates[[param]],
                                        curr.estimates[[param]])
            lower[[param]] <- rbind(lower[[param]], curr.lower[[param]])
            upper[[param]] <- rbind(upper[[param]], curr.upper[[param]])
          }
      }

    par(mfrow=c(2,2), oma=c(0,0,3,0))
    for (param in c("first", "second", "rest"))
      {

        
        curr.positions <-
          barplot(estimates[[param]], ylim=c(0,1),
                  beside=TRUE, names.arg=expected[[param]],
                  col=barplot.palette,
                  main=genotypes[param])

        ## Add some error bars.
        for (i in 1:nrow(curr.positions))
          {
            arrows(x0=curr.positions[i,],
                   y0=lower[[param]][i,],
                   y1=upper[[param]][i,],
                   code=3, angle=90, length=0.05, col="black")
          }

        if (param %in% c("first", "second"))
          {
            ## Add horizontal dashed lines for the expected amounts.
            for (i in 1:ncol(curr.positions))
              {
                arrows(x0=curr.positions[1, i],
                       x1=curr.positions[nrow(cases.to.barplot), i],
                       y0=expected[[param]][i],
                       lty="dashed", length=0)
              }
          }

      }
    plot.new()
    legend("center", legend=pipeline.setting.descs,
           fill=barplot.palette, xpd=TRUE)
    mtext(ingredient.string, outer=TRUE, cex=1.5)
  }
dev.off()


## ##
## Confidence plots based on multinomial regression, as Richard H
## requested.

confidence.plot.RH.cases <- data.frame(regions.considered="full",
  ref="HCVhuman", mapq.cutoff=10, coverage=1)

pdf("mixture_confidences_RH.pdf", width=11, height=8.5)
for (i in 1:nrow(expected.mixtures))
  {
    curr.sample.name <- expected.mixtures$sample[i]
    curr.sample.data <-
      mixed.samples[grepl(curr.sample.name, mixed.samples$sample),]

    par(oma=c(0,0,3,0))
    for (idx in 1:nrow(confidence.plot.RH.cases))
      {
        relevant.cases <- merge(confidence.plot.RH.cases[idx,],
                                curr.sample.data, sort=FALSE)
        if (nrow(relevant.cases) == 0)
          {
            plot.new()
          } else {
            analyze.trial(relevant.cases)
          }
      }
    mtext(curr.sample.name, outer=TRUE, cex=1.5)
  }
dev.off()
