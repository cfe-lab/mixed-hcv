# Makes venn diagrams of read hits between different pipelines
library(VennDiagram)
library(reshape2)

library(gridExtra)

opts_chunk$set(progress = TRUE, verbose = TRUE, width=1500, tidy = FALSE, error= TRUE, warning = FALSE, message=FALSE, echo=FALSE)


RUNS  <- c("150720_M01841_0148_000000000-AE93J")
CACHEDIR <- "/media/macdatafile/mixed-hcv/cache"


IS_FIRST <- 0x040
IS_UNMAPPED <- 0x004
IS_SECONDARY_ALIGNMENT <-   0x100
IS_CHIMERIC_ALIGNMENT <-    0x800

read_runrefsam <- function (run, ref, sample) {
  
  # Find all the samples for the same run
  sam_ref_dir <- paste0(CACHEDIR, "/", run , "/sam/", ref)
  sam_basenames <- list.files(sam_ref_dir, pattern=paste0(sample, ".*--quiet--local.sam") , full.names=FALSE)
  
  if (length(sam_basenames) > 1) {
    stop(paste0("There are multiple full genome sams for the same sample in ", sam_ref_dir))
  }
                              
  # Now compare the 4 pipelines for the same sample
  samfile <- paste0(sam_ref_dir, "/", sam_basenames[1] )
    
  #qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
  sam <- read.table(samfile, sep="\t", header=FALSE, fill=TRUE,
                    # skip the columns that we don't care about with NULL col class
                    colClasses = c("character", "integer", "character", "NULL", "integer", 
                                   rep("NULL", 16))
  )
  
  
  colnames(sam) <- c("qname", "flag", "rname", "mapq")
  
  hits <- subset(sam,
                 bitwAnd(IS_FIRST, as.integer(sam$flag)) > 0 )
  
  hits$subtype <- as.factor(sapply(hits$rname, function(x) {
    if (grepl("hg38", x) == TRUE) {
      return ("Human")
    } else if (x == "*") {
      return ("Unaligned")
    } else if (grepl("HCV-1a", x) == TRUE) {
      return ("1a")
    } else if (grepl("HCV-1b", x) == TRUE) {
      return ("1b")
    } else {
      return (substr(x, 5, 5))
    }
  }))
  
  hits$category <- as.factor(sapply(hits$rname, function(x) {
    if (grepl("hg38", x) ==  TRUE) {
      return ("Human")
    } else if (x == "*") {
      return ("Unaligned")
    } else  {
      return ("HCV")
    }
  }))
  
  return (hits)
}


venn_category <- function(hcv_human_q0_hits, hcv_human_q10_hits,
                            hcv_q0_hits, hcv_q10_hits,                            
                            category, title) {
  filt_hcv_human_q0_hits <- hcv_human_q0_hits$qname[hcv_human_q0_hits$category == category]
  #filt_hcv_human_q10_hits <- hcv_human_q10_hits$qname[hcv_human_q10_hits$category == category]
  filt_hcv_q0_hits <- hcv_q0_hits$qname[hcv_q0_hits$category == category]
  filt_hcv_q10_hits <- hcv_q10_hits$qname[hcv_q10_hits$category == category]

#   venn.plot <- venn.diagram(list(filt_hcv_human_q0_hits, filt_hcv_human_q10_hits, filt_hcv_q0_hits, filt_hcv_q10_hits),
#                filename=NULL,
#                main=title,
#                category=c("HCV+Human Mapq>=0",
#                            "HCV+Human Mapq>=10",
#                            "HCV Mapq>=0",
#                            "HCV Mapq>=10"),
#                 fill = c("red", "green", "blue", "yellow"))
  venn.plot <- venn.diagram(list(filt_hcv_human_q0_hits, filt_hcv_q0_hits, filt_hcv_q10_hits),
               filename=NULL,
               main=title,
               category=c("HCV+Human Mapq>=0",
                           "HCV Mapq>=0",
                           "HCV Mapq>=10"),
                fill = c("red", "blue", "yellow"))
  return (venn.plot)

}

#+ fig.width=10
# for (run in RUNS) {
  
  # TODO:  remove me
  run <- "150720_M01841_0148_000000000-AE93J"
  
  # Find all the samples for the same run
  sam_ref_dir <- paste0(CACHEDIR, "/", run , "/sam/gb-ref+hg38_v2")
  sam_basenames <- list.files(sam_ref_dir, pattern="*.--quiet--local.sam" , full.names=FALSE)
  samples <- sapply(sam_basenames, function(x) {unlist(strsplit(x, "_L001_R"))[1]})
  
  # TODO:  remove me
samples <- samples[1]
  # Now compare the 4 pipelines for the same sample

  for (sample in samples) {
    
    #hcv_human_q0_hits <- read_runrefsam(run, ref="gb-ref+hg38_v2", sample)
    #hcv_q0_hits <- read_runrefsam(run, ref="gb-ref", sample)
    
    
    #hcv_human_q10_hits <- hcv_human_q0_hits
    #hcv_human_q10_hits$subtype[hcv_human_q10_hits$mapq < 10] <- "Unaligned"
    #hcv_human_q10_hits$category[hcv_human_q10_hits$mapq < 10] <- "Unaligned"
    
    #hcv_q10_hits <- hcv_q0_hits
    #hcv_q10_hits$subtype[hcv_q10_hits$mapq < 10] <- "Unaligned"
    #hcv_q10_hits$category[hcv_q10_hits$mapq < 10] <- "Unaligned"
    
    # TODO:  we only want to compare categories for Full Genome + Deli
#     for (category in levels(hcv_human_q0_hits$category)) {
#       venn.plot <- venn_category (hcv_human_q0_hits, hcv_human_q10_hits,
#                                   hcv_q0_hits, hcv_q10_hits,
#                                   category=category, 
#                                   title=paste0("Run=", run, " Sample=", sample, " Category=", category))
#       plot.new()
#       grid.draw(venn.plot)
#     }
    
    
    print("Are there reads that are ambiguously human or HCV, depending on Reference?")
    
    filt_hcv_human_q0_hits <- hcv_human_q0_hits[hcv_human_q0_hits$category == "Human" &
                                                  bitwAnd(hcv_human_q0_hits$flag, IS_SECONDARY_ALIGNMENT) == 0 & 
                                                  bitwAnd(hcv_human_q0_hits$flag, IS_CHIMERIC_ALIGNMENT) == 0 ,
                                                c("qname", "mapq")]
    filt_hcv_q0_hits <- hcv_q0_hits[hcv_q0_hits$category == "HCV" &
                                      bitwAnd(hcv_q0_hits$flag, IS_SECONDARY_ALIGNMENT) == 0 & 
                                      bitwAnd(hcv_q0_hits$flag, IS_CHIMERIC_ALIGNMENT) == 0 ,
                                    c("qname", "mapq")]
    venn.plot <- venn.diagram(list(filt_hcv_human_q0_hits$qname, filt_hcv_q0_hits$qname),
                              filename=NULL,
                              height=7, width=12, units="in",
                              main=paste0("Run=", run, " Sample=", sample),
                              category=c("Ref=HCV+Human: Human Hits",
                                         "Ref=HCV: HCV Hits"),
                              cat.just=list(c(-1, 1), c(1, 1)),
                              fill = c("red", "blue"))
    
    plot.new()
    grid.draw(venn.plot)
    
    
    print("What is the Map Quality Distro of Ambiguous and Unambiguous Hits?")
    
    ambig <- merge(x=filt_hcv_human_q0_hits,
                   y=filt_hcv_q0_hits,
                   by="qname", 
                   suffixes=c("RefHCVHuman.Q0.HumanHits", "RefHCV.Q0.HCVHits"),
                   all=FALSE)
    dim(ambig)
    head(ambig)
    summary(ambig)
    
    
    ambig_melt <- reshape2:::melt(data=ambig, id.vars=c("qname"),
                                  variable.name="dfsource", value.name="mapq")


    
                      
    unambig <- rbind(data.frame(dfsource="RefHCVHuman", hcv_human_q0_hits[!(hcv_human_q0_hits$qname %in% ambig$qname) &
                                                                            bitwAnd(hcv_human_q0_hits$flag, IS_SECONDARY_ALIGNMENT) == 0 & 
                                                                            bitwAnd(hcv_human_q0_hits$flag, IS_CHIMERIC_ALIGNMENT) == 0 &
                                                                            hcv_human_q0_hits$category != "Unaligned",]),
                     data.frame(dfsource="RefHCV", hcv_q0_hits[!(hcv_q0_hits$qname %in% ambig$qname) &
                                                                 bitwAnd(hcv_q0_hits$flag, IS_SECONDARY_ALIGNMENT) == 0 & 
                                                                 bitwAnd(hcv_q0_hits$flag, IS_CHIMERIC_ALIGNMENT) == 0 &
                                                                 hcv_q0_hits$category != "Unaligned",]))

    full <- rbind(unambig[, c("qname", "dfsource", "mapq")], ambig_melt)
    
    fig <-  ggplot(full,  aes(x=mapq, color=dfsource)) + 
      geom_density(size=2) + 
      xlab("Mapping Quality") + 
      ylab("Density") + 
      scale_color_discrete(name="", labels=c("Ref=HCV+Human: Unambiguous Hits",
                                             "Ref=HCV: Unambiguous Hits",
                                             "Ref=HCV+Human: Ambiguous Hits",
                                             "Ref=HCV: Ambiguous Hits")) + 
      ggtitle(paste0("Run=", run, "\nSample=", sample, "\nMap Quality Distribution")) + 
      theme_bw(base_size = 12) + 
      theme(panel.grid.major = element_line(size = .3, color = "lightgray")) + 
      theme(strip.text.x = element_text(size=rel(1)), 
            axis.title=element_text(size=rel(2)), 
            axis.text=element_text(size=rel(1.5)),
            legend.title=element_text(size=rel(2)),
            legend.text=element_text(size=rel(1.5)))
    print(fig)
    
  }
    
  
# }

