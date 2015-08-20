# Makes venn diagrams of read hits between different pipelines
library(VennDiagram)

HITS_CSV <- "/media/macdatafile/mixed-hcv/staging/hits.csv"

# run,ref,sample,read,subtype,mapq
hits <- read.csv(HITS_CSV)

# pipelines we want to compare
# HCVHuman mapq>0  | HCV  mapq > 0  | HCV mapq > 10  | HCVHuman mapq > 10


hits$category <- as.factor(sapply(hits$subtype, 
                                  function(x) {
                                    if (x == "Human") {
                                      return ("Human")
                                    } else if (x == "Unaligned") {
                                      return ("Unaligned")
                                    } else {
                                      return ("HCV")
                                    }
                                  }))

pipeline_ref = "gb-ref+hg38_v2"
pipeline_mapq = 10

run = "150720_M01841_0148_000000000-AE93J"
subtype = "1a"
run_type_hits <- subset(hits, hits$run == run & hits$subtype == subtype)

# for each run,  subtype, and category, draw a venn diagram containing each pipeline

venn_runsubtype <- function(filterhits, title) {

  hcv_human_q0_hits <- filterhits$read[filterhits$ref == "gb-ref+hg38_v2"]
  hcv_human_q10_hits <- filterhits$read[filterhits$ref == "gb-ref+hg38_v2" & filterhits$mapq >= 10]
  hcv_q0_hits <- filterhits$read[filterhits$ref == "gb-ref"]
  hcv_q10_hits <- filterhits$read[filterhits$ref == "gb-ref+hg38_v2" & filterhits$mapq >= 10]
  
  i12=intersect(hcv_human_q0_hits, hcv_human_q10_hits)
  i13=intersect(hcv_human_q0_hits, hcv_q0_hits)
  i14=intersect(hcv_human_q0_hits, hcv_q10_hits)
  i23=intersect(hcv_human_q10_hits, hcv_q0_hits)
  i24=intersect(hcv_human_q10_hits, hcv_q10_hits)
  i34=intersect(hcv_q0_hits, hcv_q10_hits)
  
  venn.plot <- draw.quad.venn(area1=length(hcv_human_q0_hits),
                              area2=length(hcv_human_q10_hits),
                              area3=length(hcv_q0_hits),
                              area4=length(hcv_q10_hits),
                              n12=length(i12),
                              n13=length(i13),
                              n14=length(i14),
                              n23=length(i23),
                              n24=length(i24),
                              n34=length(i34),
                              n123=length(intersect(i12, hcv_q0_hits)),
                              n134=length(intersect(i13, hcv_q10_hits)),
                              n234=length(intersect(i23, hcv_q10_hits)),
                              n1234=length(intersect(i12, i34)),
                              category=c("HCV+Human Mapq>=0",
                                         "HCV+Human Mapq>=10",
                                         "HCV Mapq>=0",
                                         "HCV Mapq>=10"),
                              main=title)
  grid.draw(venn.plot)
  dev.off()
}

for run in levels(hits$run) {
  for category in levels(hit$category) {
    filterhits <- hits[hits$run == run & hits$category == category]
    venn_runsubtype(filterhits=filterhits, title=paste0(run + " " + category))
  }
}


for run in levels(hits$run) {
  for subtype in levels(hit$subtype) {
    filterhits <- hits[hits$run == run & hits$subtype == subtype]
    venn_runsubtype(filterhits=filterhits, title=paste0(run + " " + category))
  }
}
