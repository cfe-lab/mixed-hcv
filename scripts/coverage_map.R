#! /usr/bin/env Rscript --vanilla

## A Kiveified version of the original, encased-in-amber, coverage_map.R,
## which worked with a hardcoded file path.

args <- commandArgs(TRUE)

if (length(args) != 3)
  {
    cat("Error: must specify exactly three arguments",
        "(input CSV filename, title CSV filename, output PDF filename)\n",
        file=stderr())
    quit(save="no", status=1)
  }

## This has columns:
## gene,ref.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
input.csv <- args[1]

## This file simply specifies the title for the plot, and has a single column:
## title
title.csv <- args[2]

## This specifies where to write the output plot to.
output.pdf <- args[3]

## For reference:
gene.seqs <-
  c(NS3="APITAYAQQTRGLLGCIITSLTGRDKNQVEGEVQIVSTATQTFLATCINGVCWTVYHGAGTRTIASPKGPVIQMYTNVDQDLVGWPAPQGSRSLTPCTCGSSDLYLVTRHADVIPVRRRGDSRGSLLSPRPISYLKGSSGGPLLCPAGHAVGLFRAAVCTRGVAKAVDFIPVENLETTMRSPVFTDNSSPPAVPQSFQVAHLHAPTGSGKSTKVPAAYAAQGYKVLVLNPSVAATLGFGAYMSKAHGVDPNIRTGVRTITTGSPITYSTYGKFLADGGCSGGAYDIIICDECHSTDATSILGIGTVLDQAETAGARLVVLATATPPGSVTVSHPNIEEVALSTTGEIPFYGKAIPLEVIKGGRHLIFCHSKKKCDELAAKLVALGINAVAYYRGLDVSVIPTSGDVVVVSTDALMTGFTGDFDSVIDCNTCVTQTVDFSLDPTFTIETTTLPQDAVSRTQRRGRTGRGKPGIYRFVAPGERPSGMFDSSVLCECYDAGCAWYELTPAETTVRLRAYMNTPGLPVCQDHLEFWEGVFTGLTHIDAHFLSQTKQSGENFPYLVAYQATVCARAQAPPPSWDQMWKCLIRLKPTLHGPTPLLYRLGAVQNEVTLTHPITKYIMTCMSADLEVVT",
    NS5a="SGSWLRDIWDWICEVLSDFKTWLKAKLMPQLPGIPFVSCQRGYRGVWRGDGIMHTRCHCGAEITGHVKNGTMRIVGPRTCRNMWSGTFPINAYTTGPCTPLPAPNYKFALWRVSAEEYVEIRRVGDFHYVSGMTTDNLKCPCQIPSPEFFTELDGVRLHRFAPPCKPLLREEVSFRVGLHEYPVGSQLPCEPEPDVAVLTSMLTDPSHITAEAAGRRLARGSPPSMASSSASQLSAPSLKATCTANHDSPDAELIEANLLWRQEMGGNITRVESENKVVILDSFDPLVAEEDEREVSVPAEILRKSRRFARALPVWARPDYNPPLVETWKKPDYEPPVVHGCPLPPPRSPPVPPPRKKRTVVLTESTLSTALAELATKSFGSSSTSGITGDNTTTSSEPAPSGCPPDSDVESYSSMPPLEGEPGDPDLSDGSWSTVSSGADTEDVVCC",
    NS5b="SMSYSWTGALVTPCAAEEQKLPINALSNSLLRHHNLVYSTTSRSACQRQKKVTFDRLQVLDSHYQDVLKEVKAAASKVKANLLSVEEACSLTPPHSAKSKFGYGAKDVRCHARKAVAHINSVWKDLLEDSVTPIDTTIMAKNEVFCVQPEKGGRKPARLIVFPDLGVRVCEKMALYDVVSKLPLAVMGSSYGFQYSPGQRVEFLVQAWKSKKTPMGFSYDTRCFDSTVTESDIRTEEAIYQCCDLDPQARVAIKSLTERLYVGGPLTNSRGENCGYRRCRASGVLTTSCGNTLTCYIKARAACRAAGLQDCTMLVCGDDLVVICESAGVQEDAASLRAFTEAMTRYSAPPGDPPQPEYDLELITSCSSNVSVAHDGAGKRVYYLTRDPTTPLARAAWETARHTPVNSWLGNIIMFAPTLWARMILMTHFFSVLIARDQLEQALNCEIYGACYSIEPLDLPPIIQRLHGLSAFSLHSYSPGEINRVAACLRKLGVPPLRAWRHRARSVRARLLSRGGRAAICGKYLFNWAVRTKLKLTPIAAAGRLDLSGWFTAGYSGGDIYHSVSHARPRWFWFCLLLLAAGVGIYLLPNR")

coverage.plot <- function(coverage.table, overall.title)
  {
    coverage.table$total.coverage <-
      sapply(1:nrow(coverage.table),
             function (idx) sum(coverage.table[idx,3:23]))

    par(mfrow=c(3,1), mar=c(5,5,4,2) + 0.1,
        oma=c(0,0,3,0))
    for (gene in names(gene.seqs))
      {
        curr.gene.coverage <- coverage.table[coverage.table$gene == gene,]
        
        ## We'll want to plot total.coverage against ref.pos, but we first
        ## need to buttress our data with 0 coverage for any missing values
        ## of ref.pos.
        curr.ref.pos <- curr.gene.coverage$ref.pos
        curr.total.coverage <- curr.gene.coverage$total.coverage
        
        missing.positions <- setdiff(min(curr.ref.pos):max(curr.ref.pos),
                                     curr.ref.pos)
        missing.coverages <- rep(0, length(missing.positions))
        
        curr.ref.pos <- c(curr.ref.pos, missing.positions)
        curr.total.coverage <- c(curr.total.coverage, missing.coverages)
        curr.total.coverage <- curr.total.coverage[order(curr.ref.pos)]
        
        plot(min(curr.ref.pos):max(curr.ref.pos),
             curr.total.coverage, type="s",
             xlim=c(1, nchar(gene.seqs[gene])),
             ylim=c(1, 10000), log="y",
             main=gene, xlab="Position (H77 gene coordinates)",
             ylab="Count",
             cex.axis=1.5,
             cex.lab=1.5,
             cex.main=1.5)
        mtext(overall.title, outer=TRUE, cex=1.5)
      }
  }

cat("Plotting coverage maps for file ", input.csv, "...", sep="")

curr.coverage.table <- read.csv(input.csv, stringsAsFactors=FALSE)
title.table <- read.csv(title.csv, stringsAsFactors=FALSE)

pdf(output.pdf, width=8.5, height=11, bg="white")
coverage.plot(curr.coverage.table, title.table$title[1])
dev.off()

cat(" done.\n")
