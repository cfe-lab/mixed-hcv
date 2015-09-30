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

coverage.dir <- "coverage"
run.names <- list.dirs(coverage.dir, full.names=FALSE, recursive=FALSE)

sample.file.re <- "(.+)_coverage\\.csv"
for (run.name in run.names)
  {
    run.dir <- paste(coverage.dir, run.name, sep="/")
    ref.names <- list.dirs(run.dir, full.names=FALSE, recursive=FALSE)

    for (ref.name in ref.names)
      {
        ref.dir <- paste(run.dir, ref.name, sep="/")
        sample.files <- list.files(ref.dir, pattern=sample.file.re,
                                   full.names=FALSE)

        for (sample.file in sample.files)
          {
            sample.basename <- sub(sample.file.re, "\\1", sample.file)
            curr.amino.file.path <- paste(ref.dir, sample.file, sep="/")
            cat("Plotting coverage maps for file ",
                curr.amino.file.path, "...", sep="")
            
            curr.coverage.table <-
              read.csv(curr.amino.file.path,
                       stringsAsFactors=FALSE)
            
            pdf(paste(ref.dir, "/", sample.basename, "_coverage.pdf", sep=""),
                width=8.5, height=11, bg="white")
            coverage.plot(curr.coverage.table, sample.basename)
            dev.off()

            cat(" done.\n")
          }
      }
  }
