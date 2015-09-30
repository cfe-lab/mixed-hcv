# 2 runs x 8 samples x 2 treatments (--no-mixed)

hcv <- read.csv('/Volumes/MACDATAFILE/mixed-hcv/aln2coverage.csv', header=T)

# omit the 150720 run
hcv.1 <- hcv[hcv$runname != '150720_M01841_0148_000000000-AE93J' & hcv$flags=='--quiet--local', ]
hcv.1$dilution <- grepl('5LOG', hcv.1$sample)
hcv.1$sample <- substr(hcv.1$sample, 1, 6)

write.csv(hcv.1, file='~/Desktop/all-coverage.csv', quote=FALSE)

setwd('~/git/mixed-hcv/results/')

for (sample in unique(hcv.1$sample)) {
	#pdf(paste(sample, ".pdf", sep=''), width=6, height=5, bg='white')
	
	temp1 <- hcv.1[hcv.1$sample==sample & hcv.1$dilution==FALSE, ]
	temp2 <- hcv.1[hcv.1$sample==sample & hcv.1$dilution==TRUE, ]
	
	par(mar=c(5,5,4,2))
	plot(temp1$rcoord, temp1$coverage, type='l', log='y', xlab='H77 coordinates', ylab=expression('Coverage'), cex.lab=1.2, col='dodgerblue', ylim=c(1,1E4), main=sample, cex.main=1.4, lwd=2)
	lines(temp2$rcoord, temp2$coverage, col='coral2', lwd=2)
	
	legend(x=3000, y=10, legend=c('undiluted', '5LOG'), col=c('dodgerblue', 'coral2'), lty=1, bty='n', lwd=2)
	#dev.off()
}
