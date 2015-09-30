files <- Sys.glob('/Volumes/MACDATAFILE/mixed-hcv/coverage/150724_M01841_0150_000000000-AE8E0/gb-ref+hg38_v2/*coverage.csv')

alphabet <- c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X.')

all <- {}
for (f in files) {
	items <- strsplit(f, split='/')[[1]]
	filename <- items[length(items)]
	sample <- strsplit(filename, split='_')[[1]][1]
	
	temp <- read.csv(f, header=TRUE)
	ns5b <- temp[temp$gene == 'NS5b', ]
	ns5b$total <- apply(ns5b[,names(ns5b) %in% alphabet], 1, sum)
	
	temp <- data.frame(sample=sample, gene=ns5b$gene, aa.pos=ns5b$ref.pos, coverage=ns5b$total)
	all <- rbind(all, temp)
}

write.csv(all, '~/Desktop/Coverage-all-ns5b.csv', row.names=FALSE, quote=FALSE)