library(sequenza)
source("/farm/home/rabino01/github/ngs_r/sequenza_onlyRatioSNP.R")

# Rscript sequenza.for.R <cnvkit.out.cns> <sample.name> <outputfile> <outputplot>

args <- commandArgs(trailingOnly=TRUE)
file.loc <- args[1]
sample.name <- args[2]
file.out <- args[3]
plot.out <- args[4]

chrom <- paste0('chr', c(1:19,'X'))

array.logR <- read.table(file.loc, stringsAsFactors = F,
                         header = TRUE , sep = "\t")

array.logR <- array.logR[array.logR$chromosome%in%chrom,]

array.logR$ratio <- 2^(array.logR$log2)

seg.size <- round((array.logR$end - array.logR$start)/1e6,0)
seg.filt <- seg.size >= 3

avg.depth.ratio = weighted.mean(array.logR$ratio,w = sqrt(1+seg.size))

CP <- ratio.model.fit(depth.ratio = array.logR$ratio[seg.filt],
                      weight.ratio = seg.size[seg.filt] + 150, 
                      priors.table = data.frame(CN = c(0, 1, 2, 3, 4),
                                                value = c(0.5, 1, 2, 1.5, 1)),
                      avg.depth.ratio = avg.depth.ratio,
                      mc.cores = 2)

cint <- get.ci(CP)

seg.cn <- ratio.bayes(depth.ratio = array.logR$ratio,
                      weight.ratio = 300,
                      priors.table = data.frame(CN = 2,value = 1),
                      avg.depth.ratio = avg.depth.ratio,
                      cellularity = cint$max.cellularity,
                      ploidy = cint$max.ploidy)
seg.all <- cbind(array.logR,seg.cn)


# add cellularity
seg.all <- cbind(seg.all,rep(cint$max.cellularity,dim(seg.all)[1]))

# add ploidy
seg.all <- cbind(seg.all,rep(cint$max.ploidy,dim(seg.all)[1]))

# add sample name as the first column
seg.all <- cbind(rep(sample.name,dim(seg.all)[1]),seg.all)

colnames(seg.all) <- c("sample.name", "chromosome", "start.pos", "end.pos", "gene", "log2", "probes.num", "weight", "ratio", "CNt", "L", "cellularity", "ploidy")

seg.all$A <- seg.all$CNt
seg.all$B <- 0


#genome.view(seg.all,info.type="CNt")
pdf(file=plot.out)
copynumber.view(seg.all,info.type="CNt")
cp.plot(CP)
cp.plot.contours(CP,add=T)
dev.off()
#genome.view.logR(seg.all, ylim = c(-2.0,2.0))
#abline(h=log2(avg.depth.ratio), lty = 2)

write.table(seg.all,file=file.out,sep="\t",quote=FALSE,row.names=FALSE)

#ploidy.i <- w.median(seg.all$CNt, w= sqrt(1+seg.size), c= seg.all$chromosome)
#ploidy <- round(ploidy.i, 0)
#seg.all$sample <- nat2012.files[i,2]

