library(dmrseq)
setwd("/rsrch4/scratch/genomic_med/METI/Suresh/RRBS/bismark_methcall")
bismarkBSseq <- read.bismark(files = c('KO-1000nM-5-AZA_dedup.sorted.bismark.cov.gz','KO-1000nM-DAC_dedup.sorted.bismark.cov.gz','KO-150nM-NTX-301_dedup.sorted.bismark.cov.gz','KO-50nM-NTX-301_dedup.sorted.bismark.cov.gz','KO-con_dedup.sorted.bismark.cov.gz','R175H-1000nM-DAC_dedup.sorted.bismark.cov.gz','R175H-120nM-NTX-301_dedup.sorted.bismark.cov.gz','R175H-2000nM-5-AZA_dedup.sorted.bismark.cov.gz','R175H-40nM-NTX-301_dedup.sorted.bismark.cov.gz','R175H-con_dedup.sorted.bismark.cov.gz','R248Q-1000nM-DAC_dedup.sorted.bismark.cov.gz','R248Q-100nM-NTX-301_dedup.sorted.bismark.cov.gz','R248Q-2000nM-5-AZA_dedup.sorted.bismark.cov.gz','R248Q-300nM-NTX-301_dedup.sorted.bismark.cov.gz','R248Q-con_dedup.sorted.bismark.cov.gz','WT-10nM-NTX-301_dedup.sorted.bismark.cov.gz','WT-300nM-5-AZA_dedup.sorted.bismark.cov.gz','WT-300nM-DAC_dedup.sorted.bismark.cov.gz','WT-30nM-NTX-301_dedup.sorted.bismark.cov.gz','WT-con_dedup.sorted.bismark.cov.gz'),colData = DataFrame(row.names = c('KO-1000nM-5-AZA','KO-1000nM-DAC','KO-150nM-NTX-301','KO-50nM-NTX-301','KO-con','R175H-1000nM-DAC','R175H-120nM-NTX-301','R175H-2000nM-5-AZA','R175H-40nM-NTX-301','R175H-con','R248Q-1000nM-DAC','R248Q-100nM-NTX-301','R248Q-2000nM-5-AZA','R248Q-300nM-NTX-301','R248Q-con','WT-10nM-NTX-301','WT-300nM-5-AZA','WT-300nM-DAC','WT-30nM-NTX-301','WT-con')), rmZeroCov = FALSE, strandCollapse = FALSE)

bismarkBSseq

pData(bismarkBSseq)$CellType <- c("KO","KO","KO","KO","KO","R175H","R175H","R175H","R175H","R175H","R248Q","R248Q","R248Q","R248Q","R248Q","WT","WT","WT","WT","WT")

pData(bismarkBSseq)

pData(bismarkBSseq)$Replicate <- c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5)

loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bismarkBSseq, type="Cov")==0) == 0)
sample.idx <- which(pData(bismarkBSseq)$CellType %in% c("KO","R175H","R248Q","WT"))
bismarkBSseq.filtered <- bismarkBSseq[loci.idx, sample.idx]
testCovariate <- "CellType"
regions <- dmrseq(bismarkBSseq.filtered, cutoff = 0.05, testCovariate=testCovariate)
blocks <- dmrseq(bismarkBSseq.filtered, cutoff = 0.05, testCovariate=testCovariate, block = TRUE, minInSpan = 500, bpSpan = 5e4, maxGapSmooth = 1e6, maxGap = 5e3)
head(blocks)

write.csv(regions, file="regions.markers.xls")
write.csv(blocks, file="blocks.markers.xls")

sum(regions$qval < 0.05)

#regions_qval <- regions$qval < 0.05

sigRegions <- regions[regions$qval < 0.05,]

sum(sigRegions$stat > 0) / length(sigRegions)

annoTrack <- getAnnot("hg38")

plotDMRs(bs, regions=regions[1,], testCovariate="CellType", annoTrack=annoTrack)
plotDMRs(bs, regions=regions[2,], testCovariate="CellType", annoTrack=annoTrack)
plotDMRs(bs, regions=regions[3,], testCovariate="CellType", annoTrack=annoTrack)
plotDMRs(bs, regions=regions[4,], testCovariate="CellType", annoTrack=annoTrack)
plotDMRs(bs, regions=regions[5,], testCovariate="CellType", annoTrack=annoTrack)

plotDMRs(bs, regions=blocks[1,], testCovariate="CellType", annoTrack=annoTrack)
plotDMRs(bs, regions=blocks[2,], testCovariate="CellType", annoTrack=annoTrack)
plotDMRs(bs, regions=blocks[3,], testCovariate="CellType", annoTrack=annoTrack)
plotDMRs(bs, regions=blocks[4,], testCovariate="CellType", annoTrack=annoTrack)

plotDMRs(bs, regions=blocks[5,], testCovariate="CellType", annoTrack=annoTrack)


plotEmpiricalDistribution(bs, testCovariate="CellType")


plotEmpiricalDistribution(bs, testCovariate="CellType", type="Cov", bySample=TRUE)

write.csv(as.data.frame(regions), file="regions_results.csv")

write.csv(as.data.frame(blocks), file="blocks_results.csv")

rawDiff <- meanDiff(bs, dmrs=sigRegions, testCovariate="CellType")
str(rawDiff)

write.csv(as.data.frame(rawDiff), file="rawDiff_results.csv")


