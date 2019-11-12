packages <- c("Rsubread", "limma", "edgeR")
lapply(packages, library, character.only =TRUE)

out_dir <- "/gpfs/fs1/data/taylorlab/HbAS/readCounts/working/primary_allowSingle_featureCounts"
bam_dir <- "/gpfs/fs1/data/taylorlab/HbAS/readCounts/bamFiles/all"
gtf_dir <- "/gpfs/fs1/data/taylorlab/Genomes/3D7_2019/GTF/"
gtffile <- file.path(gtf_dir,"Plasmodium_falciparum.EPr1.43.gtf")
csvfile <- file.path(out_dir, "sampleTable.csv")
sampleTable <- read.csv(csvfile, row.names = 1)
filenames <- file.path(bam_dir, paste0(sampleTable$sampleName, ".sort.bam"))
file.exists(filenames)

fc <- featureCounts(filenames, annot.ext=gtffile, isGTFAnnotationFile=TRUE, nthreads=4, isPairedEnd=T, requireBothEndsMapped=F, primaryOnly=T)

colnames(fc$counts) <- sampleTable$sampleName 

x <- DGEList(counts=fc$counts, genes = fc$annotation[,c("GeneID", "Length")])

x <- calcNormFactors(x)

calc_tpm <- function(x, gene.length) {
  x <- as.matrix(x)
  len.norm.lib.size <- colSums(x / gene.length)
  return((t(t(x) / len.norm.lib.size) * 1e06) / gene.length)
}

tpm <- calc_tpm(x, gene.length=x$genes$Length)

write.csv(fc$counts, file = "counts.txt", sep = "\t")
write.csv(tpm, file = "tpm.txt", sep = "\t")

save(file="primary_allowSingle_featureCounts.RData")

