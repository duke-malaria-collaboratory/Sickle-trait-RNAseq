## Assay to detect differential transcript usage (DTU)

# Directory of directories containing "quant.sf" from salmon quant, gtf file, and sample table
data_dir <- "path_to_dir/"
cvsfile <- file.path(data_dir, "sampleTable.csv"
sampleTable <- read.csv(cvsfile, row.names = 1)
sampleTable$HPI <- factor(sampleTable$HPI)
files <- file.path(data_dir, "quant", sampleTable$sampleName, "quant.sf")
file.exists(files)

library(tximport)
txi <- tximport(files, type="salmon", txOut=TRUE, countsFromAbundance="scaledTPM")
cts <- txi$counts
cts <- cts[rowSums(cts) > 0,]
colnames(cts) <- rownames(sampleTable)

library(GenomicFeatures)
gtf <- file.path(data_dir, "gencode.v32.annotation.gtf.gz")
txdb <- makeTxDbFromGFF(gtf)

txdf <- select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")
tab <- table(txdf$GENEID)
txdf$ntx <- tab[match(txdf$GENEID, names(tab))]

txdf <- txdf[match(rownames(cts),txdf$TXNAME),]
all(rownames(cts) == txdf$TXNAME)
all(rownames(cts) %in% txdf$TXNAME)

counts <- data.frame(gene_id=txdf$GENEID, feature_id=txdf$TXNAME, cts)

library(DRIMSeq)
d <- dmDSdata(counts=counts, samples=sampleTable)
d <- dmFilter(d,
  min_samps_feature_expr=n.small, min_feature_expr=1,
  min_samps_feature_prop=n.small, min_feature_prop=0.1,
  min_samps_gene_expr=n, min_gene_expr=1)

d

design_full <- model.matrix(~Genotype, data=DRIMSeq::samples(d))

set.seed(1)
system.time({
  d <- dmPrecision(d, design=design_full)
  d <- dmFit(d, design=design_full)
  d <- dmTest(d, coef="GenotypeAS")
})

res.txp <- DRIMSeq::results(d, level="feature")
res.txp