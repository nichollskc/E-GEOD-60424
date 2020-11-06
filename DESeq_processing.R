library("DESeq2")
library("tximeta")
library("rhdf5")

sample_info_file <- snakemake@input[["sample_info"]]
tx2gene_file <- snakemake@input[["tx2gene"]]

# Read in samples csv - use convention of calling this column data (coldata)
coldata <- read.csv(sample_info_file, sep='\t')

# Read in counts from kallisto and generate summarised experiment (se)
# Condense to genes using the transcripts_to_genes file we generated before
coldata[["names"]] <- coldata[["Comment..ENA_SAMPLE."]]
coldata[["files"]] <- paste0("data/kallisto/", coldata[["names"]], "/abundance.h5")
print(names(coldata))
tx2gene <- read.table(tx2gene_file, sep="\t", header=FALSE)
gse <- tximeta(coldata, type="kallisto", tx2gene=tx2gene, txOut=FALSE)

print("Finished txi import")

gse[["condition"]] <- factor(gse[["FactorValue..disease."]])
gse[["celltype"]] <- factor(gse[["FactorValue..cell.type."]])

print("Setting up DESeq dataset")

dds <- DESeqDataSet(gse, design = ~ condition + celltype)

print(nrow(dds))
# We will manually discard genes later, for now just discard the absolute necessary
keep <- rowSums(counts(dds)) > 1
# keep <- rowSums(counts(dds) >= 10) >= 30
dds <- dds[keep,]
print(nrow(dds))

# Apply variance-stabilizing transformation and size-factor adjustment
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = FALSE)
sf_normalised <- counts(dds, normalized=TRUE)

write.table(t(assay(vsd)), snakemake@output[["vst_normalised"]], sep="\t", col.names=TRUE, row.names=coldata[["names"]])
write.table(t(sf_normalised), snakemake@output[["sf_normalised"]], sep="\t", col.names=TRUE, row.names=coldata[["names"]])

write.table(coldata, file = snakemake@output[["vst_sample_info"]], row.names = FALSE, col.names = TRUE, sep="\t")
write.table(coldata, file = snakemake@output[["sf_sample_info"]], row.names = FALSE, col.names = TRUE, sep="\t")
