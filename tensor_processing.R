library(dplyr)

# Read in counts, leaving sample name column as part of the df
counts <- read.csv(snakemake@input[['counts']], sep='\t', row.names=NULL)
sample_info <- read.csv(snakemake@input[['sample_info']], sep='\t')
print(colnames(counts)[1:10])
sample_name_col <- colnames(counts)[1]

# Sort the samples by cell type and then individual
# (The samples are already sorted like this I think, but best to be safe)
# Discard NK cells as this cell type isn't complete
sorted_samples <- sample_info %>%
    filter(FactorValue..cell.type. != 'NK') %>%
    arrange(FactorValue..cell.type., Characteristics..individual.)

# Pick out the samples in the required order
sorted_counts <- sorted_samples %>%
    select(Comment..ENA_SAMPLE.) %>%
    left_join(counts,
              by=c('Comment..ENA_SAMPLE.' = sample_name_col))

# Write out sorted sample info
write.table(sorted_samples, file = snakemake@output[["sample_info"]], row.names = FALSE, col.names = TRUE, sep="\t")

# Write out sorted counts matrix, without gene names or sample names
write.table(sorted_counts[, -1], file = snakemake@output[["Y"]], row.names = FALSE, col.names = FALSE, sep="\t")
writeLines(colnames(sorted_counts[, -1]), con = snakemake@output[["gene_names"]])

# Write out number of individuals (N), needed by tensor methods
N <- length(unique(sorted_samples$Characteristics..individual.))
writeLines(as.character(), con = snakemake@output[["N"]])
