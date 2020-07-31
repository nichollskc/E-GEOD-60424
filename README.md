# E-GEOD-60424

Pre-processing of RNA-seq dataset E-GEOD-60424 for Kath's biclustering comparison project.

# Steps

First collect sample info and process it in a form that we can easily extract fastq files for each sample.

`snakemake construct_fastq_dict`

Then run the rest of the analysis:

`snakemake data/tpm.tsv --cluster "sbatch --time=20 --account=CWALLACE-SL3-CPU" --jobs 1000`
