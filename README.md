# E-GEOD-60424

Pre-processing of RNA-seq dataset E-GEOD-60424 for Kath's biclustering comparison project.

# Steps

Install the conda environment and activate it:

`conda env create -f environment.yml && conda activate rnaseq-kallisto`

First collect sample info and process it in a form that we can easily extract fastq files for each sample.

`snakemake construct_fastq_dict`

Then run the rest of the analysis. This step will cause an error if the `fastq_dict.json` file has not been created by the first snakemake command.

`snakemake all_datasets --cluster "sbatch --time=30 --account=CWALLACE-SL3-CPU" --jobs 1000`

If using in biclustering comparison, copy the files to the relevant directory:

`BIC_COMP_DIR=~/rds/rds-cew54-wallace-ngs/kath_biclust_comp/biclustering_comparison/data/real/sortedblood/; RSYNC_OPTS="-avrn"  mkdir -p $BIC_COMP_DIR; rsync $RSYNC_OPTS data/raw $BIC_COMP_DIR; rsync $RSYNC_OPTS data/tensor $BIC_COMP_DIR; rsync $RSYNC_OPTS data/deseq $BIC_COMP_DIR`

## Installing conda environment - alternative instructions

I actually used mamba for the installation instead, but it should have the same result as the commands above. Supposedly mamba has a better method for searching for dependencies and finding the most recent version of a package.

Outside of the conda environment, install mamba.

`conda install -c conda-forge mamba`

Then use mamba to create the environment

`mamba env create -f environment.yml`
