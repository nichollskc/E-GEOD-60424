configfile: "config.yml"

localrules: generate_count_matrix

import json

def all_abundance_files(wildcards):
    return all_kallisto_output(wildcards, 'tsv')

def all_h5_files(wildcards):
    return all_kallisto_output(wildcards, 'h5')

def all_kallisto_output(wildcards, extension):
    fastq_dict = get_fastq_dict()
    return [f"data/kallisto/{SAMPLE_ID}/abundance.{extension}"
            for SAMPLE_ID in fastq_dict.keys()]

def get_fastq_files_for_ID(wildcards):
    fastq_dict = get_fastq_dict()
    sample_ID = wildcards.SAMPLE_ID
    files = [f"data/fastq/{sample_ID}/{full_id}.fastq.gz"
            for full_id in fastq_dict[sample_ID].keys()]
    return files

def get_fastq_url_from_fastq_id(wildcards):
    fastq_dict = get_fastq_dict()
    return fastq_dict[wildcards.SAMPLE_ID][wildcards.FASTQ_ID]

def get_fastq_dict():
    # Ideally this wouldn't be hard-coded, but this way is easier..
    fastq_dict_file = "fastq_dict.json"
    with open(fastq_dict_file, 'r') as f:
        fastq_dict = json.load(f)
    return fastq_dict

rule construct_fastq_dict:
    input:
        "data/raw/sample_info.txt"
    output:
        "fastq_dict.json"
    run:
        import json
        import pandas as pd

        sample_info = pd.read_csv(input[0], sep="\t")
        sample_col = 'Comment [ENA_SAMPLE]'
        url_col = 'Comment [FASTQ_URI]'
        id_col = 'FASTQ_ID'
        sample_info[id_col] = sample_info[url_col].str.extract(r'ftp://.*/(\w+)\.fastq\.gz')

        # Generate a dictionary with each key one of the samples,
        # and the value a dictionary from fastq_id to fastq_url for each
        # fastq file associated with the sample
        grouped = sample_info.groupby(sample_col)[id_col, url_col]
        fastq_dict = grouped.apply(lambda df : dict(zip(df[id_col], df[url_col]))).to_dict()

        with open(output[0], 'w') as f:
            json.dump(fastq_dict, f, indent=2)

rule fetch_sample_info:
    output:
        "data/raw/sample_info.txt"
    shell:
        "wget -O {output} --no-verbose {config[sample_info_url]}"

rule fetch_index:
    output:
        temp("transcriptome.idx.tar.gz")
    shell:
        "wget -O {output} --no-verbose {config[index_tar_url]}"

rule extract_index:
    input:
        "transcriptome.idx.tar.gz"
    output:
        "homo_sapiens/transcriptome.idx"
        "homo_sapiens/transcripts_to_genes.txt",
    shell:
        "tar -xvzf {input}"

rule fetch_fastq:
    input:
        "fastq_dict.json"
    group:
        "kallisto_sample"
    params:
        url=get_fastq_url_from_fastq_id
    log:
        "logs/data/fastq/{SAMPLE_ID}/{FASTQ_ID}.log"
    output:
        temp("data/fastq/{SAMPLE_ID}/{FASTQ_ID}.fastq.gz")
    shell:
        "wget -O {output} --no-verbose {params.url} 2>&1 | tee {log}"

rule run_kallisto:
    group:
        "kallisto_sample"
    input:
        "fastq_dict.json",
        index="homo_sapiens/transcriptome.idx",
        fastq=get_fastq_files_for_ID
    log:
        "logs/data/kallisto/{SAMPLE_ID}.log"
    output:
        "data/kallisto/{SAMPLE_ID}/abundance.tsv"
    shell:
        "kallisto quant --index={input.index} --output-dir=$(dirname {output}) {input.fastq} 2>&1 | tee {log}"

rule generate_count_matrix:
    input:
        "fastq_dict.json",
        tx_to_gene="homo_sapiens/transcripts_to_genes.txt",
        abundance=all_abundance_files
    output:
        df="data/raw/tpm.tsv"
    run:
        import os
        import pandas as pd
        tpms = []
        names = []
        print(input)
        print(output.df)

        transcript_info = pd.read_csv(input['tx_to_gene'], sep="\t", header=None, index_col=0)
        transcript_info.columns = ['GENE_ID', 'GENE_SYMBOL']

        for file in input.abundance:
            print(file)
            df = pd.read_csv(file, sep="\t")
            df['gene'] = df['target_id'].map(transcript_info['GENE_ID'])
            tpm = df.groupby('gene').sum()['tpm']
            tpms.append(tpm)
            names.append(os.path.basename(os.path.dirname(file)))
            print(len(names))

        combined = pd.DataFrame(tpms, index=names)
        combined.to_csv(output.df, sep="\t", header=True, index=True)

rule generate_deseq_normalised_matrix:
    input:
        "fastq_dict.json",
        sample_info="data/raw/sample_info.txt",
        tx2gene="homo_sapiens/transcripts_to_genes.txt",
        h5=all_h5_files
    output:
        normalised="data/deseq/raw/tpm.tsv",
        sample_info="data/deseq/raw/sample_info.txt"
    script:
        "DESeq_processing.R"

rule construct_tensor_dataset:
    input:
        counts="data/{folder}/tpm.tsv",
        sample_info="data/raw/sample_info.txt"
    output:
        Y="data/tensor/{folder}/Y.txt",
        N="data/tensor/{folder}/N.txt",
        sample_info="data/tensor/{folder}/sample_info.txt",
        gene_names="data/tensor/{folder}/gene_names.txt",
    script:
        "tensor_processing.R"

rule extract_gene_names:
    input:
        raw_tpm="data/{folder}/tpm.tsv"
    output:
        Y="data/{folder}/Y.txt",
        gene_names="data/{folder}/gene_names.txt",
        sample_names="data/{folder}/sample_names.txt"
    shell:
        "tail -n +2 {input.raw_tpm} | cut -f 2- > {output.Y} && "\
        "head -n 1 {input.raw_tpm} | sed 's/\\t/\\n/g' | tail -n +2 > {output.gene_names} && "\
        "cut -f 1 {input.raw_tpm} | tail -n +2 > {output.sample_names}"

rule all_datasets:
    input:
        expand("data/{folder}/{file}",
               folder=["raw", "tensor/raw", "deseq/raw", "tensor/deseq/raw"],
               file=["Y.txt", "sample_info.txt", "gene_names.txt"])
