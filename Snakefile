configfile: "config.yml"

import json

def all_abundance_files():
    fastq_dict = get_fastq_dict()
    return ["data/kallisto/{SAMPLE_ID}/abundance.tsv"
            for SAMPLE_ID in fastq_dict.keys()]

def get_fastq_files_for_ID(wildcards):
    fastq_dict = get_fastq_dict()
    sample_ID = wildcards.SAMPLE_ID
    return ["data/fastq/sample_ID/{full_id}.fastq.gz"
            for full_id in fastq_dict[sample_ID].keys()]

def get_fastq_url_from_fastq_id(wildcards):
    fastq_dict = get_fastq_dict()
    return {"url": fastq_dict[wildcards.SAMPLE_ID][wildcards.FASTQ_ID]}

def get_fastq_dict()
    fastq_dict_file = checkpoints.construct_fastq_dict.get().output[0]
    with open(fastq_dict_file, 'r') as f:
        fastq_dict = json.load(f)
    return fastq_dict

rule fetch_sample_info:
    output:
        "fastq_info.txt"
    shell:
        "wget -O {output} --no-verbose {config[fastq_info_url]}"

checkpoint construct_fastq_dict:
    input:
        "fastq_info.txt"
    output:
        "fastq_dict.json"
    run:
        import json
        import pandas as pd

        fastq_info = pd.read_csv(input[0], sep="\t")
        sample_col = 'Comment [ENA_SAMPLE]'
        url_col = 'Comment [FASTQ_URI]'
        id_col = 'FASTQ_ID'
        fastq_info[id_col] = fastq_info[url_col].str.extract(r'ftp://.*/(\w+)\.fastq\.gz')

        # Generate a dictionary with each key one of the samples,
        # and the value a dictionary from fastq_id to fastq_url for each
        # fastq file associated with the sample
        grouped = fastq_info.groupby(sample_col).loc[[id_col, url_col]]
        fastq_dict = grouped.apply(lambda df : dict(zip(df[id_col], df[url_col]))).to_dict()

        with open(output[0], 'w') as f:
            json.dump(fastq_dict, f)

rule fetch_fastq:
    group:
        "kallisto_sample"
    params:
        get_fastq_url_from_fastq_id
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
        get_fastq_files_for_ID
    log:
        "logs/data/kallisto/{SAMPLE_ID}.log"
    output:
        "data/kallisto/{SAMPLE_ID}/abundance.tsv"
    shell:
        "kallisto quant --index={config[INDEX_FILE]} --output-dir=$(dirname {output}) --single -l 50 -s 2 {input} 2>&1 | tee {log}"

rule generate_count_matrix:
    input:
        all_abundance_files
    output:
        df="data/tpm.tsv"
    run:
        import os
        import pandas as pd
        tpms = []
        names = []
        print(input)
        print(output.df)

        transcript_info = pd.read_csv("mus_musculus/transcripts_to_genes.txt", sep="\t", header=None, index_col=0)
        transcript_info.columns = ['GENE_ID', 'GENE_SYMBOL']

        for file in input:
            print(file)
            df = pd.read_csv(file, sep="\t")
            df['gene'] = df['target_id'].map(transcript_info['GENE_ID'])
            tpm = df.groupby('gene').sum()['tpm']
            tpms.append(tpm)
            names.append(os.path.basename(os.path.dirname(file)))
            print(len(names))

        combined = pd.DataFrame(tpms, index=names)
        combined.to_csv(output.df, sep="\t", header=True, index=True)
