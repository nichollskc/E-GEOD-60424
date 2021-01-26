configfile: "config.yml"

include: "biclustering_comparison/Snakefile"

localrules: generate_count_matrix, restrict_to_expressed_protein, find_expressed_protein_genes

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
        grouped = fastq_info.groupby(sample_col)[id_col, url_col]
        fastq_dict = grouped.apply(lambda df : dict(zip(df[id_col], df[url_col]))).to_dict()

        with open(output[0], 'w') as f:
            json.dump(fastq_dict, f, indent=2)

rule fetch_fastq_info:
    output:
        "fastq_info.txt"
    shell:
        "wget -O {output} --no-verbose {config[fastq_info]}"

rule process_sample_info:
    input:
        fastq_info="fastq_info.txt"
    output:
        sample_info="data/real/raw/sample_info.txt"
    script:
        "process_sample_info.R"

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
        df="data/real/raw/tpm.tsv"
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
        sample_info="data/real/raw/sample_info.txt",
        tx2gene="homo_sapiens/transcripts_to_genes.txt",
        h5=all_h5_files
    output:
        vst_normalised="data/real/deseq/raw/tpm.tsv",
        vst_sample_info="data/real/deseq/raw/sample_info.txt",
        sf_normalised="data/real/deseq_sf/raw/tpm.tsv",
        sf_sample_info="data/real/deseq_sf/raw/sample_info.txt",
    script:
        "DESeq_processing.R"

rule extract_gene_names:
    input:
        raw_tpm="data/real/{folder}/tpm.tsv"
    output:
        Y="data/real/{folder}/Y.txt",
        gene_names="data/real/{folder}/gene_names.txt",
        sample_names="data/real/{folder}/sample_names.txt"
    shell:
        "tail -n +2 {input.raw_tpm} | cut -f 2- > {output.Y} && "\
        "head -n 1 {input.raw_tpm} | sed 's/\\t/\\n/g' | tail -n +1 > {output.gene_names} && "\
        "cut -f 1 {input.raw_tpm} | tail -n +2 > {output.sample_names}"

rule find_expressed_protein_genes:
    input:
        tpm="data/real/deseq_sf/raw/tpm.tsv",
        sample_info="data/real/deseq_sf/raw/sample_info.txt",
    output:
        gene_names="data/real/deseq_sf/raw/expressed_protein_gene_names.txt"
    run:
        import select_genes as sg
        gene_names = sg.restrict_to_expressed_protein(input.tpm, input.sample_info)
        with open(output.gene_names, 'w') as f:
            f.write('\n'.join(gene_names))

rule restrict_to_expressed_protein:
    input:
        Y="data/real/{folder}/Y.txt",
        gene_names="data/real/{folder}/gene_names.txt",
        expressed_protein_genes="data/real/deseq_sf/raw/expressed_protein_gene_names.txt"
    output:
        Y="data/real/{folder}/expressed/Y.txt",
        gene_names="data/real/{folder}/expressed/gene_names.txt",
    run:
        import shutil
        import biclust_comp.utils as utils

        Y = utils.read_matrix_tsv(input.Y)
        all_genes = utils.read_list_from_file(input.gene_names, strip_quotes=True)
        keep_genes = utils.read_list_from_file(input.expressed_protein_genes)
        # Find the indices of these genes in our whole gene list
        keep_gene_indices = list(map(lambda x: all_genes.index(x),
                                     keep_genes))

        restricted = Y.iloc[:, keep_gene_indices]
        restricted.to_csv(output.Y, sep='\t', header=False, index=False)

        shutil.copy(input.gene_names, output.gene_names)

rule tensor_dataset:
    input:
        sample_info="data/real/raw/sample_info.txt",
        Y="data/real/{folder}/expressed/Y.txt",
        gene_names="data/real/{folder}/expressed/gene_names.txt",
    output:
        Y="data/real/{folder}/expressed/tensor/Y.txt",
        N="data/real/{folder}/expressed/tensor/N.txt",
        sample_info="data/real/{folder}/expressed/tensor/sample_info.txt",
        gene_names="data/real/{folder}/expressed/tensor/gene_names.txt",
    run:
        import shutil
        import tensor_processing as tp

        tp.convert_to_tensor_dataset(input.Y,
                                     output.Y,
                                     input.sample_info,
                                     output.sample_info,
                                     output.N)

        shutil.copy(input.gene_names, output.gene_names)

rule resave_with_np:
    input:
        Y="data/real/{folder}/Y.txt",
    output:
        Y="data/real/{folder}/Y_resaved.txt",
    run:
        import np
        Y = np.loadtxt(input.Y, delimiter='\t')
        np.savetxt(output.Y, Y, delimiter='\t')

rule all_datasets:
    input:
        expand("data/real/{folder}/expressed{tensor}/{file}",
               folder=["raw", "deseq/raw", "deseq_sf/raw",
                       "log", "deseq/log", "deseq_sf/log"],
               tensor=["/tensor", ""],
               file=["Y.txt", "Y_resaved.txt"])

rule all_runs:
    input:
        expand("logs/{method_dataset_runid}.log",
               method_dataset_runid=config['REAL_DATASET_METHOD_RUNIDS'])

rule all_runs_force:
    input:
        expand("results/{method_dataset_runid}/X.txt",
               method_dataset_runid=config['REAL_DATASET_METHOD_RUNIDS'])

