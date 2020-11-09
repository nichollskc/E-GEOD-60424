import pandas as pd

import mygene

# http://localhost:8888/notebooks/Rough_2020-11-06_Presnell.ipynb

def restrict_to_expressed_genes(gene_expression_df, sample_info_df, threshold):
    """Restrict to genes with median expression greater than threshold in at least one cell type"""
    merged = gene_expression_df.merge(sample_info_df[['sample_id', 'cell_type']],
                                      left_index=True,
                                      right_on='sample_id').set_index('sample_id')
    grouped = merged.groupby(['cell_type']).median()
    max_median = grouped.max().reset_index()
    max_median.columns = ['gene', 'max_median_expression']

    expressed_genes = max_median['gene'][(max_median['max_median_expression'] > threshold)]
    expressed = gene_expression_df[expressed_genes]
    return expressed


def restrict_to_protein_coding(gene_expression_df):
    unversioned_gene_names = [name.split('.')[0] for name in gene_expression_df.columns]
    assert len(unversioned_gene_names) == len(set(unversioned_gene_names))

    mg = mygene.MyGeneInfo()
    gene_info = mg.getgenes(unversioned_gene_names, fields=['accession', 'pathway'])

    reduced_gene_info = {}
    for gene in gene_info:
        reduced_gene_info[gene['query']] = {'has_protein': 'accession' in gene and 'protein' in gene['accession'],
                                            'in_pathway': 'pathway' in gene}
    gene_info_df = pd.DataFrame.from_dict(reduced_gene_info, orient='index')

    # So there should be one-one mapping from unversioned to versioned
    reversion_dict = dict(zip(unversioned_gene_names, gene_expression_df.columns))

    protein_coding = gene_info_df.index[gene_info_df['has_protein']]
    protein_coding_versioned = protein_coding.map(reversion_dict)
    restricted = gene_expression_df[protein_coding_versioned]
    return restricted

def restrict_to_expressed_protein(tpm_file, sample_info):
    deseq_sf = pd.read_csv(tpm_file, sep='\t', index_col=0)
    sample_info = pd.read_csv(sample_info, sep='\t')
    sample_info['sample_id'] = sample_info['Comment..ENA_SAMPLE.']
    sample_info['cell_type'] = sample_info['FactorValue..cell.type.']

    expressed = restrict_to_expressed_genes(deseq_sf, sample_info, 0)
    print(expressed.shape)
    protein_coding = restrict_to_protein_coding(expressed)
    print(protein_coding.shape)
    return protein_coding.columns
