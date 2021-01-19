import pandas as pd
import logging

def convert_to_tensor_dataset(Y_file,
                              tensor_Y_file,
                              sample_info_file,
                              tensor_sample_info_file,
                              N_file):
    sample_info = pd.read_csv(sample_info_file, sep='\t')
    Y = pd.read_csv(Y_file,
                    sep='\t',
                    header=None)
    logging.info(sample_info)
    logging.info(Y.shape)
    logging.info(sample_info.shape)

    # Sort by cell type, then individual
    sorted_sample_info = sample_info.sort_values(['FactorValue..cell.type.',
                                                  'Characteristics..individual.'])
    tensor_sample_info = sorted_sample_info[sorted_sample_info['FactorValue..cell.type.'] != 'NK']
    tensor_samples = tensor_sample_info['Comment..ENA_SAMPLE.']
    logging.info(tensor_sample_info)
    logging.info(f"Should have 120 tensor samples, have {len(tensor_samples)}")

    sample_info.index = sample_info['Comment..ENA_SAMPLE.']
    row_indices = [sample_info.index.get_loc(sample_id)
                   for sample_id in tensor_samples]


    tensor_Y = Y.iloc[row_indices, :]
    tensor_Y.to_csv(tensor_Y_file, sep='\t', header=False, index=False)

    tensor_sample_info.to_csv(tensor_sample_info_file, sep='\t', index=False)

    N = tensor_sample_info['FactorValue..cell.type.'].nunique()
    logging.info(f"Should be 6 cell types, found {N}")

    with open(N_file, 'w') as f:
        f.write(str(N))

