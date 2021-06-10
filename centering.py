import pandas as pd

def center_Y(Y_file, Y_outfile):
    """Remove any genes with 0 variance and center Y so that each gene has mean 0.
    Read the nxp Y matrix in from Y_file and write centered version out to Y_outfile"""

    Y = pd.read_csv(Y_file,
                    sep='\t',
                    header=None)
    Y_trimmed = Y.loc[:, (Y.var(axis=0) != 0)]
    Y_centered = Y_trimmed - Y_trimmed.mean(axis=0)

    Y_centered.to_csv(Y_outfile, sep='\t', index=False)
