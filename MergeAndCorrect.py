from collections import defaultdict, Counter

import pandas as pd
import numpy as np

from BarcodeErrorCorrector import error_correct_file_or_df

def merge_and_correct_batch(
    batchfile,
    outdir,
    out_corrected,
    out_merged=None,
    filtered=False,
    min_counts_for_centroid=2, 
    max_edits=3, 
    poisson_error_rate=0.1, 
    do_not_correct=[]
):
    """
    Reads in a bunch of count files, merges them, and error corrects.
    Note this can be very memory intensive.

    Args:
        batchfile (str): The path to the CSV batch file. The file should 
            contain the following columns (additional columns are allowed):
                - 'Sample': The sample name.
        outdir (str): The output directory where results and log files are saved.
        out_corrected (str): output file for corrected data.
        out_merged (str, optional): Output filename for merged (uncorrected)
            file. If None (default), no file will be written.
        filtered (boolean, optional): if True, will look for files like 
            Sample_index_filtered.csv, otherwise will look for
            Sample_counts.csv
        min_counts_for_centroid (int or dict, optional): The minimum number 
            of counts a barcode needs to be considered a cluster centroid 
            during error correction (default: 2) 
                - int: Applies the same threshold to all barcode columns.
                - dict: A dictionary mapping barcode column names to 
                  specific thresholds. 
        max_edits (int or dict, optional): The maximum number of edit 
            operations (substitutions, insertions, deletions) allowed 
            when correcting barcodes (default: 3).
                - int: Applies the same edit distance to all barcode columns.
                - dict:  A dictionary mapping barcode column names to 
                  specific edit distances.
        poisson_error_rate (float or dict, optional):  The estimated error 
            rate of the sequencing process, modeled using a Poisson 
            distribution. This value is used to calculate the probability 
            of different barcode errors (default: 0.1). Note that this is
            a VERY permissive estimate of the error rate, which works well
            in practice for error rates up to 10% (10% might occur in long
            mononucleotide runs).
                - float: Applies the same error rate to all barcode columns.
                - dict:  A dictionary mapping barcode column names to 
                  specific error rates.
        do_not_correct (list, optional): A list of barcode column names 
            to exclude from error correction (default: []).
    """
    bf = pd.read_csv(batchfile)
    comb_counts = defaultdict(Counter)
    samples = []
    for j, row in bf.iterrows():
        s = row['Sample']
        samples.append(s)
        out = outdir+'/'+s
        fname = out+'/'+s+'_index_filtered.csv' if filtered else out+'/'+s+'_counts.csv'
        td = pd.read_csv(fname)
        # (These columns have to be the same in all files)
        unknown_cols = [i for i in td if i!='Count']
        td['Combined'] = td.apply(lambda row: '_'.join([row[c] for c in unknown_cols]), axis=1)
        assert len(td) == len(set(td['Combined'])) # There should be no repeat unknown region combos
        for comb, count in np.array(td[['Combined', 'Count']]):
            comb_counts[comb][s] = count
            
    mat = []
    for comb in comb_counts:
        tmp = comb_counts[comb]
        mat.append(comb.split('_')+[tmp[s] for s in samples])
        
    merged = pd.DataFrame(mat, columns=unknown_cols+samples)
    merged['Total_Counts'] = np.sum(merged[samples], axis=1)
    merged = merged.sort_values(by='Total_Counts', ascending=False)
    
    if out_merged:
        merged.to_csv(out_merged, index=False)
    
    error_correct_file_or_df(
            merged, 
            outfile=out_corrected,
            bc_cols=unknown_cols, 
            count_col='Total_Counts', 
            all_count_cols=['Total_Counts']+samples, 
            min_counts_for_centroid=min_counts_for_centroid, 
            max_edits=max_edits,
            poisson_error_rate=poisson_error_rate, 
            do_not_correct=do_not_correct
        )
       