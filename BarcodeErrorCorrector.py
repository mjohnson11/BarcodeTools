import json
import re

import numpy as np
import pandas as pd

from DeletionCorrect import correct_bc_errors
from util import customJSONEncoder


def dictify_param(p, cols):
    if isinstance(p, dict):
        return p
    else:
        return {c: p for c in cols}
    
def record_run_stats(td, tmp, b, count_col, info_dict):
    # Just putting this up here to not clutter the main function
    # This just records a bunch of statistics about the error correction
    not_excluded = td[td[b+'_corrected'].apply(lambda bc: 'xcluded' not in str(bc))]
    excluded = td[td[b+'_corrected'].apply(lambda bc: 'xcluded' in str(bc))]
    info_dict[b]['nUnique'] = len(set(td[b]))
    info_dict[b]['nReads'] = np.sum(td[count_col])
    info_dict[b]['nUniqueAfterFiltering'] = len(set(tmp[b]))
    info_dict[b]['nReadsAfterFiltering'] = np.sum(tmp[count_col])
    info_dict[b]['nTrueBCs'] = len(set(not_excluded[b+'_corrected']))
    info_dict[b]['nCorrected'] = len(not_excluded[not_excluded[b]!=not_excluded[b+'_corrected']])
    info_dict[b]['readsUsed'] = np.sum(not_excluded[count_col])
    info_dict[b]['readsTrueBCs'] = np.sum(not_excluded[not_excluded[b]==not_excluded[b+'_corrected']][count_col])
    info_dict[b]['readsCorrected'] = np.sum(not_excluded[not_excluded[b]!=not_excluded[b+'_corrected']][count_col])
    info_dict[b]['nExcludedLowCount'] = len(excluded[excluded[b+'_corrected']=='excluded low count'])
    info_dict[b]['readsExcludedLowCount'] = np.sum(excluded[excluded[b+'_corrected']=='excluded low count'][count_col])
    info_dict[b]['nExcludedError'] = len(excluded[excluded[b+'_corrected']=='excluded error'])
    info_dict[b]['readsExcludedError'] = np.sum(excluded[excluded[b+'_corrected']=='excluded error'][count_col])

def error_correct_file_or_df(fileOrDf, outfile='return', logfile='auto', remove_off_bases=True, bc_cols='infer', count_col='Count', all_count_cols='infer', min_counts_for_centroid=2, max_edits=3, poisson_error_rate=0.1, do_not_correct=[]):
    """
    Corrects errors in barcode columns of a CSV file or pandas DataFrame.

    This function applies barcode error correction to specified columns,
    grouping similar barcodes and consolidating their counts. 

    Args:
        fileOrDf (str or pd.DataFrame): The path to the input CSV file or 
            a pandas DataFrame containing the data. 
        outfile (str, optional): The path to the output CSV file where the 
            corrected DataFrame will be saved. If set to 'return' (default), 
            the function returns the corrected DataFrame instead of saving.
        logfile (str, optional): The path to the log file (JSON format)
            containing information about the correction process. 
            If 'auto' (default), the log file will be saved in the same 
            directory as `outfile` with a '_log.json' suffix. 
            If set to None, no log file will be created.
        remove_off_bases (bool, optional): If True (default), barcode sequences 
            containing characters other than 'A', 'T', 'C', or 'G' will be 
            removed before error correction. Note that this will remove the
            "RegexFail" and "QualityFail" rows from a file produced by
            UnknownRegionParser (which is generally desirable).
        bc_cols (str or list, optional): Specifies the columns containing the 
            barcodes to correct:
                - 'infer' (default): Automatically detect columns with 
                  `object` dtype as barcode columns.
                - list of str: A list of column names to treat as barcode 
                  columns.
        count_col (str, optional): The name of the column containing the 
            counts associated with each barcode (default: 'Count').
        all_count_cols (str or list, optional): Specifies columns besides 
            `count_col` that contain count data to be aggregated during 
            the correction process.
                - 'infer' (default):  Automatically detect columns with numeric 
                  dtypes (`int`, `float`, 'int64', 'float64') as count 
                  columns.
                - list of str: A list of column names to treat as additional 
                  count columns.
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

    Returns:
        pd.DataFrame or None:  If `outfile` is set to 'return', the function 
        returns the corrected pandas DataFrame. Otherwise, it returns None 
        and saves the DataFrame to the specified `outfile`.
    """
       
    if isinstance(fileOrDf, pd.core.frame.DataFrame):
        td = fileOrDf.copy(deep=True)
    else:
        td = pd.read_csv(fileOrDf)
    if bc_cols == 'infer':
        bc_cols = [i for i in td if td[i].dtype==object]
        
    if all_count_cols == 'infer':
        all_count_cols = [i for i in td if td[i].dtype in (int, float, 'int64', 'float64')]
        
    min_counts_for_centroid = dictify_param(min_counts_for_centroid, bc_cols) 
    max_edits = dictify_param(max_edits, bc_cols) 
    poisson_error_rate = dictify_param(poisson_error_rate, bc_cols)
    
    run_info = {
        'outfile': outfile,
        'logfile': logfile,
        'remove_off_bases': remove_off_bases,
        'bc_cols': bc_cols,
        'count_col': count_col,
        'min_counts_for_centroid': min_counts_for_centroid,
        'max_edits': max_edits,
        'poisson_error_rate': poisson_error_rate,
        'do_not_correct': do_not_correct,
        'stats': {b: dict() for b in bc_cols}
    } 
        
    for b in bc_cols:
        if b in do_not_correct:
            continue
        if len(td.columns)==2: # no need to mess with groupby here
            tmp = td[[b, count_col]]
        else:
            tmp = td[[b, count_col]].groupby(b).sum().reset_index()
        if remove_off_bases:
            tmp = tmp[tmp[b].apply(lambda seq: bool(re.match(r'^[ATCG]+$', seq)))]
        corrector = correct_bc_errors(np.array(tmp[[b, count_col]]), min_counts_for_centroid=min_counts_for_centroid[b], max_edits=max_edits[b], poisson_error_rate=poisson_error_rate[b])
        td[b+'_corrected'] = td[b].map(corrector)
        record_run_stats(td, tmp, b, count_col, run_info['stats'])

    corr_cols = [b if b in do_not_correct else b+'_corrected' for b in bc_cols]
    final_df = td[corr_cols+all_count_cols].groupby(corr_cols).sum().reset_index().sort_values(by=count_col, ascending=False)
    if outfile == 'return':
        return final_df
    else:
        final_df.to_csv(outfile, index=False)
        if logfile == 'auto':
            assert '.csv' in outfile
            with open(outfile.replace('.csv', '_log.json'), 'w') as logout:
                json.dump(run_info, logout, cls=customJSONEncoder)