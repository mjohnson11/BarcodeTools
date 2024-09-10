import os
import argparse
import json

import pandas as pd
import Levenshtein

from UnknownRegionParser import parse_by_regex


def make_dir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    
def run_regex_batch(
    batchfile,
    outdir,
    row_index=None,
    row_batch_size=1,
    filter_column=None,
    construct=None,
    trim_all_reads=False,
    autodetect_barcodes=False,
    regex_flanking_len=8,
    unknown_lens=None,
    quality_thresh=False,
):
    """
    Processes a batch of sequencing data using a CSV batch file.

    This function reads a CSV batch file specifying FASTQ files, 
    construct sequences, and processing options, and then performs 
    sequence parsing and error correction for each sample in the batch.

    Args:
        batchfile (str): The path to the CSV batch file. The file should 
            contain the following columns (additional columns are allowed):
                - 'Sample': The sample name.
                - 'fastq_file': The path to the input FASTQ file.
                - 'Construct_Seq': The DNA sequence of the known construct 
                   (can be omitted if `construct` is provided).
        outdir (str): The output directory where results and log files will
            be saved.
        row_index (int, optional):  If specified, processes only the row 
            with this index from the batch file. This is useful for 
            running the function in parallel using job arrays. 0-indexed 
            Default: None (processes all rows).
        row_batch_size (int, optional):  Allows you to batch a few rows
            at once. So if row_index is 3 and row_batch_size is 10, this
            will parse rows 30-39 (0-indexed). Default: 1.
        filter_column (str, optional): Used to check for inline indices.
            Should be the same as the name of one of the unknown regions
            in the construct, and should be a column in the batchfile. If 
            specified, this function will filter the parsed dataframe for 
            only rows where the region matches the sequence in the batchfile.
        construct (str, optional): The DNA sequence of the known 
            construct. If provided, this overrides any 'Construct_Seq' 
            entries in the batch file. Default: None.
        trim_all_reads (str or False): A number of bp to trim off the
            start of all reads. This argument is overridden if there is a
            'trim_read_start' column in the batch file. Default: False.
        error_correction_on (bool, optional): If True (default), performs 
            barcode error correction on the parsed results.
        autodetect_barcodes (bool, optional):  If True, unknown regions 
            are automatically detected in the construct sequence. 
            Default: False.
        regex_flanking_len (int, optional): Length of flanking sequences 
            for regex matching (used only if `parsing_method` is 'regex'). 
            Default: 8. 
        unknown_lens (str or list, optional): Controls how lengths of 
            unknown regions are handled (see `UnknownRegionParser` 
            documentation). Default: None.
        quality_thresh (int or False, optional):  The minimum average 
            quality score required for an unknown region for a read to be 
            considered valid. If False, quality filtering is skipped. 
            Default: False.
    """
    bf = pd.read_csv(batchfile)
    for j, row in bf.iterrows():
        # if row_index is defined, only parse one row
        if isinstance(row_index, int) and (j < row_index*row_batch_size or j >= (row_index+1)*row_batch_size): 
            continue
        s = row['Sample']
        refSeq = construct or row['Construct_Seq']
        if 'trim_read_start' in row:
            trim_read_start = row['trim_read_start']
        else:
            trim_read_start = trim_all_reads
        out = outdir+'/'+s
        make_dir(out)
        count_out = out+'/'+s+'_counts.csv'
        parse_by_regex(
            refSeq, 
            row['fastq_file'], 
            outfile=count_out,
            autodetect_barcodes=autodetect_barcodes,
            regex_flanking_len=regex_flanking_len,
            unknown_lens=unknown_lens,
            trim_read_start=trim_read_start,
            quality_thresh=quality_thresh
        )
        if filter_column:
            filtered_out = out+'/'+s+'_index_filtered.csv'
            expected = row[filter_column]
            td = pd.read_csv(count_out)
            # One error allowed
            td['Matches'] = td[filter_column].apply(lambda seq: seq==expected or Levenshtein.distance(seq,expected)<=1)
            filt_fail = td[td[filter_column].apply(lambda seq: 'Fail' not in seq)]
            filt_match = filt_fail[filt_fail['Matches']]
            keep_cols = [i for i in filt_match if i not in ['Matches', filter_column, 'Count']]
            filt_match[keep_cols+['Count']].groupby(keep_cols).sum().reset_index().sort_values(by='Count', ascending=False).to_csv(filtered_out, index=False)
            # updating log
            logfile = out+'/'+s+'_counts_log.json'
            with open(logfile, 'r') as infile:
                jd = json.load(infile)
            jd['stats']['nIndexWrong'] = len(filt_fail)-len(filt_match)
            jd['stats']['readsIndexWrong'] = int(filt_fail[~filt_fail['Matches']]['Count'].sum())
            jd['stats']['nIndexRight'] = len(filt_match)
            jd['stats']['readsIndexRight'] = int(filt_match['Count'].sum())
            with open(logfile, 'w') as logout:
                json.dump(jd, logout, indent=4)
            
            
        
        

            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process a batch of sequencing data."
    )
    parser.add_argument(
        "batchfile", type=str, help="Path to the CSV batch file"
    )
    parser.add_argument(
        "outdir", type=str, help="Output directory for results"
    )
    parser.add_argument(
        "--row_index",
        type=int,
        default=None,
        help="Process only the specified row from the batch file",
    )
    parser.add_argument(
        "--row_batch_size",
        type=int,
        default=1,
        help="Makes the row index go in this length chunks of rows",
    )
    parser.add_argument(
        "--filter_column",
        type=str,
        default=None,
        help="Column to use for inline index filtering (must match an unknown region name)",
    )
    parser.add_argument(
        "--construct",
        type=str,
        default=None,
        help="DNA sequence of the known construct (overrides batch file)",
    )
    parser.add_argument(
        "--trim_all_reads",
        type=int,
        default=False,
        help="Number of bp to trim from the start of all reads",
    )
    parser.add_argument(
        "--autodetect_barcodes",
        type=bool,
        default=False,
        help="Automatically detect unknown regions in the construct sequence",
    )
    parser.add_argument(
        "--regex_flanking_len",
        type=int,
        default=8,
        help="Length of flanking sequences for regex matching",
    )
    parser.add_argument(
        "--unknown_lens",
        type=str,
        default=None,
        help="How to handle lengths of unknown regions (see documentation)",
    )
    parser.add_argument(
        "--quality_thresh",
        type=int,
        default=False,
        help="Minimum average quality score for valid unknown regions",
    )

    args = parser.parse_args()

    run_regex_batch(
        args.batchfile,
        args.outdir,
        row_index=args.row_index,
        row_batch_size=args.row_batch_size,
        filter_column=args.filter_column,
        construct=args.construct,
        trim_all_reads=args.trim_all_reads,
        autodetect_barcodes=args.autodetect_barcodes,
        regex_flanking_len=args.regex_flanking_len,
        unknown_lens=args.unknown_lens,
        quality_thresh=args.quality_thresh
    )