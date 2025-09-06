import pandas as pd
import subprocess
import tempfile
import argparse
import csv

def format_bowtie_row(raw_row):
    tmp_dict = {
        'bowtie_code': raw_row[1], 
        'contig': raw_row[2], 
        'mapq': raw_row[4], 
        'cigar_match': raw_row[5]
    }
    # changes insertion location depending on the orientation of alignment
    if raw_row[1] == '16':
        tmp_dict['strand'] = '-'
        tmp_dict['insertion_edge'] = int(raw_row[3]) + len(raw_row[0]) - 1
    else:
        tmp_dict['strand'] = '+'
        tmp_dict['insertion_edge'] = int(raw_row[3])

    return tmp_dict

def run_bowtie(seqs, bowtie_base_directory):
    with tempfile.NamedTemporaryFile(mode='w+t', delete=True) as temp_fasta_file:
        with tempfile.NamedTemporaryFile(mode='w+t', delete=True) as temp_bowtie_output:
            for seq in seqs:
                temp_fasta_file.write('>' + seq + '\n' + seq + '\n')
            temp_fasta_file.flush()  # Flush to ensure data is written
            try:
                subprocess.call(['bowtie2', '-x', bowtie_base_directory, '-U', temp_fasta_file.name, '-S', temp_bowtie_output.name, '-f', '--local'])
            except subprocess.CalledProcessError as e:
                print(f"Error running Bowtie2: {e}")
                return None
            align_info = dict()
            with open(temp_bowtie_output.name, 'r') as infile:
                reader = csv.reader(infile, delimiter='\t')
                reading = False
                for row in reader:
                    if not reading:
                        if row[0] == '@PG':
                            reading = True
                    else:
                        align_info[row[0]] = format_bowtie_row(row)
    return align_info

def align_column_to_ref(infile_or_df, bowtie_ref, outfile, alignment_column):
    
    if isinstance(infile_or_df, pd.core.frame.DataFrame):
        df = infile_or_df
    else:
        df = pd.read_csv(infile_or_df)
    seqs_to_align = set(df[pd.notnull(df[alignment_column]) & df[alignment_column].apply(lambda a: 'excluded' not in a)][alignment_column])
    align_info = run_bowtie(seqs_to_align, bowtie_ref)
    align_cols = ['bowtie_code', 'contig', 'insertion_edge', 'strand', 'mapq', 'cigar_match']
    
    for col in align_cols:
        df[col] = df[alignment_column].apply(lambda s: align_info.get(s, {}).get(col, ''))
    
    if outfile == 'return':
        return df
    else:
        df.to_csv(outfile, index=False)
             


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Align one column of a csv to a reference using bowtie2"
    )
    parser.add_argument(
        "inputfile", type=str, help="Path to the CSV"
    )
    parser.add_argument(
        "aligncol", type=str, help="Column to align"
    )
    parser.add_argument(
        "reference", type=str, help="Path to the reference bowtie build"
    )
    parser.add_argument(
        "outfile", type=str, help="Output file for results"
    )

    args = parser.parse_args()

    align_column_to_ref(
        args.inputfile,
        args.reference,
        args.outfile,
        args.aligncol
    )
