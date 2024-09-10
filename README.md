# BarcodeTools

This is a collection of scripts meant to be useful for processing sequencing data with random barcodes or other variable regions. I've tried to make the tools somewhat modular so that if you have some python experience you can jump in and out, copy and paste and use what you need. Maybe at some point I will package it up more nicely, for now I am still working on making it usable for myself and others. 

Specific dependencies:
- [regex](https://pypi.org/project/regex/) (not the built-in module)
- [Levenshtein](https://pypi.org/project/python-Levenshtein/)
- [mappy](https://pypi.org/project/mappy/)

Common dependencies: pandas, numpy, scipy

Here is a breakdown of the files:

#### UnknownRegionParser.py
This contains a class for parsing fastq files and extracting unknown regions. Extraction can be done using regex or alignment (using minimap2 via mappy). The alignment approach currently works best for medium (150 bp) to long reads, and I haven't tested it extensively. In both cases, the idea is that you provide a "construct" sequence to the parser, and then give it a fastq file to parse. It can detect regions in the construct with N/S/W bases automatically, but I have been generally giving it a formatted construct sequence, e.g.
`'NNNAACCAGAGGGGGATAGATGTCCACGAGGTCTCT(BC:NNNNNNNNNNNNNNNNNNNN)CGTACGCTGCAGGTCGACCAGCGGCCGGCCAGACCGGGGACTTATCAGCCAACCTGTTA(Edge:NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN)'`

This is the expected amplicon sequence, and now the parser will know to extract the BC and Edge sequences. Using the regex option, this will create regex patterns around both of these regions and use them to extract them (since Edge is at the end, it will only use the bases before it in the regex). I'd like to add an option to input your own regex patterns, but for now you'd need to edit the code. Currently the program uses an increasingly lenient list of regex patterns to extract these regions: 
1. exact
2. exact flanking regions with +-1 base in the region
3. exact flanking regions with +-2 bases in the region
4. flanking regions with up to one error
5. flanking regions with up to one error with +-1 base in the region
6. flanking regions with up to one error with +-2 bases in the region

There are options to change the regex flanking sequence lengths (default is up to 8 bp), or let an unknown region have any length. There is also a test_regex method that can be used to see how the extraction is working on a few reads.

For the alignment-based extraction, there is an option to reindex the read so it starts where the construct sequence starts (assuming the molecule is circular).

#### BarcodeErrorCorrector.py
This file contains the error_correct_file_or_df function, which will error correct the sequences all string columns in a csv file or pandas dataframe (or a set of columns provided by the user). The error correction is done by shared-deletion-networks and double checked by Levenshtein distance.

#### DeletionCorrect.py
This is the algorithm for correction called by BarcodeErrorCorrector.py

#### BatchRegexParser.py
This is just a wrapper to call UnknownRegionParser.py on a batch of files. I think the likely use case would be making a copy of this file for your own use.

#### MergeAndCorrect.py
This is a wrapper to merge several count files and then error correct them. I think the likely use case would be making a copy of this file for your own use.

#### util.py
A handful of simple functions used by these tools.


## More detailed documentation for the two main files

## Unknown Region Parser

This script parses FASTQ files and extracts sequences from unknown regions within DNA sequences, based on a provided known construct sequence (e.g., a plasmid). It offers two parsing methods: regular expression-based matching and alignment. 

### Functions:

**1. `parse_by_regex`**

   - **Description:** Parses a FASTQ file and extracts sequences from unknown regions using regular expressions. 
   - **Usage:**
     ```python
     parse_by_regex(construct, fastq_file, outfile='return', logfile='auto', construct_is_file=False, autodetect_barcodes=False, regex_flanking_len=8, unknown_lens=None, trim_read_start=False, quality_thresh=False, read_cutoff=None)
     ```
   - **Arguments:**
     - `construct (str)`: The DNA sequence of the known construct. This can include placeholders for unknown regions (see below).
     - `fastq_file (str)`: Path to the input FASTQ file.
     - `outfile (str, optional)`: The path to the output CSV file where the parsed results will be saved. If set to 'return' (default), the function returns the parsed DataFrame.
     - `logfile (str, optional)`: The path to the output log file. If set to 'auto' (default), the log file will be saved in the same directory as outfile with a '_log.json' suffix. If None, no log file is created.
     - `construct_is_file (bool, optional)`: If True, `construct` is interpreted as the path to a file containing the construct sequence. Default: False.
     - `autodetect_barcodes (bool, optional)`:  If True, unknown regions are automatically detected as consecutive stretches of 'N', 'W', or 'S' characters in the `construct` sequence. If False, unknown regions must be explicitly defined using parentheses with a region name (e.g., "(region1:NNNNN)") within the `construct` sequence. Default: False.
     - `regex_flanking_len (int, optional)`: The length of flanking sequences around unknown regions to use for regular expression matching. Default: 8.
     - `unknown_lens (str or list, optional)`: Controls how the lengths of unknown regions are handled: - None (default): Use the actual lengths of the unknown regions defined in `construct`. - 'All': Treat all unknown regions as having unknown lengths during parsing. - list of str: A list of unknown region names to treat as having unknown lengths.
     - `trim_read_start (int or False, optional)`: A number of bases to trim off the start of the read. If False, no trimming is performed. Default: False.
     - `quality_thresh (int or False, optional)`:  The minimum average quality score required for the flanking regions around an unknown region for a read to be considered valid. If False, quality filtering is skipped. Default: False.
     - `read_cutoff (int or None, optional)`: If not None, limits the parsing to the specified number of reads. Default: None.
   - **Returns:**
     - `pd.DataFrame or None`: If `outfile` is 'return', returns the parsed DataFrame. Otherwise, returns None and saves the DataFrame to the specified file. 

**2. `parse_by_alignment`**

   - **Description:** Parses a FASTQ file and extracts sequences from unknown regions using alignment to the known construct.
   - **Usage:**
     ```python
     parse_by_alignment(construct, fastq_file, outfile='return', logfile='auto', construct_is_file=False, autodetect_barcodes=False, unknown_lens=None, read_cutoff=None)
     ```
   - **Arguments:**
     - `construct (str)`: The DNA sequence of the known construct. This can include placeholders for unknown regions (see below).
     - `fastq_file (str)`: Path to the input FASTQ file.
     - `outfile (str, optional)`: The path to the output CSV file where the parsed results will be saved. If set to 'return' (default), the function returns the parsed DataFrame.
     - `logfile (str, optional)`:  The path to the output log file. If set to 'auto' (default), the log file will be saved in the same directory as outfile with a '_log.json' suffix. If None, no log file is created.
     - `construct_is_file (bool, optional)`: If True, `construct` is interpreted as the path to a file containing the construct sequence. Default: False.
     - `autodetect_barcodes (bool, optional)`:  If True, unknown regions are automatically detected as consecutive stretches of 'N', 'W', or 'S' characters in the `construct` sequence. If False, unknown regions must be explicitly defined using parentheses with a region name (e.g., "(region1:NNNNN)") within the `construct` sequence. Default: False.
     - `unknown_lens (str or list, optional)`: Controls how the lengths of unknown regions are handled:
        - None (default): Use the actual lengths of the unknown regions defined in `construct`.
        - 'All': Treat all unknown regions as having unknown lengths during parsing. 
        - list of str: A list of unknown region names to treat as having unknown lengths.
     - `read_cutoff (int or None, optional)`: If not None, limits the parsing to the specified number of reads. Default: None.
   - **Returns:**
     - `pd.DataFrame or None`: If `outfile` is 'return', returns the parsed DataFrame. Otherwise, returns None and saves the DataFrame to the specified file. 

### Example:
```python
parse_by_regex("ATGC(region1:NNNNN)GCTT(region2:WWWWW)A", "reads.fastq", outfile="results.csv") 
```
This will parse the FASTQ file "reads.fastq", extracting sequences from regions named "region1" and "region2" as defined in the provided construct sequence, and save the results to "results.csv".

## Barcode Error Correction

This script corrects errors in barcode columns within a CSV file or pandas DataFrame. It groups similar barcodes together, consolidating their counts and correcting likely errors.

### Function:

**`error_correct_file_or_df`**

- **Description:**  Corrects errors in barcode columns of a CSV file or pandas DataFrame.
- **Usage:**
  ```python
  error_correct_file_or_df(fileOrDf, outfile='return', logfile='auto', remove_off_bases=True, bc_cols='infer', count_col='Count', all_count_cols='infer', min_counts_for_centroid=2, max_edits=3, poisson_error_rate=0.1, do_not_correct=[])
  ```
- **Arguments:**
  - `fileOrDf (str or pd.DataFrame)`: The path to the input CSV file or a pandas DataFrame containing the data. 
  - `outfile (str, optional)`: The path to the output CSV file where the corrected DataFrame will be saved. If set to 'return' (default), the function returns the corrected DataFrame instead of saving.
  - `logfile (str, optional)`: The path to the log file (JSON format) containing information about the correction process. If 'auto' (default), the log file will be saved in the same directory as `outfile` with a '_log.json' suffix. If set to None, no log file will be created.
  - `remove_off_bases (bool, optional)`: If True (default), barcode sequences containing characters other than 'A', 'T', 'C', or 'G' will be removed before error correction. Note that this will remove the "RegexFail" and "QualityFail" rows from a file produced by UnknownRegionParser (which is generally desirable).
  - `bc_cols (str or list, optional)`: Specifies the columns containing the barcodes to correct:
    - 'infer' (default): Automatically detect columns with `object` dtype as barcode columns.
    - list of str: A list of column names to treat as barcode columns.
  - `count_col (str, optional)`: The name of the column containing the counts associated with each barcode (default: 'Count').
  - `all_count_cols (str or list, optional)`: Specifies columns besides `count_col` that contain count data to be aggregated during the correction process.
    - 'infer' (default):  Automatically detect columns with numeric dtypes (`int`, `float`, 'int64', 'float64') as count columns.
    - list of str: A list of column names to treat as additional count columns.
  - `min_counts_for_centroid (int or dict, optional)`: The minimum number of counts a barcode needs to be considered a cluster centroid during error correction (default: 2). 
    - int: Applies the same threshold to all barcode columns.
    - dict: A dictionary mapping barcode column names to specific thresholds. 
  - `max_edits (int or dict, optional)`: The maximum number of edit operations (substitutions, insertions, deletions) allowed when correcting barcodes (default: 3).
    - int: Applies the same edit distance to all barcode columns.
    - dict:  A dictionary mapping barcode column names to specific edit distances.
  - `poisson_error_rate (float or dict, optional)`:  The estimated error rate of the sequencing process, modeled using a Poisson distribution. This value is used to calculate the probability of different barcode errors (default: 0.1). Note that this is a VERY permissive estimate of the error rate, which works well in practice for error rates up to 10% (10% might occur in long mononucleotide runs).
    - float: Applies the same error rate to all barcode columns.
    - dict:  A dictionary mapping barcode column names to specific error rates.
  - `do_not_correct (list, optional)`: A list of barcode column names to exclude from error correction (default: []).
- **Returns:**
  - `pd.DataFrame or None`:  If `outfile` is set to 'return', the function returns the corrected pandas DataFrame. Otherwise, it returns None and saves the DataFrame to the specified `outfile`.

### Example:
```python
error_correct_file_or_df("barcode_counts.csv", outfile="corrected_barcodes.csv", bc_cols=['barcode1', 'barcode2'], do_not_correct=['long_unknown_region'], min_counts_for_centroid=1)
```

This will read barcodes from "barcode_counts.csv" and correct errors in columns 'barcode1' and 'barcode2' but not in the 'long_unknown_region' column, using a minimum count of 1 for true barcodes, and save the corrected results to "corrected_barcodes.csv".



