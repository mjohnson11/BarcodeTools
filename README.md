# BarcodeTools

This is a collection of scripts meant to be useful for processing sequencing data with random barcodes or other variable regions. I've tried to make the tools somewhat modular so that if you have some python experience you can jump in and out, copy and paste and use what you need. Maybe at some point I will package it up more nicely! For now, here is a breakdown of the files:

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
This is just a wrapper to call UnknownRegionParser.py on a batch of files.

#### MergeAndCorrect.py
This is a wrapper to merge several count files and then error correct them.

#### util.py
A handful of simple functions used by these tools.


