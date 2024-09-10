import gzip
import json
from collections import Counter

import mappy as mp
import numpy as np
import pandas as pd
import regex

from util import (
    FourLineFastq,
    cigar_to_ref_array,
    color_dna,
    customJSONEncoder,
    reorient_circular_read,
)


class UnknownRegionParser:
    """
    Parses and extracts sequences from unknown regions within DNA sequences.

    This class provides functionality to define unknown regions within a 
    known construct sequence (e.g., a plasmid) and then extract the 
    sequences from those regions in a set of sequencing reads. It offers 
    two parsing methods: regular expression-based matching and alignment. 
    """
    def __init__(self, construct, construct_is_file=False, autodetect_barcodes=False, regex_flanking_len=8, unknown_lens=None):
        """
        Initializes the UnknownRegionParser with construct sequence and options.

        Args:
            construct (str): The DNA sequence of the known construct. This 
                can include placeholders for unknown regions (see below).
            construct_is_file (bool, optional): If True, `construct` is 
                interpreted as the path to a file containing the construct 
                sequence. Default: False.
            autodetect_barcodes (bool, optional):  If True, unknown regions 
                are automatically detected as consecutive stretches of 'N', 
                'W', or 'S' characters in the `construct` sequence. 
                If False, unknown regions must be explicitly defined using 
                parentheses with a region name (e.g., "(region1:NNNNN)") 
                within the `construct` sequence. Default: False.
            regex_flanking_len (int, optional): The length of flanking 
                sequences around unknown regions to use for regular 
                expression matching. Default: 8.
            unknown_lens (str or list, optional): Controls how the lengths 
                of unknown regions are handled:
                - None (default): Use the actual lengths of the unknown 
                  regions defined in `construct`.
                - 'All': Treat all unknown regions as having unknown 
                  lengths during parsing. 
                - list of str: A list of unknown region names to treat as 
                  having unknown lengths.
        """
        if construct_is_file:
            with open(construct, 'r') as infile:
                self.construct = ''.join([line.strip() for line in infile if line[0]!='>'])
        else:
            self.construct = construct
        self.autodetect_barcodes = autodetect_barcodes


        self.regex_flanking_len = regex_flanking_len
        # unknown_lens is either "All" or a list of unknown region names
        self.unknown_lens = unknown_lens
        # currently defaults, may become params
        self.regex_buffer = 10
        
        
        self.refSeq = None
        self.unknown_regions = []

        self._parse_construct()

    def _parse_construct(self):
        """
        Internal method to parse the construct sequence and define unknown regions. 

        This method populates the `refSeq` (reference sequence) and 
        `unknown_regions` attributes based on the `construct` and other 
        initialization parameters.
        """
        if self.autodetect_barcodes:
            self.refSeq = self.construct
            unknowns = regex.findall(r"[NWS]+", self.construct)
            # note that split returns an empty string if the sequence starts
            # or ends with NWS
            inbetweens = regex.split(r"[NWS]+", self.construct)
            running_index = 0
            for i in range(len(unknowns)):
                running_index += len(inbetweens[i])
                self.unknown_regions.append({
                    'name': 'unknown_region_'+str(i+1),
                    'start': running_index,
                    'end': running_index + len(unknowns[i]),
                    'flanking_seq_left': inbetweens[i][-1*self.regex_flanking_len:],
                    'flanking_seq_right': inbetweens[i+1][:self.regex_flanking_len],
                    'seq': unknowns[i],
                    'len': len(unknowns[i])
                })
                running_index += len(unknowns[i])
        else:
            pattern = r"\((.*?):(.*?)\)"
            unknowns = regex.findall(pattern, self.construct)
            # note that this includes the groups from the pattern,
            # so we have to count by 3's (skip 2)
            inbetweens = regex.split(pattern, self.construct)
            running_index = 0
            self.refSeq = ''
            for i in range(len(unknowns)):
                running_index += len(inbetweens[i*3])
                self.refSeq += inbetweens[i*3]+unknowns[i][1]
                self.unknown_regions.append({
                    'name': unknowns[i][0],
                    'start': running_index,
                    'end': running_index + len(unknowns[i][1]),
                    'flanking_seq_left': inbetweens[i*3][-1*self.regex_flanking_len:],
                    'flanking_seq_right': inbetweens[(i+1)*3][:self.regex_flanking_len],
                    'seq': unknowns[i][1],
                    'len': len(unknowns[i][1])
                })
                running_index += len(unknowns[i][1])
            self.refSeq += inbetweens[-1]
        if self.unknown_lens == 'All':
            for u in self.unknown_regions:
                u['len'] = 'unknown_len'
        elif self.unknown_lens: # it must be a list of names
            for u in self.unknown_regions:
                if u['name'] in self.unknown_lens:
                    u['len'] = 'unknown_len'
        self.refSeqLen = len(self.refSeq)

    def create_regex_lists(self):
        """
        creates an increasingly lenient list of regexes 
        for each unknown region
        """
        for unknown in self.unknown_regions:
            unknown['search_start'] = max(0, unknown['start']-self.regex_flanking_len-self.regex_buffer)
            unknown['search_end'] = min(len(self.refSeq), unknown['end']+self.regex_flanking_len+self.regex_buffer)
            c1 = '('+unknown['flanking_seq_left']+')'
            c2 = '('+unknown['flanking_seq_right']+')'
            if unknown['len'] == 'unknown_len':
                bc_reg = '(.*?)'
                unknown['regexes'] = [
                    regex.compile(c1+bc_reg+c2),
                    regex.compile(c1+'{e<=1}'+bc_reg+c2+'{e<=1}'),
                ]
            else:
                bc_reg = '(.{'+str(unknown['len'])+'})'
                bc_1off = '(.{'+str(unknown['len']-1)+','+str(unknown['len']+1)+'})'
                bc_2off = '(.{'+str(unknown['len']-2)+','+str(unknown['len']+2)+'})'
                unknown['regexes'] = [
                    regex.compile(c1+bc_reg+c2),
                    regex.compile(c1+bc_1off+c2),
                    regex.compile(c1+bc_2off+c2),
                    regex.compile(c1+'{e<=1}'+bc_reg+c2+'{e<=1}'),
                    regex.compile(c1+'{e<=1}'+bc_1off+c2+'{e<=1}'),
                    regex.compile(c1+'{e<=1}'+bc_2off+c2+'{e<=1}'),
                ]
        
    def regex_example(self, seq, qual):
        # makes a nice-to-look at regex hit display for testing
        
        def set_chars(char_array, position, string):
            # utility function to make hit string
            while (position+len(string)) > len(char_array): 
                # in case the hit text exceeds the sequence length
                char_array.append(' ')
            for i in range(len(string)):
                char_array[position+i] = string[i]
            
        
        hit_string = [' '] * len(self.refSeq)
        for u in self.unknown_regions:
            s = u['search_start']
            e = u['search_end']
            if np.mean(np.frombuffer(qual[s:e].encode('ascii'), dtype=np.uint8))-33 < self.quality_cutoff:
                set_chars(hit_string, u['start'], 'QualityFail')
            else:
                regex_worked = False
                for regex_pattern in u['regexes']:
                    reghit = regex_pattern.search(seq[s:e])
                    if reghit:
                        regex_worked = True
                        set_chars(hit_string, reghit.start()+s+len(reghit.group(1)), reghit.group(2))
                        break
                if not regex_worked:
                    set_chars(hit_string, u['start'], 'RegexFail')
                
        print(color_dna('\n'.join([self.refSeq, seq, ''.join(hit_string)])))

    def test_regex(self, fastq_file, trim_read_start=False, read_lim=10, tail_from=1000):
        """
        Looks at the last 10 reads of the first 1000 reads
        or the last 10 reads if there are < 1000.
        Just used for visualization and checking regexes are working.
        """
        self.create_regex_lists()
        if fastq_file.endswith(".gz"):
            opener = gzip.open
        else:
            opener = open
        rc = 0
        recs = []
        with opener(fastq_file, 'rt') as infile:
            for title, seq, qual in FourLineFastq(infile):
                if trim_read_start:
                    seq = seq[trim_read_start:]
                    qual = qual[trim_read_start:]
                recs.append((title, seq, qual))
                rc += 1
                if rc >= tail_from:
                    break
        for title, seq, qual in recs[-1*read_lim:]:
            print(title)
            self.regex_example(seq, qual)

    def regex_parse(self, seq, qual, u, quality_thresh):
        """
        Applies regular expressions to extract the sequence from a single unknown region.
        
        Args:
            seq (str): The DNA sequence read.
            qual (str): The quality scores for the read.
            u (dict): A dictionary representing a single unknown region (from 
                the `unknown_regions` list).
            quality_thresh (int or False): The minimum average quality score 
                required for the flanking regions around the unknown region 
                to be considered a valid match. If False, quality filtering 
                is disabled.

        Returns:
            str: The extracted sequence from the unknown region if a match 
                 is found and passes quality filtering.  Otherwise, returns
                 'QualityFail' or 'RegexFail'.
        """
        s = u['search_start']
        e = u['search_end']
        # if quality_thresh if False, we skip this step, which greatly improves speed
        if quality_thresh and np.mean(np.frombuffer(qual[s:e].encode('ascii'), dtype=np.uint8))-33 < quality_thresh:
            return 'QualityFail'
        for regex_pattern in u['regexes']:
            reghit = regex_pattern.search(seq[s:e])
            if reghit:
                return reghit.group(2)
        return 'RegexFail'
        
    def parse_fastq_regex(self, fastq_file, trim_read_start=False, quality_thresh=False, progress_read_count=100000, read_cutoff=None):
        """
        Parses a FASTQ file and extracts sequences from unknown regions using regular expressions. 

        Args:
            fastq_file (str): Path to the input FASTQ file (can be gzipped).
            trim_read_start (int or False, optional): A number of
                bases to trim off the start of the read. If False, no 
                trimming is performed. Default: False.
            quality_thresh (int or False, optional):  The minimum average 
                quality score required for the flanking regions around 
                an unknown region for a read to be considered valid. If 
                False, quality filtering is skipped. Default: False.
            progress_read_count (int, optional): The number of reads processed 
                before printing a progress message. Default: 100000.
            read_cutoff (int or None, optional): If not None, limits the parsing 
                to the specified number of reads from the beginning of the 
                FASTQ file. Default: None.

        Returns:
            pd.DataFrame: A DataFrame summarizing the parsed results. Columns 
                include the unknown region names and 'Count'. Each row 
                represents a unique combination of sequences extracted from 
                the unknown regions and their observed frequency.
        """
        print('Parsing', fastq_file)
        self.create_regex_lists()
        if fastq_file.endswith(".gz"):
            opener = gzip.open
        else:
            opener = open
        bc_counter = Counter()
        c = 0
        with opener(fastq_file, 'rt') as infile:
            for title, seq, qual in FourLineFastq(infile):
                if trim_read_start:
                    seq = seq[trim_read_start:]
                    qual = qual[trim_read_start:]
                result = '_'.join([self.regex_parse(seq, qual, u, quality_thresh) for u in self.unknown_regions])
                bc_counter[result] += 1
                c += 1
                if read_cutoff and c >= read_cutoff:
                    break
                if c % progress_read_count == 0:
                    print('Parsed', c, 'reads')
                 
        self.stats = {'totalReads': c}
        mat = [b.split('_')+[bc_counter[b]] for b in bc_counter]
        names = [u['name'] for u in self.unknown_regions]
        td = pd.DataFrame(mat, columns=names+['Count']).sort_values(by='Count', ascending=False)
        for col in names:
            lens = td[col].apply(len)
            lend = dict(lens.value_counts())
            self.stats[col] = {
                'nUnique': len(set(td[col])), 
                # have to cast to int here to avoid a  json serializing bug
                'lenCounts': {int(i): lend[i] for i in lend},
                'lenReads': {int(i[0]):i[1] for i in np.array(td[[col, 'Count']].groupby(lens).sum(numeric_only=True).reset_index())},
            }
            for fail in ['RegexFail', 'QualityFail']:
                if fail in set(td[col]):
                    self.stats[col][fail] = np.sum(td[td[col]==fail]['Count'])
                else:
                    self.stats[col][fail] = 0
        return td
    
    def parse_fastq_alignment(self, fastq_file, progress_read_count=100000, read_cutoff=None):
        """
        Parses a FASTQ file and extracts sequences from unknown regions using alignment.

        This method aligns each read to the reference sequence and then
        extracts the sequences from the defined unknown regions based on
        the alignment coordinates. 

        Args:
            fastq_file (str): Path to the input FASTQ file.
            progress_read_count (int, optional):  The number of reads 
                processed before printing a progress message. 
                Default: 100000.
            read_cutoff (int or None, optional):  If not None, limits 
                the parsing to the specified number of reads from the 
                beginning of the FASTQ file. Default: None.

        Returns:
            pd.DataFrame:  A DataFrame summarizing the parsed results. 
                Columns include the unknown region names and 'Count'. 
                Each row represents a unique combination of sequences 
                extracted from the unknown regions and their observed 
                frequency.
        """
        print('Parsing', fastq_file)
        a = mp.Aligner(seq=self.refSeq)
        bc_counter = Counter()
        c = 0
        h = 0
        for title, seq, qual in mp.fastx_read(fastq_file):
            seq = reorient_circular_read(seq, self.refSeq[:200])
            hits = list(a.map(seq))
            if len(hits) > 0:
                hit = hits[0]
                h += 1
                al = cigar_to_ref_array(self.refSeqLen, hit, seq)
                result = '_'.join([''.join(al[u['start']:u['end']]) for u in self.unknown_regions])
                bc_counter[result] += 1
                
            c += 1
            if read_cutoff and c >= read_cutoff:
                break
            if c % progress_read_count == 0:
                print('Parsed', c, 'reads')
        print('Parsed', c, 'reads.', h, 'aligned.')
        self.stats = {'totalReads': c, 'alignedReads': h}
        mat = [b.split('_')+[bc_counter[b]] for b in bc_counter]
        names = [u['name'] for u in self.unknown_regions]
        td = pd.DataFrame(mat, columns=names+['Count']).sort_values(by='Count', ascending=False)
        for col in names:
            lens = td[col].apply(len)
            lend = dict(lens.value_counts())
            self.stats[col] = {
                'nUnique': len(set(td[col])), 
                # have to cast to int here to avoid a  json serializing bug
                'lenCounts': {int(i): lend[i] for i in lend},
                'lenReads': {int(i[0]):i[1] for i in np.array(td[[col, 'Count']].groupby(lens).sum(numeric_only=True).reset_index())},
            }
        return td

def parse_by_regex(construct, fastq_file, outfile='return', logfile='auto', construct_is_file=False, autodetect_barcodes=False, regex_flanking_len=8, unknown_lens=None, trim_read_start=False, quality_thresh=False, read_cutoff=None):
    """
    Parses a FASTQ file and extracts sequences from unknown regions using regular expressions.

    This function provides a higher-level interface to the 
    `UnknownRegionParser` class, simplifying the process of parsing 
    sequencing reads and extracting sequences from unknown regions
    using regular expression matching.

    Args:
        construct (str): The DNA sequence of the known construct. This 
            can include placeholders for unknown regions (see below).
        fastq_file (str): Path to the input FASTQ file.
        outfile (str, optional): The path to the output CSV file where the 
            parsed results will be saved. If set to 'return' (default), 
            the function returns the parsed DataFrame.
        logfile (str, optional):  The path to the output log file. If set 
            to 'auto' (default), the log file will be saved in the same 
            directory as outfile with a '_log.json' suffix. If None, no 
            log file is created.
        construct_is_file (bool, optional): If True, `construct` is 
            interpreted as the path to a file containing the construct 
            sequence. Default: False.
        autodetect_barcodes (bool, optional):  If True, unknown regions 
            are automatically detected as consecutive stretches of 'N', 
            'W', or 'S' characters in the `construct` sequence. 
            If False, unknown regions must be explicitly defined using 
            parentheses with a region name (e.g., "(region1:NNNNN)") 
            within the `construct` sequence. Default: False.
        regex_flanking_len (int, optional): The length of flanking 
            sequences around unknown regions to use for regular 
            expression matching. Default: 8.
        unknown_lens (str or list, optional): Controls how the lengths 
            of unknown regions are handled:
            - None (default): Use the actual lengths of the unknown 
                regions defined in `construct`.
            - 'All': Treat all unknown regions as having unknown 
                lengths during parsing. 
            - list of str: A list of unknown region names to treat as 
                having unknown lengths.
        trim_read_start (int or False, optional): A number of
            bases to trim off the start of the read. If False, no 
            trimming is performed. Default: False.
        quality_thresh (int or False, optional):  The minimum average 
            quality score required for the flanking regions around 
            an unknown region for a read to be considered valid. If 
            False, quality filtering is skipped. Default: False.
        read_cutoff (int or None, optional): If not None, limits the parsing 
            to the specified number of reads. Default: None.

    Returns:
        pd.DataFrame or None: If `outfile` is 'return', returns the parsed 
            DataFrame. Otherwise, returns None and saves the DataFrame to 
            the specified file. 
    """
    urp = UnknownRegionParser(
        construct,
        construct_is_file=construct_is_file, 
        autodetect_barcodes=autodetect_barcodes, 
        regex_flanking_len=regex_flanking_len, 
        unknown_lens=unknown_lens
        )
    result = urp.parse_fastq_regex(
        fastq_file,
        quality_thresh=quality_thresh,
        trim_read_start=trim_read_start, 
        read_cutoff=read_cutoff
        )
    if outfile == 'return':
        return result
    else:
        result.to_csv(outfile, index=False)
        if logfile == 'auto':
            print('Log:', outfile.replace('.csv', '_log.json'))
            assert '.csv' in outfile
            with open(outfile.replace('.csv', '_log.json'), 'w') as logout:
                json.dump(vars(urp), logout, cls=customJSONEncoder, indent=4)

def parse_by_alignment(construct, fastq_file, outfile='return', logfile='auto', construct_is_file=False, autodetect_barcodes=False, unknown_lens=None, read_cutoff=None):
    """
    Parses a FASTQ file and extracts sequences from unknown regions using alignment.

    This function provides a higher-level interface to the 
    `UnknownRegionParser` class, using its alignment-based parsing
    method to extract sequences from unknown regions.

    Args:
        construct (str): The DNA sequence of the known construct. This 
            can include placeholders for unknown regions (see below).
        fastq_file (str): Path to the input FASTQ file.
        outfile (str, optional): The path to the output CSV file where the 
            parsed results will be saved. If set to 'return' (default), 
            the function returns the parsed DataFrame.
        logfile (str, optional):  The path to the output log file. If set 
            to 'auto' (default), the log file will be saved in the same 
            directory as outfile with a '_log.json' suffix. If None, no 
            log file is created.
        construct_is_file (bool, optional): If True, `construct` is 
            interpreted as the path to a file containing the construct 
            sequence. Default: False.
        autodetect_barcodes (bool, optional):  If True, unknown regions 
            are automatically detected as consecutive stretches of 'N', 
            'W', or 'S' characters in the `construct` sequence. 
            If False, unknown regions must be explicitly defined using 
            parentheses with a region name (e.g., "(region1:NNNNN)") 
            within the `construct` sequence. Default: False.
        unknown_lens (str or list, optional): Controls how the lengths 
            of unknown regions are handled:
            - None (default): Use the actual lengths of the unknown 
                regions defined in `construct`.
            - 'All': Treat all unknown regions as having unknown 
                lengths during parsing. 
            - list of str: A list of unknown region names to treat as 
                having unknown lengths.
        read_cutoff (int or None, optional): If not None, limits the parsing 
            to the specified number of reads. Default: None.

    Returns:
        pd.DataFrame or None: If `outfile` is 'return', returns the parsed 
            DataFrame. Otherwise, returns None and saves the DataFrame to 
            the specified file. 
    """
    urp = UnknownRegionParser(
        construct, 
        construct_is_file=construct_is_file, 
        autodetect_barcodes=autodetect_barcodes, 
        unknown_lens=unknown_lens
        )
    result = urp.parse_fastq_alignment(fastq_file, read_cutoff=read_cutoff)
    if outfile == 'return':
        return result
    else:
        result.to_csv(outfile, index=False)
        if logfile == 'auto':
            assert '.csv' in outfile
            with open(outfile.replace('.csv', '_log.json'), 'w') as logout:
                json.dump(vars(urp), logout, cls=customJSONEncoder) 