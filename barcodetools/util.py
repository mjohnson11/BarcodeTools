import numpy as np
import regex
import json

def rc(s):
    # reverse complement
    return ''.join([{'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}.get(i,i) for i in s[::-1]])


def FourLineFastq(handle, decode=False):
    """
    Reads 4 lines from the file and returns 1, 2, and 4 with newlines stripped
    The only check for fastq format is that line 3 starts with '+'
    """
    # line reading function depends on decode flag
    parse_line = (lambda h: h.readline().decode()) if decode else (lambda h: h.readline())
        
    while True:
        line = parse_line(handle)     
        if not line:
            # end of file
            break       
        title = line.rstrip()
        seq = parse_line(handle).rstrip()
        jnk_line = parse_line(handle)
        if jnk_line[0] != "+":
            print(title, seq, jnk_line)
            raise ValueError("Looks like this isn’t a strictly 4-line fastq file")
        qual = parse_line(handle).rstrip()
        yield (title, seq, qual)
        

def color_dna(s):
    # simple way to print DNA with colors
    esc = "\033"
    colors = {'A': '[42m', 'T': '[41m', 'G': '[43m', 'C': '[44m'}
    return esc + esc.join([colors.get(i,'[0m')+i for i in s])+esc+'[0m'


def cigar_to_ref_array(refLen, hit, seq):
    """
    Takes a hit from mappy (minimap2) and creates a list of strings,
    of the length of the reference. Goes through the cigar string 
    and puts query bases at the appropriate reference position
    """
    alignment = [''] * refLen 
    query_pos = hit.q_st - 1
    ref_pos = hit.r_st - 1   

    for length, opcode in hit.cigar:
        op = "MIDNSHP=XB"[opcode]
        length = int(length)

        if op in 'M=X':  # Match or mismatch
            for _ in range(length):
                alignment[ref_pos] += seq[query_pos] 
                query_pos += 1
                ref_pos += 1
        elif op == 'I':  # Insertion
            for _ in range(length):
                alignment[ref_pos] += seq[query_pos] 
                query_pos += 1
        elif op == 'D':  # Deletion
            ref_pos += length
        elif op == 'S': # pre or post alignment
            query_pos += length

    return alignment

def reorient_circular_read(read, orientation_seq):
    """
    reorients a circular read by looking for 20mers from 
    an orientation sequence (which has to be at least 200 bp long)
    Note that this will cause a small deletion where the read actually started
    """
    def get_offsets(r, kmer20s):
        offsets = []
        for i in range(len(kmer20s)):
            if kmer20s[i] in r:
                offsets.append((r.index(kmer20s[i])-i*20+len(r)) % len(r))
        return offsets
    
    assert len(orientation_seq) >= 200 # orientation_seq must be at least 200 bp
    orient_kmers = [orientation_seq[i*20:(i+1)*20] for i in range(10)]
    forward_offsets = get_offsets(read, orient_kmers)
    if len(forward_offsets) >= 5:
        # rough median
        reindex = sorted(forward_offsets)[len(forward_offsets) // 2]
        return read[reindex:]+read[:reindex]
    else:
        rc_read = rc(read)
        back_offsets = get_offsets(rc_read, orient_kmers)
        if len(back_offsets) >= 5:
            # rough median
            reindex = sorted(back_offsets)[len(back_offsets) // 2]
            return rc_read[reindex:]+rc_read[:reindex]
        else:
            # no reindexing, didn't find the kmers
            return read
    
    
class customJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, regex.Pattern):
            return obj.pattern
        elif isinstance(obj, np.int64):
            return int(obj)
        else:
            return super().default(obj)
        
    # never got this to work, ended up converting keys
    # to int in the main script
    """
    def encode(self, obj):
        # Convert int64 keys to int before encoding
        if isinstance(obj, dict):
            return super().encode({int(k) if isinstance(k, np.int64) else k: v for k, v in obj.items()})
        return super().encode(obj)
    """   
              