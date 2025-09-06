import mappy as mp
import pandas as pd

from .util import rc

def format_mappy_hit(hit):
    # returns all params separated by ; except cigar
    # (includes cigar_str which has the same info)
    return f"is_primary={hit.is_primary};ctg={hit.ctg};ctg_len={hit.ctg_len};r_st={hit.r_st};r_en={hit.r_en};q_st={hit.q_st};q_en={hit.q_en};strand={hit.strand};mapq={hit.mapq};blen={hit.blen};mlen={hit.mlen};NM={hit.NM};trans_strand={hit.trans_strand};read_num={hit.read_num};cigar_str={hit.cigar_str};MD={hit.MD};cs={hit.cs}"

def mappy_hits_processor(hits):
    if len(hits) == 0:
        return {'nHits': 0}
    elif len(hits) == 1:
        result = (hits[0], 'SingleHit')
    else:
        extra_hits = '|'.join([format_mappy_hit(h) for h in hits if not h.is_primary])
        primary_hit = [h for h in hits if h.is_primary][0]
        result = (primary_hit, extra_hits)
    return {
        'nHits': len(hits),
        'extraHits': result[1],
        'primaryHit': format_mappy_hit(result[0]),
        'contig': result[0].ctg,
        'ref_start': result[0].r_st,
        'ref_end': result[0].r_en,
        'strand': result[0].strand,
        'mapq': result[0].mapq
    }

def get_align_info(seqs_to_align_info, s, col):
    if pd.isnull(s):
        return ''
    else:
        return seqs_to_align_info[s].get(col, '') 

def align_column_to_ref(infile_or_df, outfile, fasta_ref, alignment_column, alignment_method='infer'):
    
    if isinstance(infile_or_df, pd.core.frame.DataFrame):
        df = infile_or_df
    else:
        df = pd.read_csv(infile_or_df)
    seqs_to_align = set(df[pd.notnull(df[alignment_column])][alignment_column])
    if alignment_method == 'infer':
        median_len = sorted([len(i) for i in seqs_to_align])[len(seqs_to_align)//2]
        alignment_method = 'mappy' if median_len > 60 else 'exact'
     
    seqs_to_align_info = dict()   
    if alignment_method == 'mappy':
        print('Aligning with mappy (minimap2)')
        aligner = mp.Aligner(fasta_ref)
        for s in seqs_to_align:
            hits = list(aligner.map(s))
            seqs_to_align_info[s] = mappy_hits_processor(hits)
        align_cols = ['nHits', 'extraHits', 'primaryHit', 'contig', 'ref_start', 'ref_end', 'strand', 'mapq']
    else:
        print('Aligning by exact matches (short sequences (still not ideal))')
        fasta_contigs = {name: seq for name, seq, qual in mp.fastx_read(fasta_ref)}
        rc_contigs = {name: rc(fasta_contigs[name]) for name in fasta_contigs}
        for s in seqs_to_align:
            hits = []
            for n in fasta_contigs:
                refseq = fasta_contigs[n]
                while s in refseq:
                    index = refseq.index(s)
                    hits.append({
                        'contig': n,
                        'strand': 1,
                        'ref_start': index,
                        'ref_end': index+len(s),
                    })
                    refseq = refseq[index+1:]
                refseq_rc = rc_contigs[n]
                while s in refseq_rc:
                    index = refseq_rc.index(s)
                    hits.append({
                        'contig': n,
                        'strand': -1,
                        'ref_start': len(refseq_rc)-index-len(s),
                        'ref_end': len(refseq_rc)-index,
                    })
                    refseq_rc = refseq_rc[index+1:]
            if len(hits) == 1:
                seqs_to_align_info[s] = hits[0]
                seqs_to_align_info[s]['nHits'] = 1
            else:
                seqs_to_align_info[s] = {'nHits': len(hits)}
            align_cols = ['nHits', 'contig', 'ref_start', 'ref_end', 'strand']

    for col in align_cols:
        df[col] = df[alignment_column].apply(lambda s: get_align_info(seqs_to_align_info, s, col))
    
    if outfile == 'return':
        return df
    else:
        df.to_csv(outfile, index=False)
             
                