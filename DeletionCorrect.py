import Levenshtein
import numpy as np
import scipy.stats as sci_stats


def get_deletion_neighborhood(stringer):
    """
    Generates a set of all single-base deletions of a given string.

    For example, for the input "ACTG", the function will return:
    {'ACTG', 'CTG', 'ATG', 'ACG', 'ACT'}.

    Args:
        stringer (str): The input string.

    Returns:
        set: A set containing the input string and all its single-base
             deletion variants.
    """
    return set([stringer] + [stringer[:x] + stringer[x+1:] for x in range(len(stringer))])


def correct_bc_errors(arr, min_counts_for_centroid=2, max_edits=3, poisson_error_rate=0.1):
    """
    Corrects errors in barcode sequences using a network of sequences
    linked by shared single base deletions.

    Args:
        arr (list or numpy.ndarray): A 2D array or list where each sublist 
            or row contains a barcode sequence and its corresponding read 
            count (e.g., [['ATCGG', 100], ['ATCGT', 2], ...]).
        min_counts_for_centroid (int, optional): The minimum number of reads 
            a barcode needs to have to be considered a potential "true" 
            barcode. Default: 2.
        max_edits (int, optional): The maximum edit distance (number of 
            substitutions, insertions, or deletions) allowed between a 
            barcode and a centroid for it to be considered an error of 
            that centroid. Default: 3.
        poisson_error_rate (float, optional): The estimated error rate of 
            the sequencing process, modeled using a Poisson distribution. 
            This parameter influences how conservative the error correction 
            is. A higher value indicates a higher expected error rate. 
            Default: 0.1. Note that this is a VERY permissive estimate of 
            the error rate, which works well in practice for error rates up
            to 10% (10% might occur in long mononucleotide runs).

    Returns:
        dict: A dictionary mapping each input barcode sequence to either: 
              - Its corrected form (a "centroid" barcode) if it was 
                successfully corrected.
              - 'excluded error' if it appears to be an error but is 
                outside the allowed edit distance.
              - 'excluded low count' if it has too few reads to be 
                considered a centroid or corrected.
    """
    
    # Sort array
    bc_counts = sorted(arr, key=lambda x: -1*x[1])
    error_correction_results = dict() # like {BC: [BC, "centroid" or BC_parent, count, dist_from_centroid]}
    deletion_dict = dict() # points from all single bp deletions to the highest count BC with that deletion
    seen_deletions = set() # set of all seen deletions
    corrector = dict()

    for bc, count in bc_counts:
        del_net = get_deletion_neighborhood(bc) # get the deletion neighborhood
        bc_parent = None                        # bc_parent will record the barcode we think this barcode is an error of
        dist_from_centroid = 0                  # dist_from_centroid is the number of edits away from a "centroid" or "true" barcode this bc is
        bc_parents = set()                      # bc_parents is the set of potential bc_parent bcs based on a deletion neighborhood overlap
        for d in del_net:
            if d in deletion_dict:
                bc_parents.add(deletion_dict[d])
        if bc_parents:
            if len(bc_parents) == 1:
                top_hit = bc_parents.pop()
            else:
                # multiple hits - could in theory do something more complicated here, 
                # but I am just taking the hit with the most reads
                top_hit = sorted(bc_parents, key=lambda b: error_correction_results[b][2])[-1]
            # the default is to be very generous and say we expect errors to occur with a rate of 0.1
            # if the probability of getting this many error counts is below 0.05, we will reject this error
            # assuming poisson probabilities
            # Note: this is slow, could definitely improve
            if sci_stats.poisson.sf(count-1, error_correction_results[top_hit][2]*poisson_error_rate) > 0.05:
                bc_parent = top_hit # the top hit is the bc we think this bc is an error off of
                dist_from_centroid = error_correction_results[top_hit][3]+1 # the distance from a centroid is 1 + the distance of this top hit
                # if we are within range based on the edit distance threshold
                # (note here that edits are defined as having a deletion-neighborhood overlap)
                if dist_from_centroid <= max_edits:
                    # now we'll search up the chain of parents to find the true/centroid bc
                    while error_correction_results[top_hit][1] != 'centroid':
                        top_hit = error_correction_results[top_hit][1]
                    # final check using Levenshtein
                    if Levenshtein.distance(top_hit, bc) <= max_edits:
                        corrector[bc] = top_hit
                    else:
                        # this is the case where this bc looks like an error, but it is out of the max edit distance threshold
                        corrector[bc] = 'excluded error'
                else:
                    # this is the case where this bc looks like an error, but it is out of the max edit distance threshold
                    corrector[bc] = 'excluded error'
        if not bc_parent:
            if count >= min_counts_for_centroid:
                bc_parent = 'centroid'
                corrector[bc] = bc
            else:
                bc_parent = 'excluded low count'
                corrector[bc] = 'excluded low count'
        if 'exclude' not in bc_parent:
            error_correction_results[bc] = [bc, bc_parent, count, dist_from_centroid]
            # adding to deletion dictionary
            for d in del_net:
                if d not in seen_deletions:
                    deletion_dict[d] = bc
            seen_deletions.update(del_net)
    
    trues = [i[2] for i in error_correction_results.values() if i[1]=='centroid']
    errors = [i[2] for i in error_correction_results.values() if i[1]!='centroid']
    print(f'Dataset: {len(bc_counts)} bcs, {np.sum([i[1] for i in bc_counts])} reads.')
    print(f'"True" barcodes: {len(trues)}, totaling {sum(trues)} reads')
    print(f'Error barcodes: {len(errors)}, totaling {sum(errors)} reads')
    return corrector
