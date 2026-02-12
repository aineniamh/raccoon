#!/usr/bin/env python3
import sys
import os
import logging
from raccoon.utils.misc import green, cyan

import select
from Bio import SeqIO
from Bio import AlignIO
import collections
import csv
import math

from .constants import *

def sliding_window(elements, window_size):
    """Yield sliding windows over elements."""
    if len(elements) <= window_size:
        yield elements
        return
    for i in range(len(elements) - window_size + 1):
        yield elements[i:i+window_size]


def read_alignment(path, fmt='fasta'):
    """Read an alignment file using AlignIO and return a MultipleSeqAlignment."""
    with open(path, 'r') as fh:
        aln = AlignIO.read(fh, fmt)
    return aln


def find_high_N_sequences(aln, threshold):
    """Find sequences with N content greater than threshold.

    Returns a list of tuples (record_id, n_content_float).
    """
    flagged = []
    for rec in aln:
        seq = str(rec.seq).upper()
        n_count = seq.count('N')
        n_content = n_count / len(seq) if len(seq) else 0.0
        if n_content > threshold:
            flagged.append((rec.id, round(n_content, 3)))
    return flagged


def find_frame_breaking_indels(aln, genbank_path, reference_id=None):
    """Find frame-breaking indels using coding coordinates from a GenBank file.

    This function expects the alignment to contain the reference sequence with
    id matching `reference_id`. It returns a dict of sites (1-based positions)
    flagged due to frame-breaking indels in coding sequences.
    """
    if not genbank_path:
        return {}

    # load genbank
    try:
        gb = None
        for rec in SeqIO.parse(genbank_path, 'genbank'):
            if reference_id is None or rec.id == reference_id or rec.name == reference_id:
                gb = rec
                break
        if gb is None:
            raise FileNotFoundError(f"Reference {reference_id} not found in GenBank file")
    except Exception as exc:
        logging.error("Error reading genbank file: %s", exc)
        raise

    # build mapping from reference ungapped index -> alignment column index
    ref_rec = None
    for rec in aln:
        if rec.id == gb.id or rec.id == gb.name or (reference_id and rec.id == reference_id):
            ref_rec = rec
            break
    if ref_rec is None:
        raise ValueError("Reference sequence not found in alignment")

    ref_to_col = {}
    ungapped_pos = 0
    for col_idx in range(aln.get_alignment_length()):
        base = ref_rec.seq[col_idx]
        if base != '-':
            ungapped_pos += 1
            ref_to_col[ungapped_pos] = col_idx

    sites_to_mask = {}
    # iterate CDS features
    for feat in gb.features:
        if feat.type != 'CDS':
            continue
        # location using feature.location: start/end are 0-based on unaligned ref
        start = int(feat.location.start) + 1
        end = int(feat.location.end)
        # map to alignment columns
        cols = [ref_to_col.get(pos) for pos in range(start, end + 1)]
        cols = [c for c in cols if c is not None]
        if not cols:
            continue
        # for each sequence, if CDS length is not multiple of 3, flag indel sites only
        for rec in aln:
            non_gaps = sum(1 for c in cols if rec.seq[c] != '-')
            if non_gaps % 3 != 0:
                for c in cols:
                    if rec.seq[c] != '-':
                        continue
                    site = c + 1
                    if site not in sites_to_mask:
                        sites_to_mask[site] = {
                            'Name': site,
                            'Minimum': site,
                            'Maximum': site,
                            'Length': 1,
                            'present_in': [rec.id],
                            'note': set(['frame_break'])
                        }
                    else:
                        sites_to_mask[site]['present_in'].append(rec.id)
                        sites_to_mask[site]['note'].add('frame_break')
    return sites_to_mask


def analyze_alignment(alignment, n_window=2, gap_window=1, snp_window=10, snp_count=3):
    """Analyze an alignment and return helper structures used by QC checks.

    Returns a tuple (aln, aln_len, majority, unique_mutations, snps_near_n, snps_near_gap)
    """
    
    if isinstance(alignment, str):
        aln = AlignIO.read(alignment, 'fasta')
    else:
        aln = alignment
    aln_len = aln.get_alignment_length()

    bases = set(["A", "T", "G", "C"]) 
    # # get majority base per site
    # majority = [None] * aln_len
    # for i in range(aln_len):
    #     col = [str(rec.seq[i]).upper() for rec in aln if str(rec.seq[i]).upper() in bases]
    #     if col:
    #         c = collections.Counter(col)
    #         majority[i] = c.most_common(1)[0][0]

    # set up dicts keyed by sequence id and values a set of indexes
    unique_mutations = collections.defaultdict(set)
    snps_near_n = collections.defaultdict(set)
    snps_near_gap = collections.defaultdict(set)
    clustered_snps = collections.defaultdict(set)
    
    # the assumption here is that there are relatively few variable sites, 
    # so we can iterate column-wise and skip most sites quickly, then do more detailed checks only on the variable sites.
    
    snp_cols = set()
    for i in range(aln_len):
        
        col = set()
        
        for s in aln:

            #find the columns with variable sites
            base = s.seq[i].upper()
            if base in bases:
                col.add(base)
                
        if col:
            if len(col)>1:
                #do this for only the variable sites to save time
                snp_cols.add(i)
                    
                col_dict = collections.defaultdict(list) # key is base, value is list of sequence ids with that base at this position
                col_counter = collections.Counter() # key is base, value is count of sequences with that base at this position
                
                for s in aln:
                    base = s.seq[i].upper()
                    col_dict[base].append(s.id)
                    if base in bases:
                        col_counter[base]+=1
                    
                majority = col_counter.most_common(1)[0][0] # majority base at this position

                for base in col_dict:
                    if len(col_dict[base]) == 1:
                        unique_mutations[col_dict[base][0]].add(i) # if only one sequence has this base, it's a unique mutation
                
                # if the snp is within a couple bases of an N, may be an issue with coverage/ alignment
                for s in aln:
                    # if the variant itself isn't N or a gap and isn't the majority base
                    if s.seq[i] != "N" and s.seq[i] != majority and s.seq[i] != "-":
                        if "N" in s.seq[i-n_window:i+n_window+1]: # check a window around the snp for Ns 
                            snps_near_n[s.id].add(i)

                for s in aln:
                    if s.seq[i] != "N" and s.seq[i] != majority and s.seq[i] != "-":
                        if "-" in s.seq[i-gap_window:i+gap_window+1]:# check a window around the snp for gaps
                            snps_near_gap[s.id].add(i)

        # identify clustered snps (e.g. x or more unique snps within y bases)
        clustered_sites = set()
        for s in aln:
            unique = unique_mutations[s.id]
            s_unique = sorted(unique)
            for i,val in enumerate(s_unique):

                if len(s_unique) > i+snp_count-1: # if there are enough snps remaining to meet the threshold
                    
                    snps_in_window = s_unique[i:i+snp_count] # get the next snp_count snps
                    # if x snps are within y bases
                    if snps_in_window[-1] < val+snp_window: 
                        for snp in snps_in_window:
                            clustered_snps[s.id].add(snp)

    return unique_mutations, snps_near_n, snps_near_gap, clustered_snps


def run_alignment_qc(
    alignment_path,
    outdir,
    genbank_path=None,
    reference_id=None,
    n_threshold=DEFAULT_N_THRESHOLD,
    n_window=2,
    gap_window=1,
    cluster_window=DEFAULT_CLUSTER_WINDOW,
    cluster_count=DEFAULT_CLUSTER_COUNT,
    mask_clustered=True,
    mask_n_adjacent=True,
    mask_gap_adjacent=True,
    mask_frame_break=True,
):
    """Run a series of alignment QC checks and write summary outputs.

    Returns a dict summary with details.
    """
    aln = read_alignment(alignment_path)

    summary = {}
    flagged_n = find_high_N_sequences(aln, threshold=n_threshold)
    summary['high_n_sequences'] = flagged_n

    unique_mutations, snps_near_n, snps_near_gap, clustered_snps = analyze_alignment(
        aln, n_window=n_window, gap_window=gap_window, snp_window=cluster_window, snp_count=cluster_count
    )
    
    frame_sites = {}
    if mask_frame_break and genbank_path:
        try:
            frame_sites = find_frame_breaking_indels(aln, genbank_path, reference_id=reference_id)
        except Exception as exc:
            logging.warning("Frame-breaking indel check skipped due to error: %s", exc)

    def _add_site(merged_sites, pos, seq_id, note):
        if pos not in merged_sites:
            merged_sites[pos] = {
                KEY_NAME: pos,
                KEY_MINIMUM: pos,
                KEY_MAXIMUM: pos,
                KEY_LENGTH: 1,
                KEY_PRESENT_IN: [seq_id],
                KEY_NOTE: set([note])
            }
            return
        if seq_id not in merged_sites[pos][KEY_PRESENT_IN]:
            merged_sites[pos][KEY_PRESENT_IN].append(seq_id)
        merged_sites[pos][KEY_NOTE].add(note)

    # build mask-site dict
    merged = {}
    if mask_n_adjacent:
        for seq_id, sites in snps_near_n.items():
            for site in sites:
                _add_site(merged, site + 1, seq_id, NOTE_N_ADJACENT)

    if mask_gap_adjacent:
        for seq_id, sites in snps_near_gap.items():
            for site in sites:
                _add_site(merged, site + 1, seq_id, NOTE_GAP_ADJACENT)

    if mask_clustered:
        for seq_id, sites in clustered_snps.items():
            for site in sites:
                _add_site(merged, site + 1, seq_id, NOTE_CLUSTERED_SNPS)

    # merge frame-breaking sites
    for site, info in frame_sites.items():
        if site not in merged:
            merged[site] = info
        else:
            merged[site][KEY_PRESENT_IN].extend(info[KEY_PRESENT_IN])
            merged[site][KEY_NOTE].update(info[KEY_NOTE])

    # write mask report
    mask_file = os.path.join(outdir, 'mask_sites.csv')
    with open(mask_file, 'w') as fw:
        writer = csv.DictWriter(fw, lineterminator='\n', fieldnames=["Name","Minimum","Maximum","Length","present_in","note"])
        writer.writeheader()
        for site in sorted(merged):
            row = merged[site]
            new_row = dict(row)
            if len(new_row[KEY_PRESENT_IN]) > 10:
                new_row[KEY_PRESENT_IN] = 'many'
            else:
                new_row[KEY_PRESENT_IN] = ';'.join(new_row[KEY_PRESENT_IN])
            new_row[KEY_NOTE] = ';'.join(sorted(new_row[KEY_NOTE]))
            writer.writerow(new_row)

    summary[KEY_SITES_TO_MASK] = merged
    summary[KEY_MASK_FILE] = mask_file
    # determine overall alignment health
    summary[KEY_ISSUES_FOUND] = bool(flagged_n or merged)
    return summary


