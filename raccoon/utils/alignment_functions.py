#!/usr/bin/env python3
import sys
import os
import logging
try:
    from raccoon.utils.misc import green, cyan
except Exception:
    # fallback simple helpers
    def green(s):
        return s
    def cyan(s):
        return s
import select
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.AlignInfo import SummaryInfo
import collections
import csv
import math
import baltic as bt
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

plt.switch_backend('Agg') 


mpl.rcParams.update({'font.size': 18})

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}
mpl.rcParams.update(new_rc_params)
plt.rcParams['font.family'] = 'Helvetica'


def check_flag_N_content(input_fasta,exclude_file,config):
    """Write a CSV of sequences with high N content. Uses find_high_N_sequences under the hood."""
    flagged = []
    for rec in AlignIO.read(input_fasta, 'fasta'):
        seq = str(rec.seq).upper()
        n_count = seq.count('N')
        n_content = n_count / len(seq) if len(seq) else 0.0
        if n_content > 0.2:
            flagged.append((rec.description, round(n_content, 3)))

    with open(exclude_file, "w") as fw:
        writer = csv.DictWriter(fw, fieldnames=["name", "note"], delimiter=",", lineterminator="\n")
        writer.writeheader()
        for name, n_content in flagged:
            writer.writerow({"name": name, "note": f"N content is {n_content}"})

    logging.info(f"{green(f'{len(flagged)} sequences flagged as high N content (>0.2): ')} {exclude_file}")


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


def find_high_N_sequences(aln, threshold=0.2):
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


def find_clustered_snps(aln=None, window=10, min_snps=3, outgroup_ids=None, alignment_path=None):
    """Identify clustered SNPs per sequence using shared helper logic.

    Accepts either an alignment object (`aln`) or an `alignment_path` (fasta).
    Returns a dict mapping 1-based site -> info dict using constants for keys.
    """
    if outgroup_ids is None:
        outgroup_ids = []

    if aln is None and alignment_path is not None:
        aln, _, _, unique_mutations, _, _ = analyze_alignment(alignment_path)
    elif aln is not None:
        _, _, _, unique_mutations, _, _ = analyze_alignment(aln)
    else:
        raise ValueError("Either 'aln' or 'alignment_path' must be provided")

    from .constants import KEY_NAME, KEY_MINIMUM, KEY_MAXIMUM, KEY_LENGTH, KEY_PRESENT_IN, KEY_NOTE, NOTE_CLUSTERED_SNPS

    sites_to_mask = {}
    for rec_id, uniq_list in unique_mutations.items():
        if rec_id in outgroup_ids:
            continue
        uniq = sorted(uniq_list)
        for idx, val in enumerate(uniq):
            if idx + 1 < len(uniq) and uniq[idx + 1] < val + 2:
                window_sites = [val, uniq[idx + 1]]
            elif idx + 2 < len(uniq) and uniq[idx + 2] < val + 10:
                window_sites = [val, uniq[idx + 1], uniq[idx + 2]]
            else:
                continue
            for pos in window_sites:
                site = pos + 1
                if site not in sites_to_mask:
                    sites_to_mask[site] = {
                        KEY_NAME: site,
                        KEY_MINIMUM: site,
                        KEY_MAXIMUM: site,
                        KEY_LENGTH: 1,
                        KEY_PRESENT_IN: [rec_id],
                        KEY_NOTE: set([NOTE_CLUSTERED_SNPS]),
                    }
                else:
                    sites_to_mask[site][KEY_PRESENT_IN].append(rec_id)
    return sites_to_mask


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
        # for each sequence, compute length of non-gap within cols
        for rec in aln:
            non_gaps = sum(1 for c in cols if rec.seq[c] != '-')
            if non_gaps % 3 != 0:
                for c in cols:
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


def run_alignment_qc(alignment_path, outdir='.', genbank_path=None, reference_id=None, n_threshold=0.2, cluster_window=10, cluster_count=3):
    """Run a series of alignment QC checks and write summary outputs.

    Returns a dict summary with details.
    """
    aln = read_alignment(alignment_path)

    summary = {}
    flagged_n = find_high_N_sequences(aln, threshold=n_threshold)
    summary['high_n_sequences'] = flagged_n

    clustered_sites = find_clustered_snps(aln, window=cluster_window, min_snps=cluster_count)

    frame_sites = {}
    if genbank_path:
        try:
            frame_sites = find_frame_breaking_indels(aln, genbank_path, reference_id=reference_id)
        except Exception as exc:
            logging.warning("Frame-breaking indel check skipped due to error: %s", exc)

    # merge sites
    merged = clustered_sites
    for site, info in frame_sites.items():
        if site not in merged:
            merged[site] = info
        else:
            merged[site]['present_in'].extend(info['present_in'])
            merged[site]['note'].update(info['note'])

    # write mask report
    mask_file = os.path.join(outdir, 'mask_sites.csv')
    with open(mask_file, 'w') as fw:
        writer = csv.DictWriter(fw, lineterminator='\n', fieldnames=["Name","Minimum","Maximum","Length","present_in","note"])
        writer.writeheader()
        for site in sorted(merged):
            row = merged[site]
            new_row = dict(row)
            if len(new_row['present_in']) > 10:
                new_row['present_in'] = 'many'
            else:
                new_row['present_in'] = ';'.join(new_row['present_in'])
            new_row['note'] = ';'.join(sorted(new_row['note']))
            writer.writerow(new_row)

    summary['sites_to_mask'] = merged
    summary['mask_file'] = mask_file
    # determine overall alignment health
    summary['issues_found'] = bool(flagged_n or merged)
    return summary

def analyze_alignment(alignment):
    """Analyze an alignment and return helper structures used by QC checks.

    Returns a tuple (aln, aln_len, majority, unique_mutations, snps_near_n, snps_near_gap)
    """
    bases = set(["A", "T", "G", "C"]) 
    if isinstance(alignment, str):
        aln = AlignIO.read(alignment, 'fasta')
    else:
        aln = alignment
    aln_len = aln.get_alignment_length()

    majority = [None] * aln_len
    for i in range(aln_len):
        col = [str(rec.seq[i]).upper() for rec in aln if str(rec.seq[i]).upper() in bases]
        if col:
            c = collections.Counter(col)
            majority[i] = c.most_common(1)[0][0]

    unique_mutations = collections.defaultdict(list)
    snps_near_n = collections.defaultdict(list)
    snps_near_gap = collections.defaultdict(list)

    for i in range(aln_len):
        for rec in aln:
            b = str(rec.seq[i]).upper()
            cns = majority[i]
            if b in bases and cns and b != cns:
                # unique if only one sequence carries this non-major base
                col_vals = [str(r.seq[i]).upper() for r in aln]
                if col_vals.count(b) == 1:
                    unique_mutations[rec.id].append(i)

            # check near Ns/gaps for non-major variants
            if b != 'N' and b != '-' and cns and b != cns:
                window_n = ''.join([str(rec.seq[j]).upper() for rec in [rec] for j in range(max(0, i-2), min(aln_len, i+3))])
                # note: above is a small optimization; alternatively check column-wise
                if 'N' in window_n:
                    snps_near_n[rec.id].append(i)
                window_gap = ''.join([str(rec.seq[j]).upper() for rec in [rec] for j in range(max(0, i-1), min(aln_len, i+2))])
                if '-' in window_gap:
                    snps_near_gap[rec.id].append(i)

    return aln, aln_len, majority, unique_mutations, snps_near_n, snps_near_gap


def check_for_alignment_issues(alignment, outgroup_ids=None):
    """High-level check that returns sites to mask based on multiple heuristics."""
    if outgroup_ids is None:
        outgroup_ids = []

    aln, aln_len, majority, unique_mutations, snps_near_n, snps_near_gap = analyze_alignment(alignment)

    # clustered snps detection
    clustered_snps = collections.defaultdict(set)
    for rec in aln:
        uniq = sorted(unique_mutations.get(rec.id, []))
        for idx, val in enumerate(uniq):
            if idx + 1 < len(uniq) and uniq[idx + 1] < val + 2:
                clustered_snps[rec.id].add(val)
                clustered_snps[rec.id].add(uniq[idx + 1])
            if idx + 2 < len(uniq) and uniq[idx + 2] < val + 10:
                clustered_snps[rec.id].add(val)
                clustered_snps[rec.id].add(uniq[idx + 1])
                clustered_snps[rec.id].add(uniq[idx + 2])

    sites_to_mask = {}

    from .constants import KEY_NAME, KEY_MINIMUM, KEY_MAXIMUM, KEY_LENGTH, KEY_PRESENT_IN, KEY_NOTE, NOTE_CLUSTERED_SNPS, NOTE_N_ADJACENT, NOTE_GAP_ADJACENT

    for rec in aln:
        if rec.id in outgroup_ids:
            continue

        if rec.id in clustered_snps:
            sites = [i + 1 for i in sorted(clustered_snps[rec.id])]
            for site in sites:
                if site not in sites_to_mask:
                    sites_to_mask[site] = {
                        KEY_NAME: site,
                        KEY_MINIMUM: site,
                        KEY_MAXIMUM: site,
                        KEY_LENGTH: 1,
                        KEY_PRESENT_IN: [rec.id],
                        KEY_NOTE: set([NOTE_CLUSTERED_SNPS]),
                    }
                else:
                    sites_to_mask[site][KEY_PRESENT_IN].append(rec.id)

        if rec.id in snps_near_n:
            sites = [i + 1 for i in sorted(snps_near_n[rec.id])]
            for site in sites:
                if site not in sites_to_mask:
                    sites_to_mask[site] = {
                        KEY_NAME: site,
                        KEY_MINIMUM: site,
                        KEY_MAXIMUM: site,
                        KEY_LENGTH: 1,
                        KEY_PRESENT_IN: [rec.id],
                        KEY_NOTE: set([NOTE_N_ADJACENT]),
                    }
                else:
                    sites_to_mask[site][KEY_PRESENT_IN].append(rec.id)
                    sites_to_mask[site][KEY_NOTE].add(NOTE_N_ADJACENT)

        if rec.id in snps_near_gap:
            sites = [i + 1 for i in sorted(snps_near_gap[rec.id])]
            for site in sites:
                if site not in sites_to_mask:
                    sites_to_mask[site] = {
                        KEY_NAME: site,
                        KEY_MINIMUM: site,
                        KEY_MAXIMUM: site,
                        KEY_LENGTH: 1,
                        KEY_PRESENT_IN: [rec.id],
                        KEY_NOTE: set([NOTE_GAP_ADJACENT]),
                    }
                else:
                    sites_to_mask[site][KEY_PRESENT_IN].append(rec.id)
                    sites_to_mask[site][KEY_NOTE].add(NOTE_GAP_ADJACENT)

    logging.info(f"{green('Number of possibly problematic SNPs: ')}{len(sites_to_mask)}")
    return sites_to_mask

def merge_flagged_sites(sites_to_mask,branch_reversions,branch_convergence,out_report):

    for branch in branch_reversions:
        for j in branch_reversions[branch]:
            
            site = int(j[:-1])
            allele = j[-1]
            if site not in sites_to_mask:
                sites_to_mask[site] = {
                            "Name": site,
                            "Minimum": site,
                            "Maximum": site,
                            "Length": 1,
                            "present_in": [f"{allele}|{branch}"],
                            "note": {"reversion"}
                        }
            else:
                sites_to_mask[site]["present_in"].append(f"{j}|{branch}")
                sites_to_mask[site]["note"].add("reversion")
    
    for branch in branch_convergence:
        for snp in branch_convergence[branch]:
            site = int(snp[1:-1])
            if site not in sites_to_mask:
                sites_to_mask[site] = {
                            "Name": site,
                            "Minimum": site,
                            "Maximum": site,
                            "Length": 1,
                            "present_in": [f"{snp}|{branch}"],
                            "note": {"convergent_snp"}
                        }
            else:
                sites_to_mask[site]["present_in"].append(f"{snp}|{branch}")
                sites_to_mask[site]["note"].add("convergent_snp")

    with open(out_report,"w") as fw:
        writer = csv.DictWriter(fw,lineterminator="\n",fieldnames=["Name","Minimum","Maximum","Length","present_in","note"])
        writer.writeheader()

        for site in sorted(sites_to_mask):
            row = sites_to_mask[site]
            new_row = row
            if len(row["present_in"]) > 10:
                new_row["present_in"] = "many"
            else:
                new_row["present_in"] = ";".join(row["present_in"])
            new_row["note"] = ";".join(row["note"])
            writer.writerow(new_row)

