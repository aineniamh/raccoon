#!/usr/bin/env python3
import sys
import os
import logging
try:
    from raccoon.utils.misc import green, cyan
except Exception:
    def green(s):
        return s
    def cyan(s):
        return s
import statistics
import select
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.AlignInfo import SummaryInfo
import collections
import csv
from .constants import KEY_ASSEMBLY_REFERENCES, KEY_OUTDIR, KEY_OUTFILENAME, KEY_PHYLOGENY, KEY_RUN_APOBEC3_PHYLO
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

from raccoon.utils.reconstruction_functions import load_tree, ensure_node_label


def find_assembly_refs(cwd,assembly_refs,config):
    refs = []
    ref_ids = []
    if not assembly_refs:
        logging.info(cyan(f'Note: no assembly references supplied.\nDefaulting to installed assembly references:'))
        for record in SeqIO.parse(config[KEY_ASSEMBLY_REFERENCES],"fasta"):
            refs.append(record)
            ref_ids.append(record.id)
        for i in ref_ids:
            logging.info(f"- {i}")
    else:
        path_to_try = os.path.join(cwd,assembly_refs)
        try:
            for record in SeqIO.parse(path_to_try,"fasta"):
                refs.append(record)
                ref_ids.append(record.id)
                
        except Exception as exc:
            msg = cyan(f'Error: cannot find/parse reference fasta file at: ') + f'{path_to_try}\n' + cyan('Please check file path and format.\n')
            logging.error(msg)
            raise FileNotFoundError(msg) from exc

        logging.info(green(f'Assembly references supplied:'))
        for i in ref_ids:
            logging.info(f"- {i}")

    config[KEY_ASSEMBLY_REFERENCES] = ref_ids

    return refs

def recurse_back(node,root_node):
    path = [node.traits["label"]]
    p = node.parent
    root_name = root_node.traits["label"]
    while p.traits["label"] != root_name:
        path.append(p.traits["label"])
        p = p.parent

    path.append(root_name)
    return path

def get_path_to_root(treefile, tree_format="auto"):
    
    """
    traces back from a given tip through each parent 
    and returns a list of branches that make up the 
    phylo path from root to tip

    """
    my_tree = load_tree(treefile, tree_format=tree_format)
    root_node = ""
    for k in my_tree.Objects:
        current_node = k
        if not root_node:
            root_node = k
            ensure_node_label(root_node)

        ensure_node_label(current_node)
        if getattr(current_node, "parent", None):
            ensure_node_label(current_node.parent)
            
    
    branch_paths = {}
    
    for k in my_tree.Objects:
        current_node = k
        if current_node.branchType == 'leaf':
            ensure_node_label(current_node)
            path = recurse_back(current_node,root_node)
            branch_path = []
            for i in range(len(path)-1, 0, -1):
                try:
                    branch = f"{path[i]}_{path[i-1]}"
                    branch_path.append(branch)
                except Exception:
                    logging.debug(f"breaks {path}")
#             print(current_node.traits["label"])
#             print(branch_path)
            branch_paths[current_node.traits["label"]] = branch_path
    return branch_paths


def read_in_branch_snps(branch_snps):
    branch_snps_dict = collections.defaultdict(list)
    with open(branch_snps,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            branch = f"{row['parent']}_{row['child']}"
            branch_snps_dict[branch].append((row['site'],row['snp'],row['dimer'])) 
    return branch_snps_dict

def get_seq_at_node(state_file,nodename):
    
    """
    returns a dict keys off 1-based positions in the genome
    and the value is a list of tuples (id, base) at a given node
    allows you to look up for a given site what the base is for a
    given internal node or tip
    """
    seq = ""
    with open(f"{state_file}","r") as f:
        for l in f:

            if not l.startswith("#"):
                try:
                    node,site,state,probA,probC,probG,probT = l.rstrip("\n").split("\t")
                except Exception:
                    logging.debug(l)
                    break
                if node == nodename:
                    seq+=state

    return seq

def load_assembly_refs(assembly_refs):

    refs = {}
    if isinstance(assembly_refs, str):
        for record in SeqIO.parse(assembly_refs, "fasta"):
            refs[record.id] = record.seq
        return refs
    for record in assembly_refs:
        refs[record.id] = record.seq
    return refs


def flag_reversions(branch_paths, branch_snp_dict,state_file, refs):
    root_node = get_seq_at_node(state_file,"Node1")
    possible_reversions = []
    branch_reversions = collections.defaultdict(set)
    snp_to_branch = {}
    will_be_reverted = collections.defaultdict(set)
    for tip in branch_paths:
        path_snps = []
        for branch in branch_paths[tip]:
            snps = branch_snp_dict[branch]
            for i in snps:
                
                if i[0] in [j[0] for j in path_snps]:
                    base = int(i[0])
                    
                    allele = i[1][-1]
                    
                    branch_reversions[branch].add(f"{base}{allele}")
                    will_be_reverted[snp_to_branch[i[0]]].add(f"{base}{i[1][0]}")
                    reversion_to = []
                    
                    ref_alleles = []
                    for ref in refs:
                        var = refs[ref][base-1]
                        ref_alleles.append(f"{ref}:{var}")
                        if allele == var:
                            reversion_to.append(ref)
                        
                    root_var = root_node[base-1]
                    if allele == root_var:
                        reversion_to.append("Root")
                    
                    original_snp = [k for k in path_snps if k[0] == i[0]][0]
                    
                    possible_reversions.append({
                        "taxon":tip,
                        "site":i[0],
                        "original_snp": original_snp[1],
                        "original_branch":snp_to_branch[i[0]],
                        "reversion_branch":branch,
                        "dinucleotide_context": original_snp[2],
                        "reversion_snp": i[1],
                        "reference_alleles": ";".join(ref_alleles),
                        "root_allele":root_node[base-1],
                        "reversion_to":";".join(reversion_to)
                    })

                
                snp_to_branch[i[0]] = branch
                path_snps.append(i)
    if branch_reversions:
        logging.info(green("Reversions flagged:"))
        for i in branch_reversions:
            for j in branch_reversions[i]:
                logging.info(f"- {j} ({i})")
    return possible_reversions,branch_reversions,will_be_reverted


def flag_convergence(treefile, branch_snp_dict, tree_format="auto"):
    snp_to_branch = collections.defaultdict(set)
    branch_convergence = collections.defaultdict(set)

    my_tree = load_tree(treefile, tree_format=tree_format)

    for k in my_tree.Objects:
        current_node = k
        ensure_node_label(current_node)
        if getattr(current_node, "parent", None):
            ensure_node_label(current_node.parent)
        if k.branchType == 'leaf':
            node_name = current_node.traits.get("label")
            
        try:
            parent_name = current_node.parent.traits.get("label")
            node_name = current_node.traits.get("label")
        except Exception:
            continue
        if not parent_name or not node_name:
            continue
        
        branch = f"{parent_name}_{node_name}"

        snps = branch_snp_dict[branch]
        for i in snps:
            snp_to_branch[i].add(branch)
    convergent_snps= []
    for snp in snp_to_branch:
        if len(snp_to_branch[snp])> 1:
            report_snp = f"{snp[1][0]}{snp[0]}{snp[1][-1]}"
            convergent_snps.append(report_snp)
            for branch in snp_to_branch[snp]:
                branch_convergence[branch].add(report_snp)

    if convergent_snps:
        logging.info("Convergent snps flagged:")
        for i in convergent_snps:
            logging.info(f"- {i}")
    return branch_convergence


def make_reversion_tree_figure(outfile,branch_snps,branch_reversions,will_be_reverted,treefile,w,h):
    """Wrapper that delegates to plotting.make_reversion_tree_figure"""
    from .plotting import make_reversion_tree_figure as _plot
    return _plot(outfile, branch_snps, branch_reversions, will_be_reverted, treefile, w, h)
    # plt.show()

def make_convergence_tree_figure(outfile,branch_snps,branch_convergence,treefile,w,h):
    """Wrapper that delegates to plotting.make_convergence_tree_figure"""
    from .plotting import make_convergence_tree_figure as _plot
    return _plot(outfile, branch_snps, branch_convergence, treefile, w, h)


def run_phylo_snp_checks(assembly_references,config,h):

    state_file = os.path.join(config[KEY_OUTDIR],f"{config[KEY_PHYLOGENY]}.state")
    treefile = os.path.join(config[KEY_OUTDIR],f"{config[KEY_PHYLOGENY]}")

    branch_snps = os.path.join(config[KEY_OUTDIR],f"{config[KEY_PHYLOGENY]}.branch_snps.reconstruction.csv")
    reversion_figure_out = os.path.join(config[KEY_OUTDIR],f"{config[KEY_OUTFILENAME]}.reversions_fig")
    convergence_figure_out = os.path.join(config[KEY_OUTDIR],f"{config[KEY_OUTFILENAME]}.convergence_fig")

    refs = load_assembly_refs(assembly_references)
    

    branch_snp_dict = read_in_branch_snps(branch_snps)
    branch_paths= get_path_to_root(treefile)

    possible_reversions,branch_reversions,will_be_reverted = flag_reversions(branch_paths, branch_snp_dict,state_file, refs)  
    branch_convergence = flag_convergence(treefile, branch_snp_dict)
    make_reversion_tree_figure(reversion_figure_out,branch_snps,branch_reversions,will_be_reverted,treefile,25,h)
    make_convergence_tree_figure(convergence_figure_out,branch_snps,branch_convergence,treefile,25,h)

    return branch_reversions, branch_convergence

def check_for_snp_anomalies(assembly_references,outgroup_ids,mask_file,config,h):

    alignment = os.path.join(config[KEY_OUTDIR],config[KEY_OUTFILENAME])
    branch_reversions, branch_convergence = {},{}

    if config[KEY_RUN_APOBEC3_PHYLO]:
        branch_reversions, branch_convergence = run_phylo_snp_checks(assembly_references,config,h)

    sites_to_mask = check_for_alignment_issues(alignment,outgroup_ids)

    merge_flagged_sites(sites_to_mask,branch_reversions,branch_convergence,mask_file)

def flag_long_branches(treefile, tree_format="auto", sd_threshold=3.0):
    my_tree = load_tree(treefile, tree_format=tree_format)
    tip_heights = []
    tip_names = []
    for node in my_tree.Objects:
        if node.branchType == 'leaf':
            ensure_node_label(node)
            tip_names.append(node.traits.get("label"))
            tip_heights.append(node.height)

    if not tip_heights:
        return {}

    mean_val = statistics.mean(tip_heights)
    std_val = statistics.pstdev(tip_heights)
    if std_val == 0:
        return {}

    flags = {}
    for name, height in zip(tip_names, tip_heights):
        if height > mean_val + (sd_threshold * std_val):
            flags[name] = height
    return flags

def classify_snp_for_edit(snp, dimer):
    if snp == "G->A" and dimer == "GA":
        return "apobec3"
    if snp == "C->T" and dimer == "TC":
        return "apobec3"
    if snp == "A->G" or snp == "T->C":
        return "adar"
    return None

def _has_adr_cluster(sites, window, min_count):
    if len(sites) < min_count:
        return False
    sites = sorted(sites)
    for i in range(len(sites) - min_count + 1):
        if sites[i + min_count - 1] - sites[i] <= window:
            return True
    return False

def parse_site_from_report(report_snp):
    digits = "".join([c for c in report_snp if c.isdigit()])
    return digits or ""

def build_branch_snps(treefile, state_file, alignment, outdir, tree_format, phylogeny_base):
    from raccoon.utils import reconstruction_functions as rf

    state_differences = os.path.join(outdir, f"{phylogeny_base}.state_differences.csv")
    branch_snps_out = os.path.join(outdir, f"{phylogeny_base}.branch_snps.reconstruction.csv")

    node_states = rf.generate_reconstruction_files(alignment, state_file, state_differences)
    node_states_diff = rf.load_unambiguous_varying_sites(state_differences)
    rf.map_site_changes_to_branches(treefile, branch_snps_out, node_states, node_states_diff, tree_format=tree_format)

    return branch_snps_out

def write_phylo_flags_csv(outfile, rows):
    with open(outfile, "w") as handle:
        writer = csv.DictWriter(handle, fieldnames=["site", "mutation_type", "present_in", "mask_boolean"], lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

def run_phylo_qc(
    treefile,
    tree_format,
    outdir,
    alignment=None,
    state_file=None,
    assembly_refs=None,
    long_branch_sd=3.0,
    include_apobec=False,
    include_adar=False,
    adar_window=300,
    adar_min_count=3,
    output_csv=None,
):
    rows = []

    long_branch_flags = flag_long_branches(treefile, tree_format=tree_format, sd_threshold=long_branch_sd)
    for tip_name in sorted(long_branch_flags):
        rows.append({
            "site": "",
            "mutation_type": "long_branch",
            "present_in": tip_name,
            "mask_boolean": True,
        })

    if state_file:
        if not alignment:
            raise ValueError("Alignment is required when using an ASR state file")

        phylogeny_base = os.path.splitext(os.path.basename(treefile))[0]
        branch_snps = build_branch_snps(treefile, state_file, alignment, outdir, tree_format, phylogeny_base)
        branch_snp_dict = read_in_branch_snps(branch_snps)

        if include_apobec or include_adar:
            adar_sites_by_branch = collections.defaultdict(list)
            adar_all_by_branch = collections.defaultdict(list)
            for branch, snps in branch_snp_dict.items():
                for site, snp, dimer in snps:
                    edit_type = classify_snp_for_edit(snp, dimer)
                    if edit_type == "apobec3" and include_apobec:
                        rows.append({
                            "site": site,
                            "mutation_type": "apobec3",
                            "present_in": branch,
                            "mask_boolean": False,
                        })
                    elif edit_type == "adar" and include_adar:
                        try:
                            adar_sites_by_branch[branch].append(int(site))
                        except Exception:
                            continue
                        adar_all_by_branch[branch].append(site)

            if include_adar:
                adar_first_allowed = False
                for branch, sites in adar_sites_by_branch.items():
                    if _has_adr_cluster(sites, adar_window, adar_min_count):
                        ordered_sites = sorted(adar_all_by_branch[branch], key=lambda s: int(s))
                        for site in ordered_sites:
                            mask_value = True
                            if not adar_first_allowed:
                                mask_value = False
                                adar_first_allowed = True
                            rows.append({
                                "site": site,
                                "mutation_type": "adar",
                                "present_in": branch,
                                "mask_boolean": mask_value,
                            })

        branch_paths = get_path_to_root(treefile, tree_format=tree_format)
        refs = load_assembly_refs(assembly_refs) if assembly_refs else {}
        possible_reversions, branch_reversions, will_be_reverted = flag_reversions(
            branch_paths, branch_snp_dict, state_file, refs
        )

        for item in possible_reversions:
            rows.append({
                "site": item.get("site", ""),
                "mutation_type": "reversion",
                "present_in": item.get("taxon", ""),
                "mask_boolean": True,
            })

        branch_convergence = flag_convergence(treefile, branch_snp_dict, tree_format=tree_format)
        for branch, snps in branch_convergence.items():
            for report_snp in snps:
                rows.append({
                    "site": parse_site_from_report(report_snp),
                    "mutation_type": "convergent_snp",
                    "present_in": branch,
                    "mask_boolean": True,
                })

    output_csv = output_csv or os.path.join(outdir, f"{os.path.splitext(os.path.basename(treefile))[0]}.phylo_flags.csv")
    write_phylo_flags_csv(output_csv, rows)
    return output_csv

