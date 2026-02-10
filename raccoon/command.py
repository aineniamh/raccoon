#!/usr/bin/env python3
"""Command-line interface for raccoon toolkit."""

import sys
import os
import argparse
import logging
from raccoon import __version__, _program
from raccoon.commands import alignment as alignment_cmd, phylo as phylo_cmd


def build_parser():
    parser = argparse.ArgumentParser(prog=_program, description="raccoon toolkit")
    parser.add_argument('-v','--version', action='version', version=f'raccoon {__version__}', help='Show version and exit')
    parser.add_argument('-V','--verbose', action='count', default=0, help='Increase verbosity (can specify multiple times)')
    sub = parser.add_subparsers(dest='command')
    sub.required = True

    a = sub.add_parser('alignment', help='alignment QC')
    a.add_argument('alignment', help='Input alignment fasta file')
    a.add_argument('-t','--sequence-type', choices=['nt','aa'], default='nt', dest='sequence_type', help='Sequence type (default: nt)')
    a.add_argument('-d','--output-dir', default='.', dest='output_dir', help='Output directory')
    a.add_argument('--genbank', help='GenBank file for frame-breaking indel checks', dest='genbank', default=None)
    a.add_argument('--reference-id', help='Reference sequence ID in alignment (for GenBank features)', dest='reference_id', default=None)
    a.add_argument('--n-threshold', type=float, default=0.2, dest='n_threshold', help='N content threshold for flagging (default: 0.2)')
    a.add_argument('--cluster-window', type=int, default=10, dest='cluster_window', help='Window size for clustered SNP detection (default: 10bp)')
    a.add_argument('--cluster-count', type=int, default=3, dest='cluster_count', help='Min SNPs in window to flag as clustered (default: 3)')
    a.set_defaults(func=alignment_cmd.main)

    p = sub.add_parser('phylo', help='phylogenetic QC')
    p.add_argument('--assembly-refs', help='Assembly references fasta', dest='assembly_refs', default=None)
    p.add_argument('--phylogeny', help='Phylogeny base name (used to locate .state etc)', dest='phylogeny', required=True)
    p.add_argument('--outdir','-d', help='Output directory', dest='outdir', default='.')
    p.add_argument('--mask-file', help='Mask output file', dest='mask_file', default=None)
    p.add_argument('--outgroup-ids', help='Comma-separated outgroup ids', dest='outgroup_ids', default=None)
    p.add_argument('--height', help='Figure height', dest='height', default=None)
    p.add_argument('--run-apobec', action='store_true', dest='run_apobec', help='Run APOBEC3 phylo checks')
    p.add_argument('--run-adar', action='store_true', dest='run_adar', help='Run ADAR phylo checks')
    p.set_defaults(func=phylo_cmd.main)

    return parser


def main(sysargs=None):
    if sysargs is None:
        sysargs = sys.argv[1:]
    parser = build_parser()
    args = parser.parse_args(sysargs)

    # configure logging
    level = logging.WARNING
    if getattr(args,'verbose',0) >=2:
        level = logging.DEBUG
    elif getattr(args,'verbose',0) == 1:
        level = logging.INFO
    logging.basicConfig(level=level, format='%(levelname)s: %(message)s')

    try:
        return args.func(args)
    except Exception as e:
        logging.exception("Error running raccoon")
        return 2


if __name__ == '__main__':
    sys.exit(main())