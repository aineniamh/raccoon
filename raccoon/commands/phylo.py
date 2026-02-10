#!/usr/bin/env python3
"""Phylogenetic QC subcommand wrapper."""
import os
import logging


def main(args):
    """Run phylogenetic QC.

    The implementation defers heavy imports until called so tests can import the package without needing all scientific deps.
    """
    logging.info("Starting phylogenetic QC")
    if not hasattr(args, 'phylogeny'):
        raise ValueError("Expected argparse Namespace from raccoon.command.build_parser")

    try:
        # Lazy import
        from raccoon.utils import phylo_functions as pf
        from raccoon.utils import io
        from raccoon.utils.constants import KEY_OUTDIR, KEY_OUTFILENAME, KEY_PHYLOGENY, KEY_RUN_APOBEC3_PHYLO

        outdir = args.outdir or os.getcwd()
        
        # validate output directory
        if not io.ensure_output_directory(outdir):
            return 1
        
        # validate input files
        if not io.validate_alignment_file(args.assembly_refs):
            return 1
        
        phylogeny_file = getattr(args, 'phylogeny', None)
        if phylogeny_file and not io.validate_input_file(phylogeny_file, "Phylogeny file"):
            return 1
        
        mask_file = getattr(args, 'mask_file', None)
        if mask_file and not io.validate_input_file(mask_file, "Mask file"):
            return 1
        
        config = {}
        config[KEY_OUTDIR] = outdir
        config[KEY_OUTFILENAME] = args.phylogeny
        config[KEY_PHYLOGENY] = args.phylogeny
        config[KEY_RUN_APOBEC3_PHYLO] = args.run_apobec

        outgroup_ids = []
        if args.outgroup_ids:
            outgroup_ids = [x.strip() for x in args.outgroup_ids.split(',') if x.strip()]

        mask_file = mask_file or os.path.join(outdir, f"{args.phylogeny}.mask.csv")

        pf.check_for_snp_anomalies(args.assembly_refs, outgroup_ids, mask_file, config, args.height)
        logging.info("Phylogenetic QC finished")
        return 0
    except Exception:
        logging.exception("Phylo QC failed")
        return 2
