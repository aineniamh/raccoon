#!/usr/bin/env python3
"""Combine FASTA files and optionally harmonise headers from metadata."""
import csv
import logging
import sys
from typing import Dict, Iterable, Optional

from Bio import SeqIO


def get_field(row: Dict[str, str], field: str) -> str:
    if field in row:
        return row.get(field, "") or ""
    bom_field = f"\ufeff{field}"
    return row.get(bom_field, "") or ""


def load_metadata_map(metadata_path: str, id_field: str, delimiter: str) -> Dict[str, Dict[str, str]]:
    metadata = {}
    with open(metadata_path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        for row in reader:
            key = get_field(row, id_field)
            if not key:
                continue
            metadata[key] = row
    return metadata


def load_metadata_maps(metadata_paths: Iterable[str], id_field: str, delimiter: str) -> Dict[str, Dict[str, str]]:
    merged = {}
    for path in metadata_paths:
        merged.update(load_metadata_map(path, id_field, delimiter))
    return merged


def format_header(
    record_id: str,
    row: Dict[str, str],
    location_field: str,
    date_field: str,
    sep: str,
) -> str:
    location = get_field(row, location_field)
    date = get_field(row, date_field)
    return f"{record_id}{sep}{location}{sep}{date}"


def write_fasta_record(handle, header: str, sequence: str) -> None:
    sequence = "".join(sequence.split()).upper()
    handle.write(f">{header}\n{sequence}\n")


def main(args):
    """Combine fasta files into a single upper-case, unwrapped FASTA."""
    if not hasattr(args, "inputs"):
        raise ValueError("Expected argparse Namespace from raccoon.command.build_parser")

    try:
        from raccoon.utils import io

        inputs = args.inputs or []
        if not inputs:
            logging.error("No input FASTA files provided")
            return 1

        for path in inputs:
            if not io.validate_input_file(path, "FASTA file"):
                return 1

        metadata_map: Optional[Dict[str, Dict[str, str]]] = None
        if args.metadata:
            metadata_paths = list(args.metadata)
            for path in metadata_paths:
                if not io.validate_input_file(path, "Metadata file"):
                    return 1
            metadata_map = load_metadata_maps(metadata_paths, args.metadata_id_field, args.metadata_delimiter)

        output_path = args.output or "combined.fasta"
        if output_path == "-":
            out_handle = sys.stdout
            close_handle = False
        else:
            out_handle = open(output_path, "w")
            close_handle = True

        try:
            for path in inputs:
                for rec in SeqIO.parse(path, "fasta"):
                    header = rec.id
                    if metadata_map is not None:
                        row = metadata_map.get(rec.id)
                        if row:
                            header = format_header(
                                rec.id,
                                row,
                                args.metadata_location_field,
                                args.metadata_date_field,
                                args.header_separator,
                            )
                        else:
                            logging.warning("No metadata row found for %s", rec.id)
                    write_fasta_record(out_handle, header, str(rec.seq))
        finally:
            if close_handle:
                out_handle.close()

        logging.info("Combined %d input FASTA files", len(inputs))
        return 0
    except Exception:
        logging.exception("Combine failed")
        return 2
