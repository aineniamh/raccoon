"""I/O validation and file handling utilities."""
import os
import logging
from typing import Iterable

from Bio import SeqIO
import sys
import tempfile
from typing import Tuple, Optional


def ensure_output_directory(outdir: str) -> bool:
    """Ensure output directory exists and is writable.
    
    Creates the directory if it doesn't exist. Returns True on success,
    logs error and returns False on failure.
    """
    if not os.path.isdir(outdir):
        try:
            os.makedirs(outdir, exist_ok=True)
            logging.info(f"Created output directory: {outdir}")
        except OSError as exc:
            logging.error(f"Failed to create output directory '{outdir}': {exc}")
            return False
    
    if not os.access(outdir, os.W_OK):
        logging.error(f"Output directory '{outdir}' is not writable")
        return False
    
    return True


def validate_input_file(filepath: str, name: str = "input file") -> bool:
    """Validate that an input file exists and is readable.
    
    Returns True if file is valid, logs error and returns False otherwise.
    """
    if not filepath:
        return True  # optional files are OK
    
    if not os.path.isfile(filepath):
        logging.error(f"{name} '{filepath}' does not exist")
        return False
    
    if not os.access(filepath, os.R_OK):
        logging.error(f"{name} '{filepath}' is not readable")
        return False
    
    return True


def validate_alignment_file(filepath: str) -> bool:
    """Validate that an alignment file exists and is readable."""
    # support stdin (use '-' to indicate reading from STDIN)
    if filepath == '-' or filepath == '':
        # cannot fully validate stdin here (consuming it would prevent later reads)
        # assume caller will handle reading from stdin into a temp file
        return True

    if not os.path.isfile(filepath):
        logging.error(f"Alignment file '{filepath}' does not exist")
        return False
    
    if not os.access(filepath, os.R_OK):
        logging.error(f"Alignment file '{filepath}' is not readable")
        return False
    
    # ensure sequences in the alignment are consistent lengths (requires Biopython)
    try:
        lengths = set()
        for rec in SeqIO.parse(filepath, "fasta"):
            lengths.add(len(rec.seq))

        if not lengths:
            logging.error(f"Alignment file '{filepath}' contains no sequences")
            return False

        if len(lengths) != 1:
            logging.error(f"Alignment file '{filepath}' contains sequences of differing lengths: {sorted(lengths)}")
            return False

    except Exception as exc:  # pragma: no cover - defensive
        logging.error(f"Failed to read alignment file '{filepath}': {exc}")
        return False

    return True


def validate_genbank_file(filepath: str) -> bool:
    """Validate that a GenBank file exists and is readable (if provided)."""
    if not filepath:
        return True  # optional
    return validate_input_file(filepath, "GenBank file")


def validate_reference_file(filepath: str, name: str = "Reference file") -> bool:
    """Validate that a reference file exists and is readable (if provided)."""
    if not filepath:
        return True  # optional
    return validate_input_file(filepath, name)


def prepare_alignment_input(filepath: str) -> Tuple[str, Optional[str]]:
    """Prepare an alignment input path for downstream processing.

    If `filepath` is '-' the function reads stdin and writes it into a temporary
    FASTA file, returning the temp path and the temp path as the second element.
    For normal files it returns the original path and None.
    Caller is responsible for calling `cleanup_temp_file` on the returned
    temp path when finished.
    """
    if filepath == '-' or filepath == '':
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix='.fasta')
        try:
            data = sys.stdin.buffer.read()
            tmp.write(data)
        finally:
            tmp.close()
        return tmp.name, tmp.name
    return filepath, None


def cleanup_temp_file(temp_path: Optional[str]) -> None:
    """Remove a temporary file created by `prepare_alignment_input`.

    Silently ignores missing files and errors.
    """
    if not temp_path:
        return
    try:
        if os.path.exists(temp_path):
            os.unlink(temp_path)
    except Exception:
        logging.debug("Failed to remove temp file %s", temp_path)
