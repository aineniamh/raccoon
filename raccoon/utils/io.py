"""I/O validation and file handling utilities."""
import os
import logging


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
    if not os.path.isfile(filepath):
        logging.error(f"Alignment file '{filepath}' does not exist")
        return False
    
    if not os.access(filepath, os.R_OK):
        logging.error(f"Alignment file '{filepath}' is not readable")
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
