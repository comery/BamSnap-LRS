"""Indexed access to reference FASTA sequences using pyfaidx."""

import atexit
import os
from typing import Dict, Any


_OPEN_REFERENCES: Dict[str, Any] = {}


def open_reference(fa_path: str):
    """Open and cache a pyfaidx FASTA object for the lifetime of this process.

    pyfaidx creates a .fai index automatically when one is not present. The
    FASTA is expected to be uncompressed, as documented by BamSnap-LRS.
    """
    if not fa_path:
        raise ValueError("A reference FASTA path is required.")

    path = os.path.abspath(os.path.expanduser(os.fspath(fa_path)))
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Reference FASTA file was not found: '{path}'")

    cached = _OPEN_REFERENCES.get(path)
    if cached is not None:
        return cached

    try:
        from pyfaidx import Fasta
    except ImportError as exc:
        raise RuntimeError(
            "pyfaidx is required for reference FASTA access. Install dependencies "
            "with 'pip install -r requirements.txt'."
        ) from exc

    try:
        # Fasta creates a .fai index when needed. Cache this object so it is not
        # reopened for each alignment record or each region in batch mode.
        fasta = Fasta(path, as_raw=True, sequence_always_upper=True)
    except Exception as exc:
        raise RuntimeError(
            f"Could not open/index reference FASTA '{path}' with pyfaidx: {exc}. "
            "Check that the FASTA is uncompressed and its directory is writable "
            "if a .fai index needs to be created."
        ) from exc

    _OPEN_REFERENCES[path] = fasta
    return fasta


def validate_reference_contig(fa_path: str, chrom: str) -> None:
    """Raise a clear error if the requested contig is absent from the FASTA."""
    fasta = open_reference(fa_path)
    if chrom not in fasta.keys():
        contigs = list(fasta.keys())
        preview = ", ".join(contigs[:10])
        if len(contigs) > 10:
            preview += ", ..."
        available = f" Available FASTA contigs include: {preview}." if preview else ""
        raise ValueError(
            f"Reference contig '{chrom}' was not found in FASTA '{fa_path}'. "
            "Contig names must match exactly between the alignment and reference files."
            f"{available}"
        )


def get_ref_subseq(fa_path: str, chrom: str, start: int, end: int):
    """Fetch a 0-based, half-open reference interval from the cached FASTA."""
    fasta = open_reference(fa_path)
    validate_reference_contig(fa_path, chrom)
    try:
        return str(fasta[chrom][start:end]).upper()
    except Exception as exc:
        raise ValueError(
            f"Could not fetch reference interval {chrom}:{start}-{end} from "
            f"'{fa_path}': {exc}"
        ) from exc


def close_references() -> None:
    """Close cached FASTA handles at process exit."""
    for fasta in _OPEN_REFERENCES.values():
        try:
            fasta.close()
        except Exception:
            pass
    _OPEN_REFERENCES.clear()


atexit.register(close_references)
