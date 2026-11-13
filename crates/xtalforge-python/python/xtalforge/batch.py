"""Batch helpers that preserve input order.

Results are deterministic and indexed like the input list (AGENTS.md rule 6).
"""

from . import detect_symmetry as _detect_symmetry
from . import pxrd as _pxrd
from . import read as _read


def read(paths):
    """Reads a list of structure paths."""
    return [_read(path) for path in paths]


def detect_symmetry(texts, symprec=0.001):
    """Detects symmetry for a list of structure texts."""
    return [
        _detect_symmetry(text, symprec=symprec) for text in texts
    ]


def pxrd(texts, radiation="cu-ka"):
    """Calculates PXRD for a list of structure texts."""
    return [_pxrd(text, radiation=radiation) for text in texts]
