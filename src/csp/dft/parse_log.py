"""Extract energies from Gaussian16 .log output.

`read_counterpoise_energies` ports Ono's `get_E`
(legacy/ono_scripts/stepwise_optimization/utils.py), now that the log format
is confirmed: a Counterpoise=2 job prints five "SCF Done: E(R...)" lines per
pair (supermolecule; fragment 1 and 2 in the full (ghost) basis; fragment 1
and 2 in their own basis), and the BSSE-corrected interaction energy is
E[0] - E[1] - E[2]. Legacy inputs chain several pairs with --Link1, giving
one value per group of five.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import List, Optional

HARTREE_TO_KCAL = 627.510

_SCF_DONE_RE = re.compile(r"SCF Done:\s*E\([^)]*\)\s*=\s*(-?\d+\.\d+)")


def read_scf_energy(log_path: str | Path) -> Optional[float]:
    """Return the last "SCF Done" energy (Hartree) in the log, or None."""
    energy = None
    with open(log_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            m = _SCF_DONE_RE.search(line)
            if m:
                energy = float(m.group(1))
    return energy


def read_counterpoise_energies(log_path: str | Path) -> List[float]:
    """BSSE-corrected interaction energies (kcal/mol), one per CP=2 pair.

    Incomplete trailing groups (job still running or aborted mid-pair) are
    ignored, mirroring the legacy behaviour of only consuming full groups
    of five SCF energies.
    """
    energies = []
    with open(log_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            m = _SCF_DONE_RE.search(line)
            if m:
                energies.append(float(m.group(1)) * HARTREE_TO_KCAL)
    return [
        energies[5 * i] - energies[5 * i + 1] - energies[5 * i + 2]
        for i in range(len(energies) // 5)
    ]
