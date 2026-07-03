"""Extract energies from Gaussian16 .log output.

Only the generic, well-defined "SCF Done" energy line is parsed here. Reading
the Counterpoise=2 BSSE-corrected complexation energy (spec.md §Step 1b's
`E(接触) = E(ダイマー) - 2×E(モノマー)` with BSSE) needs a real Gaussian log
to confirm the exact "Counterpoise corrected energy" line format, so that
parser is left as a stub — fill it in once the first Step 1b jobs come back
(or once Ono's legacy scripts, which presumably already do this, land in
legacy/ono_scripts/).
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

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


def read_counterpoise_energy(log_path: str | Path) -> Optional[float]:
    """Return the BSSE-corrected (Counterpoise) complexation energy — not yet implemented."""
    raise NotImplementedError(
        "read_counterpoise_energy: needs a real Counterpoise=2 log to confirm "
        "the output line format before parsing it."
    )
