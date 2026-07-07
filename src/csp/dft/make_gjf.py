"""Gaussian16 input (.inp/.gjf) generation, in Ono's legacy format.

Ports get_xyzR_lines / make_gjf_xyz from
legacy/ono_scripts/stepwise_optimization/make_step1.py: counterpoise dimers
use the Fragment=N atom markers with a "0 1 0 1 0 1" charge/multiplicity
line, and several pairs are chained into one input file with --Link1.

The route line is `EmpiricalDispersion=GD3` (zero damping) — confirmed
against the legacy input generators (see spec.md).
"""
from __future__ import annotations

from typing import List, Sequence

import numpy as np

from csp.structure.molecule import Molecule
from csp.structure.intralayer import dimer, to_layer_frame

DEFAULT_ROUTE = "#P TEST b3lyp/6-311G** EmpiricalDispersion=GD3 counterpoise=2"


def fragment_dimer_lines(
    symbols: Sequence[str],
    coords_1: np.ndarray,
    coords_2: np.ndarray,
    description: str,
    *,
    route: str = DEFAULT_ROUTE,
    mem: str = "24GB",
    nproc: int = 42,
) -> List[str]:
    """One counterpoise dimer job as a list of input-file lines.

    `symbols` is the per-monomer atom list (both fragments share it)."""
    lines = [
        f"%mem={mem}\n",
        f"%nproc={nproc}\n",
        f"{route}\n",
        "\n",
        f"{description}\n",
        "\n",
        "0 1 0 1 0 1\n",
    ]
    for frag_idx, coords in ((1, coords_1), (2, coords_2)):
        for sym, (x, y, z) in zip(symbols, coords):
            lines.append(f"{sym}(Fragment={frag_idx}) {x} {y} {z}\n")
    return lines


def make_step1_input(
    mol: Molecule,
    a: float,
    b: float,
    theta: float,
    *,
    A2: float = 0.0,
    route: str = DEFAULT_ROUTE,
    mem: str = "24GB",
    nproc: int = 42,
) -> str:
    """Step 1 input: T-shaped, SP-along-b and SP-along-a dimers via --Link1.

    Matches legacy make_gjf_xyz; the resulting log parses with
    parse_log.read_counterpoise_energies into [E_t, E_p1, E_p2], and
    E_intra(8) = 4*E_t + 2*E_p1 + 2*E_p2.
    """
    mol_l = to_layer_frame(mol)
    desc = f"{mol_l.name}_A2={int(A2)}_A3={round(theta, 2)}"

    blocks = []
    for kind, tag in (("t", "_t1"), ("p1", "_p1"), ("p2", "_p2")):
        c_i, c_j = dimer(mol_l, kind, a, b, theta, A2, in_layer_frame=True)
        blocks.append(fragment_dimer_lines(
            mol_l.symbols, c_i, c_j, desc + tag,
            route=route, mem=mem, nproc=nproc,
        ))

    lines = ["$ RunGauss\n"] + blocks[0]
    for block in blocks[1:]:
        lines += ["\n\n--Link1--\n"] + block
    lines += ["\n\n\n"]
    return "".join(lines)
