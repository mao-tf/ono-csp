"""Gaussian16 (.gjf) input file generation.

Implements the Step 1b dimer input template from spec.md ("Gaussian16 入力
テンプレート"): B3LYP-D3(BJ)/6-311G** with Counterpoise=2 BSSE correction,
one fragment per monomer.
"""
from __future__ import annotations

from typing import List, Sequence

DEFAULT_KEYWORDS = "# B3LYP empiricaldispersion=gd3bj 6-311g** counterpoise=2 nosymm"


def _atom_block(symbols: Sequence[str], coords) -> str:
    return "\n".join(
        f"{s} {x:.6f} {y:.6f} {z:.6f}" for s, (x, y, z) in zip(symbols, coords)
    )


def make_dimer_gjf(
    title: str,
    symbols_1: Sequence[str], coords_1,
    symbols_2: Sequence[str], coords_2,
    *,
    charge: int = 0,
    mult: int = 1,
    keywords: str = DEFAULT_KEYWORDS,
) -> str:
    """Build a Gaussian16 counterpoise-corrected dimer input (spec.md §Step 1b).

    `coords_1` / `coords_2`: (N, 3) arrays of Cartesian coordinates (Å),
    already placed in the lab frame (e.g. via
    `structure.molecule.place_monomer`).
    """
    block_1 = _atom_block(symbols_1, coords_1)
    block_2 = _atom_block(symbols_2, coords_2)
    return (
        f"{keywords}\n"
        f"\n"
        f"{title}\n"
        f"\n"
        f"{charge} {mult}\n"
        f"{block_1}\n"
        f"--\n"
        f"{charge} {mult}\n"
        f"{block_2}\n"
        f"\n"
    )


def make_monomer_gjf(
    title: str,
    symbols: Sequence[str], coords,
    *,
    charge: int = 0,
    mult: int = 1,
    keywords: str = DEFAULT_KEYWORDS,
) -> str:
    """Single-monomer counterpart of `make_dimer_gjf`, for the E(monomer) term."""
    block = _atom_block(symbols, coords)
    return (
        f"{keywords}\n"
        f"\n"
        f"{title}\n"
        f"\n"
        f"{charge} {mult}\n"
        f"{block}\n"
        f"\n"
    )
