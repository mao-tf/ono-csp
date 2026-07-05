#!/usr/bin/env python3
"""Plot S = a×b vs alpha from the Step 1a vdW sweep candidates.

The vdW analogue of the paper's Fig. 2(b): instead of the DFT-optimized
E_intra(8), the y-axis is the minimum VdW-contact unit-cell area, one curve
per candidate structure type from extract_init.py (a-stack / b-stack
endpoints, and the interior local minima of S that appear at larger alpha).

Usage:
  python -m csp.plot.vdw_scan \\
      --init-csv runs/pentacene/step1_init_params.csv \\
      --out runs/pentacene/vdw_scan_pentacene.png
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


_STYLE = {
    # structure_type: (color, marker, label)
    'a-stack':   ('tab:blue',   'o', 'a-stack endpoint (a minimized)'),
    'b-stack':   ('tab:red',    's', 'b-stack endpoint (b minimized)'),
    'local_min': ('tab:green',  '^', 'local minimum of S'),
}


def plot_vdw_scan(
    init_csv: str | Path,
    out_path: str | Path | None = None,
    title: str | None = None,
) -> None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    df = pd.read_csv(init_csv)
    required = {'alpha', 'S', 'structure_type'}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in {init_csv}: {missing}")

    fig, ax = plt.subplots(figsize=(7, 5))

    for stype, (color, marker, label) in _STYLE.items():
        sub = df[df['structure_type'] == stype].sort_values('alpha')
        if sub.empty:
            continue
        if stype == 'local_min':
            # possibly several minima per alpha -> scatter, no connecting line
            ax.scatter(sub['alpha'], sub['S'], color=color, marker=marker,
                       s=30, label=label, zorder=3)
        else:
            ax.plot(sub['alpha'], sub['S'], color=color, marker=marker,
                    markersize=4, linewidth=1.2, label=label)

    # star the global minimum over all candidates
    i_min = df['S'].idxmin()
    ax.scatter(df.loc[i_min, 'alpha'], df.loc[i_min, 'S'],
               marker='*', s=280, color='gold', edgecolors='k',
               linewidths=0.6, zorder=5,
               label=(f"global min: alpha={df.loc[i_min, 'alpha']:.0f}°, "
                      f"S={df.loc[i_min, 'S']:.1f} Å²"))

    ax.set_xlabel('alpha (herringbone angle, °)', fontsize=12)
    ax.set_ylabel('S = a×b (unit cell area, Å²)', fontsize=12)
    if title is None:
        title = f"Step 1a vdW scan — {Path(init_csv).parent.name}"
    ax.set_title(title)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    plt.tight_layout()

    if out_path:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_path, dpi=150, bbox_inches='tight')
        print(f"Wrote {out_path}")
    else:
        plt.show()
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Plot S=a×b vs alpha from Step 1a vdW sweep candidates"
    )
    ap.add_argument('--init-csv', required=True,
                    help='step1_init_params.csv from extract_init.py')
    ap.add_argument('--out', default=None,
                    help='Output image path (.png/.pdf); omit to display')
    ap.add_argument('--title', default=None)
    args = ap.parse_args()
    plot_vdw_scan(args.init_csv, args.out, args.title)


if __name__ == '__main__':
    main()
