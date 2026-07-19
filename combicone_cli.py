#!/usr/bin/env python
"""Command-line interface for CombiCone: triage and certify combinatorial screens.

One command from a real screen file to a ranked / certified table:

    # Rank all unmeasured pairs of the measured singles by predicted emergence
    combicone triage screen.h5ad --condition-key condition --order 2 -o triage.csv

    # Certify every measured combination against the single-gene cone
    combicone certify screen.h5ad --condition-key condition -o certs.csv

    # Just parse a screen and report its structure (no scoring)
    combicone ingest screen.h5ad --condition-key condition

Inputs
------
An AnnData ``.h5ad`` (``screen.h5ad``), a dense expression CSV plus a metadata
CSV (``expr.csv --conditions-csv meta.csv``), or an ``.npz`` substrate with the
keys ``atoms, single_genes`` (and, for certify, ``means/conditions/ctrl`` or a
precomputed ``noise_sd``). See ``screen_ingest`` for the grammar.

Everything is training-free and deterministic. Certificates carry the two-bar
verdict (significance AND a noise-floor ratio); triage scores are rank-only.
This CLI never fabricates data: if a combination lacks the cells to estimate
measurement noise, its certificate reports that honestly rather than guessing.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import combicone as cc
import screen_ingest as si


# --------------------------------------------------------------------------- #
# Loading
# --------------------------------------------------------------------------- #
def _load_substrate(args) -> si.ScreenSubstrate:
    """Dispatch on the input path/extension to a ScreenSubstrate."""
    src = args.screen
    suffix = Path(src).suffix.lower()
    if suffix == ".npz":
        return _substrate_from_npz(src, args)
    if args.conditions_csv:
        return si.ingest_screen(
            src,
            conditions_csv=args.conditions_csv,
            condition_key=args.condition_key,
            control_label=args.control_label,
            separator=args.separator,
            arm_handling=args.arm_handling,
            compute_noise=not args.no_noise,
            noise_seed=args.noise_seed,
            min_cells_per_half=args.min_cells_per_half,
        )
    # default: AnnData
    return si.ingest_screen(
        src,
        condition_key=args.condition_key,
        gene_name_key=args.gene_name_key,
        layer=args.layer,
        control_label=args.control_label,
        separator=args.separator,
        arm_handling=args.arm_handling,
        compute_noise=not args.no_noise,
        noise_seed=args.noise_seed,
        min_cells_per_half=args.min_cells_per_half,
    )


def _substrate_from_npz(path: str, args) -> si.ScreenSubstrate:
    """Build a ScreenSubstrate from a precomputed .npz (e.g. combicone_substrate.npz).

    Expects at least ``atoms`` (n_atoms, n_genes) and ``single_genes`` (n_atoms,).
    Measured combinations come from ``doubles`` + ``means`` + ``conditions`` +
    ``ctrl`` when present; otherwise only triage is available.
    """
    z = np.load(path, allow_pickle=True)
    keys = set(z.files)
    if not {"atoms", "single_genes"} <= keys:
        raise SystemExit(
            f"{path}: an .npz substrate needs at least 'atoms' and 'single_genes' "
            f"(found {sorted(keys)})"
        )
    atoms = np.asarray(z["atoms"], dtype=float)
    atom_names = np.asarray([str(g) for g in z["single_genes"]], dtype=str)
    genes = np.asarray([str(g) for g in z["genes"]], dtype=str) if "genes" in keys else np.asarray(
        [f"gene_{i}" for i in range(atoms.shape[1])], dtype=str
    )
    combos: list[si.ComboRecord] = []
    if {"doubles", "means", "conditions", "ctrl"} <= keys:
        cond_list = [str(c) for c in z["conditions"]]
        means = np.asarray(z["means"], dtype=float)
        ctrl = np.asarray(z["ctrl"], dtype=float)
        atom_set = set(atom_names.tolist())
        # Split-half SE = |m1 - m2| / 2 when the substrate ships the two halves.
        has_halves = {"means1", "means2"} <= keys
        if has_halves:
            m1 = np.asarray(z["means1"], dtype=float)
            m2 = np.asarray(z["means2"], dtype=float)
        for d in z["doubles"]:
            d = str(d)
            g = tuple(sorted(si.parse_condition(d, control_label=args.control_label, separator=args.separator)))
            if d not in cond_list:
                continue
            row = cond_list.index(d)
            eff = means[row] - ctrl
            noise_sd = np.abs(m1[row] - m2[row]) / 2.0 if has_halves else None
            constituents = tuple(x for x in g if x in atom_set)
            combos.append(
                si.ComboRecord(
                    name=d, genes=g, effect=eff, noise_sd=noise_sd,
                    constituent_singles=constituents,
                    has_all_singles=(len(constituents) == len(g)),
                )
            )
    return si.ScreenSubstrate(
        atoms=atoms, atom_names=atom_names, genes=genes, combos=combos,
        control_label=args.control_label, separator=args.separator,
        provenance={"source": path, "from_npz": True},
    )


# --------------------------------------------------------------------------- #
# Subcommands
# --------------------------------------------------------------------------- #
def cmd_ingest(args) -> int:
    sub = _load_substrate(args)
    summary = sub.summary()
    print(json.dumps(summary, indent=2, default=str))
    if args.output:
        Path(args.output).write_text(json.dumps(summary, indent=2, default=str))
        print(f"\nwrote {args.output}", file=sys.stderr)
    return 0


def cmd_triage(args) -> int:
    sub = _load_substrate(args)
    atoms, names = sub.triage_ready()
    candidates = None
    if args.measured_only:
        candidates = [c.genes for c in sub.combos if c.has_all_singles]
        if not candidates:
            raise SystemExit("--measured-only given but no measured combos with all singles")
    tr = cc.triage_combinations(
        atoms, names, candidates,
        order=args.order, pairwise=args.pairwise, use_gap=args.use_gap,
    )
    order = np.argsort(-tr.score, kind="stable")
    rows = []
    for r, i in enumerate(order, start=1):
        rows.append({
            "rank": r,
            "combination": args.separator.join(tr.pairs[i]),
            "score": float(tr.score[i]),
            "agg_cosine": float(tr.cos_ab[i]),
            "gap": float(tr.gap[i]) if tr.gap is not None else float("nan"),
        })
    _emit_table(rows, args.output, ["rank", "combination", "score", "agg_cosine", "gap"])
    print(f"\nscope: {tr.scope}", file=sys.stderr)
    if args.top:
        print(f"\nTop {args.top} to run first:", file=sys.stderr)
        for row in rows[: args.top]:
            print(f"  {row['rank']:>3}. {row['combination']:20s} score={row['score']:+.4f}", file=sys.stderr)
    return 0


def cmd_certify(args) -> int:
    sub = _load_substrate(args)
    targets = sub.combos
    if args.combo:
        wanted = set(args.combo)
        targets = [c for c in sub.combos if c.name in wanted]
        if not targets:
            raise SystemExit(f"none of {args.combo} are measured combinations")
    if not targets:
        raise SystemExit("no measured combinations to certify in this screen")
    rows = []
    for rec in targets:
        if rec.noise_sd is None or not np.all(np.isfinite(rec.noise_sd)):
            rows.append({
                "combination": rec.name,
                "verdict": "insufficient cells for noise estimate",
                "unreachable_fraction": float("nan"), "z": float("nan"),
                "floor_ratio": float("nan"), "p_value": float("nan"),
                "geometry_status": "not_tested",
            })
            continue
        cert = cc.certify_emergence(
            cone_atoms=sub.atoms, measured_combo=rec.effect, noise_sd=rec.noise_sd,
            n_boot=args.n_boot, floor_threshold=args.floor_threshold,
            alpha=args.alpha, seed=args.seed,
        )
        rows.append({
            "combination": rec.name,
            "verdict": cert.verdict,
            "unreachable_fraction": float(cert.unreachable_fraction),
            "z": float(cert.z),
            "floor_ratio": float(cert.floor_ratio),
            "p_value": float(cert.p_value),
            "geometry_status": cert.geometry_status,
        })
    rows.sort(key=lambda r: (-(r["z"] if np.isfinite(r["z"]) else -1e9)))
    _emit_table(
        rows, args.output,
        ["combination", "verdict", "unreachable_fraction", "z", "floor_ratio", "p_value", "geometry_status"],
    )
    n_cert = sum(1 for r in rows if r["verdict"].startswith("certified"))
    print(f"\n{n_cert}/{len(rows)} certified emergent (both bars)", file=sys.stderr)
    return 0


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #
def _emit_table(rows: list[dict], output: str | None, columns: list[str]) -> None:
    import csv
    if output:
        with open(output, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=columns)
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {output} ({len(rows)} rows)", file=sys.stderr)
    else:
        # pretty print to stdout
        widths = {c: max(len(c), *(len(f"{r[c]:.4g}" if isinstance(r[c], float) else str(r[c])) for r in rows)) for c in columns} if rows else {c: len(c) for c in columns}
        print("  ".join(c.ljust(widths[c]) for c in columns))
        for r in rows:
            print("  ".join(
                (f"{r[c]:.4g}" if isinstance(r[c], float) else str(r[c])).ljust(widths[c])
                for c in columns
            ))


# --------------------------------------------------------------------------- #
# Arg parsing
# --------------------------------------------------------------------------- #
def _add_common(p: argparse.ArgumentParser) -> None:
    p.add_argument("screen", help="path to .h5ad, expression .csv, or .npz substrate")
    p.add_argument("--conditions-csv", default=None, help="metadata CSV (when screen is an expression CSV)")
    p.add_argument("--condition-key", default="condition", help="obs/column with condition labels")
    p.add_argument("--gene-name-key", default=None, help="var column for gene names (AnnData)")
    p.add_argument("--layer", default=None, help="AnnData layer to use instead of X")
    p.add_argument("--control-label", default="ctrl", help="control token (default: ctrl)")
    p.add_argument("--separator", default="+", help="condition-label separator (default: +)")
    p.add_argument("--arm-handling", default="merge", choices=["merge", "control_left", "control_right"],
                   help="how to pool a single spelled on both control arms")
    p.add_argument("-o", "--output", default=None, help="write CSV/JSON here (else stdout)")


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="combicone", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)

    pi = sub.add_parser("ingest", help="parse a screen and report its structure")
    _add_common(pi)
    pi.add_argument("--no-noise", action="store_true", help="skip split-half noise estimation")
    pi.add_argument("--noise-seed", type=int, default=7)
    pi.add_argument("--min-cells-per-half", type=int, default=6)
    pi.set_defaults(func=cmd_ingest)

    pt = sub.add_parser("triage", help="rank unmeasured combinations by predicted emergence")
    _add_common(pt)
    pt.add_argument("--order", type=int, default=2, help="combination order k (2=pairs, 3=triples)")
    pt.add_argument("--pairwise", default="mean", choices=["mean", "max"])
    pt.add_argument("--use-gap", action="store_true", help="add leave-combo-out gap feature (slower)")
    pt.add_argument("--measured-only", action="store_true", help="only score combos actually measured in the screen")
    pt.add_argument("--top", type=int, default=10, help="print the top-N to run first")
    pt.add_argument("--no-noise", action="store_true")
    pt.add_argument("--noise-seed", type=int, default=7)
    pt.add_argument("--min-cells-per-half", type=int, default=6)
    pt.set_defaults(func=cmd_triage)

    pc = sub.add_parser("certify", help="certify measured combinations against the singles cone")
    _add_common(pc)
    pc.add_argument("--combo", nargs="+", default=None, help="only certify these named combinations")
    pc.add_argument("--n-boot", type=int, default=200)
    pc.add_argument("--floor-threshold", type=float, default=1.9)
    pc.add_argument("--alpha", type=float, default=0.05)
    pc.add_argument("--seed", type=int, default=0)
    pc.add_argument("--no-noise", action="store_true")
    pc.add_argument("--noise-seed", type=int, default=7)
    pc.add_argument("--min-cells-per-half", type=int, default=6)
    pc.set_defaults(func=cmd_certify)

    return ap


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
