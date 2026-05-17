#!/usr/bin/env python3
"""
=============================================================================
MMPBSA Decomposition Analyser — Wild Type vs Mutant Comparison
=============================================================================
Usage (single run):
    python analyze_mmpbsa.py \
        --wt-decomp   APC1/WT/FINAL_DECOMP_MMPBSA.dat \
        --wt-results  APC1/WT/FINAL_RESULTS_MMPBSA.dat \
        --mut-decomp  APC1/MUT/FINAL_DECOMP_MMPBSA.dat \
        --mut-results APC1/MUT/FINAL_RESULTS_MMPBSA.dat \
        --name APC1 \
        --out  results/

Usage (batch — all mutations at once):
    python analyze_mmpbsa.py --batch batch_list.csv --out results/

    batch_list.csv columns:
        name, wt_decomp, wt_results, mut_decomp, mut_results
=============================================================================
"""

import argparse
import sys
import warnings
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

warnings.filterwarnings("ignore")

# ── Colour palette ────────────────────────────────────────────────────────────
C_WT      = "#2196F3"   # blue
C_MUT     = "#F44336"   # red
C_FAV     = "#4CAF50"   # green  (stabilising)
C_UNFAV   = "#FF5722"   # orange (destabilising)
C_NEUTRAL = "#9E9E9E"   # grey

# =============================================================================
# 1.  PARSE FINAL_RESULTS_MMPBSA.dat  — true total binding energy
# =============================================================================

def parse_results(filepath):
    """
    Read FINAL_RESULTS_MMPBSA.dat and extract the DELTA total binding energy
    and its SEM from the Delta section.

    Returns a dict with keys:
        delta_total, delta_total_sem,
        delta_vdw, delta_elec, delta_polar, delta_nonpolar
    """
    result   = {}
    in_delta = False

    # Key mappings — handle both Unicode delta and plain ASCII
    key_map = {
        "VDWAALS" : ("delta_vdw",      "delta_vdw_sem"),
        "EEL"     : ("delta_elec",     "delta_elec_sem"),
        "EGB"     : ("delta_polar",    "delta_polar_sem"),
        "ESURF"   : ("delta_nonpolar", "delta_nonpolar_sem"),
        "TOTAL"   : ("delta_total",    "delta_total_sem"),
    }

    with open(filepath) as fh:
        for line in fh:
            # Enter delta section
            if re.search(r"Delta\s*\(Complex", line, re.IGNORECASE):
                in_delta = True
                continue

            if not in_delta:
                continue

            # Exit after the closing dashes
            if re.match(r"^-{10,}", line) and result:
                break

            # Strip any leading Greek delta character or whitespace
            clean = line.strip().lstrip("Δ").lstrip("δ").strip()

            # Match:  KEYWORD   avg_value   sd   sd   sem_value
            m = re.match(
                r"^(VDWAALS|EEL|EGB|ESURF|TOTAL)\s+"
                r"([-\d.]+)\s+\S+\s+\S+\s+([-\d.]+)",
                clean
            )
            if m:
                key = m.group(1)
                avg = float(m.group(2))
                sem = float(m.group(3))
                if key in key_map:
                    result[key_map[key][0]] = avg
                    result[key_map[key][1]] = sem

    if "delta_total" not in result:
        sys.exit(
            f"\n[ERROR] Could not parse DELTA TOTAL from:\n  {filepath}\n"
            "Make sure this is a valid FINAL_RESULTS_MMPBSA.dat file and "
            "that it contains a 'Delta (Complex - Receptor - Ligand)' section."
        )

    return result


# =============================================================================
# 2.  PARSE FINAL_DECOMP_MMPBSA.dat  — per-residue breakdown
# =============================================================================

def parse_decomp(filepath):
    """
    Parse FINAL_DECOMP_MMPBSA.dat into a DataFrame.
    Only reads lines from the DELTAS section.
    """
    records  = []
    in_delta = False

    with open(filepath) as fh:
        for line in fh:
            stripped = line.strip()

            if stripped.startswith("DELTAS:"):
                in_delta = True
                continue

            if not in_delta:
                continue

            if not (stripped.startswith("R:") or stripped.startswith("L:")):
                continue

            parts = stripped.split(",")
            if len(parts) < 19:
                continue

            try:
                res_id     = parts[0]
                seg        = res_id.split(":")
                chain_type = seg[0]
                chain      = seg[1] if len(seg) > 1 else "?"
                resname    = seg[2] if len(seg) > 2 else "?"
                resnum     = int(seg[3]) if len(seg) > 3 else 0

                # Column layout (0-indexed):
                # 0:res  1:int_avg 2:int_sd 3:int_sem
                # 4:vdw_avg 5:vdw_sd 6:vdw_sem
                # 7:elec_avg 8:elec_sd 9:elec_sem
                # 10:polar_avg 11:polar_sd 12:polar_sem
                # 13:nonpolar_avg 14:nonpolar_sd 15:nonpolar_sem
                # 16:total_avg 17:total_sd 18:total_sem

                records.append({
                    "residue"     : res_id,
                    "chain_type"  : chain_type,
                    "chain"       : chain,
                    "resname"     : resname,
                    "resnum"      : resnum,
                    "vdw"         : float(parts[4]),
                    "vdw_sem"     : float(parts[6]),
                    "elec"        : float(parts[7]),
                    "elec_sem"    : float(parts[9]),
                    "polar"       : float(parts[10]),
                    "polar_sem"   : float(parts[12]),
                    "nonpolar"    : float(parts[13]),
                    "nonpolar_sem": float(parts[15]),
                    "total"       : float(parts[16]),
                    "total_sem"   : float(parts[18]),
                })
            except (ValueError, IndexError):
                continue

    if not records:
        sys.exit(
            f"\n[ERROR] No residue data found in:\n  {filepath}\n"
            "Check the file is a valid FINAL_DECOMP_MMPBSA.dat."
        )

    return pd.DataFrame(records)


# =============================================================================
# 3.  COMPARISON
# =============================================================================

def compare(wt_df, mut_df):
    """
    Merge WT and mutant decomp data, compute per-residue ΔΔG.
    ΔΔG = mutant − WT  (negative = stabilising, positive = destabilising)
    """
    merged = pd.merge(
        wt_df, mut_df,
        on=["residue", "chain_type", "chain", "resname", "resnum"],
        suffixes=("_wt", "_mut")
    )

    merged["delta_total"]    = merged["total_mut"]    - merged["total_wt"]
    merged["delta_vdw"]      = merged["vdw_mut"]      - merged["vdw_wt"]
    merged["delta_elec"]     = merged["elec_mut"]     - merged["elec_wt"]
    merged["delta_polar"]    = merged["polar_mut"]    - merged["polar_wt"]
    merged["delta_nonpolar"] = merged["nonpolar_mut"] - merged["nonpolar_wt"]

    # Propagated SEM
    merged["delta_sem"] = np.sqrt(
        merged["total_sem_wt"]**2 + merged["total_sem_mut"]**2
    )

    # Significant: |ΔΔG| > 2 × propagated SEM
    merged["significant"] = (
        np.abs(merged["delta_total"]) > 2 * merged["delta_sem"]
    )

    return merged.sort_values("delta_total")


# =============================================================================
# 4.  SUMMARY
# =============================================================================

def summarise(merged, wt_res, mut_res, name):
    """
    Print full summary using true ΔG from FINAL_RESULTS_MMPBSA.dat
    and per-residue breakdown from FINAL_DECOMP_MMPBSA.dat.
    """
    rec = merged[merged["chain_type"] == "R"]
    lig = merged[merged["chain_type"] == "L"]

    # True total ΔΔG from FINAL_RESULTS_MMPBSA.dat
    true_ddg     = mut_res["delta_total"]     - wt_res["delta_total"]
    true_ddg_sem = np.sqrt(
        wt_res["delta_total_sem"]**2 + mut_res["delta_total_sem"]**2
    )

    # Component ΔΔGs
    ddg_vdw      = mut_res.get("delta_vdw",      0) - wt_res.get("delta_vdw",      0)
    ddg_elec     = mut_res.get("delta_elec",     0) - wt_res.get("delta_elec",     0)
    ddg_polar    = mut_res.get("delta_polar",    0) - wt_res.get("delta_polar",    0)
    ddg_nonpolar = mut_res.get("delta_nonpolar", 0) - wt_res.get("delta_nonpolar", 0)

    # Interface contribution from decomp (partial — within 6 Å only)
    interface_ddg     = merged["delta_total"].sum()
    interface_rec_ddg = rec["delta_total"].sum()
    interface_lig_ddg = lig["delta_total"].sum()

    n_sig    = merged["significant"].sum()
    n_stab   = ((merged["delta_total"] < -1) & merged["significant"]).sum()
    n_destab = ((merged["delta_total"] >  1) & merged["significant"]).sum()

    top_stab   = merged.nsmallest(5, "delta_total")[
        ["residue", "delta_total", "delta_vdw", "delta_elec", "significant"]
    ]
    top_destab = merged.nlargest(5, "delta_total")[
        ["residue", "delta_total", "delta_vdw", "delta_elec", "significant"]
    ]

    print(f"\n{'='*65}")
    print(f"  MUTATION: {name}")
    print(f"{'='*65}")
    print(f"\n  ── True Binding Energies (from FINAL_RESULTS_MMPBSA.dat) ──")
    print(f"  WT  total ΔG          : {wt_res['delta_total']:+.2f} kcal/mol "
          f"(±{wt_res['delta_total_sem']:.2f})")
    print(f"  MUT total ΔG          : {mut_res['delta_total']:+.2f} kcal/mol "
          f"(±{mut_res['delta_total_sem']:.2f})")
    print(f"  ΔΔG  (MUT − WT)       : {true_ddg:+.2f} kcal/mol "
          f"(±{true_ddg_sem:.2f})")
    print(f"    ΔΔG VdW             : {ddg_vdw:+.2f} kcal/mol")
    print(f"    ΔΔG Electrostatic   : {ddg_elec:+.2f} kcal/mol")
    print(f"    ΔΔG Polar solvation : {ddg_polar:+.2f} kcal/mol")
    print(f"    ΔΔG Non-polar solv. : {ddg_nonpolar:+.2f} kcal/mol")
    stab = "STABILISING" if true_ddg < 0 else "DESTABILISING"
    print(f"\n  → Mutation is {stab} ({true_ddg:+.2f} kcal/mol)")

    print(f"\n  ── Interface Residues (from FINAL_DECOMP_MMPBSA.dat) ──")
    print(f"  Interface ΔΔG sum     : {interface_ddg:+.2f} kcal/mol  "
          f"(within 6 Å only — not the full system)")
    print(f"    Protein contribution: {interface_rec_ddg:+.2f} kcal/mol")
    print(f"    DNA contribution    : {interface_lig_ddg:+.2f} kcal/mol")
    print(f"  Significant residues  : {n_sig}")
    print(f"    Stabilising  (< -1) : {n_stab}")
    print(f"    Destabilising (> +1): {n_destab}")
    print(f"\n  Top 5 stabilising residues:")
    print(top_stab.to_string(index=False))
    print(f"\n  Top 5 destabilising residues:")
    print(top_destab.to_string(index=False))

    return {
        "mutation"          : name,
        "wt_dg"             : round(wt_res["delta_total"], 3),
        "wt_dg_sem"         : round(wt_res["delta_total_sem"], 3),
        "mut_dg"            : round(mut_res["delta_total"], 3),
        "mut_dg_sem"        : round(mut_res["delta_total_sem"], 3),
        "true_ddg"          : round(true_ddg, 3),
        "true_ddg_sem"      : round(true_ddg_sem, 3),
        "ddg_vdw"           : round(ddg_vdw, 3),
        "ddg_elec"          : round(ddg_elec, 3),
        "ddg_polar"         : round(ddg_polar, 3),
        "ddg_nonpolar"      : round(ddg_nonpolar, 3),
        "interface_ddg"     : round(interface_ddg, 3),
        "interface_prot_ddg": round(interface_rec_ddg, 3),
        "interface_dna_ddg" : round(interface_lig_ddg, 3),
        "n_significant"     : int(n_sig),
        "n_stabilising"     : int(n_stab),
        "n_destabilising"   : int(n_destab),
    }


# =============================================================================
# 5.  PLOTS
# =============================================================================

def plot_per_residue(merged, name, outdir):
    """Bar plot of ΔΔG per residue, coloured by significance."""
    fig, axes = plt.subplots(2, 1, figsize=(16, 10), sharex=False)

    for ax, mol, label in zip(
        axes, ["R", "L"], ["Protein (Receptor)", "DNA (Ligand)"]
    ):
        sub = merged[merged["chain_type"] == mol].copy()
        sub = sub.sort_values(["chain", "resnum"])
        sub["label"] = sub["resname"] + sub["resnum"].astype(str)

        colours = [
            C_FAV    if (d < -1 and s) else
            C_UNFAV  if (d >  1 and s) else
            C_NEUTRAL
            for d, s in zip(sub["delta_total"], sub["significant"])
        ]

        ax.bar(range(len(sub)), sub["delta_total"],
               color=colours, edgecolor="none", alpha=0.85)
        ax.errorbar(range(len(sub)), sub["delta_total"],
                    yerr=sub["delta_sem"], fmt="none",
                    ecolor="black", elinewidth=0.5, capsize=2)
        ax.axhline(0,  color="black",  linewidth=0.8)
        ax.axhline( 1, color=C_UNFAV, linewidth=0.6, linestyle="--", alpha=0.5)
        ax.axhline(-1, color=C_FAV,   linewidth=0.6, linestyle="--", alpha=0.5)

        for i, row in sub[sub["significant"]].iterrows():
            idx = sub.index.get_loc(i)
            ax.text(idx, row["delta_total"],
                    f"  {row['label']}", fontsize=6,
                    va="bottom" if row["delta_total"] > 0 else "top",
                    rotation=90)

        ax.set_title(f"{label} — Per-Residue ΔΔG  [{name}]", fontsize=11)
        ax.set_ylabel("ΔΔG (kcal/mol)", fontsize=9)
        ax.set_xticks(range(len(sub)))
        ax.set_xticklabels(sub["label"], rotation=90, fontsize=5)

    patches = [
        mpatches.Patch(color=C_FAV,    label="Stabilising  (< −1, significant)"),
        mpatches.Patch(color=C_UNFAV,  label="Destabilising (> +1, significant)"),
        mpatches.Patch(color=C_NEUTRAL,label="Not significant"),
    ]
    fig.legend(handles=patches, loc="upper right", fontsize=8, framealpha=0.8)
    fig.tight_layout()
    path = outdir / f"{name}_per_residue.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  [saved] {path}")


def plot_energy_breakdown(merged, name, outdir):
    """Stacked bar of ΔΔG components for top significant residues."""
    sig = merged[merged["significant"]].copy()
    sig["abs_delta"] = np.abs(sig["delta_total"])
    sig = sig.nlargest(20, "abs_delta").sort_values("delta_total")

    if sig.empty:
        print(f"  [skip] No significant residues to plot for {name}")
        return

    sig["label"] = sig["resname"] + sig["resnum"].astype(str)
    x = np.arange(len(sig))
    w = 0.6

    fig, ax = plt.subplots(figsize=(14, 6))
    ax.bar(x, sig["delta_vdw"],      w, label="ΔvdW",            color="#1565C0")
    ax.bar(x, sig["delta_elec"],     w, label="ΔElec",            color="#E53935",
           bottom=sig["delta_vdw"])
    ax.bar(x, sig["delta_polar"],    w, label="ΔPolar solv.",     color="#43A047",
           bottom=sig["delta_vdw"] + sig["delta_elec"])
    ax.bar(x, sig["delta_nonpolar"], w, label="ΔNon-polar solv.", color="#FB8C00",
           bottom=sig["delta_vdw"] + sig["delta_elec"] + sig["delta_polar"])

    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(sig["label"], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("ΔΔG component (kcal/mol)", fontsize=9)
    ax.set_title(
        f"Energy Component Breakdown — Top Significant Residues [{name}]",
        fontsize=11
    )
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    path = outdir / f"{name}_energy_breakdown.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  [saved] {path}")


def plot_wt_vs_mut_scatter(merged, name, outdir):
    """Scatter: WT total vs Mutant total per residue."""
    fig, ax = plt.subplots(figsize=(8, 8))

    rec = merged[merged["chain_type"] == "R"]
    lig = merged[merged["chain_type"] == "L"]

    ax.scatter(rec["total_wt"], rec["total_mut"],
               alpha=0.6, s=30, color=C_WT, label="Protein residues")
    ax.scatter(lig["total_wt"], lig["total_mut"],
               alpha=0.6, s=30, color=C_MUT, label="DNA residues")

    for _, row in merged[merged["significant"]].iterrows():
        ax.annotate(f"{row['resname']}{row['resnum']}",
                    (row["total_wt"], row["total_mut"]),
                    fontsize=6, alpha=0.8)

    lims = [
        min(merged["total_wt"].min(), merged["total_mut"].min()) - 1,
        max(merged["total_wt"].max(), merged["total_mut"].max()) + 1,
    ]
    ax.plot(lims, lims, "k--", linewidth=0.8, alpha=0.5, label="No change")
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_xlabel("WT residue energy (kcal/mol)", fontsize=9)
    ax.set_ylabel("Mutant residue energy (kcal/mol)", fontsize=9)
    ax.set_title(f"WT vs Mutant Per-Residue Energies [{name}]", fontsize=11)
    ax.legend(fontsize=8)
    fig.tight_layout()
    path = outdir / f"{name}_wt_vs_mut.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  [saved] {path}")


def plot_batch_summary(summaries, outdir):
    """Ranked bar chart of all mutations by true ΔΔG."""
    df = pd.DataFrame(summaries).sort_values("true_ddg")
    colours = [C_FAV if v < 0 else C_UNFAV for v in df["true_ddg"]]

    fig, ax = plt.subplots(figsize=(max(10, len(df) * 0.6), 7))
    ax.bar(range(len(df)), df["true_ddg"],
           color=colours, edgecolor="none")
    ax.errorbar(range(len(df)), df["true_ddg"],
                yerr=df["true_ddg_sem"], fmt="none",
                ecolor="black", elinewidth=1, capsize=4)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df["mutation"], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("True ΔΔG (kcal/mol)", fontsize=10)
    ax.set_title(
        "Mutation Ranking — True ΔΔG vs Wild Type\n"
        "(from FINAL_RESULTS_MMPBSA.dat)",
        fontsize=12
    )
    patches = [
        mpatches.Patch(color=C_FAV,   label="Stabilising (ΔΔG < 0)"),
        mpatches.Patch(color=C_UNFAV, label="Destabilising (ΔΔG > 0)"),
    ]
    ax.legend(handles=patches, fontsize=9)
    fig.tight_layout()
    path = outdir / "batch_mutation_ranking.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  [saved] {path}")


# =============================================================================
# 6.  SAVE OUTPUTS
# =============================================================================

def save_csvs(merged, name, outdir):
    cols = [
        "residue", "chain_type", "chain", "resname", "resnum",
        "total_wt", "total_sem_wt", "total_mut", "total_sem_mut",
        "delta_total", "delta_sem", "significant",
        "delta_vdw", "delta_elec", "delta_polar", "delta_nonpolar",
    ]
    full_path = outdir / f"{name}_full_comparison.csv"
    sig_path  = outdir / f"{name}_significant_residues.csv"
    merged[cols].to_csv(full_path, index=False)
    merged[merged["significant"]][cols].to_csv(sig_path, index=False)
    print(f"  [saved] {full_path}")
    print(f"  [saved] {sig_path}")


# =============================================================================
# 7.  SINGLE-RUN PIPELINE
# =============================================================================

def run_single(wt_decomp, wt_results, mut_decomp, mut_results, name, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"\nParsing WT  decomp  : {wt_decomp}")
    wt_df  = parse_decomp(wt_decomp)
    print(f"Parsing WT  results : {wt_results}")
    wt_res = parse_results(wt_results)

    print(f"Parsing MUT decomp  : {mut_decomp}")
    mut_df  = parse_decomp(mut_decomp)
    print(f"Parsing MUT results : {mut_results}")
    mut_res = parse_results(mut_results)

    print(f"WT residues : {len(wt_df)}  |  Mutant residues : {len(mut_df)}")

    merged  = compare(wt_df, mut_df)
    summary = summarise(merged, wt_res, mut_res, name)

    print(f"\nGenerating plots...")
    plot_per_residue(merged, name, outdir)
    plot_energy_breakdown(merged, name, outdir)
    plot_wt_vs_mut_scatter(merged, name, outdir)
    save_csvs(merged, name, outdir)

    return summary


# =============================================================================
# 8.  BATCH PIPELINE
# =============================================================================

def run_batch(batch_csv, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    batch = pd.read_csv(batch_csv)
    required = {"name", "wt_decomp", "wt_results", "mut_decomp", "mut_results"}
    if not required.issubset(batch.columns):
        sys.exit(
            f"\n[ERROR] batch CSV must have columns:\n"
            f"  name, wt_decomp, wt_results, mut_decomp, mut_results"
        )

    summaries = []
    for _, row in batch.iterrows():
        sub_dir = outdir / row["name"]
        summary = run_single(
            row["wt_decomp"], row["wt_results"],
            row["mut_decomp"], row["mut_results"],
            row["name"], sub_dir
        )
        summaries.append(summary)

    summary_df   = pd.DataFrame(summaries).sort_values("true_ddg")
    summary_path = outdir / "all_mutations_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    print(f"\n  [saved] {summary_path}")

    print(f"\n{'='*65}")
    print("  MUTATION RANKING (most stabilising → most destabilising)")
    print(f"{'='*65}")
    cols = ["mutation", "wt_dg", "mut_dg", "true_ddg", "true_ddg_sem",
            "n_significant", "n_stabilising", "n_destabilising"]
    print(summary_df[cols].to_string(index=False))

    plot_batch_summary(summaries, outdir)


# =============================================================================
# 9.  ARGUMENT PARSING
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Compare WT vs Mutant MMPBSA results and decomposition."
    )
    parser.add_argument("--wt-decomp",   help="WT  FINAL_DECOMP_MMPBSA.dat")
    parser.add_argument("--wt-results",  help="WT  FINAL_RESULTS_MMPBSA.dat")
    parser.add_argument("--mut-decomp",  help="MUT FINAL_DECOMP_MMPBSA.dat")
    parser.add_argument("--mut-results", help="MUT FINAL_RESULTS_MMPBSA.dat")
    parser.add_argument("--name", default="mutation",
                        help="Label for this mutation (e.g. APC1)")
    parser.add_argument("--batch",
                        help="CSV for batch mode "
                             "(columns: name, wt_decomp, wt_results, "
                             "mut_decomp, mut_results)")
    parser.add_argument("--out", default="mmpbsa_analysis",
                        help="Output directory (default: mmpbsa_analysis/)")

    args = parser.parse_args()

    if args.batch:
        run_batch(args.batch, args.out)
    elif all([args.wt_decomp, args.wt_results,
              args.mut_decomp, args.mut_results]):
        run_single(
            args.wt_decomp, args.wt_results,
            args.mut_decomp, args.mut_results,
            args.name, args.out
        )
    else:
        parser.print_help()
        sys.exit(
            "\n[ERROR] Provide either:\n"
            "  --wt-decomp + --wt-results + --mut-decomp + --mut-results\n"
            "  or --batch"
        )


if __name__ == "__main__":
    main()
