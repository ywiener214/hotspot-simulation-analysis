#!/usr/bin/env python3
"""
Generate non-hotspot control sequences from gene FASTA files.

The script:
1) Reads hotspot rows from an Excel file (expects Gene + 9/15mer columns).
2) Maps hotspot sequences to each gene FASTA.
3) Samples nearby, non-overlapping control centers from the same genes.
4) Writes controls with Gene, Sequence_15_Long, Sequence_9_Long.
"""

from __future__ import annotations

import argparse
import math
import random
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import pandas as pd


DNA = set("ACGTN")


def clean_seq(value: object) -> Optional[str]:
    if value is None:
        return None
    if isinstance(value, float) and math.isnan(value):
        return None
    seq = str(value).strip().upper()
    if not seq:
        return None
    return "".join(ch for ch in seq if ch in DNA)


def read_fasta_sequence(path: Path) -> str:
    if not path.exists():
        raise FileNotFoundError(f"Missing FASTA file: {path}")
    chunks: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            chunks.append(line.upper())
    if not chunks:
        raise ValueError(f"No sequence content in FASTA: {path}")
    return "".join(chunks)


def matches_at(target: str, query: str, start: int) -> bool:
    # Treat N as wildcard in either query or target.
    for i, q in enumerate(query):
        t = target[start + i]
        if q != "N" and t != "N" and q != t:
            return False
    return True


def find_all_matches(target: str, query: str) -> List[int]:
    if not query or len(query) > len(target):
        return []
    out: List[int] = []
    max_start = len(target) - len(query)
    for s in range(max_start + 1):
        if matches_at(target, query, s):
            out.append(s)
    return out


def intervals_overlap(a_start: int, a_end: int, b_start: int, b_end: int, min_gap: int) -> bool:
    # Closed intervals [start, end], with optional required gap between them.
    return not (a_end + min_gap < b_start or b_end + min_gap < a_start)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate control non-hotspot sequences.")
    parser.add_argument("--hotspots-xlsx", required=True, help="Path to Excel file with hotspot rows.")
    parser.add_argument("--sheet", default=0, help="Excel sheet name or index (default: 0).")
    parser.add_argument("--fasta-dir", required=True, help="Directory containing <Gene>.fasta files.")
    parser.add_argument("--output", default="controls_generated.xlsx", help="Output Excel file path.")
    parser.add_argument("--target-controls", type=int, default=2000, help="Desired number of controls.")
    parser.add_argument(
        "--nearby-max-distance",
        type=int,
        default=500,
        help="Maximum distance in bp from any hotspot center for control candidates.",
    )
    parser.add_argument(
        "--min-gap",
        type=int,
        default=0,
        help="Minimum required gap between control interval and hotspot interval (default: 0 means no overlap).",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed.")
    parser.add_argument("--gene-col", default="Gene", help="Gene column name in the hotspot file.")
    parser.add_argument("--seq9-col", default="Sequence_9_Long", help="9-mer column name in the hotspot file.")
    parser.add_argument("--seq15-col", default="Sequence_15_Long", help="15-mer column name in the hotspot file.")
    parser.add_argument(
        "--skip-n-bases",
        action="store_true",
        help="If set, candidate controls containing 'N' are excluded.",
    )
    return parser.parse_args()


def allocate_gene_quotas(
    hotspot_counts: Dict[str, int],
    candidate_counts: Dict[str, int],
    target_controls: int,
) -> Dict[str, int]:
    genes = [g for g in hotspot_counts if candidate_counts.get(g, 0) > 0]
    if not genes:
        return {}
    total_hotspots = sum(hotspot_counts[g] for g in genes)
    quotas: Dict[str, int] = {}
    for g in genes:
        raw = target_controls * hotspot_counts[g] / total_hotspots
        quotas[g] = int(round(raw))

    # Normalize rounding so sum == target_controls before cap by availability.
    diff = target_controls - sum(quotas.values())
    if diff != 0:
        ranked = sorted(
            genes,
            key=lambda x: (
                target_controls * hotspot_counts[x] / total_hotspots - quotas[x],
                hotspot_counts[x],
            ),
            reverse=(diff > 0),
        )
        i = 0
        while diff != 0 and ranked:
            g = ranked[i % len(ranked)]
            if diff > 0:
                quotas[g] += 1
                diff -= 1
            else:
                if quotas[g] > 0:
                    quotas[g] -= 1
                    diff += 1
            i += 1

    # Cap by candidate availability.
    for g in list(quotas):
        quotas[g] = min(quotas[g], candidate_counts.get(g, 0))

    # Redistribute leftover demand to genes with spare candidates.
    current = sum(quotas.values())
    needed = target_controls - current
    if needed > 0:
        spare = {g: candidate_counts[g] - quotas[g] for g in quotas if candidate_counts[g] > quotas[g]}
        if spare:
            spare_genes = sorted(spare, key=lambda x: spare[x], reverse=True)
            idx = 0
            while needed > 0 and spare_genes:
                g = spare_genes[idx % len(spare_genes)]
                if quotas[g] < candidate_counts[g]:
                    quotas[g] += 1
                    needed -= 1
                idx += 1
                if idx > 10_000_000:
                    break
                spare_genes = [k for k in spare_genes if quotas[k] < candidate_counts[k]]

    return quotas


def main() -> None:
    args = parse_args()
    rng = random.Random(args.seed)

    hotspots_path = Path(args.hotspots_xlsx)
    fasta_dir = Path(args.fasta_dir)
    output_path = Path(args.output)

    df = pd.read_excel(hotspots_path, sheet_name=args.sheet)
    required = {args.gene_col, args.seq9_col, args.seq15_col}
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in hotspot file: {missing}")

    hotspot_rows = df[[args.gene_col, args.seq9_col, args.seq15_col]].copy()
    hotspot_rows[args.gene_col] = hotspot_rows[args.gene_col].astype(str).str.strip()
    hotspot_rows[args.seq9_col] = hotspot_rows[args.seq9_col].map(clean_seq)
    hotspot_rows[args.seq15_col] = hotspot_rows[args.seq15_col].map(clean_seq)
    hotspot_rows = hotspot_rows[hotspot_rows[args.gene_col] != ""]

    genes = sorted(set(hotspot_rows[args.gene_col]))
    gene_sequences: Dict[str, str] = {}
    for gene in genes:
        fasta_path = fasta_dir / f"{gene}.fasta"
        gene_sequences[gene] = read_fasta_sequence(fasta_path)

    hotspot_counts = Counter(hotspot_rows[args.gene_col])
    hotspot_centers_by_gene: Dict[str, Set[int]] = defaultdict(set)
    hotspot_intervals_by_gene: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    hotspot_seq9_set: Set[str] = set()
    hotspot_seq15_set: Set[str] = set()
    unmatched_rows = 0

    for _, row in hotspot_rows.iterrows():
        gene = row[args.gene_col]
        seq9 = row[args.seq9_col]
        seq15 = row[args.seq15_col]
        gseq = gene_sequences[gene]

        if seq9:
            hotspot_seq9_set.add(seq9)
        if seq15:
            hotspot_seq15_set.add(seq15)

        centers: Set[int] = set()
        if seq15 and len(seq15) == 15:
            starts15 = find_all_matches(gseq, seq15)
            for s in starts15:
                centers.add(s + 7)
        if not centers and seq9 and len(seq9) == 9:
            starts9 = find_all_matches(gseq, seq9)
            for s in starts9:
                centers.add(s + 4)

        if not centers:
            unmatched_rows += 1
            continue

        for c in centers:
            hotspot_centers_by_gene[gene].add(c)
            # Protect 15bp around hotspot center so controls do not overlap.
            left = max(0, c - 7)
            right = min(len(gseq) - 1, c + 7)
            hotspot_intervals_by_gene[gene].append((left, right))

    control_candidates: Dict[str, List[int]] = {}
    for gene in genes:
        gseq = gene_sequences[gene]
        centers = hotspot_centers_by_gene.get(gene, set())
        if not centers:
            control_candidates[gene] = []
            continue

        protected = hotspot_intervals_by_gene.get(gene, [])
        candidate_set: Set[int] = set()
        for hc in centers:
            low = max(7, hc - args.nearby_max_distance)
            high = min(len(gseq) - 8, hc + args.nearby_max_distance)
            for c in range(low, high + 1):
                c15_start, c15_end = c - 7, c + 7
                c9_start, c9_end = c - 4, c + 4
                blocked = False
                for hs, he in protected:
                    if intervals_overlap(c15_start, c15_end, hs, he, args.min_gap):
                        blocked = True
                        break
                if blocked:
                    continue

                seq15 = gseq[c15_start : c15_end + 1]
                seq9 = gseq[c9_start : c9_end + 1]
                if len(seq15) != 15 or len(seq9) != 9:
                    continue
                if args.skip_n_bases and ("N" in seq15 or "N" in seq9):
                    continue
                if seq9 in hotspot_seq9_set or seq15 in hotspot_seq15_set:
                    continue
                candidate_set.add(c)
        control_candidates[gene] = sorted(candidate_set)

    candidate_counts = {g: len(v) for g, v in control_candidates.items()}
    quotas = allocate_gene_quotas(hotspot_counts, candidate_counts, args.target_controls)

    sampled_rows: List[dict] = []
    for gene, quota in quotas.items():
        centers = control_candidates[gene]
        if not centers or quota <= 0:
            continue
        if quota >= len(centers):
            chosen = centers
        else:
            chosen = rng.sample(centers, quota)
        gseq = gene_sequences[gene]
        hotspot_centers = sorted(hotspot_centers_by_gene.get(gene, []))
        for c in chosen:
            seq15 = gseq[c - 7 : c + 8]
            seq9 = gseq[c - 4 : c + 5]
            nearest = min(abs(c - h) for h in hotspot_centers) if hotspot_centers else None
            sampled_rows.append(
                {
                    "Gene": gene,
                    "Sequence_15_Long": seq15,
                    "Sequence_9_Long": seq9,
                    "Control_Center_1Based": c + 1,
                    "Nearest_Hotspot_Distance_bp": nearest,
                    "Label": "control",
                }
            )

    controls_df = pd.DataFrame(sampled_rows)
    controls_df = controls_df.sample(frac=1.0, random_state=args.seed).reset_index(drop=True)
    controls_df.to_excel(output_path, index=False)

    requested = args.target_controls
    made = len(controls_df)
    print(f"Generated {made} controls (requested {requested}).")
    print(f"Unique genes with controls: {controls_df['Gene'].nunique() if made else 0}")
    print(f"Hotspot rows not mapped to FASTA sequence: {unmatched_rows}")
    if made < requested:
        print("Warning: Fewer controls were available under current constraints.")


if __name__ == "__main__":
    main()
