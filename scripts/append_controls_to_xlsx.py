#!/usr/bin/env python3
"""
Append non-hotspot control rows to an existing XLSX sheet without external packages.

This script reads hotspot rows from an XLSX (OOXML), samples control sequences from
gene FASTA files, and writes a new XLSX with controls appended to the same sheet.
"""

from __future__ import annotations

import argparse
import math
import random
import re
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple
import xml.etree.ElementTree as ET
from xml.sax.saxutils import escape


NS_MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
NS_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
NS = {"a": NS_MAIN, "r": NS_REL}
ET.register_namespace("", NS_MAIN)
ET.register_namespace("r", NS_REL)

DNA = set("ACGTN")


def clean_seq(value: str) -> Optional[str]:
    if value is None:
        return None
    seq = str(value).strip().upper()
    if not seq:
        return None
    seq = "".join(ch for ch in seq if ch in DNA)
    return seq or None


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
    return not (a_end + min_gap < b_start or b_end + min_gap < a_start)


def col_to_index(col: str) -> int:
    result = 0
    for ch in col:
        result = result * 26 + (ord(ch) - ord("A") + 1)
    return result


def index_to_col(index: int) -> str:
    out = []
    i = index
    while i > 0:
        i, rem = divmod(i - 1, 26)
        out.append(chr(ord("A") + rem))
    return "".join(reversed(out))


CELL_REF_RE = re.compile(r"^([A-Z]+)(\d+)$")


def parse_cell_ref(cell_ref: str) -> Tuple[str, int]:
    m = CELL_REF_RE.match(cell_ref)
    if not m:
        raise ValueError(f"Bad cell ref: {cell_ref}")
    return m.group(1), int(m.group(2))


def read_shared_strings(zf: zipfile.ZipFile) -> List[str]:
    try:
        sst_root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
    except KeyError:
        return []
    strings: List[str] = []
    for si in sst_root.findall("a:si", NS):
        text = "".join(t.text or "" for t in si.iterfind(".//a:t", NS))
        strings.append(text)
    return strings


def cell_value(cell: ET.Element, shared_strings: Sequence[str]) -> str:
    t = cell.attrib.get("t")
    if t == "inlineStr":
        text = "".join(x.text or "" for x in cell.iterfind(".//a:t", NS))
        return text
    v = cell.find("a:v", NS)
    if v is None or v.text is None:
        return ""
    if t == "s":
        idx = int(v.text)
        if 0 <= idx < len(shared_strings):
            return shared_strings[idx]
        return ""
    return v.text


def parse_sheet_rows(sheet_root: ET.Element, shared_strings: Sequence[str]) -> List[Dict[str, str]]:
    sheet_data = sheet_root.find("a:sheetData", NS)
    if sheet_data is None:
        return []
    rows_out: List[Dict[str, str]] = []
    for row in sheet_data.findall("a:row", NS):
        row_map: Dict[str, str] = {}
        for cell in row.findall("a:c", NS):
            cref = cell.attrib.get("r")
            if not cref:
                continue
            col, _ = parse_cell_ref(cref)
            row_map[col] = cell_value(cell, shared_strings)
        rows_out.append(row_map)
    return rows_out


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
        quotas[g] = int(round(target_controls * hotspot_counts[g] / total_hotspots))

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

    for g in list(quotas):
        quotas[g] = min(quotas[g], candidate_counts.get(g, 0))

    current = sum(quotas.values())
    needed = target_controls - current
    if needed > 0:
        genes_by_spare = sorted(
            (g for g in quotas if candidate_counts[g] > quotas[g]),
            key=lambda g: candidate_counts[g] - quotas[g],
            reverse=True,
        )
        i = 0
        while needed > 0 and genes_by_spare:
            g = genes_by_spare[i % len(genes_by_spare)]
            if quotas[g] < candidate_counts[g]:
                quotas[g] += 1
                needed -= 1
            i += 1
            genes_by_spare = [x for x in genes_by_spare if quotas[x] < candidate_counts[x]]

    return quotas


def build_inline_cell_xml(ref: str, value: str) -> str:
    esc = escape(value)
    return f'<c r="{ref}" t="inlineStr"><is><t>{esc}</t></is></c>'


def patch_sheet_xml_text(
    original_xml: bytes,
    rows_xml: Sequence[str],
    max_col_idx: int,
    max_row_idx: int,
) -> bytes:
    text = original_xml.decode("utf-8")
    row_blob = "".join(rows_xml)

    # Inject new rows right before the closing sheetData tag.
    if "</sheetData>" not in text:
        raise ValueError("Could not locate </sheetData> in worksheet XML.")
    text = text.replace("</sheetData>", f"{row_blob}</sheetData>", 1)

    new_dim = f'A1:{index_to_col(max_col_idx)}{max_row_idx}'
    dim_re = re.compile(r'(<dimension[^>]*\bref=")([^"]*)(")')
    if dim_re.search(text):
        text = dim_re.sub(rf"\1{new_dim}\3", text, count=1)
    else:
        # Insert dimension after opening worksheet tag if missing.
        ws_open = re.search(r"<worksheet[^>]*>", text)
        if not ws_open:
            raise ValueError("Could not locate worksheet root tag.")
        insert_at = ws_open.end()
        text = text[:insert_at] + f'<dimension ref="{new_dim}"/>' + text[insert_at:]

    return text.encode("utf-8")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Append sampled controls to XLSX Hotspots sheet.")
    p.add_argument("--input-xlsx", required=True, help="Input XLSX path.")
    p.add_argument("--fasta-dir", required=True, help="Directory containing <Gene>.fasta files.")
    p.add_argument("--output-xlsx", required=True, help="Output XLSX path.")
    p.add_argument("--sheet-name", default="Hotspots", help="Sheet containing hotspot rows.")
    p.add_argument("--target-controls", type=int, default=2000, help="Desired control rows.")
    p.add_argument("--nearby-max-distance", type=int, default=500, help="Max bp distance from hotspots.")
    p.add_argument("--min-gap", type=int, default=0, help="Min bp gap from hotspot 15bp intervals.")
    p.add_argument("--seed", type=int, default=42, help="Random seed.")
    p.add_argument("--gene-col", default="Gene", help="Gene column header.")
    p.add_argument("--seq15-col", default="Sequence_15_Long", help="15-mer column header.")
    p.add_argument("--seq9-col", default="Sequence_9_Long", help="9-mer column header.")
    p.add_argument("--hotspot-col", default="Hotspot", help="Hotspot label column header.")
    p.add_argument("--skip-n-bases", action="store_true", help="Exclude controls containing N.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    rng = random.Random(args.seed)

    input_xlsx = Path(args.input_xlsx)
    output_xlsx = Path(args.output_xlsx)
    fasta_dir = Path(args.fasta_dir)

    with zipfile.ZipFile(input_xlsx, "r") as zin:
        shared_strings = read_shared_strings(zin)
        wb_root = ET.fromstring(zin.read("xl/workbook.xml"))
        rels_root = ET.fromstring(zin.read("xl/_rels/workbook.xml.rels"))
        rel_map = {r.attrib["Id"]: r.attrib["Target"] for r in rels_root}

        target_sheet_rel = None
        for sheet in wb_root.findall("a:sheets/a:sheet", NS):
            if sheet.attrib.get("name") == args.sheet_name:
                target_sheet_rel = sheet.attrib.get(f"{{{NS_REL}}}id")
                break
        if not target_sheet_rel:
            raise ValueError(f"Sheet '{args.sheet_name}' not found")

        rel_target = rel_map[target_sheet_rel].lstrip("/")
        target_sheet_path = "xl/" + rel_target if not rel_target.startswith("xl/") else rel_target
        original_sheet_xml = zin.read(target_sheet_path)
        sheet_root = ET.fromstring(original_sheet_xml)
        rows = parse_sheet_rows(sheet_root, shared_strings)
        if not rows:
            raise ValueError("No rows found in sheet")

        header = rows[0]
        col_by_name: Dict[str, str] = {}
        for col, val in header.items():
            col_by_name[val.strip()] = col

        required = [args.gene_col, args.seq15_col, args.seq9_col, args.hotspot_col]
        missing = [h for h in required if h not in col_by_name]
        if missing:
            raise ValueError(f"Missing required columns in header: {missing}")

        gene_col = col_by_name[args.gene_col]
        seq15_col = col_by_name[args.seq15_col]
        seq9_col = col_by_name[args.seq9_col]
        hotspot_col = col_by_name[args.hotspot_col]

        hotspot_rows: List[Tuple[str, Optional[str], Optional[str]]] = []
        for r in rows[1:]:
            gene = (r.get(gene_col) or "").strip()
            if not gene:
                continue
            seq15 = clean_seq(r.get(seq15_col, ""))
            seq9 = clean_seq(r.get(seq9_col, ""))
            hotspot_rows.append((gene, seq15, seq9))

        genes = sorted({g for g, _, _ in hotspot_rows})
        gene_sequences: Dict[str, str] = {}
        for gene in genes:
            gene_sequences[gene] = read_fasta_sequence(fasta_dir / f"{gene}.fasta")

        hotspot_counts = Counter(g for g, _, _ in hotspot_rows)
        hotspot_centers_by_gene: Dict[str, Set[int]] = defaultdict(set)
        hotspot_intervals_by_gene: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
        hotspot_seq9_set: Set[str] = set()
        hotspot_seq15_set: Set[str] = set()
        unmatched_rows = 0

        for gene, seq15, seq9 in hotspot_rows:
            gseq = gene_sequences[gene]
            if seq9:
                hotspot_seq9_set.add(seq9)
            if seq15:
                hotspot_seq15_set.add(seq15)

            centers: Set[int] = set()
            if seq15 and len(seq15) == 15:
                for s in find_all_matches(gseq, seq15):
                    centers.add(s + 7)
            if not centers and seq9 and len(seq9) == 9:
                for s in find_all_matches(gseq, seq9):
                    centers.add(s + 4)

            if not centers:
                unmatched_rows += 1
                continue

            for c in centers:
                hotspot_centers_by_gene[gene].add(c)
                hs_start = max(0, c - 7)
                hs_end = min(len(gseq) - 1, c + 7)
                hotspot_intervals_by_gene[gene].append((hs_start, hs_end))

        control_candidates: Dict[str, List[int]] = {}
        for gene in genes:
            gseq = gene_sequences[gene]
            hotspots = hotspot_centers_by_gene.get(gene, set())
            if not hotspots:
                control_candidates[gene] = []
                continue
            protected = hotspot_intervals_by_gene.get(gene, [])
            candidate_set: Set[int] = set()
            for hc in hotspots:
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
                    if seq15 in hotspot_seq15_set or seq9 in hotspot_seq9_set:
                        continue
                    candidate_set.add(c)
            control_candidates[gene] = sorted(candidate_set)

        candidate_counts = {g: len(v) for g, v in control_candidates.items()}
        quotas = allocate_gene_quotas(hotspot_counts, candidate_counts, args.target_controls)

        controls: List[Tuple[str, str, str]] = []
        for gene, quota in quotas.items():
            centers = control_candidates[gene]
            if not centers or quota <= 0:
                continue
            chosen = centers if quota >= len(centers) else rng.sample(centers, quota)
            gseq = gene_sequences[gene]
            for c in chosen:
                seq15 = gseq[c - 7 : c + 8]
                seq9 = gseq[c - 4 : c + 5]
                controls.append((gene, seq15, seq9))

        rng.shuffle(controls)

        sheet_data = sheet_root.find("a:sheetData", NS)
        if sheet_data is None:
            raise ValueError("sheetData missing")
        existing_rows = sheet_data.findall("a:row", NS)
        last_row_idx = int(existing_rows[-1].attrib["r"]) if existing_rows else 0

        rows_xml: List[str] = []
        for i, (gene, seq15, seq9) in enumerate(controls, start=1):
            row_idx = last_row_idx + i
            row_xml = (
                f'<row r="{row_idx}">'
                f'{build_inline_cell_xml(f"{gene_col}{row_idx}", gene)}'
                f'{build_inline_cell_xml(f"{seq15_col}{row_idx}", seq15)}'
                f'{build_inline_cell_xml(f"{seq9_col}{row_idx}", seq9)}'
                f'{build_inline_cell_xml(f"{hotspot_col}{row_idx}", "0")}'
                "</row>"
            )
            rows_xml.append(row_xml)

        max_col_idx = max(col_to_index(c) for c in header.keys())
        new_sheet_xml = patch_sheet_xml_text(
            original_xml=original_sheet_xml,
            rows_xml=rows_xml,
            max_col_idx=max_col_idx,
            max_row_idx=last_row_idx + len(controls),
        )

        with zipfile.ZipFile(output_xlsx, "w", compression=zipfile.ZIP_DEFLATED) as zout:
            for info in zin.infolist():
                data = zin.read(info.filename)
                if info.filename == target_sheet_path:
                    data = new_sheet_xml
                zout.writestr(info.filename, data)

    print(f"Input rows (excluding header): {len(rows) - 1}")
    print(f"Unmatched hotspot rows: {unmatched_rows}")
    print(f"Controls appended: {len(controls)}")
    print(f"Output workbook: {output_xlsx}")
    if len(controls) < args.target_controls:
        print("Warning: fewer controls than requested were available under constraints.")


if __name__ == "__main__":
    main()
