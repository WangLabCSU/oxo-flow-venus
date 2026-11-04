#!/usr/bin/env python3
"""Screen STAR Chimeric.out.junction into per-sample fusion candidates (N12).

SCREENING LEVEL ONLY — parses split-read junction evidence STAR already
produced (--chimSegmentMin 12 --chimOutType WithinBAM SoftClip Junctions, see
rules/rna.oxoflow) and annotates gene names from the GENCODE GTF. NOT a fusion
caller: no spanning-fragment support, no transcript reconstruction, no
RT-artifact filtering. Every output row carries the screening disclosure;
positives require orthogonal RNA confirmation (FISH/RT-PCR/panel) before
clinical use. Full reconstruction is STAR-Fusion's job (P2, see README).

Chimeric.out.junction (STAR >= 2.5.1b, ONE ROW PER SPLIT READ, no header when
--chimMultimapNmax=0; 15 columns with --outSAMattrRGline, 14 without):
  1 donor_chr            6 acceptor_strand
  2 donor junction bp*   7 motif (-1 none / 0 noncanonical / 1 GT/AG / 2 CT/AC)
  3 donor_strand         8 homology left
  4 acceptor_chr         9 homology right
  5 acceptor junction bp*   10 read name
  11 donor first aligned base (1-based)   12 donor CIGAR (p-extended)
  13 acceptor first aligned base          14 acceptor CIGAR
  15 [optional RG]
  * junction bp = first INTRONIC base (1-based) = last aligned base + 1.

"donor/acceptor" here follow STAR's read-segment order, NOT transcript 5'/3'
direction (panel matching is order-insensitive, clinical_sv semantics).

Anchors are recomputed from the CIGARs (sum of M/I/=/X) so per-side support is
explicit in the output instead of trusting the aligner threshold blindly.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys

# --- thresholds (disclosed constants; keep in sync with star_align flags) ----
# STAR's --chimJunctionOverhangMin 12 already guarantees >=12bp anchors; the
# explicit floor guards against junction files produced by looser runs.
MIN_ANCHOR = 12
MIN_SPLIT_READS = 2       # generic candidate floor (unique split reads)
SCREENING_DISCLOSURE = ("screening_level: STAR split-read junction only; "
                        "no spanning-fragment support / transcript "
                        "reconstruction / RT-artifact filtering; direction "
                        "is STAR read-order not 5'-3'; requires orthogonal "
                        "RNA confirmation")

MOTIFS = {"-1": "none", "0": "non_canonical", "1": "GT/AG", "2": "CT/AC"}
CIGAR_OP = re.compile(r"(\d+)([MIDNSHPX=])")
READ_CONSUMING = {"M", "I", "=", "X"}

COLUMNS = ["sample", "fusion_pair", "gene_a", "gene_b",
           "breakpoint_a", "strand_a", "breakpoint_b", "strand_b",
           "motif", "anchor_a", "anchor_b", "split_reads",
           "panel_hit", "panel_label", "evidence_class", "disclosure"]


def anchor_len(cigar: str) -> int:
    """Read-aligned bases of a STAR p-extended CIGAR (M/I/=/X; p = mate gap)."""
    return sum(int(n) for n, op in CIGAR_OP.findall(cigar or "")
               if op in READ_CONSUMING)


def gene_spans_from_gtf(gtf_path):
    """All gene spans: {chrom: [(start, end, name), ...]} (1-based, inclusive)."""
    spans: dict[str, list[tuple[int, int, str]]] = {}
    with open(gtf_path) as fh:
        for line in fh:
            if "\tgene\t" not in line or 'gene_name "' not in line:
                continue
            fl = line.rstrip("\n").split("\t")
            name = line.split('gene_name "')[1].split('"')[0]
            spans.setdefault(fl[0], []).append((int(fl[3]), int(fl[4]), name))
    return spans


def gene_at(spans, chrom, pos):
    """Gene whose span contains pos, else nearest within 10kb (promoter-
    proximal breakpoints), else ''."""
    best, best_dist = "", 10_000
    for start, end, name in spans.get(chrom, ()):
        if start <= pos <= end:
            return name
        dist = start - pos if pos < start else pos - end
        if dist < best_dist:
            best, best_dist = name, dist
    return best


def load_panel(path):
    with open(path) as fh:
        return [dict(r) for r in csv.DictReader(fh, delimiter="\t")]


def panel_match(panel, g1, g2):
    """clinical_sv semantics: unordered exact pair match, or single-gene match
    against a wildcard ('*') side only. Returns (hit, label)."""
    for p in panel:
        a = (p.get("gene_5") or "").upper()
        b = (p.get("gene_3") or "").upper()
        if g1 and g2 and a != "*" and b != "*" and {g1, g2} == {a, b}:
            return True, f"{p.get('gene_5')}--{p.get('gene_3')}"
        if b == "*" and a != "*" and (g1 == a or g2 == a):
            other = g2 if g1 == a else g1
            if other:
                return True, (f"{p.get('gene_5')}--{other} "
                              "(wildcard single-gene rule)")
            return True, (f"{p.get('gene_5')} (single breakpoint; "
                          "partner unresolved)")
        if a == "*" and b != "*" and (g1 == b or g2 == b):
            other = g2 if g1 == b else g1
            if other:
                return True, (f"{other}--{p.get('gene_3')} "
                              "(wildcard single-gene rule)")
            return True, (f"{p.get('gene_3')} (single breakpoint; "
                          "partner unresolved)")
    return False, ""


def read_junctions(path):
    """Yield parsed per-read rows (skips a STAR-Fusion-format header)."""
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("chr_donorA"):
                continue  # header (chimOutJunctionFormat 1 runs)
            f = line.split("\t")
            if len(f) < 14:
                continue
            try:
                yield {"d_chr": f[0], "d_junc": int(f[1]), "d_strand": f[2],
                       "a_chr": f[3], "a_junc": int(f[4]), "a_strand": f[5],
                       "motif": f[6], "read": f[9],
                       "d_cigar": f[11], "a_cigar": f[13]}
            except ValueError:
                continue  # malformed line: skipped


def screen(junction_path, gtf_path, panel_path, sample,
           min_anchor=MIN_ANCHOR, min_split=MIN_SPLIT_READS):
    spans = gene_spans_from_gtf(gtf_path)
    panel = load_panel(panel_path)

    # canonical junction key: position-sorted ends + motif. STAR may report
    # the same junction with donor/acceptor swapped across reads — collapse.
    groups: dict[tuple, dict] = {}
    for j in read_junctions(junction_path):
        ends = sorted([(j["d_chr"], j["d_junc"]), (j["a_chr"], j["a_junc"])])
        key = (ends[0][0], ends[0][1], ends[1][0], ends[1][1], j["motif"])
        g = groups.setdefault(key, {
            "meta": (j["d_chr"], j["d_junc"], j["d_strand"],
                     j["a_chr"], j["a_junc"], j["a_strand"]),
            "reads": {}, "anchors": {},
        })
        prev = g["reads"].get(j["read"])
        da, aa = anchor_len(j["d_cigar"]), anchor_len(j["a_cigar"])
        if prev is None or da + aa > prev:
            g["reads"][j["read"]] = da + aa
        for c, p, anc in ((j["d_chr"], j["d_junc"], da),
                          (j["a_chr"], j["a_junc"], aa)):
            k = f"{c}:{p}"
            g["anchors"][k] = max(g["anchors"].get(k, 0), anc)

    rows = []
    for key, g in groups.items():
        d_chr, d_pos, d_strand, a_chr, a_pos, a_strand = g["meta"]
        split_reads = len(g["reads"])
        anc_a = g["anchors"].get(f"{d_chr}:{d_pos}", 0)
        anc_b = g["anchors"].get(f"{a_chr}:{a_pos}", 0)
        if min(anc_a, anc_b) < min_anchor:
            continue
        gene_a = gene_at(spans, d_chr, d_pos - 1) or gene_at(spans, d_chr, d_pos)
        gene_b = gene_at(spans, a_chr, a_pos - 1) or gene_at(spans, a_chr, a_pos)
        hit, label = panel_match(panel, gene_a.upper(), gene_b.upper())
        if hit:
            # panel hits are ALWAYS kept (min_split governs non-panel
            # candidates only — a stricter floor must not silently discard
            # actionable single-read hits; the evidence class discloses it)
            evidence = ("panel_hit" if split_reads >= 2
                        else "panel_hit_single_read")
        elif split_reads >= max(3, min_split):
            evidence = f"supported_{split_reads}reads"
        elif split_reads >= min_split:
            evidence = f"low_support_{split_reads}reads"
        else:
            continue
        if gene_a and gene_b:
            fusion = f"{gene_a}--{gene_b}"
        elif gene_a or gene_b:
            fusion = f"{gene_a or gene_b} (partner gene unresolved)"
        else:
            fusion = f"intergenic_{d_chr}:{d_pos}--{a_chr}:{a_pos}"
        rows.append({
            "sample": sample, "fusion_pair": fusion,
            "gene_a": gene_a, "gene_b": gene_b,
            "breakpoint_a": f"{d_chr}:{d_pos}", "strand_a": d_strand,
            "breakpoint_b": f"{a_chr}:{a_pos}", "strand_b": a_strand,
            "motif": MOTIFS.get(key[4], f"code_{key[4]}"),
            "anchor_a": anc_a, "anchor_b": anc_b,
            "split_reads": split_reads,
            "panel_hit": "yes" if hit else "no", "panel_label": label,
            "evidence_class": evidence, "disclosure": SCREENING_DISCLOSURE,
        })
    rows.sort(key=lambda r: (r["panel_hit"] != "yes",
                             -int(r["split_reads"]), r["fusion_pair"]))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--junction", required=True,
                    help="Chimeric.out.junction path (per-read rows)")
    ap.add_argument("--gtf", required=True, help="GENCODE GTF (gene spans)")
    ap.add_argument("--panel", required=True, help="actionable_fusion_pairs.tsv")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--min-anchor", type=int, default=MIN_ANCHOR)
    ap.add_argument("--min-split-reads", type=int, default=MIN_SPLIT_READS)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    rows = screen(args.junction, args.gtf, args.panel, args.sample,
                  min_anchor=args.min_anchor, min_split=args.min_split_reads)
    with open(args.output, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=COLUMNS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    hits = sum(1 for r in rows if r["panel_hit"] == "yes")
    print(f"[fusion_screen] {args.sample}: {len(rows)} candidates "
          f"({hits} panel hits) -> {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
