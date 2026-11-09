#!/usr/bin/env python3
"""Convert Manta BND somatic SVs into Arriba-format fusions for pVACfuse.

venus#31 P1 item 1 (pvacfuse_run). pVACtools 7.1.5 consumes Arriba TSVs via
``InputConverter.parse_arriba_file`` (lib/input_file_converter.py), which
reads these named columns (DictReader, so order is free; extra columns ok):

  breakpoint1 ("chr:start"), breakpoint2, peptide_sequence (MUST contain the
  '|' junction marker -- ``determine_fusion_sequence`` sys.exits without it),
  #gene1, gene2, transcript_id1, transcript_id2, reading_frame (exactly
  'in-frame' or 'out-of-frame'; any other value drops the row), split_reads1,
  split_reads2, discordant_mates (pVACfuse read support = split1 + split2 +
  discordant; rows below --read-support (default 5) are hard-filtered later;
  Expression 'NA' passes, since we have no StarFusion expression).

Rewrite rules for Manta BND records (deliver/{pair}.sv.vcf.gz = PASS paired
deliverable, or deliver/{pair}.sv.candidate.vcf.gz = tumor-only):

  - MATEID pairs are deduplicated (Manta double-emits every BND); solo BNDs
    (mate record missing) and unparseable ALTs are skipped with a reason.
  - The BND ALT bracket position says which genomic side each side keeps:
      N[c:p[  -> record keeps LEFT,  mate keeps RIGHT   (simple 5'->3')
      N]c:p]  -> record keeps LEFT,  mate keeps LEFT    (head-to-head)
      [c:p[N  -> record keeps RIGHT, mate keeps RIGHT   (tail-to-tail)
      ]c:p]N  -> record keeps RIGHT, mate keeps LEFT
    A side is a usable 5' piece only if it keeps the transcript-UPSTREAM side
    of its gene (LEFT for + strand, RIGHT for - strand) and a usable 3' piece
    keeps the transcript-DOWNSTREAM side (RIGHT for +, LEFT for -). Pairs
    that do not yield exactly one 5' and one 3' piece are skipped with a
    reason -- no guessed orientations.
  - Genes are assigned by overlap (nearest gene span within --window bp
    otherwise, logged); each gene uses its longest-CDS transcript from the
    GENCODE GTF.
  - peptide_sequence = translate(donor CDS piece + acceptor CDS piece),
    '|' inserted at the junction codon boundary (len(donor)//3), capped at
    150 codons per side (codon-aligned, junction frame preserved). Both
    breakpoint anchor bases are kept (Manta's REF-base convention); the
    strict VCF-spec reading would drop the mate's anchor and shift the
    junction by <=1 base -- immaterial downstream (pVACfuse re-scores the
    8-11mers around the junction). A stop codon BEFORE the junction
    drops the row (the fusion ORF ends there); stops after it are stripped
    (pVACtools strips '*' itself). DNA-level approximation: a BND carries no
    transcript/splice evidence -- the FUSION_ID trace column and the
    header-only-when-empty honesty live in the rule, and RNA evidence from
    the fusion_screen module is the corroboration path.
  - reading_frame = len(donor CDS retained) % 3 (0 -> 'in-frame', else
    'out-of-frame'); pVACfuse uses it only for the inframe_fusion vs
    frameshift_fusion label.
  - read support from Manta FORMAT PR/SR alt counts of the TUMOR sample
    (--tumor-sample; auto-detected for single-sample VCFs = the tumor-only
    deliverable, else the highest-PR+SR sample, logged). split_reads1/2 =
    donor/acceptor SR, discordant_mates = donor PR.

Zero fusions is a valid outcome: header-only TSV, exit 0 -- the RULE owns
the honest disclosure. Sequence access uses pysam (envs/pvac); tests inject
a fetch callable, so the conversion core is stdlib-only.
"""
import argparse
import csv
import gzip
import itertools
import re
import sys

ARRIBA_COLUMNS = [
    "#gene1", "gene2", "breakpoint1", "breakpoint2",
    "transcript_id1", "transcript_id2", "reading_frame",
    "peptide_sequence", "split_reads1", "split_reads2",
    "discordant_mates", "FUSION_ID",
]

# Standard genetic code, TCAG base order (NCBI translation table 1); kept
# local so the conversion core stays stdlib-only for tests.
_AA64 = ("FFLLSSSSYY**CC*W"
         "LLLLPPPPHHQQRRRR"
         "IIIMTTTTNNKKSSRR"
         "VVVVAAAADDEEGGGG")
CODON_TABLE = {"".join(c): _AA64[i]
               for i, c in enumerate(itertools.product("TCAG", repeat=3))}

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# (record retained side, mate retained side); LEFT = [.. pos], RIGHT = [pos ..]
_ALT_SIDES = (
    (re.compile(r"[ACGTNacgtn]\[([^:\]]+):(\d+)\["), "LEFT", "RIGHT"),
    (re.compile(r"[ACGTNacgtn]\]([^:\]]+):(\d+)\]"), "LEFT", "LEFT"),
    (re.compile(r"\[([^:\]]+):(\d+)\][ACGTNacgtn]"), "RIGHT", "RIGHT"),
    (re.compile(r"\]([^:\]]+):(\d+)\]\[?"), "RIGHT", "LEFT"),
)

FLANK_AA = 150  # codons kept per side around the junction (codon-aligned)


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def translate(seq):
    """Translate a nucleotide fragment; unknown codons -> X, stops kept '*'
    (pVACtools strips '*' from fusion sequences itself)."""
    aa = []
    for i in range(0, len(seq) - len(seq) % 3, 3):
        aa.append(CODON_TABLE.get(seq[i:i + 3].upper(), "X"))
    return "".join(aa)


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def parse_gtf_cds(path):
    """GENCODE GTF -> (tx_by_id, genes).

    tx_by_id: {transcript_id: {tid, gene, chrom, strand, cds:[(s,e) sorted]}}
    genes:    {gene_name: (chrom, min_start, max_end, strand)}
    """
    tx_by_id, genes = {}, {}
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            attrs = {}
            for kv in f[8].split(";"):
                k, _, v = kv.strip().partition(" ")
                if k:
                    attrs[k] = v.strip('"')
            tid = attrs.get("transcript_id")
            gname = attrs.get("gene_name") or attrs.get("gene_id")
            if not tid or not gname:
                continue
            s, e = int(f[3]), int(f[4])
            rec = tx_by_id.setdefault(tid, {"tid": tid, "gene": gname,
                                            "chrom": f[0], "strand": f[6],
                                            "cds": []})
            rec["cds"].append((s, e))
            g = genes.setdefault(gname, (f[0], s, e, f[6]))
            genes[gname] = (g[0], min(g[1], s), max(g[2], e), g[3])
    for rec in tx_by_id.values():
        rec["cds"].sort()
    return tx_by_id, genes


def longest_tx_for_gene(tx_by_id, gene):
    cands = [t for t in tx_by_id.values() if t["gene"] == gene]
    if not cands:
        return None
    return max(cands, key=lambda r: sum(e - s + 1 for s, e in r["cds"]))


def gene_at(genes, chrom, pos, window):
    """Gene overlapping pos, else the nearest gene on the chromosome within
    `window` bp of its span (None -> caller skips with a reason)."""
    best, best_gap = None, None
    for name, (c, s, e, _strand) in genes.items():
        if c != chrom:
            continue
        gap = 0 if s <= pos <= e else (s - pos if pos < s else pos - e)
        if best is None or gap < best_gap:
            best, best_gap = name, gap
    return best if best is not None and best_gap <= window else None


def piece_kind(strand, retained):
    """'five' if the retained genomic side is transcript-UPSTREAM for this
    gene's strand, 'three' if transcript-DOWNSTREAM, else None."""
    if retained == "LEFT":
        return "five" if strand == "+" else "three"
    return "three" if strand == "+" else "five"


def retained_piece_seq(tx, pos, retained, fetch):
    """CDS bases of `tx` on the retained side of `pos`, in TRANSCRIPT
    orientation (5'->3'). Non-CDS overlap is fine; empty -> not coding."""
    chunks = []
    for s, e in tx["cds"]:
        if retained == "LEFT":
            lo, hi = s, min(e, pos)
        else:
            lo, hi = max(s, pos), e
        if lo > hi:
            continue
        chunks.append(fetch(tx["chrom"], lo - 1, hi).upper())
    if tx["strand"] == "-":
        chunks = [revcomp(c) for c in reversed(chunks)]
    return "".join(chunks)


def parse_alt_side(alt):
    """BND ALT -> (mate_chrom, mate_pos, rec_retained, mate_retained)."""
    for rx, rec_side, mate_side in _ALT_SIDES:
        m = rx.search(alt)
        if m:
            return m.group(1), int(m.group(2)), rec_side, mate_side
    return None, None, None, None


def read_bnd_pairs(vcf_path):
    """Deduplicated MATEID pairs [(rec_a, rec_b)] plus skip counters.

    Each rec = {chrom, pos, id, alt, fmt: {sample: {tag: value}}}.
    Returns (pairs, sample_names, skipped) where skipped counts
    solo/unparseable records that never became a pair.
    """
    samples, by_id = [], {}
    with open_text(vcf_path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = f[9:]
                continue
            if len(f) < 8:
                continue
            info = {}
            for kv in f[7].split(";"):
                k, _, v = kv.partition("=")
                info[k] = v
            if info.get("SVTYPE") != "BND":
                continue
            fmt_keys = f[8].split(":") if len(f) > 8 and f[8] else []
            per_sample = {name: dict(zip(fmt_keys, val.split(":")))
                          for name, val in zip(samples, f[9:])}
            by_id[f[2]] = {"chrom": f[0], "pos": int(f[1]), "id": f[2],
                           "alt": f[4], "fmt": per_sample,
                           "mate": info.get("MATEID", "")}
    pairs, seen = [], set()
    skipped = {}
    for rid, rec in by_id.items():
        mate_id = rec["mate"]
        if not mate_id or mate_id not in by_id:
            skipped["solo_bnd_no_mate"] = \
                skipped.get("solo_bnd_no_mate", 0) + 1
            continue
        key = tuple(sorted((rid, mate_id)))
        if key in seen:
            continue
        seen.add(key)
        mc, mp, rec_side, mate_side = parse_alt_side(rec["alt"])
        if mc is None:
            skipped["unparseable_bnd_alt"] = \
                skipped.get("unparseable_bnd_alt", 0) + 1
            continue
        pairs.append((rec, by_id[mate_id], mc, mp, rec_side, mate_side))
    return pairs, samples, skipped


def alt_count(fmt_map, tag):
    """Alt-allele count from a Manta FORMAT field (``ref,alt`` comma list)."""
    v = (fmt_map or {}).get(tag, "")
    part = v.split(",")[-1].strip() if "," in v else ""
    return int(part) if part.isdigit() else 0


def pick_tumor_sample(pairs, sample_names, preferred):
    """Tumor sample column: --tumor-sample if present; a single-sample VCF
    (the tumor-only deliverable) auto-resolves; otherwise the sample with
    the highest total PR+SR alt support across the calls (logged)."""
    if preferred and preferred in sample_names:
        return preferred
    if len(sample_names) == 1:
        return sample_names[0]

    def total(name):
        t = 0
        for rec, _, _, _, _, _ in pairs:
            fmt = rec["fmt"].get(name, {})
            t += alt_count(fmt, "PR") + alt_count(fmt, "SR")
        return t
    ranked = sorted(sample_names, key=lambda n: -total(n))
    print(f"[manta_to_arriba] --tumor-sample not found in VCF; using "
          f"'{ranked[0]}' (highest PR+SR support) of {sorted(sample_names)}",
          file=sys.stderr)
    return ranked[0]


def build_fusion_row(rec, mate, mate_chrom, mate_pos, rec_side, mate_side,
                     tx_by_id, genes, fetch, window, tumor, skipped):
    """One MATEID pair -> one Arriba-format row (or None + reason)."""
    def note(reason):
        skipped[reason] = skipped.get(reason, 0) + 1

    gene1 = gene_at(genes, rec["chrom"], rec["pos"], window)
    gene2 = gene_at(genes, mate_chrom, mate_pos, window)
    if gene1 is None or gene2 is None:
        note("breakpoint outside annotated gene spans (window)")
        return None
    if gene1 == gene2:
        note("intragenic BND (same-gene; not a gene-pair fusion)")
        return None
    tx1 = longest_tx_for_gene(tx_by_id, gene1)
    tx2 = longest_tx_for_gene(tx_by_id, gene2)
    if tx1 is None or tx2 is None:
        note("assigned gene has no CDS transcript in GTF")
        return None
    strand1, strand2 = tx1["strand"], tx2["strand"]
    kind1 = piece_kind(strand1, rec_side)
    kind2 = piece_kind(strand2, mate_side)
    if {kind1, kind2} != {"five", "three"}:
        note("no 5'/3' transcript model (head/tail-to-tail or strand clash)")
        return None
    # donor = the five' side; map back to (record, pos, retained, tx)
    if kind1 == "five":
        dtx, dpos, dret, dtid = tx1, rec["pos"], rec_side, tx1["tid"]
        atx, apos, aret, atid = tx2, mate_pos, mate_side, tx2["tid"]
        dchrom, achrom = rec["chrom"], mate_chrom
        drec, arec = rec, mate
        dgene, agene = gene1, gene2
    else:
        dtx, dpos, dret, dtid = tx2, mate_pos, mate_side, tx2["tid"]
        atx, apos, aret, atid = tx1, rec["pos"], rec_side, tx1["tid"]
        achrom, dchrom = rec["chrom"], mate_chrom
        arec, drec = rec, mate
        dgene, agene = gene2, gene1
    dseq = retained_piece_seq(dtx, dpos, dret, fetch)
    aseq = retained_piece_seq(atx, apos, aret, fetch)
    if len(dseq) < 3 or len(aseq) < 3:
        note("retained CDS < 3 nt on one side")
        return None
    # Truncate the donor 5' end WITHOUT shifting the junction reading frame:
    # keep K <= FLANK_AA*3 nucleotides with K == full donor length (mod 3).
    full_len = len(dseq)
    keep = min(full_len, FLANK_AA * 3)
    keep -= (keep - full_len) % 3
    dseq = dseq[-keep:]
    aseq = aseq[:FLANK_AA * 3]
    cut = keep // 3                      # complete donor codons before '|'
    aa = translate(dseq + aseq)
    if len(aa) <= cut:
        note("empty fused translation")
        return None
    if "*" in aa[:cut]:
        note("stop codon before junction (ORF ends pre-junction)")
        return None
    peptide = aa[:cut] + "|" + aa[cut:].replace("*", "")
    pr1 = alt_count(drec["fmt"].get(tumor, {}), "PR")
    sr1 = alt_count(drec["fmt"].get(tumor, {}), "SR")
    sr2 = alt_count(arec["fmt"].get(tumor, {}), "SR")
    return {
        "#gene1": dgene, "gene2": agene,
        "breakpoint1": f"{dchrom}:{dpos}",
        "breakpoint2": f"{achrom}:{apos}",
        "transcript_id1": dtid, "transcript_id2": atid,
        "reading_frame": "in-frame" if full_len % 3 == 0
                         else "out-of-frame",
        "peptide_sequence": peptide,
        "split_reads1": sr1, "split_reads2": sr2,
        "discordant_mates": pr1,
        "FUSION_ID": f"{drec['id']}|{arec['id']}",
    }


def convert(vcf_path, gtf_path, fetch, tumor_sample=None, window=5000):
    """VCF + GTF + fetch -> (rows, skipped{reason: count})."""
    tx_by_id, genes = parse_gtf_cds(gtf_path)
    return convert_with(vcf_path, tx_by_id, genes, fetch, tumor_sample,
                        window)


def convert_with(vcf_path, tx_by_id, genes, fetch, tumor_sample=None,
                 window=5000):
    """Core loop on a pre-parsed GTF (shared by MantaToArriba.run)."""
    pairs, samples, skipped = read_bnd_pairs(vcf_path)
    rows = []
    tumor = None
    if pairs:
        tumor = pick_tumor_sample(pairs, samples, tumor_sample)
        for rec, mate, mc, mp, rside, mside in pairs:
            row = build_fusion_row(rec, mate, mc, mp, rside, mside,
                                   tx_by_id, genes, fetch, window, tumor,
                                   skipped)
            if row is not None:
                rows.append(row)
    return rows, skipped


class MantaToArriba:
    """Convert one Manta BND VCF into one Arriba-format TSV (venus#31)."""

    def __init__(self, gtf_path, fetch_fn, tumor_sample=None, window=5000):
        self.tx_by_id, self.genes = parse_gtf_cds(gtf_path)
        self.fetch = fetch_fn
        self.tumor_sample = tumor_sample
        self.window = window
        self.skipped = {}

    def run(self, vcf_path, output_path):
        tx_by_id, genes = self.tx_by_id, self.genes
        pairs, samples, skipped = read_bnd_pairs(vcf_path)
        rows = []
        tumor = None
        if pairs:
            tumor = pick_tumor_sample(pairs, samples, self.tumor_sample)
        for rec, mate, mc, mp, rside, mside in pairs:
            row = build_fusion_row(rec, mate, mc, mp, rside, mside,
                                   tx_by_id, genes, self.fetch, self.window,
                                   tumor, self.skipped)
            if row is not None:
                rows.append(row)
        with open(output_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=ARRIBA_COLUMNS,
                                    delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        total = f"{len(rows)} fusion(s) -> {output_path}; skipped: " \
                f"{skipped or self.skipped or 'none'}"
        print(f"[manta_to_arriba] {total}", file=sys.stderr)
        return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf", required=True,
                    help="Manta BND VCF (.vcf.gz or .vcf): the PASS "
                         "deliverable (paired) or candidateSV (tumor-only)")
    ap.add_argument("--gencode", required=True, help="GENCODE GTF (gz ok)")
    ap.add_argument("--genome-fa", required=True,
                    help="reference FASTA for peptide translation")
    ap.add_argument("--tumor-sample", default=None,
                    help="tumor sample column (auto-detected for "
                         "single-sample VCFs, else highest-support sample)")
    ap.add_argument("--window", type=int, default=5000,
                    help="max distance from a gene span for breakpoint "
                         "assignment (default 5000)")
    ap.add_argument("--output", required=True, help="Arriba-format TSV")
    args = ap.parse_args()
    try:
        import pysam  # envs/pvac ships pysam >= 0.22
    except ImportError:
        sys.exit("pysam is required to read the reference FASTA "
                 "(envs/pvac provides it); refusing to fake peptides")
    fa = pysam.FastaFile(args.genome_fa)
    fetch = lambda c, s, e: fa.fetch(c, s, e).upper()
    converter = MantaToArriba(args.gencode, fetch, args.tumor_sample,
                              args.window)
    return converter.run(args.vcf, args.output)


if __name__ == "__main__":
    sys.exit(main())
