#!/usr/bin/env python3
"""Convert a PASS-filtered somatic VCF (VEP-annotated) into a minimal MAF.

Reads the CSQ INFO field written by vep-rs (--vcf output). The CSQ subfield
names and their order are read from the ##INFO=<ID=CSQ,...,Format=...> header
line, so the parser stays correct regardless of which annotation flags were
enabled. With --pick there is one picked consequence per variant; the first
CSQ entry is used. Chosen deliberately over the full vcf2maf stack so the
pipeline has zero Perl dependencies while retaining the fields that matter:
gene symbol, consequence, HGVS c/p, depth, VAF, and variant classification.
"""
from __future__ import annotations

import argparse
import gzip
import re
import sys

# MAF column order (2.4.1 subset). Reference_Allele appears once; Tumor_Seq_Allele1
# carries REF and Tumor_Seq_Allele2 the alternate allele.
MAF_HEADER = [
    "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
    "Reference_Allele", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode",
    "Variant_Classification", "Variant_Type", "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2", "dbSNP_RS", "t_depth", "t_ref_count", "t_alt_count",
    "n_depth", "n_ref_count", "n_alt_count", "t_vaf", "HGVSc", "HGVSp",
    "Transcript_ID", "Consequence",
]

CSQ_HEADER_RE = re.compile(
    r'##INFO=<ID=CSQ,.*?Format=([^">]+)', re.IGNORECASE
)


def parse_csq_format(line: str) -> list[str]:
    """Extract the ordered subfield names from the CSQ INFO header line."""
    m = CSQ_HEADER_RE.search(line)
    if not m:
        return []
    return [f.strip() for f in m.group(1).split("|")]


def parse_csq_entries(csq_value: str, fields: list[str]) -> list[dict[str, str]]:
    """Split a CSQ value into per-transcript dicts keyed by subfield name."""
    entries: list[dict[str, str]] = []
    for chunk in csq_value.split(","):
        parts = chunk.split("|")
        if fields and len(parts) == len(fields):
            entries.append(dict(zip(fields, parts)))
    return entries


def classify(conseq: str) -> str:
    """Map VEP consequence term to a MAF Variant_Classification."""
    first = conseq.split(",")[0]
    table = {
        "missense_variant": "Missense_Mutation",
        "synonymous_variant": "Silent",
        "stop_gained": "Nonsense_Mutation",
        "stop_lost": "Nonstop_Mutation",
        "frameshift_variant": "Frame_Shift_Del",
        "inframe_insertion": "In_Frame_Ins",
        "inframe_deletion": "In_Frame_Del",
        "start_lost": "Translation_Start_Site",
        "splice_acceptor_variant": "Splice_Site",
        "splice_donor_variant": "Splice_Site",
        "stop_retained_variant": "Silent",
        "coding_sequence_variant": "Coding_Variant",
        "intron_variant": "Intron",
        "5_prime_UTR_variant": "5'UTR",
        "3_prime_UTR_variant": "3'UTR",
        "upstream_gene_variant": "5'Flank",
        "downstream_gene_variant": "3'Flank",
        "intergenic_variant": "IGR",
    }
    return table.get(first, "Other")


def snp_or_indel(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "SNP"
    if len(ref) > len(alt):
        return "DEL"
    if len(alt) > len(ref):
        return "INS"
    return "ONINV"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--tumor-id", required=True)
    ap.add_argument("--normal-id", default="")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    opener = gzip.open if args.vcf.endswith(".gz") else open
    csq_fields: list[str] = []
    rows: list[list[str]] = []
    with opener(args.vcf, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                if not csq_fields and "ID=CSQ" in line:
                    csq_fields = parse_csq_format(line)
                continue
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, _id, ref, alt = f[0], f[1], f[2], f[3], f[4]
            info: dict[str, str] = {}
            for kv in f[7].split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    info[k] = v
            entries = parse_csq_entries(info.get("CSQ", ""), csq_fields)
            csq = entries[0] if entries else {}
            fmt_keys = f[8].split(":") if len(f) > 9 else []
            tvals = f[9].split(":") if len(f) > 9 else []
            t = dict(zip(fmt_keys, tvals))
            n = {}
            if len(f) > 10:
                n = dict(zip(fmt_keys, f[10].split(":")))

            def ad_of(d: dict) -> tuple[str, str]:
                ad = d.get("AD", "").split(",")
                if len(ad) >= 2:
                    return ad[0], ad[1]
                return "", ""

            t_ref, t_alt = ad_of(t)
            n_ref, n_alt = ad_of(n)
            t_depth = t.get("DP", "")
            n_depth = n.get("DP", "")
            vaf = ""
            if t_depth and t_alt:
                try:
                    vaf = f"{int(t_alt) / int(t_depth):.4f}"
                except (ValueError, ZeroDivisionError):
                    vaf = ""

            conseq = csq.get("Consequence", "")
            rows.append([
                csq.get("SYMBOL", ""),
                chrom, pos, str(int(pos) + len(ref) - 1),
                ref, args.tumor_id, args.normal_id,
                classify(conseq),
                snp_or_indel(ref, alt),
                ref, alt,
                _id if _id.startswith("rs") else "",
                t_depth, t_ref, t_alt,
                n_depth, n_ref, n_alt,
                vaf,
                csq.get("HGVSc", ""), csq.get("HGVSp", ""),
                csq.get("Feature", ""),
                conseq,
            ])

    with open(args.output, "w") as out:
        out.write("## MAF derived from {} (vep-rs CSQ) by venus\n".format(args.vcf))
        out.write("#version 2.4.1\n")
        out.write("\t".join(MAF_HEADER) + "\n")
        for r in rows:
            out.write("\t".join(r) + "\n")
    print(f"[vcf2maf_light] wrote {len(rows)} variants -> {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
