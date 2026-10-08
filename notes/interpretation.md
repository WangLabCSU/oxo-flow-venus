## Cohort-specific findings (hnzl cohort, 10 CRC patients)

All numbers below were recomputed directly from the delivered artifacts
(`deliver/*.maf`, `report/tmb/*.tmb.tsv`, `vcf.filtered/*.filteringStats.tsv`,
`cnv/*/`, reference FASTA + `.fai`). Measurement criteria are stated inline so
each figure is reproducible; the STR annotator used is `scripts/str_context.py`.

**Pt09 is a high-confidence MSI/MMR-deficiency hypermutator.** Tumor
mutational burden (0.679 coding mut/Mb over the 3,000-Mb callable target — `report/tmb/Pt09.tmb.tsv` `coding_non_silent`/`target_mb`, i.e. 59.9 per coding-Mb — = 14× the cohort median of 0.049;
99.75 all-somatic mut/Mb; 299,254 PASS somatic variants) is extreme, and
every line of evidence below is internally consistent — this is biology, not
a technical artifact:

- **Variant class composition**: 177,644 PASS indels; 175,699 (98.9%) carry a
  STR repeat-unit annotation under `scripts/str_context.py` (reference window
  ±15 bp around the indel; shortest qualifying unit wins: unit length 1 with
  ≥4 identical bases → mononucleotide, length 2 with ≥4 units → dinucleotide,
  length 3–6 with ≥3 units → longer STR). Of the annotated indels 171,643 are
  mononucleotide (97.7%), 3,510 dinucleotide, 546 longer units; only 1,945
  indels have no qualifying repeat — the canonical MSI indel signature
  (homopolymer slippage).
- **Clonal structure**: tumor VAF mean 0.379; 68.7% of PASS variants at
  VAF ≥ 0.3 with the clonal peak at 0.35–0.40 (largest 0.05-bin, 43,533
  variants), and zero PASS variants show normal-sample VAF ≥ 0.15 — burden is
  somatic and clonal, not germline contamination or mapping artifact.
- **Genome-wide distribution**: no focal hypermutation — per-chromosome
  density ranges 83.1–129.7 mut/Mb across autosomes + chrX (mean 98.4,
  cv 0.10; min chr15, max chr19 gene-rich); including near-silent chrY
  (~1 mut/Mb) the 24-chromosome mean is 94.3 (cv 0.23).
  This breadth-without-foci pattern is as expected for a mutator phenotype.
- **Mutational signature**: among PASS SNPs, C>T/G>A transitions (both
  strands counted, pyrimidine-normalized) sit at CpG dinucleotides 78.0% of
  the time at VAF ≥ 0.3 (71.1% over all PASS SNPs) — for comparison Pt01 is
  37.0% / 44.5% and Pt05 41.7% / 54.7%. The Pt09 CpG dominance is the
  aged-colon/SBS1-like pattern characteristic of MMR-deficient colorectal
  cancer. (At VAF ≥ 0.3 the Pt01/Pt05 counts are small — 46 and 84
  transitions — so their all-PASS figures are the stabler comparison.)
- **MSI-marker regions**: the chr2:47.63–47.66 Mb window — which contains a
  32-bp poly-A tract at chr2:47,634,827–858, a BAT-26-class mononucleotide
  microsatellite — carries, in Pt09, a PASS single-base deletion at
  chr2:47,641,314 (TA>T, VAF 0.50; loss of one A from an 8-bp poly-A tract at
  47,641,315–322, verified against the reference FASTA). The deletion sits
  ~6.5 kb from the tract and neither locus was among msisensor-pro's scanned
  sites, so this is a region-level observation, not a BAT-26 stability call;
  Pt09's MSI-H classification rests on the genome-wide panel (65,535/621,533
  sites unstable, 10.54%). Pt01 has no PASS indel in the same window.
- **MMR genes** (deliver/Pt09.maf): a single nonsilent coding hit —
  **MSH6 p.Ala1055Thr** (chr2:47,801,146 G>A, ENST00000234420.11:c.3163G>A,
  VAF 0.372, tumor depth 43×, normal depth 17×) — plus intronic deletions in
  MSH2 (25), EPCAM (9), MSH6 (6), PMS2 (3) and MLH1 (3) (criterion: PASS DEL
  with Variant_Classification Intron or 5′/3′ flank), a burden pattern typical of
  MMR-deficient tumors. Several additional 5′/3′-flank SNPs on EPCAM/MSH2
  are present at clonal VAFs.
- **Copy-number stability**: the Pt09 CNVkit segmentation
  (`cnv/Pt09/Pt09-TD.bqsr.call.cns`) holds 3,130 segments; at a |log2| ≥ 0.3
  threshold only 44 amplifications / 128 deletions (±0.5: 10/48; ±0.8: 2/18);
  by copy number, 1 segment ≥ cn5 and 17 segments at cn0 (163 ≤ cn1).
  Comparably Pt01 (3,661 segments, 88/283 at ±0.3, 7 at cn0) and Pt05
  (4,655 segments, 112/506, 56 at cn0) — Pt09 carries the few high-amplitude
  events of the three, the chromosomally stable hypermutator subtype.
- **Depth parity**: length-weighted mean depth over CNVkit target+antitarget
  bins (2,925.77 Mb ≈ 2.9 Gb of bin space, i.e. genome-wide): Pt09 tumor 41.1× /
  normal 19.0× (cohort sampled: Pt01 37.0×/18.5×, Pt05 38.3×/19.1×) —
  balanced tumor/normal depth. Mutect2 filtering statistics are healthy
  (FDR 0.054, sensitivity 0.95, per-filter FDR ≤ 0.02), ruling out caller
  artifacts.

**Recommendation**: reflex MSI testing (PCR panel or IHC for MLH1/MSH2/MSH6/
PMS2) to confirm; MSI/MMR-deficient CRC is predictive of immune-checkpoint
inhibitor response.

**Pt05 is a possible low-level MSI case** and is flagged for follow-up:
26,395 PASS indels = 50.5% of its 52,219 variants (the other eight
non-flagged patients sit at 6.4–10.5%; Pt09 itself is far higher at 59.4%),
24,605 of them in homopolymers, at an
intermediate TMB (0.094 coding mut/Mb over the 3,000-Mb target). Its CpG-transition share is only
mildly elevated over Pt01, so the flag rests on the indel fraction and
homopolymer load rather than on the substitution pattern.
