## Cohort-specific findings (hnzl cohort, 10 CRC patients)

**Pt09 is a high-confidence MSI/MMR-deficiency hypermutator.** Tumor
mutational burden (0.679 coding mut/Mb = 14× the cohort median of 0.049;
99.75 all-somatic mut/Mb; 299,254 PASS somatic variants) is extreme, and
every line of evidence below is internally
consistent — this is biology, not a technical artifact:

- **Variant class composition**: 156,615 PASS indels carry STR repeat-unit
  annotations; 152,939 are mononucleotide (97%) and 2,833 dinucleotide
  repeats — the canonical MSI indel signature (homopolymer slippage).
- **Clonal structure**: corrected tumor VAF mean 0.385; 71% of PASS variants
  at VAF ≥ 0.3 with the clonal peak at 0.35–0.40, and all normal-sample VAFs
  at these sites < 0.15 — burden is somatic and clonal, not germline
  contamination or mapping artifact.
- **Genome-wide distribution**: variant density is uniform across
  chromosomes (no focal hypermutation), as expected for a mutator phenotype.
- **Mutational signature**: somatic C>T/G>A transitions are 77.8% at CpG
  dinucleotides (VAF ≥ 0.3), versus 48.9% (Pt01) and 60.8% (Pt05) —
  an aged-colon/SBS1-like pattern characteristic of MMR-deficient colorectal
  cancer.
- **MSI-marker regions**: the BAT-26 mononucleotide marker region
  (chr2:47.63–47.66 Mb) carries a PASS homopolymer indel in Pt09 (none in
  Pt01).
- **MMR genes** (deliver/Pt09.maf): a single nonsilent coding hit —
  **MSH6 p.Ala1055Thr** (chr2:47,801,146 G>A, ENST00000234420.11:c.3163G>A,
  VAF 0.372, depth 43×) — plus multiple intronic deletions in MSH2 (25),
  EPCAM (9), MSH6 (6), PMS2 (8) and MLH1 (7), a burden pattern typical of
  MMR-deficient tumors.
- **Copy-number stability**: only 51 amplifications / 126 deletions (Pt01:
  75/322, Pt05: 78/634) — the chromosomally stable hypermutator subtype.
- **Depth parity**: tumor mean depth 41× (cohort 36.7–38.1×), normal ~19×;
  Mutect2 filtering statistics are healthy (FDR 0.054, sensitivity 0.95,
  per-filter FDR ≤ 0.02), ruling out caller artifacts.

**Recommendation**: reflex MSI testing (PCR panel or IHC for MLH1/MSH2/MSH6/
PMS2) to confirm; MSI/MMR-deficient CRC is predictive of immune-checkpoint
inhibitor response.

**Pt05 is a possible low-level MSI case** and is flagged for follow-up:
26,395 PASS indels (~50% of its variants; cohort norm ~20%), 22,195 in
homopolymers, elevated indel/SNP ratio and intermediate CpG bias (60.8%),
at an intermediate TMB (0.094 coding mut/Mb).
