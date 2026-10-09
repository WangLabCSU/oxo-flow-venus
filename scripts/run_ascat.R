#!/usr/bin/env Rscript
# Venus ASCAT driver: allelic copy number, purity and ploidy for paired
# tumour/normal WGS. Implements the ASCAT >= 3.1.0 HTS chain exactly as
# documented by VanLoo-lab (ascat.prepareHTS / ascat.runAscat man pages) and
# as orchestrated by gcap::gcap.runASCAT (Wang et al. 2024,
# Nat Commun 15:3231):
#   prepareHTS -> loadData -> correctLogR -> plotRawData ->
#   aspcf(gg = NULL) -> plotSegmentedData -> runAscat(gamma = 1) -> metrics.
# WGS mode: no BED (whole genome), penalty = 70, GC + replication-timing
# correction from the official G1000 hg38 reference set, gamma = 1 (required
# for HTS data; 0.55 is the SNP-array legacy value).
# Autosomes 1:22 only: patient sex is unannotated for this cohort and ASCAT's
# sex-chromosome handling needs a reported karyotype, so gender is fixed to
# "XX" and X/Y are excluded from the fit (limitation noted in the report).

suppressPackageStartupMessages(library(ASCAT))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA_character_) {
  i <- match(flag, args)
  if (!is.na(i) && length(args) >= i + 1) args[i + 1] else default
}
need_arg <- function(flag) {
  v <- get_arg(flag)
  if (is.na(v)) stop("missing required argument ", flag, call. = FALSE)
  v
}

bam_tumour   <- need_arg("--tumor")
bam_normal   <- need_arg("--normal")
tumourname   <- need_arg("--tumor-name")
normalname   <- need_arg("--normal-name")
jobname      <- get_arg("--job", tumourname)
outdir       <- need_arg("--outdir")
refdir       <- need_arg("--refdir")
gc_file      <- need_arg("--gc")
rt_file      <- need_arg("--rt")
loci_base    <- get_arg("--loci-prefix", "G1000_loci_hg38_")
alleles_base <- get_arg("--alleles-prefix", "G1000_alleles_hg38_")
nthreads     <- as.integer(get_arg("--threads", "16"))
penalty      <- as.numeric(get_arg("--penalty", "70"))
gender       <- get_arg("--gender", "XX")
genome_build <- get_arg("--genome", "hg38")

stopifnot(gender %in% c("XX", "XY"))

check_bam <- function(path) {
  if (!file.exists(path)) stop("BAM not found: ", path, call. = FALSE)
  idx1 <- paste0(path, ".bai")           # <bam>.bam.bai
  idx2 <- sub("\\.bam$", ".bai", path)   # <bam>.bai
  if (!file.exists(idx1) && !file.exists(idx2))
    stop("BAM index not found for ", path,
         " (looked for ", basename(idx1), " and ", basename(idx2), ")",
         call. = FALSE)
}
check_bam(bam_tumour)
check_bam(bam_normal)
for (f in c(gc_file, rt_file))
  if (!file.exists(f)) stop("reference file not found: ", f, call. = FALSE)

# Absolutize file arguments BEFORE setwd(outdir) below. Rule inputs arrive
# relative to the workdir; prepareHTS forks workers whose working directory
# is outdir, so a relative seq.file would silently fail file.exists() inside
# ascat.getAlleleCounts (which warns and returns per chromosome) and the run
# would die later with "length(files) > 0 is not TRUE".
abs_path <- function(p) normalizePath(p, mustWork = TRUE)
bam_tumour <- abs_path(bam_tumour)
bam_normal <- abs_path(bam_normal)
refdir     <- abs_path(refdir)
gc_file    <- abs_path(gc_file)
rt_file    <- abs_path(rt_file)

ac_exe <- Sys.which("alleleCounter")
if (!nzchar(ac_exe)) stop("alleleCounter not on PATH", call. = FALSE)

# The G1000 hg38 reference ships one loci/alleles file per chromosome named
# with PLAIN numbers (<base>1.txt) -- ASCAT's chrom_names must stay unprefixed
# because readAllelesFiles recovers the chromosome from the filename and
# readAlleleCountFiles strips a leading "chr" from the content column, so
# plain names are what loadData's SNPpos %in% chrs filter sees. File CONTENT
# stays chr-prefixed ("chr22\t15513541"): alleleCounter matches it against
# chr-prefixed BAM contigs, and ASCAT un-prefixes it on read. Probe both
# layouts so either naming works.
loci_prefix <- file.path(refdir, loci_base)
has_loci <- function(pfx) file.exists(paste0(pfx, "1.txt"))
if (has_loci(loci_prefix)) {
  chrom_names <- as.character(1:22)
} else if (has_loci(paste0(loci_prefix, "chr"))) {
  chrom_names <- paste0("chr", 1:22)
} else {
  stop("no loci files found at ", loci_prefix, "<1|chr1>.txt", call. = FALSE)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
oldwd <- setwd(outdir)
on.exit(setwd(oldwd), add = TRUE)

message("[ascat] allele counting (", length(chrom_names), " chromosomes, ",
        nthreads, " threads)")
ascat.bc <- ascat.prepareHTS(
  tumourseqfile = bam_tumour,
  normalseqfile = bam_normal,
  tumourname = tumourname,
  normalname = normalname,
  allelecounter_exe = ac_exe,
  alleles.prefix = file.path(refdir, alleles_base),
  loci.prefix = loci_prefix,
  nthreads = nthreads,
  gender = gender,
  genomeVersion = genome_build,
  chrom_names = chrom_names,
  min_base_qual = 20,
  min_map_qual = 35
)

# prepareHTS names all four LogR/BAF tables after the TUMOUR (default
# filenames) even for the normal sample; feed the normal pair back as the
# germline slots of loadData (paired mode). gg is left NULL so germline
# genotypes derive from BAF extremes -- ascat.predictGermlineGenotypes is
# unreliable for targeted/HTS data (ASCAT issue #73) and gcap keeps NULL.
ascat.bc <- ascat.loadData(
  Tumor_LogR_file = paste0(tumourname, "_tumourLogR.txt"),
  Tumor_BAF_file = paste0(tumourname, "_tumourBAF.txt"),
  Germline_LogR_file = paste0(tumourname, "_normalLogR.txt"),
  Germline_BAF_file = paste0(tumourname, "_normalBAF.txt"),
  chrs = chrom_names,
  gender = gender,
  genomeVersion = genome_build
)
ascat.bc <- ascat.correctLogR(ascat.bc, GCcontentfile = gc_file,
                              replictimingfile = rt_file)

img_dir <- "plots"
img_prefix <- paste0(jobname, ".")
dir.create(img_dir, showWarnings = FALSE)
ascat.plotRawData(ascat.bc, img.dir = img_dir, img.prefix = img_prefix)

message("[ascat] ASPCF segmentation (penalty = ", penalty, ")")
ascat.bc <- ascat.aspcf(ascat.bc, ascat.gg = NULL, penalty = penalty)
ascat.plotSegmentedData(ascat.bc, img.dir = img_dir, img.prefix = img_prefix)

message("[ascat] purity/ploidy fit (gamma = 1)")
ascat.output <- ascat.runAscat(ascat.bc, gamma = 1, pdfPlot = TRUE,
                               img.dir = img_dir, img.prefix = img_prefix)
qc <- ascat.metrics(ascat.bc, ascat.output)
ascat.output <- c(ascat.output, list(QC = qc))

num <- function(x) as.numeric(as.character(x))
in_list <- function(l, x) length(l) > 0 && x %in% l
purity_ploidy <- data.frame(
  pair = jobname,
  sample = tumourname,
  purity = num(ascat.output$purity[1]),
  ploidy = num(ascat.output$ploidy[1]),
  psi = num(ascat.output$psi[1]),
  aberrant_cell_fraction = num(ascat.output$aberrantcellfraction[1]),
  goodness_of_fit = num(ascat.output$goodnessOfFit[1]),
  non_aberrant = in_list(ascat.output$nonaberrantarrays, tumourname),
  failed = in_list(ascat.output$failedarrays, tumourname),
  n_segments = nrow(ascat.output$segments),
  stringsAsFactors = FALSE
)
write.table(purity_ploidy, "purity_ploidy.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)

seg <- ascat.output$segments
seg <- seg[seg$sample == tumourname, ]
write.table(seg, "segments.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# ascat.metrics() returns ONE data.frame (rownames = sample names), not a
# per-sample list; qc[[1]] would grab the first COLUMN (sex). Coerce both
# possible shapes before adding the sample column.
qc_df <- if (is.data.frame(qc)) qc else do.call(rbind, qc)
qc_df <- data.frame(sample = rownames(qc_df), qc_df, row.names = NULL,
                    stringsAsFactors = FALSE)
write.table(qc_df, "metrics.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Temp allele-count files are large and fully consumed by the LogR/BAF
# tables; remove (as gcap does).
tmp <- list.files(pattern = "alleleFrequencies", full.names = TRUE)
if (length(tmp)) unlink(tmp)

saveRDS(ascat.output, file = paste0(jobname, ".ASCAT.rds"))

# The rule declares these outputs; fail loudly here rather than leaving the
# engine to discover a missing file after a silently degraded run.
expected <- c(
  purity_ploidy.tsv = "purity_ploidy.tsv",
  segments.tsv = "segments.tsv",
  metrics.tsv = "metrics.tsv",
  rds = paste0(jobname, ".ASCAT.rds"),
  logr = paste0(tumourname, "_tumourLogR.txt"),
  baf = paste0(tumourname, "_tumourBAF.txt"),
  tumour_png = file.path(img_dir, paste0(img_prefix, tumourname, ".tumour.png")),
  germline_png = file.path(img_dir, paste0(img_prefix, tumourname, ".germline.png")),
  aspcf_png = file.path(img_dir, paste0(img_prefix, tumourname, ".ASPCF.png")),
  profile_pdf = file.path(img_dir, paste0(img_prefix, tumourname,
                                          ".ASCATprofile.pdf"))
)
missing <- expected[!file.exists(expected)]
if (length(missing))
  stop("expected outputs missing: ", paste(names(missing), collapse = ", "),
       call. = FALSE)

message("[ascat] done: purity=", purity_ploidy$purity,
        " ploidy=", purity_ploidy$ploidy,
        " goodnessOfFit=", purity_ploidy$goodness_of_fit,
        " segments=", purity_ploidy$n_segments)
