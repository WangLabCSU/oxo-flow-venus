#!/usr/bin/env Rscript
#' Calculate Tumor Mutation Burden (TMB)

suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description='Calculate TMB from VCF')
parser$add_argument('--input', required=TRUE, help='Input annotated VCF')
parser$add_argument('--output', required=TRUE, help='Output TMB file')
parser$add_argument('--target-size', type='double', default=1.0, help='Target region size in Mb')
args <- parser$parse_args()

mutations <- 0
con <- file(args$input, 'r')
while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    if (!startsWith(line, '#')) {
        mutations <- mutations + 1
    }
}
close(con)

tmb <- mutations / args$target_size

write.table(
    data.frame(
        sample = basename(args$input),
        mutations = mutations,
        target_mb = args$target_size,
        tmb = round(tmb, 2)
    ),
    file = args$output,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE
)
