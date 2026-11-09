#!/usr/bin/env python3
"""scRNA corroboration for neoantigen candidate genes (venus#31 P1 item 5).

For each candidate gene of a pair's pVACseq run, compute what fraction of the
pair's single-cell profile expresses that gene, from the pair's OWN
scrna/count/{pair_id}/outs/filter_feature.h5ad (dnbc4tools), optionally split
per cluster via outs/analysis/cluster.csv (a declared scrna_count output).

HONEST DESIGN CONSTRAINT (read before extending):
  dnbc4tools' filter_feature.h5ad carries NO malignant-cell annotation —
  there is no tumor/normal cell label to restrict the denominator to. The
  percentages this script emits are therefore over ALL filtered cells of the
  sample. A low value can mean "gene not expressed by tumor cells" OR
  "tumor cells are a small fraction of the capture". Consumers must treat
  this as a coarse corroboration signal (is the mutated gene expressed at
  all in this sample?), never as tumor-cell expression prevalence. The
  caveat travels with the data: neoantigen_report.py embeds it in the
  summary detail and markdown whenever the corroboration is merged, and the
  rule comment + README restate it. The per-cluster breakdown exists so an
  operator can EYEBALL whether expression concentrates in one cluster
  (putative malignant) — no cluster is ever labeled malignant here.

Reader design: h5py-only (no scanpy/anndata import), verified against the
anndata 0.1.0 on-disk schema:
  - gene names    : var/_index — a Group {values: ...} (anndata nullable-
                    string encoding, what dnbc4tools writes) or a plain
                    string Dataset (older writers)
  - cell barcodes : obs/_index, same two encodings  - X: either a plain Dataset (dense, cells x genes) or a Group with attrs
    encoding-type in {csr_matrix, csc_matrix} + shape + datasets
    data/indices/indptr. CSC is normalized to cell-major (CSR) orientation
    once; per-gene stats are then O(nnz-of-gene) via searchsorted, so the
    full matrix is never densified in the sparse path.
Genes absent from the matrix get an explicit note
("gene_not_in_scRNA_matrix") — var namespaces differ between references and
variant annotation, and a silent zero would be indistinguishable from
"genuinely not expressed".

Runs in envs/scrna (the only venus env with h5py); no new env, no new lock
gap. Stdlib + h5py + numpy only.
"""
import argparse
import csv
import sys

import h5py
import numpy as np

OUT_HEADER = ["gene", "n_cells_total", "n_expressing", "pct_expressing",
              "mean_expr_positive", "per_cluster_pct", "note"]
GENE_COL_ALIASES = ("gene name", "gene")
CLUSTER_COLS = ("cluster", "louvain", "leiden", "seurat_clusters")
BARCODE_COLS = ("barcode", "cell", "cell_id", "cell_barcode")


def read_string_index(f, path):
    """<path>/_index -> list of names, both anndata encodings: a Group
    {values: ...} (nullable-string encoding, what dnbc4tools writes) or a
    flat string Dataset (older writers)."""
    frame = f.get(path)
    if frame is None:
        sys.exit(f"[scrna_corroboration] h5ad has no {path} — not an anndata "
                 "matrix this script understands")
    idx = frame.get("_index")
    if idx is None:
        sys.exit(f"[scrna_corroboration] h5ad has no {path}/_index — not an "
                 "anndata dataframe this script understands")
    vals = idx["values"][:] if isinstance(idx, h5py.Group) else idx[:]
    return [x.decode() if isinstance(x, bytes) else str(x) for x in vals]


class MatrixStats:
    """Cell-major view over X with per-gene positive-cell statistics.

    gene_cells(j) returns (n_expressing, positive_sum, cell_indices):
    cells = barcodes with expression > 0 for gene j; positive_sum adds those
    cells' expression only. Explicit zeros are never counted as expressing.
    """

    def __init__(self, f):
        if "X" not in f:
            sys.exit("[scrna_corroboration] h5ad has no X matrix")
        x = f["X"]
        if isinstance(x, h5py.Dataset):
            m = np.asarray(x[:])
            if m.ndim != 2:
                sys.exit(f"[scrna_corroboration] dense X has shape {m.shape}, "
                         "expected 2-D (cells x genes)")
            self.n_cells, self.n_genes = m.shape
            self._dense = m
            self._csr = None
            return
        enc = x.attrs.get("encoding-type", "")
        enc = enc.decode() if isinstance(enc, bytes) else str(enc)
        if enc not in ("csr_matrix", "csc_matrix"):
            sys.exit(f"[scrna_corroboration] unsupported X encoding {enc!r} "
                     "(expected csr_matrix or csc_matrix)")
        shape = tuple(int(s) for s in x.attrs["shape"])
        data = x["data"][:].astype(np.float64)
        if enc == "csc_matrix":
            # normalize to cell-major (CSR): expand gene ids from indptr
            # ranges, sort all nnz by cell id, rebuild the cell indptr
            gene_ids = np.repeat(np.arange(shape[1], dtype=np.int64),
                                 np.diff(x["indptr"][:].astype(np.int64)))
            cells = x["indices"][:].astype(np.int64)
            order = np.argsort(cells, kind="stable")
            data = data[order]
            indices = gene_ids[order]
            counts = np.bincount(cells, minlength=shape[0])
            indptr = np.concatenate(([0], np.cumsum(counts))).astype(np.int64)
        else:
            indices = x["indices"][:].astype(np.int64)
            indptr = x["indptr"][:].astype(np.int64)
        self.n_cells, self.n_genes = shape
        self._dense = None
        self._csr = (data, indices, indptr)

    def gene_cells(self, j):
        """(n_expressing, positive_sum, cell_indices) for gene column j."""
        if self._dense is not None:
            col = self._dense[:, j]
            cells = np.flatnonzero(col > 0)
            return len(cells), float(col[cells].sum(dtype=np.float64)), cells
        data, indices, indptr = self._csr
        nnz = np.flatnonzero(indices == j)
        pos = data[nnz] > 0
        nnz = nnz[pos]
        cells = np.searchsorted(indptr, nnz, side="right") - 1
        return len(nnz), float(data[nnz].sum()), cells


def read_cluster_map(path, barcodes):
    """cluster.csv -> (per-cell cluster codes, cluster labels in code order).

    Barcodes are matched exactly first; if that matches under half the h5ad
    cells, a trailing '-<int>' suffix is stripped on BOTH sides (10x '-1'
    convention vs dnbc4tools bare barcodes). Unmatched cells keep code -1
    and are excluded from per-cluster counts (their expression still counts
    in the overall columns).
    """
    try:
        with open(path) as fh:
            reader = csv.reader(fh)
            header = next(reader)
            lower = [h.strip().lower() for h in header]
            c_bar = next((i for i, h in enumerate(lower)
                          if h in BARCODE_COLS), 0)
            c_cl = next((i for i, h in enumerate(lower)
                         if h in CLUSTER_COLS), len(header) - 1)
            clus_of = {}
            for row in reader:
                if len(row) > max(c_bar, c_cl):
                    clus_of[row[c_bar].strip()] = row[c_cl].strip()
    except (OSError, csv.Error) as exc:
        sys.exit(f"[scrna_corroboration] cluster file unreadable ({exc})")

    def build(clus_map, keyfn):
        labels = {}
        codes = np.full(len(barcodes), -1, dtype=np.int64)
        for i, bc in enumerate(barcodes):
            cl = clus_map.get(keyfn(bc))
            if cl:
                if cl not in labels:
                    labels[cl] = len(labels)
                codes[i] = labels[cl]
        return codes, labels

    def strip(bc):
        head, _, tail = bc.rpartition("-")
        return head if tail.isdigit() else bc

    codes, labels = build(clus_of, lambda bc: bc)
    if (codes >= 0).sum() < 0.5 * len(barcodes):
        # retry with trailing '-<int>' stripped on BOTH sides (10x '-1'
        # convention vs dnbc4tools bare barcodes — either side may carry it)
        codes2, labels2 = build({strip(k): v for k, v in clus_of.items()},
                                strip)
        if (codes2 >= 0).sum() > (codes >= 0).sum():
            codes, labels = codes2, labels2
    n_matched = int((codes >= 0).sum())
    if n_matched == 0:
        sys.stderr.write(
            "[scrna_corroboration] WARNING: no h5ad barcode matched "
            f"{path} — per-cluster column will be empty\n")
    return codes, labels


def per_cluster_str(codes, labels, cells, n_per_clus):
    """'label:n_pos/n_total; ...' for the expressing cells of one gene."""
    cc = codes[cells]
    cc = cc[cc >= 0]
    counts = np.bincount(cc, minlength=len(labels))
    return "; ".join(f"{lab}:{int(counts[k])}/{int(n_per_clus[k])}"
                     for k, lab in enumerate(labels))


def read_candidate_genes(path):
    """pVACtools filtered.tsv -> unique candidate genes in file order."""
    genes = []
    seen = set()
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            sys.exit(f"[scrna_corroboration] {path} is empty")
        lower = [h.strip().lower() for h in header]
        try:
            i_gene = next(i for i, h in enumerate(lower)
                          if h in GENE_COL_ALIASES)
        except StopIteration:
            sys.exit(f"[scrna_corroboration] {path}: no gene column (looked "
                     f"for {GENE_COL_ALIASES}); header was {header}")
        for row in reader:
            if len(row) > i_gene:
                g = row[i_gene].strip()
                if g and g not in seen:
                    seen.add(g)
                    genes.append(g)
    return genes


def main():
    ap = argparse.ArgumentParser(description="scRNA corroboration percentages")
    ap.add_argument("--h5ad", required=True,
                    help="scrna/count/{pair}/outs/filter_feature.h5ad")
    ap.add_argument("--pvac-tsv", required=True,
                    help="neoantigen/{pair}/filtered.tsv (Gene Name column)")
    ap.add_argument("--cluster", default="",
                    help="optional outs/analysis/cluster.csv "
                         "(barcode -> cluster)")
    ap.add_argument("--out", required=True, help="corroboration TSV")
    args = ap.parse_args()

    genes = read_candidate_genes(args.pvac_tsv)

    with h5py.File(args.h5ad, "r") as f:
        var = read_string_index(f, "var")
        barcodes = read_string_index(f, "obs")
        mat = MatrixStats(f)
    if mat.n_cells != len(barcodes):
        sys.exit(f"[scrna_corroboration] X has {mat.n_cells} cells but obs "
                 f"has {len(barcodes)} barcodes — inconsistent h5ad")
    idx_of = {}
    for i, g in enumerate(var):
        idx_of.setdefault(g, i)

    cell_clus = clus_labels = n_per_clus = None
    if args.cluster:
        cell_clus, clus_labels = read_cluster_map(args.cluster, barcodes)
        n_per_clus = np.bincount(cell_clus[cell_clus >= 0],
                                 minlength=len(clus_labels))

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(OUT_HEADER)
        for g in genes:
            j = idx_of.get(g)
            if j is None:
                w.writerow([g, "", "", "", "", "", "gene_not_in_scRNA_matrix"])
                continue
            n_e, pos_sum, cells = mat.gene_cells(j)
            pc = (per_cluster_str(cell_clus, clus_labels, cells, n_per_clus)
                  if cell_clus is not None else "")
            mean_pos = f"{pos_sum / n_e:.2f}" if n_e else ""
            w.writerow([g, mat.n_cells, n_e,
                        f"{100.0 * n_e / mat.n_cells:.1f}"
                        if mat.n_cells else "",
                        mean_pos, pc, ""])
    n_found = sum(1 for g in genes if g in idx_of)
    sys.stderr.write(
        f"[scrna_corroboration] {mat.n_cells} cells x {len(var)} genes in "
        f"h5ad; {n_found}/{len(genes)} candidate genes found -> {args.out}\n"
        f"[scrna_corroboration] CAVEAT: pct_expressing is over ALL filtered "
        f"cells — dnbc4tools h5ad has no malignant-cell annotation; coarse "
        f"corroboration only, never tumor-cell prevalence\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
