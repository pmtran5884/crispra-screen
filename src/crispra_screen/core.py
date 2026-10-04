"""Counting and gene-scoring core for the LayV-G CRISPRa pipeline.

Everything here is pure compute with no plotting and no network, so it is easy
to test and reuse. The amplicon model and the hypergeometric statistic are the
ones validated against the Broad GPP (Apron) reference output; see
``crispra_screen.apron`` for the portal client and ``docs`` in the repo for the
validation.
"""
from __future__ import annotations

import collections
import gzip
import math
import os
from importlib import resources

import numpy as np
import pandas as pd
from scipy.stats import hypergeom, rankdata

# ---- amplicon model (px52 forward primers, U6 + tracrRNA vector) -----------
ANCHOR = "GGAAAGGACGAAACACCG"          # constant U6 3' sequence, 5' of the sgRNA
PRIMER_3P = "CTTGTGGAAAGGACGAAACACCG"  # full px52 genome-annealing part
SCAFFOLD = "GTTTAAGAGCTATGCTG"         # tracrRNA, immediately 3' of the sgRNA
STAGGERS = {0: "px52_F1", 1: "px52_F2", 2: "px52_F3", 3: "px52_F4"}
SG_LEN = 20
SCAN_WINDOW = 40                        # anchor must start within the first 40 nt

FATE_LABELS = {
    "mapped_scaffold_ok": "on-structure amplicon, guide in key",
    "mapped_scaffold_bad": "guide in key, scaffold mismatch",
    "unmapped_scaffold_ok": "on-structure amplicon, guide not in key",
    "unmapped_scaffold_bad": "anchor only, scaffold mismatch",
    "no_anchor_scaffold_present": "no U6 anchor, scaffold present",
    "no_anchor_no_scaffold": "neither anchor nor scaffold (contamination)",
}


# ---- bundled reference data -------------------------------------------------
def bundled_path(name: str) -> str:
    """Absolute path to a file shipped inside the package's data/ folder."""
    return str(resources.files("crispra_screen.data").joinpath(name))


def default_key_path() -> str:
    """The CRISPRa sgRNA key shipped with the package."""
    return bundled_path("formatted_sgrna_key.csv")


def load_key(path: str | None = None) -> pd.DataFrame:
    """Load an sgRNA key (Guide_Seq, Gene, [Guide_Num]); default = bundled key."""
    df = pd.read_csv(path or default_key_path())
    df["Guide_Seq"] = df["Guide_Seq"].str.upper().str.strip()
    if df["Guide_Seq"].duplicated().any():
        raise ValueError("duplicate guide sequences in key")
    return df


def load_bundled_mock_counts() -> pd.DataFrame:
    """Pre-counted Mock reference (Guide_Seq, Gene, count) shipped with the package.

    This stands in for the 379 MB raw Mock FASTQ, which is too large to ship; it
    is the exact per-guide Mock count table from the published analysis.
    """
    return pd.read_csv(bundled_path("mock_reference_counts.csv"))


# ---- FASTQ counting ---------------------------------------------------------
def fastq_seqs(path: str):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                yield line.rstrip("\n")


def count_fastq(path: str, key_index: dict[str, int], n_guides: int,
                progress_every: int = 2_000_000, log=print):
    """Count sgRNAs in one R1 FASTQ against a key index.

    Returns (counts ndarray aligned to key order, qc dict, unmatched Counter).
    """
    fate = collections.Counter()
    stagger = collections.Counter()
    counts = np.zeros(n_guides, dtype=np.int64)
    unmatched = collections.Counter()
    lens = collections.Counter()
    n = 0
    for seq in fastq_seqs(path):
        n += 1
        lens[len(seq)] += 1
        p = seq.find(ANCHOR, 0, SCAN_WINDOW)
        if p < 0:
            fate["no_anchor_scaffold_present" if SCAFFOLD[:12] in seq
                 else "no_anchor_no_scaffold"] += 1
            continue
        off = p - (len(PRIMER_3P) - len(ANCHOR))
        stagger[off if 0 <= off <= 3 else -1] += 1
        s = p + len(ANCHOR)
        guide = seq[s:s + SG_LEN]
        if len(guide) < SG_LEN:
            fate["anchor_truncated_read"] += 1
            continue
        scaffold_ok = seq[s + SG_LEN:s + SG_LEN + 10] == SCAFFOLD[:10]
        idx = key_index.get(guide)
        if idx is not None:
            counts[idx] += 1
            fate["mapped_scaffold_ok" if scaffold_ok else "mapped_scaffold_bad"] += 1
        else:
            unmatched[guide] += 1
            fate["unmapped_scaffold_ok" if scaffold_ok else "unmapped_scaffold_bad"] += 1
        if log and progress_every and n % progress_every == 0:
            log(f"  {n/1e6:.0f}M reads")
    qc = {
        "total_reads": n,
        "fate": dict(fate),
        "stagger": {STAGGERS.get(k, "other/unassigned"): v for k, v in stagger.items()},
        "read_len": dict(lens),
        "mapped_reads": int(counts.sum()),
        "distinct_guides_detected": int((counts > 0).sum()),
        "unmatched_20mers_distinct": len(unmatched),
        "unmatched_top": unmatched.most_common(50),
    }
    return counts, qc, unmatched


def count_to_frame(key: pd.DataFrame, counts: np.ndarray, sample: str) -> pd.DataFrame:
    out = key[["Guide_Seq", "Gene"]].copy()
    out[sample] = counts
    return out


# ---- normalisation and gene scoring (validated against Apron) --------------
def lognorm(counts: np.ndarray) -> np.ndarray:
    counts = np.asarray(counts, dtype=float)
    return np.log2(counts / counts.sum() * 1e6 + 1)


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    o = np.argsort(p)
    q = np.empty_like(p)
    m = len(p)
    prev = 1.0
    for rank, idx in enumerate(o[::-1]):
        prev = min(prev, p[idx] * m / (m - rank))
        q[idx] = prev
    return q


def hypergeom_gene_scores(lfc: pd.Series, gene: pd.Series, top_frac: float = 1.0,
                          p_method: str = "pmf") -> pd.DataFrame:
    """Broad GPP (Apron) rank hypergeometric gene score.

    Rank all N guides by LFC (mid-ranks for ties). For a gene with n guides at
    ranks r_1<...<r_n the i-th best guide gets p_i; with p_method='pmf' this is
    the hypergeometric mass P(X=i) that Apron uses, with 'sf' the upper tail
    P(X>=i). Gene score = mean -log10(p_i) over the best top_frac of guides,
    run both directions, keep the more significant.
    """
    N = len(lfc)
    out = {}
    for direction, sign in (("enriched", -1.0), ("depleted", 1.0)):
        rank = np.rint(rankdata(sign * lfc.values, method="average")).astype(np.int64)
        rank = np.clip(rank, 1, N)
        df = pd.DataFrame({"gene": gene.values, "rank": rank})
        rows = {}
        for g, sub in df.groupby("gene", sort=False):
            r = np.sort(sub["rank"].values)
            n = len(r)
            k = max(1, math.ceil(n * top_frac))
            i = np.arange(1, n + 1)
            if p_method == "pmf":
                p = np.clip(hypergeom.pmf(i, N, r, n), 1e-300, 1.0)
            else:
                p = np.clip(hypergeom.sf(i - 1, N, r, n), 1e-300, 1.0)
            rows[g] = float(np.mean(-np.log10(p)[:k]))
        out[direction] = pd.Series(rows)
    res = pd.DataFrame(out)
    res["score"] = res[["enriched", "depleted"]].max(axis=1)
    res["direction"] = res[["enriched", "depleted"]].idxmax(axis=1)
    res["p_value"] = np.power(10.0, -res["score"])
    return res


def score_screen(key: pd.DataFrame, treat_counts: np.ndarray, mock_counts: np.ndarray,
                 treat_name: str = "treatment", min_reference_reads: int = 1,
                 min_guides: int = 2, top_frac: float = 1.0, p_method: str = "pmf",
                 seed: int = 0) -> pd.DataFrame:
    """Full local scoring: LFC vs mock, hypergeometric, BH-FDR, control pseudo-genes.

    Returns a per-gene table sorted by p-value.
    """
    gene = key["Gene"].astype(str).copy()
    nt = gene.str.startswith(("NONTARGET", "NO_SITE"))
    rng = np.random.default_rng(seed)
    idx = np.flatnonzero(nt.values)
    rng.shuffle(idx)
    gene.iloc[idx] = [f"CTRL_pseudo_{i // 4:04d}" for i in range(len(idx))]

    ln_t, ln_m = lognorm(treat_counts), lognorm(mock_counts)
    usable = np.asarray(mock_counts) >= min_reference_reads
    lfc = pd.Series((ln_t - ln_m)[usable])
    gu = gene[usable]

    res = hypergeom_gene_scores(lfc, gu, top_frac=top_frac, p_method=p_method)
    agg = pd.DataFrame({
        "n_guides": gu.groupby(gu).size(),
        "n_guides_in_library": gene.groupby(gene).size(),
        "mean_lfc": lfc.groupby(gu.values).mean(),
        "median_lfc": lfc.groupby(gu.values).median(),
        f"reads_mock": pd.Series(np.asarray(mock_counts)[usable]).groupby(gu.values).sum(),
        f"reads_{treat_name}": pd.Series(np.asarray(treat_counts)[usable]).groupby(gu.values).sum(),
    })
    res = res.join(agg)
    res["is_control"] = res.index.str.startswith("CTRL_pseudo_")
    res = res[res.n_guides >= min_guides].copy()
    res["fdr"] = np.nan
    for d in ("enriched", "depleted"):
        m = res.direction == d
        res.loc[m, "fdr"] = bh(res.loc[m, "p_value"].values)
    res = res.sort_values("p_value")
    res.index.name = "gene"
    return res
