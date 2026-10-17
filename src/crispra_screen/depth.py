"""Sequencing-depth subsampling: resample a count vector without replacement
(multivariate hypergeometric) and track how detection and count-correlation
degrade. Answers 'was the depth right?'.

For shallow runs (< SHALLOW_READS mapped reads, i.e. a tightly sorted bin
rather than a whole library) three further questions decide whether to pay
for another lane, and `qc`/`score` answer them automatically:

  saturation_report  -- is there anything left in the tube to sequence?
  detection_power    -- how small a binder could this depth still see?
  hit_stability      -- do the calls move as depth changes?

The three are not interchangeable. A run can be read-saturated (nothing left
to recover) and still underpowered, in which case the lever is cells sorted,
not reads: sequencing deeper samples the same exhausted material.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import binom

from .core import lognorm, score_screen

#: Above this many reads a run is deep enough that the shallow-run diagnostics
#: are noise; below it they are run automatically.
SHALLOW_READS = 100_000


def subsample_curve(counts: np.ndarray, fracs=None, reps: int = 3, seed: int = 0) -> pd.DataFrame:
    counts = np.asarray(counts, dtype=np.int64)
    if fracs is None:
        fracs = [0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 0.85, 1.0]
    full_ln = lognorm(counts)
    rng = np.random.default_rng(seed)
    rows = []
    for f in fracs:
        for rep in range(1 if f == 1.0 else reps):
            if f == 1.0:
                c = counts.copy()
            else:
                n = int(round(counts.sum() * f))
                c = rng.multivariate_hypergeometric(counts, n, method="marginals")
            ln = lognorm(c)
            cpm = c / c.sum() * 1e6
            med = np.median(cpm[c > 0]) if (c > 0).any() else 0.0
            rows.append(dict(
                sample="library", frac=f, rep=rep, reads=int(c.sum()),
                guides_detected=int((c > 0).sum()),
                frac_of_fulldepth_guides=float((c > 0).sum() / max((counts > 0).sum(), 1)),
                r_vs_full=float(np.corrcoef(ln, full_ln)[0, 1]),
                pct_within_10x=float(((cpm > med / 10) & (cpm < med * 10)).mean() * 100),
            ))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# shallow-run diagnostics
# --------------------------------------------------------------------------

def chao1(counts) -> float:
    """Bias-corrected Chao1 estimate of the number of constructs present.

    Uses only singletons and doubletons, so it is an estimate of richness in
    the material that was sequenced -- i.e. how many distinct sgRNAs were in
    the sorted bin, whether or not a read was drawn from each.
    """
    c = np.asarray(counts, dtype=np.int64)
    c = c[c > 0]
    f1, f2 = int((c == 1).sum()), int((c == 2).sum())
    return float(len(c) + (f1 * (f1 - 1)) / (2 * (f2 + 1)))


def rarefaction(counts, fracs=None, reps: int = 3, seed: int = 0) -> pd.DataFrame:
    """Guides / singletons / Chao1 as a function of reads sampled."""
    c = np.asarray(counts, dtype=np.int64)
    if fracs is None:
        fracs = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    rng = np.random.default_rng(seed)
    rows = []
    for f in fracs:
        for rep in range(1 if f == 1.0 else reps):
            s = c.copy() if f == 1.0 else rng.multivariate_hypergeometric(
                c, int(round(c.sum() * f)), method="marginals")
            rows.append(dict(frac=f, rep=rep, reads=int(s.sum()),
                             guides=int((s > 0).sum()),
                             singletons=int((s == 1).sum()),
                             chao1=chao1(s)))
    return pd.DataFrame(rows)


def saturation_report(counts, rare: pd.DataFrame | None = None) -> dict:
    """Is there anything left in this sample to sequence?

    Returns the observed/estimated richness, the marginal yield over the last
    10% of reads, and a verdict keyed to whether more reads can still buy
    constructs.
    """
    c = np.asarray(counts, dtype=np.int64)
    if rare is None:
        rare = rarefaction(c)
    obs, est, total = int((c > 0).sum()), chao1(c), int(c.sum())
    agg = rare.groupby("frac").agg(reads=("reads", "mean"), guides=("guides", "mean"))
    lo = agg.index[agg.index < 1.0].max() if (agg.index < 1.0).any() else None
    if lo is None or total == agg.loc[lo, "reads"]:
        yield_per_1k = float("nan")
    else:
        yield_per_1k = float((obs - agg.loc[lo, "guides"]) /
                             ((total - agg.loc[lo, "reads"]) / 1000))
    frac_rec = obs / est if est else float("nan")
    if frac_rec >= 0.95 and yield_per_1k < 50:
        verdict, advice = "saturated", (
            "Sequencing is not the bottleneck: more reads cannot recover "
            "constructs that are not in the sorted material. To gain power, "
            "sort more cells.")
    elif frac_rec >= 0.85:
        verdict, advice = "near-saturated", (
            "Most of the material has been read out; extra depth would add "
            "mainly singletons, which carry no scoring weight.")
    else:
        verdict, advice = "not saturated", (
            "Constructs remain unsampled: deeper sequencing of this same "
            "library would recover more of them.")
    return {"mapped_reads": total, "guides_observed": obs,
            "chao1_estimate": round(est, 1),
            "pct_of_estimated_recovered": round(frac_rec * 100, 1),
            "singletons": int((c == 1).sum()), "doubletons": int((c == 2).sum()),
            "reads_per_observed_guide": round(total / obs, 2) if obs else float("nan"),
            "new_guides_per_1k_extra_reads": round(yield_per_1k, 1),
            "verdict": verdict, "advice": advice}


def detection_power(n_reads: int, guides_per_gene: int = 4, min_guides: int = 2,
                    freqs=None, multipliers=(1, 4)) -> pd.DataFrame:
    """P(a true binder is detectable) vs its abundance in the sorted bin.

    A gene is scoreable only if at least `min_guides` of its sgRNAs draw a
    read. If a binder's cells are a fraction p of the sorted population and
    its sgRNAs share that fraction evenly, each sgRNA is seen with
    probability 1-(1-p/g)^N; the gene is called with the binomial tail.

    `multipliers` adds hypothetical deeper runs. Those curves are only
    attainable if the sample is NOT saturated -- otherwise extra reads
    resample exhausted material and the real curve stays put.
    """
    if freqs is None:
        freqs = np.logspace(np.log10(2e-4), np.log10(3e-2), 80)
    rows = []
    for m in multipliers:
        N = int(n_reads * m)
        pg = 1 - (1 - np.asarray(freqs) / guides_per_gene) ** N
        rows.append(pd.DataFrame(dict(
            depth_multiplier=m, reads=N, gene_freq=freqs,
            p_guide_seen=pg,
            p_detected=1 - binom.cdf(min_guides - 1, guides_per_gene, pg))))
    return pd.concat(rows, ignore_index=True)


def hit_stability(key, treat_counts, mock_counts, genes, fracs=None,
                  reps: int = 3, seed: int = 0, **score_kw) -> pd.DataFrame:
    """Re-score at reduced depth and track each gene's rank and FDR.

    A call whose FDR is flat across a 4-10x span of depth is limited by the
    rank statistic over its sgRNA log-fold-changes, not by counting
    precision: more reads will not move it.
    """
    tc = np.asarray(treat_counts, dtype=np.int64)
    if fracs is None:
        fracs = [0.1, 0.25, 0.5, 0.75, 1.0]
    rng = np.random.default_rng(seed)
    genes = list(genes)
    rows = []
    for f in fracs:
        for rep in range(1 if f == 1.0 else reps):
            s = tc.copy() if f == 1.0 else rng.multivariate_hypergeometric(
                tc, int(round(tc.sum() * f)), method="marginals")
            res = score_screen(key, s, mock_counts, **score_kw)
            real = res[~res.is_control]
            order = real.score.rank(ascending=False, method="min")
            for g in genes:
                if g not in real.index:
                    continue
                rows.append(dict(gene=g, frac=f, rep=rep, reads=int(s.sum()),
                                 gene_reads=int(s[np.asarray(key.Gene) == g].sum()),
                                 n_guides=int(real.loc[g, "n_guides"]),
                                 rank=int(order.loc[g]),
                                 mean_lfc=float(real.loc[g, "mean_lfc"]),
                                 p_value=float(real.loc[g, "p_value"]),
                                 fdr=float(real.loc[g, "fdr"])))
    return pd.DataFrame(rows)
