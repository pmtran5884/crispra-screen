"""Two-stage sort-then-sequence simulator for pooled CRISPRa receptor screens.

The experiment has two independent sampling steps and they are not
interchangeable:

  1. SORT.  n_sorted cells are drawn from the ~12.5M-cell input pool. This
     fixes the DIVERSITY available downstream -- a construct absent from the
     sorted cells cannot be recovered at any read depth.
  2. SEQUENCE.  n_reads are drawn (with replacement, as PCR amplification
     implies) from the sorted cells' integrated cassettes. This fixes how
     completely that diversity is observed.

Library abundance is inherited from the empirical Mock count distribution
rather than assumed uniform, so library skew, zero-read guides and within-gene
share imbalance all carry through.

`fast_scores` reproduces crispra_screen.core.score_screen exactly for
(top_frac=1.0, p_method='pmf', min_guides=2) but vectorises the per-gene loop,
which is what makes a 20,000-replicate sweep tractable.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.special import gammaln
from scipy.stats import rankdata


# --------------------------------------------------------------------------
# library
# --------------------------------------------------------------------------

@dataclass
class Library:
    share: np.ndarray          # per-guide abundance in the input pool, sums to 1
    gene_id: np.ndarray        # int label per guide (controls already pseudo-genised)
    mock_counts: np.ndarray    # reference read counts
    is_control: np.ndarray     # per-GENE, aligned to the gene_id ordering
    gene_names: list
    min_reference_reads: int = 1
    # precomputed grouping over USABLE guides -------------------------------
    usable: np.ndarray = field(init=False)
    order0: np.ndarray = field(init=False)
    starts: np.ndarray = field(init=False)
    sizes: np.ndarray = field(init=False)
    i_flat: np.ndarray = field(init=False)
    n_flat: np.ndarray = field(init=False)
    keep: np.ndarray = field(init=False)
    N: int = field(init=False)

    def __post_init__(self):
        self.usable = np.asarray(self.mock_counts) >= self.min_reference_reads
        gu = self.gene_id[self.usable]
        self.N = int(self.usable.sum())
        # contiguous gene blocks over the usable guides
        self.order0 = np.argsort(gu, kind="stable")
        g_sorted = gu[self.order0]
        self.starts = np.flatnonzero(np.r_[True, g_sorted[1:] != g_sorted[:-1]])
        self.sizes = np.diff(np.r_[self.starts, len(g_sorted)])
        self.block_gene = g_sorted[self.starts]
        self.i_flat = (np.arange(len(g_sorted)) - np.repeat(self.starts, self.sizes) + 1)
        self.n_flat = np.repeat(self.sizes, self.sizes)
        self.keep = self.sizes >= 2          # min_guides
        self.block_names = [self.gene_names[g] for g in self.block_gene]
        self.block_is_control = np.array([self.is_control[g] for g in self.block_gene])


def _hypergeom_pmf(i, N, r, n):
    """P(X=i) for the GPP rank hypergeometric, in log space.

    scipy.stats.hypergeom.pmf is ~260x slower on arrays this size and that
    cost dominates the whole simulation, so the mass is evaluated directly
    from log-gammas. Agreement with scipy is ~1e-11 absolute, far below any
    decision boundary once scores are averaged and compared.
    """
    ok = (i <= r) & ((n - i) <= (N - r))
    out = np.zeros(len(i))
    ii, rr, nn = i[ok], r[ok], n[ok]
    out[ok] = np.exp(
        gammaln(rr + 1) - gammaln(ii + 1) - gammaln(rr - ii + 1)
        + gammaln(N - rr + 1) - gammaln(nn - ii + 1) - gammaln(N - rr - nn + ii + 1)
        - (gammaln(N + 1) - gammaln(nn + 1) - gammaln(N - nn + 1)))
    return out


def _lognorm(c):
    c = np.asarray(c, dtype=float)
    return np.log2(c / c.sum() * 1e6 + 1)


def _bh(p):
    o = np.argsort(p); q = np.empty_like(p); m = len(p); prev = 1.0
    for rank, idx in enumerate(o[::-1]):
        prev = min(prev, p[idx] * m / (m - rank))
        q[idx] = prev
    return q


def fast_scores(lib: Library, treat_counts: np.ndarray):
    """Vectorised equivalent of score_screen. Returns per-BLOCK arrays."""
    lfc = (_lognorm(treat_counts) - _lognorm(lib.mock_counts))[lib.usable]
    gu_sorted_gene = lib.gene_id[lib.usable][lib.order0]
    sc = {}
    for direction, sign in (("enriched", -1.0), ("depleted", 1.0)):
        rank = np.rint(rankdata(sign * lfc, method="average")).astype(np.int64)
        np.clip(rank, 1, lib.N, out=rank)
        r = rank[lib.order0]
        # sort ranks WITHIN each gene block (lexsort: gene major, rank minor)
        o = np.lexsort((r, gu_sorted_gene))
        r = r[o]
        p = np.clip(_hypergeom_pmf(lib.i_flat, lib.N, r, lib.n_flat), 1e-300, 1.0)
        s = -np.log10(p)
        sc[direction] = np.add.reduceat(s, lib.starts) / lib.sizes
    score = np.maximum(sc["enriched"], sc["depleted"])
    direction = np.where(sc["enriched"] >= sc["depleted"], 0, 1)   # 0=enriched
    pval = np.power(10.0, -score)
    fdr = np.full(len(score), np.nan)
    for d in (0, 1):
        msk = lib.keep & (direction == d)
        if msk.any():
            fdr[msk] = _bh(pval[msk])
    return score, pval, fdr, direction


# --------------------------------------------------------------------------
# simulation
# --------------------------------------------------------------------------

def simulate_counts(lib: Library, hit_guides: np.ndarray, E: float,
                    n_sorted: int, n_reads: int, rng) -> np.ndarray:
    """Draw n_sorted cells from the E-enriched pool, then n_reads from those cells."""
    p = lib.share.copy()
    p[hit_guides] *= E
    p /= p.sum()
    cells = rng.multinomial(n_sorted, p)                 # stage 1: the sort
    tot = cells.sum()
    if tot == 0:
        return np.zeros(len(p), dtype=np.int64)
    return rng.multinomial(n_reads, cells / tot)         # stage 2: sequencing


def run_one(lib: Library, hit_block: int, hit_guides: np.ndarray, E: float,
            n_sorted: int, n_reads: int, rng, fdr_cut: float = 0.1) -> dict:
    counts = simulate_counts(lib, hit_guides, E, n_sorted, n_reads, rng)
    score, pval, fdr, direction = fast_scores(lib, counts)
    real = lib.keep & (~lib.block_is_control)
    hit_ok = real[hit_block]
    if hit_ok:
        rank = int((score[real] > score[hit_block]).sum()) + 1
        hit_fdr = float(fdr[hit_block]); hit_dir = int(direction[hit_block])
    else:
        rank, hit_fdr, hit_dir = -1, np.nan, -1
    fp = int(((fdr < fdr_cut) & lib.keep & lib.block_is_control & (direction == 0)).sum())
    return dict(E=E, n_sorted=n_sorted, n_reads=n_reads,
                constructs_sorted=int((counts > 0).sum()),
                hit_reads=int(counts[hit_guides].sum()),
                hit_guides_seen=int((counts[hit_guides] > 0).sum()),
                hit_rank=rank, hit_fdr=hit_fdr,
                detected=bool(hit_ok and hit_dir == 0 and hit_fdr < fdr_cut),
                rank1=bool(hit_ok and rank == 1),
                n_sig_enriched=int(((fdr < fdr_cut) & lib.keep & (direction == 0)).sum()),
                n_false_pos_controls=fp)


# --------------------------------------------------------------------------
# library design: synthetic libraries at a fixed input-cell budget
# --------------------------------------------------------------------------

def build_synthetic(n_genes: int, guides_per_gene: int, n_controls: int,
                    abundance_pool: np.ndarray, rng, ref_reads: int = 12_440_964):
    """A library of arbitrary shape whose guide abundances are bootstrapped
    from the empirical Mock distribution, so skew and the zero-read fraction
    carry over to designs that were never actually built.

    Controls are pre-grouped into 4-guide pseudo-genes, matching what the
    production scorer does to non-targeting guides.
    """
    n_t = n_genes * guides_per_gene
    n_ctrl_blocks = n_controls // 4
    n = n_t + n_ctrl_blocks * 4
    counts = rng.choice(abundance_pool, size=n, replace=True)
    if counts.sum() == 0:
        counts = counts + 1
    share = counts / counts.sum()
    gene_id = np.r_[np.repeat(np.arange(n_genes), guides_per_gene),
                    n_genes + np.repeat(np.arange(n_ctrl_blocks), 4)]
    names = ([f"GENE_{i:05d}" for i in range(n_genes)]
             + [f"CTRL_pseudo_{i:04d}" for i in range(n_ctrl_blocks)])
    is_control = np.r_[np.zeros(n_genes, bool), np.ones(n_ctrl_blocks, bool)]
    mock = rng.multinomial(ref_reads, share)
    return Library(share=share, gene_id=gene_id, mock_counts=mock.astype(np.int64),
                   is_control=is_control, gene_names=names)


def simulate_counts_pooled(lib: Library, hit_guides: np.ndarray, E: float,
                           n_input_cells: int, n_sorted: int, n_reads: int, rng):
    """As simulate_counts, but the input cell pool is drawn explicitly first.

    At a fixed cell budget, library size sets coverage: a smaller library is
    represented by more cells per guide, so its input pool is a less noisy
    realisation of the intended library. That is the effect this stage adds --
    it is invisible if the sort is drawn straight from the design shares.
    """
    pool = rng.multinomial(n_input_cells, lib.share).astype(float)
    p = pool.copy()
    p[hit_guides] *= E
    if p.sum() == 0:
        return np.zeros(len(p), dtype=np.int64)
    p /= p.sum()
    cells = rng.multinomial(n_sorted, p)
    tot = cells.sum()
    if tot == 0:
        return np.zeros(len(p), dtype=np.int64)
    return rng.multinomial(n_reads, cells / tot)


# --------------------------------------------------------------------------
# the `design` command
# --------------------------------------------------------------------------

def power(cells_sorted, reads, enrichment, n_genes=1500, guides_per_gene=6,
          n_controls=500, input_cells=12_500_000, reps=200, seed=0,
          abundance_pool=None, log=print):
    """Detection probability for one design, at one sort size and read depth.

    Answers 'is this experiment worth running, and do I need deep sequencing?'
    before any cells are spent. Enrichment is the one parameter that cannot be
    known in advance -- pass the range you consider plausible and read the
    curve, rather than trusting any single number.
    """
    import pandas as pd
    if abundance_pool is None:
        from .core import load_bundled_mock_counts
        abundance_pool = load_bundled_mock_counts()["count"].values
    rng = np.random.default_rng(seed)
    lib = build_synthetic(n_genes, guides_per_gene, n_controls, abundance_pool,
                          np.random.default_rng(seed + 1))
    guides = {b: np.flatnonzero(lib.gene_id == g) for b, g in enumerate(lib.block_gene)}
    elig = np.array([b for b in range(len(lib.block_gene))
                     if lib.sizes[b] == guides_per_gene and not lib.block_is_control[b]
                     and (lib.mock_counts[guides[b]] >= 1).all()])
    rows = []
    for E in np.atleast_1d(enrichment):
        det = 0
        for _ in range(reps):
            b = int(elig[rng.integers(len(elig))])
            c = simulate_counts_pooled(lib, guides[b], float(E), input_cells,
                                       cells_sorted, reads, rng)
            score, pval, fdr, direction = fast_scores(lib, c)
            real = lib.keep & (~lib.block_is_control)
            det += bool(real[b] and direction[b] == 0 and fdr[b] < 0.1)
        rows.append(dict(n_genes=n_genes, guides_per_gene=guides_per_gene,
                         n_guides=len(lib.share),
                         coverage=round(input_cells / len(lib.share)),
                         cells_sorted=cells_sorted, reads=reads, enrichment=float(E),
                         reps=reps, detected=det / reps))
        log(f"[design] {E:>7,.0f}x enrichment -> {det/reps:5.1%} detected "
            f"({cells_sorted:,} cells, {reads:,} reads)")
    return pd.DataFrame(rows)
