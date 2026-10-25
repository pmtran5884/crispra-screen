"""Main power sweep: cells sorted x read depth x enrichment.

The hit gene is drawn fresh each replicate from the pool of 4-guide targeting
genes whose guides all clear the reference floor, so the result describes a
typical surfaceome gene rather than the particular library representation of
any one named receptor.
"""
import itertools, os, sys
import numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import simulate as S

ab = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "mock_abundance.csv"))
names = pd.unique(ab.Gene); gmap = {g: i for i, g in enumerate(names)}
LIB = S.Library(share=(ab.mock_counts / ab.mock_counts.sum()).values,
                gene_id=ab.Gene.map(gmap).values,
                mock_counts=ab.mock_counts.values.astype(np.int64),
                is_control=np.array([str(n).startswith("CTRL_pseudo_") for n in names]),
                gene_names=list(names))
GUIDES = {}
for b, gid in enumerate(LIB.block_gene):
    GUIDES[b] = np.flatnonzero(LIB.gene_id == gid)
POOL = np.array([b for b in range(len(LIB.block_gene))
                 if LIB.sizes[b] == 4 and not LIB.block_is_control[b]
                 and (LIB.mock_counts[GUIDES[b]] >= 1).all()])

CELLS = [100, 300, 1000, 3000, 10000]
READS = [3_000, 10_000, 30_000, 1_000_000]
ENRICH = [30, 100, 300, 1000, 3000]
REPS = 200

def job(args):
    ci, ns, nr, E = args
    rng = np.random.default_rng(90000 + ci)
    out = []
    for i in range(REPS):
        b = int(POOL[rng.integers(len(POOL))])
        out.append(S.run_one(LIB, b, GUIDES[b], E, ns, nr, rng))
    d = pd.DataFrame(out)
    return dict(n_sorted=ns, n_reads=nr, E=E, reps=REPS,
                detected=d.detected.mean(), rank1=d.rank1.mean(),
                constructs_med=d.constructs_sorted.median(),
                hit_reads_med=d.hit_reads.median(),
                guides_seen_med=d.hit_guides_seen.median(),
                guides_seen_ge2=(d.hit_guides_seen >= 2).mean(),
                fdr_med=d.hit_fdr.median(),
                rank_med=d.hit_rank.replace(-1, np.nan).median(),
                n_sig_med=d.n_sig_enriched.median(),
                false_pos_rate=(d.n_false_pos_controls > 0).mean())

if __name__ == "__main__":
    grid = [(i, ns, nr, E) for i, (ns, nr, E) in
            enumerate(itertools.product(CELLS, READS, ENRICH))]
    print(f"{len(grid)} conditions x {REPS} reps = {len(grid)*REPS:,} simulated screens",
          flush=True)
    with Pool(10) as p:
        res = p.map(job, grid, chunksize=1)
    pd.DataFrame(res).to_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                          "power_sweep.csv"), index=False)
    print("done", flush=True)
