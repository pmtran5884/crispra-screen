#!/usr/bin/env python
"""Generate the SYNTHETIC example FASTQs shipped in example_data/.

These files contain no real sequencing reads. Each read is assembled from the
px52 amplicon model (random adapter prefix, one of the four forward-primer
staggers, the constant U6 anchor, a 20 nt sgRNA from the library key, the
tracrRNA scaffold, then filler), so the QC output looks like a real run, but
the sequence content is generated here.

Guide abundances are drawn from the real Mock count distribution's *shape* and
then shuffled across guides, so no per-guide information from the real library
survives. The treatment file additionally plants a strong enrichment on a few
arbitrarily chosen genes (see PLANTED) so the `score` demo produces a
realistic-looking volcano with hits.

    python tools/make_example_data.py
"""
from __future__ import annotations

import argparse
import gzip
import os

import numpy as np

from crispra_screen.core import (PRIMER_3P, SCAFFOLD, load_bundled_mock_counts,
                                 load_key)

STAGGER_SEQ = ["", "A", "TC", "ACT"]
# arbitrary genes, chosen only because they are ordinary 4-guide entries;
# these are NOT hits of the real screen
PLANTED = ["ACIN1", "APLP2", "B4GALT1", "CNTN1"]


def synth_read(guide: str, rng: np.random.Generator, kind: str = "ok") -> str:
    """kind: 'ok' = in-key guide, 'offkey' = on-structure but 20-mer not in the
    key, 'offstructure' = neither anchor nor scaffold (mimics contamination)."""
    stag = STAGGER_SEQ[rng.integers(0, 4)]
    prefix = "".join(rng.choice(list("ACGT"), size=rng.integers(0, 3)))
    if kind == "ok":
        read = prefix + stag + PRIMER_3P + guide + SCAFFOLD
    elif kind == "offkey":
        rand20 = "".join(rng.choice(list("ACGT"), size=20))
        read = prefix + stag + PRIMER_3P + rand20 + SCAFFOLD
    else:
        read = "".join(rng.choice(list("ACGT"), size=60))
    filler = "".join(rng.choice(list("ACGT"), size=max(0, 151 - len(read))))
    return (read + filler)[:151]


def write_fastq(path, guides, weights, n_reads, rng,
                offkey_frac=0.07, offstructure_frac=0.02):
    """Read-fate mix is chosen to resemble a real run (~91 % assigned to key,
    ~7 % on-structure amplicons whose 20-mer is absent from the key, ~2 %
    neither anchor nor scaffold) so the QC report is instructive."""
    w = np.asarray(weights, dtype=float)
    w = w / w.sum()
    picks = rng.choice(len(guides), size=n_reads, p=w)
    u = rng.random(n_reads)
    with gzip.open(path, "wt") as fh:
        for i, (gi, r) in enumerate(zip(picks, u)):
            kind = ("offstructure" if r < offstructure_frac
                    else "offkey" if r < offstructure_frac + offkey_frac else "ok")
            seq = synth_read(guides[gi], rng, kind=kind)
            fh.write(f"@synthetic_read_{i}\n{seq}\n+\n{'I' * len(seq)}\n")
    return n_reads


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="example_data")
    ap.add_argument("--reads", type=int, default=40000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--boost", type=float, default=120.0,
                    help="abundance multiplier applied to the planted genes")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    key = load_key()
    guides = key["Guide_Seq"].tolist()

    # abundance shape from the real Mock, shuffled so it carries no per-guide info
    shape = load_bundled_mock_counts()["count"].values.astype(float)
    shape = shape + 1.0
    rng.shuffle(shape)

    write_fastq(os.path.join(a.outdir, "example_control_R1.fastq.gz"),
                guides, shape, a.reads, rng)

    boosted = shape.copy()
    planted_mask = key["Gene"].isin(PLANTED).values
    boosted[planted_mask] *= a.boost
    write_fastq(os.path.join(a.outdir, "example_treated_R1.fastq.gz"),
                guides, boosted, a.reads, rng)

    print(f"wrote synthetic examples to {a.outdir}/ "
          f"({a.reads:,} reads each; planted genes: {', '.join(PLANTED)})")


if __name__ == "__main__":
    main()
