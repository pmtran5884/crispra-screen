"""Sequencing-depth subsampling: resample a count vector without replacement
(multivariate hypergeometric) and track how detection and count-correlation
degrade. Answers 'was the depth right?'.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .core import lognorm


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
