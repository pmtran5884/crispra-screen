"""The `score` command: fold change + hypergeometric for a treatment arm vs
Mock, by the local scorer (default, no network) or by the Broad GPP portal
(Apron). Emits the full gene list and a bare vector volcano.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd

from . import apron as A
from . import figures as F
from .core import (count_fastq, load_bundled_mock_counts, load_key, lognorm,
                   score_screen)


def _counts_for(key, key_index, fastq, name, log):
    log(f"[score] counting {name} from {os.path.basename(fastq)} ...")
    counts, qc, _ = count_fastq(fastq, key_index, len(key), log=log)
    log(f"[score]   {qc['total_reads']:,} reads, {qc['mapped_reads']:,} mapped")
    return counts


def run_score(treatment, outdir, mock=None, use_bundled_mock=False, key_path=None,
              method="local", treat_name="treatment", font="Arial", font_size=10.0,
              n_pos=3, n_neg=3, bare=True, dot_scale=2.0,
              width_mm=210.0, height_mm=130.0, log=print):
    os.makedirs(outdir, exist_ok=True)
    key = load_key(key_path)
    key_index = {g: i for i, g in enumerate(key["Guide_Seq"])}

    treat_counts = _counts_for(key, key_index, treatment, treat_name, log)

    if use_bundled_mock:
        log("[score] using the bundled Mock reference counts")
        ref = load_bundled_mock_counts().set_index("Guide_Seq")["count"]
        mock_counts = ref.reindex(key["Guide_Seq"]).fillna(0).values.astype(np.int64)
    elif mock:
        mock_counts = _counts_for(key, key_index, mock, "mock", log)
    else:
        raise ValueError("provide --mock FASTQ or --use-bundled-mock")

    stem = os.path.join("", outdir, f"volcano_{treat_name}")
    if method == "local":
        res = score_screen(key, treat_counts, mock_counts, treat_name=treat_name)
        res.to_csv(os.path.join(outdir, f"gene_scores_{treat_name}.csv"))
        score_col = "score"
        title = f"{treat_name} vs Mock (local hypergeometric)"
    elif method == "apron":
        ln_t, ln_m = lognorm(treat_counts), lognorm(mock_counts)
        usable = mock_counts >= 1
        lfc = pd.Series(np.where(usable, ln_t - ln_m, np.nan), index=key.index)
        chip_path, data_path, meta = A.write_inputs(key, {treat_name: lfc}, outdir)
        log("[score] submitting to the GPP portal (Apron) ...")
        try:
            session, reqid = A.submit(chip_path, data_path, meta)
            saved = A.poll(session, reqid, outdir, log=log)
            arm_file = [p for p in saved if treat_name in os.path.basename(p)] or saved
            res = A.parse_output(arm_file[0])
            res.to_csv(os.path.join(outdir, f"apron_gene_scores_{treat_name}.csv"))
        except A.ApronUnavailable as e:
            log(f"[score] Apron unavailable ({e}); falling back to the local scorer")
            res = score_screen(key, treat_counts, mock_counts, treat_name=treat_name)
            res.to_csv(os.path.join(outdir, f"gene_scores_{treat_name}.csv"))
            method = "local"
        score_col = "score"
        title = f"{treat_name} vs Mock (Apron)" if method == "apron" else \
                f"{treat_name} vs Mock (local, Apron unavailable)"
    else:
        raise ValueError("method must be 'local' or 'apron'")

    F.volcano(res, title, stem, fdr_cut=0.1, n_pos=n_pos, n_neg=n_neg, bare=bare,
              dot_scale=dot_scale, width_mm=width_mm, height_mm=height_mm,
              font=font, font_size=font_size, score_col=score_col)
    top = res[(res.direction == "enriched") & (~res.is_control)].nsmallest(5, "p_value")
    log("[score] top enriched genes:")
    for g, r in top.iterrows():
        log(f"    {g:12s} LFC {r.mean_lfc:+.2f}  p={r.p_value:.2g}  FDR={r.fdr:.2g}")
    log(f"[score] wrote gene table and volcano_{treat_name}.(pdf|svg|png) in {outdir}")
    return res
