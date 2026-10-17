"""The `score` command: fold change + hypergeometric for a treatment arm vs
Mock, by the local scorer (default, no network) or by the Broad GPP portal
(Apron). Emits the full gene list and a bare vector volcano.
"""
from __future__ import annotations

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import apron as A
from . import depth as D
from . import figures as F
from .core import (count_fastq, load_bundled_mock_counts, load_key, lognorm,
                   score_screen)


def _counts_for(key, key_index, fastq, name, log, scan_window=0, both_strands=True):
    log(f"[score] counting {name} from {os.path.basename(fastq)} ...")
    counts, qc, _ = count_fastq(fastq, key_index, len(key), log=log,
                                scan_window=scan_window,
                                search_both_strands=both_strands)
    if qc.get("reverse_strand_reads"):
        log(f"[score]   {qc['reverse_strand_reads']:,} reads recovered from the reverse strand")
    log(f"[score]   {qc['total_reads']:,} reads, {qc['mapped_reads']:,} mapped")
    return counts, qc


def run_score(treatment, outdir, mock=None, use_bundled_mock=False, key_path=None,
              method="local", treat_name="treatment", font="Arial", font_size=10.0,
              n_pos=3, n_neg=3, bare=True, dot_scale=2.0,
              width_mm=210.0, height_mm=130.0, scan_window=0,
              label_style="column", both_strands=True,
              depth_check="auto", depth_threshold=D.SHALLOW_READS,
              n_stability_genes=3, log=print):
    os.makedirs(outdir, exist_ok=True)
    key = load_key(key_path)
    key_index = {g: i for i, g in enumerate(key["Guide_Seq"])}

    treat_counts, treat_qc = _counts_for(key, key_index, treatment, treat_name, log,
                                         scan_window=scan_window, both_strands=both_strands)

    if use_bundled_mock:
        log("[score] using the bundled Mock reference counts")
        ref = load_bundled_mock_counts().set_index("Guide_Seq")["count"]
        mock_counts = ref.reindex(key["Guide_Seq"]).fillna(0).values.astype(np.int64)
    elif mock:
        mock_counts, _ = _counts_for(key, key_index, mock, "mock", log,
                                     scan_window=scan_window, both_strands=both_strands)
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
            log(f"[score] portal job id: {reqid}")
            saved = A.poll(session, reqid, outdir, log=log)
            arm_file = [p for p in saved if treat_name in os.path.basename(p)] or saved
            res = A.parse_output(arm_file[0])
            res.to_csv(os.path.join(outdir, f"apron_gene_scores_{treat_name}.csv"))
        except (A.ApronUnavailable, OSError) as e:
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
              font=font, font_size=font_size, score_col=score_col,
              label_style=label_style)
    top = res[(res.direction == "enriched") & (~res.is_control)].nsmallest(5, "p_value")
    log("[score] top enriched genes:")
    for g, r in top.iterrows():
        log(f"    {g:12s} LFC {r.mean_lfc:+.2f}  p={r.p_value:.2g}  FDR={r.fdr:.2g}")

    # ---- shallow-run depth diagnostics ----
    stab = rep = None
    shallow = (depth_check == "always" or
               (depth_check == "auto" and treat_qc["total_reads"] < depth_threshold))
    if shallow and treat_qc["mapped_reads"] > 0:
        log(f"[score] {treat_qc['total_reads']:,} reads (< {depth_threshold:,}): "
            "checking whether depth limits these calls ...")
        rep = D.saturation_report(treat_counts)
        genes = list(top.head(n_stability_genes).index)
        stab = D.hit_stability(key, treat_counts, mock_counts, genes)
        stab.to_csv(os.path.join(outdir, f"depth_hit_stability_{treat_name}.csv"),
                    index=False)
        pd.DataFrame([rep]).to_csv(
            os.path.join(outdir, f"depth_saturation_{treat_name}.csv"), index=False)
        gpg = int(np.median(key.groupby("Gene").size())) or 4
        power = D.detection_power(treat_qc["mapped_reads"], guides_per_gene=gpg)
        rare = D.rarefaction(treat_counts)

        F.style(font, font_size)
        fig, axes = plt.subplots(2, 2, figsize=(10.0, 6.7))
        F.panel_rarefaction(axes[0, 0], rare, rep)
        F.panel_singletons(axes[0, 1], rare)
        F.panel_hit_stability(axes[1, 0], stab)
        F.panel_power(axes[1, 1], power, saturated=rep["verdict"] == "saturated")
        fig.suptitle(f"Is deeper sequencing worth it? — {treat_name} "
                     f"({rep['verdict']})", x=0.02, ha="left")
        fig.tight_layout()
        F.save(fig, os.path.join(outdir, f"depth_check_{treat_name}"),
               ("png", "pdf", "svg"))
        plt.close(fig)

        log(f"[score] richness {rep['guides_observed']:,} observed / "
            f"{rep['chao1_estimate']:,.0f} estimated "
            f"({rep['pct_of_estimated_recovered']:.0f}%) — {rep['verdict']}")
        log(f"[score] {rep['advice']}")
        # The decision is whether the CALL changes, not whether the number does.
        for g, s in stab.groupby("gene"):
            a = s.groupby("frac").fdr.mean()
            lo, hi = a.loc[a.index.min()], a.loc[a.index.max()]
            if hi < 0.1 and lo < 0.1:
                note = "called at every depth tested — depth is not limiting"
            elif hi < 0.1:
                note = "only called at full depth — more reads would firm it up"
            elif a.max() - a.min() < 0.25 * a.max():
                note = ("never called, FDR flat across depth — limited by the rank "
                        "statistic over its sgRNAs, not by read count")
            else:
                note = "never called, but FDR improving with depth — more reads may help"
            log(f"    {g:12s} FDR {lo:.2g} at {a.index.min():.0%} depth "
                f"-> {hi:.2g} at full depth: {note}")
        if rep["verdict"] == "saturated":
            log("[score] note: where a call above would benefit from more reads, "
                "this sample cannot supply them — it is read-saturated, so the "
                "extra information has to come from more cells, not another lane.")

    log(f"[score] wrote gene table and volcano_{treat_name}.(pdf|svg|png) in {outdir}")
    return res
