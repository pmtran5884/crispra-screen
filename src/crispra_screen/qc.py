"""The `qc` command: everything you need to judge a library from one FASTQ.

Produces a tidy CSV bundle plus a multi-page PDF report (matplotlib PdfPages,
no extra dependency) covering read fate, primer-stagger usage, guide detection,
library skew, within-gene guide evenness, depth-subsampling saturation, and a
table of up to N genes with no or only one detected guide.
"""
from __future__ import annotations

import json
import os
import textwrap

import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

from . import figures as F
from .core import FATE_LABELS, count_fastq, load_key
from .depth import (SHALLOW_READS, detection_power, rarefaction,
                    saturation_report, subsample_curve)


def run_qc(fastq, outdir, key_path=None, sample_name="library",
           font="Arial", font_size=10.0, n_missing=10, scan_window=0,
           both_strands=True, depth_check="auto",
           depth_threshold=SHALLOW_READS, log=print):
    os.makedirs(outdir, exist_ok=True)
    key = load_key(key_path)
    key_index = {g: i for i, g in enumerate(key["Guide_Seq"])}

    log(f"[qc] counting {sample_name} from {os.path.basename(fastq)} ...")
    counts, qc, unmatched = count_fastq(fastq, key_index, len(key), log=log,
                                       scan_window=scan_window,
                                       search_both_strands=both_strands)
    tot = qc["total_reads"]
    log(f"[qc] {tot:,} reads, {qc['mapped_reads']:,} mapped "
        f"({qc['mapped_reads']/tot:.1%}), {qc['distinct_guides_detected']:,} guides")

    # ---- read-fate / summary CSV ----
    summ = {"sample": sample_name, "total_reads": tot,
            "reads_assigned_to_key": qc["mapped_reads"],
            "reads_assigned_pct": round(qc["mapped_reads"] / tot * 100, 3),
            "distinct_guides_detected": qc["distinct_guides_detected"],
            "guides_in_key": len(key),
            "anchor_scan_window": qc["scan_window"]["window"],
            "reverse_strand_reads": qc.get("reverse_strand_reads", 0)}
    for k, lab in FATE_LABELS.items():
        summ[f"{lab} (n)"] = qc["fate"].get(k, 0)
        summ[f"{lab} (%)"] = round(qc["fate"].get(k, 0) / tot * 100, 3)
    for p in ["px52_F1", "px52_F2", "px52_F3", "px52_F4", "other/unassigned"]:
        summ[f"{p} (%)"] = round(qc["stagger"].get(p, 0) / tot * 100, 3)
    pd.DataFrame([summ]).to_csv(os.path.join(outdir, "qc_read_summary.csv"), index=False)
    json.dump(qc, open(os.path.join(outdir, "qc_raw.json"), "w"), indent=1)

    # ---- per-guide counts ----
    gl = key[["Guide_Seq", "Gene"]].copy()
    gl["count"] = counts
    gl.to_csv(os.path.join(outdir, "qc_guide_counts.csv"), index=False)

    # ---- genes with no / one detected guide ----
    real = gl[~gl.Gene.astype(str).str.startswith(("NONTARGET", "NO_SITE"))]
    per = real.groupby("Gene").agg(n_guides_in_library=("count", "size"),
                                   n_guides_detected=("count", lambda v: int((v > 0).sum())),
                                   reads_total=("count", "sum"))
    per["n_guides_missing"] = per.n_guides_in_library - per.n_guides_detected
    missing = per[per.n_guides_detected <= 1].sort_values(
        ["n_guides_detected", "reads_total"]).reset_index()
    missing.to_csv(os.path.join(outdir, "qc_missing_guide_genes.csv"), index=False)

    # ---- depth subsampling ----
    log("[qc] depth subsampling ...")
    depth = subsample_curve(counts)
    depth.to_csv(os.path.join(outdir, "qc_depth_subsampling.csv"), index=False)

    # ---- shallow-run diagnostics: is another lane worth paying for? ----
    shallow = (depth_check == "always" or
               (depth_check == "auto" and tot < depth_threshold))
    rep = rare = power = None
    if shallow and qc["mapped_reads"] > 0:
        log(f"[qc] {tot:,} reads (< {depth_threshold:,}): adding saturation "
            "and detection-power diagnostics ...")
        rare = rarefaction(counts)
        rep = saturation_report(counts, rare)
        gpg = int(np.median(key.groupby("Gene").size())) or 4
        power = detection_power(qc["mapped_reads"], guides_per_gene=gpg)
        rare.to_csv(os.path.join(outdir, "qc_rarefaction.csv"), index=False)
        power.to_csv(os.path.join(outdir, "qc_detection_power.csv"), index=False)
        pd.DataFrame([rep]).to_csv(os.path.join(outdir, "qc_saturation.csv"), index=False)
        log(f"[qc] richness {rep['guides_observed']:,} observed / "
            f"{rep['chao1_estimate']:,.0f} estimated "
            f"({rep['pct_of_estimated_recovered']:.0f}%) — {rep['verdict']}")
        log(f"[qc] {rep['advice']}")
        at90 = power[(power.depth_multiplier == 1) & (power.p_detected >= 0.9)]
        if len(at90):
            log(f"[qc] 90% detection down to {at90.gene_freq.min()*100:.2g}% "
                "of sorted cells at this depth")

    # ---- multi-page PDF report ----
    F.style(font, font_size)
    pdf_path = os.path.join(outdir, "qc_report.pdf")
    with PdfPages(pdf_path) as pdf:
        # page 1: header + key numbers
        fig = plt.figure(figsize=(8.27, 11.69))  # A4 portrait
        fig.text(0.07, 0.95, f"CRISPRa library QC — {sample_name}", fontsize=14, weight="bold")
        lines = [
            f"FASTQ: {os.path.basename(fastq)}",
            f"Total reads: {tot:,}",
            f"Reads assigned to key: {qc['mapped_reads']:,} ({qc['mapped_reads']/tot:.1%})",
            f"On-structure amplicon (anchor+scaffold): "
            f"{(qc['fate'].get('mapped_scaffold_ok',0)+qc['fate'].get('unmapped_scaffold_ok',0))/tot:.1%}",
            f"Contamination (no anchor, no scaffold): "
            f"{qc['fate'].get('no_anchor_no_scaffold',0)/tot:.2%}",
            f"Distinct guides detected: {qc['distinct_guides_detected']:,} of {len(key):,}",
            f"Genes with 0 detected guides: {int((per.n_guides_detected==0).sum()):,}",
            f"Genes with exactly 1 detected guide: {int((per.n_guides_detected==1).sum()):,}",
        ]
        if rep is not None:
            lines += [
                "",
                f"Richness (Chao1): {rep['guides_observed']:,} of "
                f"{rep['chao1_estimate']:,.0f} estimated "
                f"({rep['pct_of_estimated_recovered']:.0f}%) — {rep['verdict'].upper()}",
                f"Marginal yield: {rep['new_guides_per_1k_extra_reads']:,.1f} new "
                f"sgRNAs per 1,000 extra reads "
                f"({rep['reads_per_observed_guide']:,.1f} reads per sgRNA)",
            ]
            lines += ["  " + s for s in textwrap.wrap(rep["advice"], 92)]
        for i, ln in enumerate(lines):
            fig.text(0.09, 0.89 - i * 0.028, ln, fontsize=10)
        # missing-guide table
        fig.text(0.07, 0.60, f"Up to {n_missing} genes with no or only one detected guide:",
                 fontsize=10, weight="bold")
        tbl = missing.head(n_missing)[["Gene", "n_guides_in_library",
                                       "n_guides_detected", "reads_total"]]
        ax = fig.add_axes([0.07, 0.30, 0.86, 0.27]); ax.axis("off")
        t = ax.table(cellText=tbl.values, colLabels=tbl.columns, loc="upper left",
                     cellLoc="left", colLoc="left")
        t.auto_set_font_size(False); t.set_fontsize(9); t.scale(1, 1.3)
        pdf.savefig(fig); plt.close(fig)

        # page 2: read fate + stagger + histogram + lorenz
        fig, axes = plt.subplots(2, 2, figsize=(8.27, 8.0))
        F.panel_read_fate(axes[0, 0], {sample_name: qc})
        axes[0, 1].bar(range(5), [qc["stagger"].get(p, 0) / tot * 100 for p in
                                  ["px52_F1", "px52_F2", "px52_F3", "px52_F4", "other/unassigned"]],
                       color=F.BLUE)
        axes[0, 1].set_xticks(range(5), ["F1", "F2", "F3", "F4", "other"])
        axes[0, 1].set_ylabel("reads (% of total)")
        axes[0, 1].set_title("Forward-primer usage", loc="left")
        F.panel_count_hist(axes[1, 0], counts, sample_name)
        F.panel_lorenz(axes[1, 1], counts, sample_name)
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

        # page 3: evenness + saturation
        fig, axes = plt.subplots(1, 2, figsize=(8.27, 4.0))
        F.panel_evenness(axes[0], key, counts)
        F.panel_saturation(axes[1], depth, "library")
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

        # page 4 (shallow runs only): is another lane worth paying for?
        if rep is not None:
            fig, axes = plt.subplots(2, 2, figsize=(8.27, 7.6))
            F.panel_rarefaction(axes[0, 0], rare, rep)
            F.panel_singletons(axes[0, 1], rare)
            F.panel_power(axes[1, 0], power, saturated=rep["verdict"] == "saturated")
            axes[1, 1].axis("off")
            axes[1, 1].text(0, 1, "Verdict\n\n" + "\n".join(textwrap.wrap(
                rep["advice"], 44)) + "\n\n" + "\n".join(
                f"{k.replace('_', ' ')}: {v:,}" if isinstance(v, (int, float))
                else f"{k.replace('_', ' ')}: {v}"
                for k, v in rep.items() if k != "advice"),
                va="top", ha="left", fontsize=8.5, transform=axes[1, 1].transAxes)
            fig.suptitle(f"Is deeper sequencing worth it? — {sample_name}",
                         x=0.02, ha="left")
            fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

    log(f"[qc] wrote {pdf_path} and qc_*.csv in {outdir}")
    return {"qc": qc, "report": pdf_path, "missing": missing, "depth": depth,
            "saturation": rep, "rarefaction": rare, "power": power}
