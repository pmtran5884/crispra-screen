"""The `qc` command: everything you need to judge a library from one FASTQ.

Produces a tidy CSV bundle plus a multi-page PDF report (matplotlib PdfPages,
no extra dependency) covering read fate, primer-stagger usage, guide detection,
library skew, within-gene guide evenness, depth-subsampling saturation, and a
table of up to N genes with no or only one detected guide.
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

from . import figures as F
from .core import FATE_LABELS, count_fastq, load_key
from .depth import subsample_curve


def run_qc(fastq, outdir, key_path=None, sample_name="library",
           font="Arial", font_size=10.0, n_missing=10, scan_window=0, log=print):
    os.makedirs(outdir, exist_ok=True)
    key = load_key(key_path)
    key_index = {g: i for i, g in enumerate(key["Guide_Seq"])}

    log(f"[qc] counting {sample_name} from {os.path.basename(fastq)} ...")
    counts, qc, unmatched = count_fastq(fastq, key_index, len(key), log=log,
                                       scan_window=scan_window)
    tot = qc["total_reads"]
    log(f"[qc] {tot:,} reads, {qc['mapped_reads']:,} mapped "
        f"({qc['mapped_reads']/tot:.1%}), {qc['distinct_guides_detected']:,} guides")

    # ---- read-fate / summary CSV ----
    summ = {"sample": sample_name, "total_reads": tot,
            "reads_assigned_to_key": qc["mapped_reads"],
            "reads_assigned_pct": round(qc["mapped_reads"] / tot * 100, 3),
            "distinct_guides_detected": qc["distinct_guides_detected"],
            "guides_in_key": len(key),
            "anchor_scan_window": qc["scan_window"]["window"]}
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

    log(f"[qc] wrote {pdf_path} and qc_*.csv in {outdir}")
    return {"qc": qc, "report": pdf_path, "missing": missing, "depth": depth}
