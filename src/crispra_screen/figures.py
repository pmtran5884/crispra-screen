"""Plotting for the pipeline: the QC panels and the volcano.

All figures are emitted as fully vector PDF/SVG (plus PNG) with a single Arial
type size and zero raster image XObjects, so they drop straight into a figure
layout and every element stays editable in Illustrator.
"""
from __future__ import annotations

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox

from .figure_style import apply_figure_style, panel_letter

RED, BLUE, ORANGE, GREY, DARK = "#c8102e", "#1b5e8c", "#b3581c", "#9a9a9a", "#2d2d2d"
MM = 25.4


def style(font: str = "Arial", size: float = 10.0):
    apply_figure_style(frame="open", sizes=(size, size, size), font=font)
    mpl.rcParams.update({
        "mathtext.fontset": "custom", "mathtext.rm": font,
        "mathtext.it": f"{font}:italic", "mathtext.bf": f"{font}:bold",
        "mathtext.default": "regular", "font.size": size,
        "axes.titlesize": size, "axes.labelsize": size, "legend.fontsize": size,
        "xtick.labelsize": size, "ytick.labelsize": size,
    })


def save(fig, stem, formats, width_mm=None, height_mm=None):
    box = None
    if width_mm and height_mm:
        box = Bbox([[0, 0], [width_mm / MM, height_mm / MM]])
    for ext in formats:
        if box is not None:
            fig.savefig(f"{stem}.{ext}", dpi=300, bbox_inches=box, pad_inches=0)
        else:
            fig.savefig(f"{stem}.{ext}", dpi=300, bbox_inches="tight")


def _spread(items, x_col, y_lo, y_hi):
    items = sorted(items, key=lambda t: t[2])
    n = len(items)
    ys = [0.5 * (y_lo + y_hi)] if n == 1 else \
        [y_lo + i * (y_hi - y_lo) / (n - 1) for i in range(n)]
    return [(it, x_col, y) for it, y in zip(items, ys)]


def volcano(res, title, stem, formats=("png", "pdf", "svg"), fdr_cut=0.1,
            n_pos=3, n_neg=3, bare=False, dot_scale=2.0,
            width_mm=210.0, height_mm=130.0, font="Arial", font_size=10.0,
            score_col="score", label_span=(0.42, 0.97)):
    """One-arm volcano from a gene-score table (local or Apron-derived).

    `res` needs columns: mean_lfc, <score_col>, direction, fdr, is_control.
    bare=True drops legend/axis-titles/title (tick labels kept). Labels: the
    top n_pos/n_neg genes by score within each fold-change sign, spread evenly
    down each side so they stay apart when the figure is scaled down.
    """
    style(font, font_size)
    fig, ax = plt.subplots(figsize=(width_mm / MM, height_mm / MM), layout="constrained")
    y = res[score_col].values
    x = res["mean_lfc"].values
    ctrl = res["is_control"].values
    sig = res["fdr"].values < fdr_cut
    up = (~ctrl) & sig & (res["direction"].values == "enriched")
    dn = (~ctrl) & sig & (res["direction"].values == "depleted")
    ax.scatter(x[~ctrl & ~up & ~dn], y[~ctrl & ~up & ~dn], s=7 * dot_scale,
               c="#c9c9c9", linewidths=0)
    ax.scatter(x[ctrl], y[ctrl], s=14 * dot_scale, facecolors="none",
               edgecolors="#8a8a8a", linewidths=0.6, label="non-targeting control")
    # significance is highlighted in BOTH directions, each in its own colour, so
    # the key is satisfied by every plotted point and a labelled gene can never
    # be coloured as something its point is not
    if dn.any():
        ax.scatter(x[dn], y[dn], s=34 * dot_scale, c=BLUE, linewidths=0,
                   label=f"depleted, FDR < {fdr_cut:g}")
    if up.any():
        ax.scatter(x[up], y[up], s=34 * dot_scale, c=RED, linewidths=0,
                   label=f"enriched, FDR < {fdr_cut:g}")

    real = res[~res.is_control]
    pos = real[real.mean_lfc > 0].nlargest(n_pos, score_col)
    neg = real[real.mean_lfc < 0].nlargest(n_neg, score_col)
    ax.margins(0.08)
    x0, x1 = ax.get_xlim()
    ax.set_xlim(x0 - 0.10 * (x1 - x0), x1 + 0.16 * (x1 - x0))
    y0, y1 = ax.get_ylim()
    ax.set_ylim(y0, y1 + ((0.04 if bare else 0.14) * (y1 - y0)))
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    lo, hi = y0 + label_span[0] * (y1 - y0), y0 + label_span[1] * (y1 - y0)
    placed = []
    for sub, side in ((pos, "right"), (neg, "left")):
        items = [(g, r.mean_lfc, r[score_col], r.fdr) for g, r in sub.iterrows()]
        if not items:
            continue
        xc = (x1 - 0.015 * (x1 - x0)) if side == "right" else (x0 + 0.015 * (x1 - x0))
        placed += [(it, c, yy, side) for it, c, yy in _spread(items, xc, lo, hi)]
    for (g, gx, gy, fdr), xc, yy, side in placed:
        col = ("#444444" if fdr >= fdr_cut
               else (RED if res.loc[g, "direction"] == "enriched" else BLUE))
        ax.annotate(g, (gx, gy), xytext=(xc, yy), textcoords="data", va="center",
                    ha="right" if side == "right" else "left",
                    color=col,
                    arrowprops=dict(arrowstyle="-", lw=0.4, color="#b0b0b0",
                                    shrinkA=0, shrinkB=3))
    if not bare:
        ax.set_xlabel("mean sgRNA log$_2$ fold change vs Mock")
        ax.set_ylabel("average $-$log$_{10}$ $P$")
        ax.set_title(title, loc="left")
        ax.legend(frameon=False, loc="upper left", handletextpad=0.4,
                  borderpad=0.2, labelspacing=0.4)
    save(fig, stem, formats, width_mm, height_mm)
    plt.close(fig)
    return f"{stem}.{formats[0]}"


# ---- QC panels (rendered onto a caller-supplied axis) ----------------------
def panel_read_fate(ax, qc_by_sample):
    from .core import FATE_LABELS
    samples = list(qc_by_sample)
    order = list(FATE_LABELS)
    colors = ["#1b5e8c", "#7fa8c4", "#9e9e9e", "#c4c4c4", "#dcdcdc", "#f0f0f0"]
    yy = np.arange(len(samples))[::-1]
    left = np.zeros(len(samples))
    for k, c in zip(order, colors):
        v = np.array([qc_by_sample[s]["fate"].get(k, 0) / qc_by_sample[s]["total_reads"] * 100
                      for s in samples])
        ax.barh(yy, v, left=left, color=c, height=0.6, edgecolor="white",
                linewidth=0.4, label=FATE_LABELS[k])
        left += v
    ax.set_yticks(yy, samples)
    ax.set_xlabel("reads (% of total)")
    ax.set_xlim(0, 100)
    ax.set_title("Read fate", loc="left")
    ax.legend(frameon=False, fontsize=6, loc="upper center",
              bbox_to_anchor=(0.5, -0.30), ncol=1, handlelength=1.0, labelspacing=0.3)


def panel_count_hist(ax, counts, label):
    c = np.asarray(counts)
    bins = np.linspace(0, np.log10(c.max() + 1), 50)
    ax.hist(np.log10(c + 1), bins=bins, color=BLUE, alpha=0.85)
    med = np.median(c[c > 0])
    ax.axvline(np.log10(med + 1), color=RED, lw=0.8, ls="--")
    ax.set_xticks([0, 1, 2, 3, 4], ["0", "10", "100", "1k", "10k"])
    ax.set_xlabel("reads per sgRNA")
    ax.set_ylabel("sgRNAs")
    ax.set_title(f"{label}: median {med:.0f}/guide", loc="left")


def panel_lorenz(ax, counts, label):
    v = np.sort(np.asarray(counts))[::-1]
    cum = np.cumsum(v) / v.sum() * 100
    ax.plot(np.arange(1, len(v) + 1) / len(v) * 100, cum, lw=1.1, color=DARK)
    ax.plot([0, 100], [0, 100], color="#bbbbbb", lw=0.6, ls=":")
    ax.set_xlabel("ranked sgRNAs (% of library)")
    ax.set_ylabel("cumulative reads (%)")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 103)
    ax.set_title(f"{label}: {int((v>0).sum()):,} guides detected", loc="left")


def panel_evenness(ax, key, counts):
    g = key[["Gene"]].copy()
    g["count"] = counts
    g = g[~g.Gene.astype(str).str.startswith(("NONTARGET", "NO_SITE"))]
    n = g.groupby("Gene")["count"].transform("size")
    g = g[n == 4].copy()
    tot = g.groupby("Gene")["count"].transform("sum")
    g["share"] = np.where(tot > 0, g["count"] / tot, np.nan)
    g["rank"] = g.groupby("Gene")["share"].rank(ascending=False, method="first")
    data = [g.loc[g["rank"] == r, "share"].dropna().values for r in (1, 2, 3, 4)]
    bp = ax.boxplot(data, positions=[1, 2, 3, 4], widths=0.55, showfliers=False,
                    patch_artist=True, medianprops=dict(color="white", lw=1.0))
    for p in bp["boxes"]:
        p.set_facecolor(BLUE)
        p.set_edgecolor(BLUE)
    for key_ in ("whiskers", "caps"):
        for p in bp[key_]:
            p.set_color(BLUE)
    ax.axhline(0.25, color=RED, lw=0.8, ls="--")
    ax.set_xticks([1, 2, 3, 4], ["1st", "2nd", "3rd", "4th"])
    ax.set_xlabel("guide rank within its gene")
    ax.set_ylabel("share of gene's reads")
    ax.set_ylim(0, 1)
    ax.set_title("Within-gene evenness (even = 0.25)", loc="left")


def panel_saturation(ax, depth_df, sample):
    d = depth_df[depth_df["sample"] == sample].groupby("frac").mean(numeric_only=True).reset_index()
    ax.plot(d.reads / 1e6, d.guides_detected, "o-", ms=3, lw=1.0, color=DARK)
    ax.set_xscale("log")
    ax.set_xlabel("sequenced reads (millions)")
    ax.set_ylabel("sgRNAs detected")
    ax.set_title("Detection saturation (subsampling)", loc="left")
