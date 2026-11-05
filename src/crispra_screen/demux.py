"""The `demux` command: split one FASTQ into per-sample files by an inline barcode.

Pooling two sorts into one sequencing submission is free whenever the read
budget already exceeds what the bins need (an NGS lane against a few thousand
sorted cells is ~100x past saturation), so the only question is whether the
reads can be assigned back cleanly.

The barcode is expected inline in the amplicon, either side of the constant
U6 region:

    before:   [stagger] BARCODE CTTGTGGAAAGGACGAAACACCG [sgRNA]
    between:  [stagger] CTTGTGGAAAGGACGAAACACCG BARCODE [sgRNA]

Which one is in use is detected from the data rather than declared, because
getting it wrong silently shifts every guide by the barcode length. In the
`between` layout the barcode is EXCISED on write, so the emitted FASTQ has
the standard architecture and every downstream command works unchanged.

Reads carrying the anchor in neither orientation cannot be assigned and are
written to a separate file rather than silently dropped.
"""
from __future__ import annotations

import collections
import os

from .core import ANCHOR, PRIMER_3P, SG_LEN, fastq_seqs, revcomp

#: offset of the constant region's 5' end relative to the anchor hit
_CONST_LEAD = len(PRIMER_3P) - len(ANCHOR)


def _hamming_le(a: str, b: str, k: int) -> bool:
    if len(a) != len(b):
        return False
    d = 0
    for x, y in zip(a, b):
        if x != y:
            d += 1
            if d > k:
                return False
    return True


def fastq_records(path):
    """Yield (header, seq, plus, qual) so quality strings survive the split."""
    import gzip
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            s = fh.readline().strip()
            p = fh.readline()
            q = fh.readline().strip()
            yield h.rstrip("\n"), s, p.rstrip("\n"), q


def detect_layout(path, barcode, window=0, sample=20000, both_strands=True,
                  min_frac=0.05, specificity=3.0):
    """Work out where the barcode sits relative to the anchor, from the data.

    A real inline barcode occupies ONE offset. A string that merely resembles
    background sequence (a homopolymer against a poly-A adapter, say) matches
    at many offsets about equally, so the candidate offset is required to beat
    the surrounding offsets by `specificity`x before it is believed. Refusing
    is the right failure here: mistaking background for a barcode splits the
    pool on noise, and in the 'between' layout it also excises 6 real bases
    from every read it touches.
    """
    bl = len(barcode)
    w = window if window and window > 0 else None
    prof = collections.Counter()
    n_anchor = 0
    offs = list(range(-_CONST_LEAD - bl - 4, len(ANCHOR) + bl + 5))
    for i, seq in enumerate(fastq_seqs(path)):
        if i >= sample:
            break
        s = seq
        p = s.find(ANCHOR, 0, w) if w else s.find(ANCHOR)
        if p < 0 and both_strands:
            s = revcomp(seq)
            p = s.find(ANCHOR, 0, w) if w else s.find(ANCHOR)
        if p < 0:
            continue
        n_anchor += 1
        for o in offs:
            b0 = p + o
            if b0 >= 0 and s[b0:b0 + bl] == barcode:
                prof[o] += 1
    meta = {"anchored_reads": n_anchor}
    if n_anchor == 0:
        return None, meta
    canon = {"before": -_CONST_LEAD - bl, "between": len(ANCHOR)}
    fr = {k: prof.get(v, 0) / n_anchor for k, v in canon.items()}
    meta.update({f"frac_{k}": v for k, v in fr.items()})
    layout = max(fr, key=fr.get)
    best = fr[layout]
    meta["layout"] = layout
    meta["frac_barcoded"] = best
    # background level: the best NON-canonical offset
    bg = max([c / n_anchor for o, c in prof.items() if o not in canon.values()],
             default=0.0)
    meta["background_frac"] = bg
    if best < min_frac:
        meta["reason"] = (f"only {best:.2%} of anchored reads carry the barcode at "
                          f"either expected offset (need {min_frac:.0%})")
        return None, meta
    if bg > 0 and best < specificity * bg:
        meta["reason"] = (f"barcode matches background sequence: {best:.2%} at the "
                          f"'{layout}' offset vs {bg:.2%} at an unrelated offset, "
                          f"not the >{specificity:.0f}x separation a real inline "
                          f"barcode shows")
        return None, meta
    return layout, meta


def run_demux(fastq, barcode, outdir, names=("barcoded", "plain"), mismatches=0,
              layout=None, scan_window=0, both_strands=True, key_path=None,
              log=print):
    os.makedirs(outdir, exist_ok=True)
    barcode = barcode.strip().upper()
    bl = len(barcode)
    if layout in (None, "auto"):
        layout, meta = detect_layout(fastq, barcode, window=scan_window,
                                     both_strands=both_strands)
        if layout is None:
            raise ValueError(
                f"cannot place barcode {barcode}: {meta.get('reason', meta)}. "
                "Check the barcode sequence, or pass --layout explicitly if you "
                "are certain of the architecture.")
        log(f"[demux] barcode {barcode} detected in the '{layout}' position "
            f"({meta['frac_barcoded']:.1%} of anchored reads carry it)")
    else:
        meta = {"layout": layout, "mode": "forced"}

    # guides that would masquerade as barcoded reads in the 'between' layout
    ambiguous = []
    if layout == "between" and key_path is not False:
        from .core import load_key
        k = load_key(key_path)
        ambiguous = [g for g in k["Guide_Seq"]
                     if _hamming_le(g[:bl], barcode, mismatches)]
        if ambiguous:
            log(f"[demux] note: {len(ambiguous)} library guides start with the "
                f"barcode within {mismatches} mismatch(es); unbarcoded reads "
                "carrying them are indistinguishable and are counted as such")

    paths = {names[0]: os.path.join(outdir, f"{names[0]}.fastq"),
             names[1]: os.path.join(outdir, f"{names[1]}.fastq"),
             "unassigned": os.path.join(outdir, "unassigned.fastq")}
    fh = {k: open(v, "w") for k, v in paths.items()}
    n = collections.Counter()
    w = scan_window if scan_window and scan_window > 0 else None
    try:
        for h, seq, plus, qual in fastq_records(fastq):
            n["total"] += 1
            s, q = seq, qual
            p = s.find(ANCHOR, 0, w) if w else s.find(ANCHOR)
            if p < 0 and both_strands:
                s, q = revcomp(seq), qual[::-1]
                p = s.find(ANCHOR, 0, w) if w else s.find(ANCHOR)
            if p < 0:
                n["no_anchor"] += 1
                fh["unassigned"].write(f"{h}\n{seq}\n{plus}\n{qual}\n")
                continue
            b0 = (p - _CONST_LEAD - bl) if layout == "before" else (p + len(ANCHOR))
            obs = s[b0:b0 + bl] if b0 >= 0 else ""
            is_bc = len(obs) == bl and _hamming_le(obs, barcode, mismatches)
            if is_bc and layout == "between":
                s = s[:b0] + s[b0 + bl:]          # excise so the guide sits where
                q = q[:b0] + q[b0 + bl:]          # every downstream step expects it
            tgt = names[0] if is_bc else names[1]
            n[tgt] += 1
            fh[tgt].write(f"{h}\n{s}\n{plus}\n{q}\n")
    finally:
        for f in fh.values():
            f.close()

    tot = max(n["total"], 1)
    log(f"[demux] {n['total']:,} reads -> {names[0]} {n[names[0]]:,} "
        f"({n[names[0]]/tot:.1%}), {names[1]} {n[names[1]]:,} "
        f"({n[names[1]]/tot:.1%}), unassigned {n['no_anchor']:,} "
        f"({n['no_anchor']/tot:.1%})")
    report = {"fastq": os.path.basename(fastq), "barcode": barcode,
              "layout": layout, "mismatches_allowed": mismatches,
              "total_reads": n["total"], f"{names[0]}_reads": n[names[0]],
              f"{names[1]}_reads": n[names[1]], "unassigned_reads": n["no_anchor"],
              "ambiguous_guides": len(ambiguous), **{f"detect_{k}": v for k, v in meta.items()}}
    import pandas as pd
    pd.DataFrame([report]).to_csv(os.path.join(outdir, "demux_report.csv"), index=False)
    if ambiguous:
        pd.DataFrame({"guide": ambiguous}).to_csv(
            os.path.join(outdir, "demux_ambiguous_guides.csv"), index=False)
    return {"paths": paths, "counts": dict(n), "layout": layout, "report": report}
