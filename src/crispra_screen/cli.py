"""Command-line entry point: `crispra-screen qc|score ...`."""
from __future__ import annotations

import argparse
import sys

from . import __version__


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="crispra-screen",
        description="QC and hypergeometric gene scoring for a px52 CRISPRa screen.")
    p.add_argument("--version", action="version", version=f"crispra-screen {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    q = sub.add_parser("qc", help="full QC of one library FASTQ")
    q.add_argument("--fastq", required=True, help="library R1 FASTQ(.gz)")
    q.add_argument("--key", default=None, help="sgRNA key CSV (default: bundled CRISPRa key)")
    q.add_argument("--name", default="library", help="sample name for labels/outputs")
    q.add_argument("--outdir", default="qc_out")
    q.add_argument("--n-missing", type=int, default=10,
                   help="max genes with 0/1 guide to table in the report")
    q.add_argument("--scan-window", type=int, default=0,
                   help="how far into read 1 the U6 anchor may start; "
                        "0 (default) sizes it from the data, which handles "
                        "libraries that retain the full 5' adapter")
    q.add_argument("--font", default="Arial")
    q.add_argument("--font-size", type=float, default=10.0)

    s = sub.add_parser("score", help="fold change + hypergeometric for a treatment arm")
    s.add_argument("--treatment", required=True, help="treatment R1 FASTQ(.gz)")
    g = s.add_mutually_exclusive_group(required=True)
    g.add_argument("--mock", help="mock R1 FASTQ(.gz) to count")
    g.add_argument("--use-bundled-mock", action="store_true",
                   help="use the pre-counted Mock reference shipped with the package")
    s.add_argument("--method", choices=["local", "apron"], default="local",
                   help="local (default, no network) or the Broad GPP portal")
    s.add_argument("--key", default=None, help="sgRNA key CSV (default: bundled CRISPRa key)")
    s.add_argument("--name", default="treatment", help="treatment sample name")
    s.add_argument("--outdir", default="score_out")
    s.add_argument("--decorated", action="store_true",
                   help="keep legend/axis-titles/title on the volcano (default: bare)")
    s.add_argument("--dot-scale", type=float, default=2.0)
    s.add_argument("--label-style", choices=["column", "near"], default="column",
                   help="column: labels pinned to the axis edges, evenly spaced "
                        "(survives shrinking into a figure panel). near: labels "
                        "beside their own points (better on a slide)")
    s.add_argument("--n-label-pos", type=int, default=3)
    s.add_argument("--n-label-neg", type=int, default=3)
    s.add_argument("--width-mm", type=float, default=210.0)
    s.add_argument("--height-mm", type=float, default=130.0)
    s.add_argument("--scan-window", type=int, default=0,
                   help="how far into read 1 the U6 anchor may start; "
                        "0 (default) sizes it from the data, which handles "
                        "libraries that retain the full 5' adapter")
    s.add_argument("--font", default="Arial")
    s.add_argument("--font-size", type=float, default=10.0)

    a = p.parse_args(argv)
    if a.cmd == "qc":
        from .qc import run_qc
        run_qc(a.fastq, a.outdir, key_path=a.key, sample_name=a.name,
               font=a.font, font_size=a.font_size, n_missing=a.n_missing,
               scan_window=a.scan_window)
    elif a.cmd == "score":
        from .score_cmd import run_score
        run_score(a.treatment, a.outdir, mock=a.mock, use_bundled_mock=a.use_bundled_mock,
                  key_path=a.key, method=a.method, treat_name=a.name,
                  font=a.font, font_size=a.font_size, n_pos=a.n_label_pos,
                  n_neg=a.n_label_neg, bare=not a.decorated, dot_scale=a.dot_scale,
                  width_mm=a.width_mm, height_mm=a.height_mm,
                  scan_window=a.scan_window, label_style=a.label_style)
    return 0


if __name__ == "__main__":
    sys.exit(main())
