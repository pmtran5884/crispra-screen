"""Client for the Broad GPP pooled-screen portal (Apron).

Builds the .chip and data (log-fold-change) files, submits the Hypergeometric
analysis, polls for completion and parses the result. Network access to
portals.broadinstitute.org is required; callers should fall back to the local
scorer when it is unreachable.
"""
from __future__ import annotations

import html
import os
import re
import time

import numpy as np
import pandas as pd

URL = "https://portals.broadinstitute.org/gpp/public/analysis-tools/crispr-gene-scoring"
RESULTS = URL + "-results"
BASE = "https://portals.broadinstitute.org"


class ApronUnavailable(RuntimeError):
    pass


def write_inputs(key: pd.DataFrame, lfc_by_arm: dict[str, pd.Series], outdir: str):
    """Write the .chip (barcode->gene, controls as NO_SITE_) and the LFC data file."""
    os.makedirs(outdir, exist_ok=True)
    sym = key["Gene"].astype(str).str.replace(r"^NONTARGET_", "NO_SITE_", regex=True)
    chip = pd.DataFrame({"Barcode Sequence": key["Guide_Seq"], "Gene Symbol": sym})
    chip_path = os.path.join(outdir, "screen.chip")
    chip.to_csv(chip_path, sep="\t", index=False)
    data = pd.DataFrame({"Construct Barcode": key["Guide_Seq"]})
    for arm, lfc in lfc_by_arm.items():
        data[arm] = np.asarray(lfc)
    data = data.dropna()
    data_path = os.path.join(outdir, "screen_input.txt")
    data.to_csv(data_path, sep="\t", index=False, float_format="%.9f")
    n_ctrl = int(sym.str.startswith("NO_SITE").sum())
    avg = chip[~sym.str.startswith("NO_SITE")].groupby("Gene Symbol").size().mean()
    return chip_path, data_path, dict(n_barcodes=len(chip), n_nosite=n_ctrl, avg_guides=avg)


def submit(chip_path, data_path, meta, directionality="both", top_pct=100,
           dummy_size=4, min_guides=2, max_guides=4, genes_to_label=5):
    import requests
    s = requests.Session(); s.headers.update({"Referer": URL})
    files = {"chipfile-input": (os.path.basename(chip_path), open(chip_path, "rb"), "text/plain"),
             "datafile-input": (os.path.basename(data_path), open(data_path, "rb"), "text/plain")}
    form = {"clonepool_type": "strict", "clonepool_json_strict": "", "clonepool_json_lax": "",
            "analysistype": "Hypergeometric", "directionality": directionality,
            "percentincluded": str(top_pct), "dummygenesetsize": str(dummy_size),
            "barcodespergenemin": str(min_guides), "barcodespergenemax": str(max_guides),
            "genestohighlight": str(genes_to_label), "displaynositecontrols": "on",
            "numbarcodes": str(meta["n_barcodes"]), "numbarcodespergene": f"{meta['avg_guides']:.2f}",
            "numnositebarcodes": str(meta["n_nosite"]), "numoneintergenicsitebarcodes": "0",
            "threshold": "10", "includefirst": "Yes", "userclickedsubmit": "true"}
    try:
        r = s.post(URL, data=form, files=files, timeout=600)
    except Exception as e:  # noqa
        raise ApronUnavailable(f"could not reach the GPP portal: {e}")
    m = re.search(r"reqid=([^&\"\s]+)", r.url)
    if not m:
        raise ApronUnavailable("portal did not return a job id")
    reqid = m.group(1).replace("%3A", ":")
    return s, reqid


def poll(session, reqid, outdir, timeout_s=5400, interval_s=30, max_transient=5,
         log=print):
    os.makedirs(outdir, exist_ok=True)
    deadline = time.time() + timeout_s
    transient = 0
    while time.time() < deadline:
        # A network failure mid-poll must surface as ApronUnavailable so the
        # caller can fall back, not escape as a raw requests exception. The
        # portal can also blip, so tolerate a few consecutive failures first.
        try:
            r = session.get(RESULTS, params={"reqid": reqid}, timeout=120)
        except Exception as e:  # noqa: BLE001 - any transport error
            transient += 1
            if transient > max_transient:
                raise ApronUnavailable(
                    f"lost contact with the GPP portal while waiting for job "
                    f"{reqid} ({e})")
            log(f"[apron] portal unreachable ({transient}/{max_transient}), retrying ...")
            time.sleep(interval_s)
            continue
        transient = 0
        # Exclude only the portal's own documentation/example downloads, by
        # their actual names. A bare "sample" substring test would also throw
        # away real results whenever a user names an arm "sample_1".
        links = [l for l in sorted(set(re.findall(r'href="([^"]*download[^"]*\.txt[^"]*)"', r.text)))
                 if "README" not in l and "crispr_gene_scoring_sample" not in l]
        if links:
            saved = []
            for l in links:
                u = l if l.startswith("http") else BASE + l
                fn = re.search(r"filename=([^&]+)", u)
                fn = fn.group(1) if fn else os.path.basename(u.split("?")[0])
                try:
                    g = session.get(u, timeout=600)
                except Exception as e:  # noqa: BLE001
                    raise ApronUnavailable(f"could not download {fn} ({e})")
                if g.ok and len(g.content) > 100:
                    path = os.path.join(outdir, fn)
                    open(path, "wb").write(g.content)
                    saved.append(path)
            if not saved:
                raise ApronUnavailable("portal returned result links but no usable files")
            return saved
        log(f"[apron] job {reqid} processing ...")
        time.sleep(interval_s)
    raise ApronUnavailable("timed out waiting for the portal")


def _mean_of_list(series):
    return np.array([np.mean([float(x) for x in str(v).split(";") if x not in ("", "nan")])
                     if isinstance(v, str) else np.nan for v in series])


def parse_output(path: str) -> pd.DataFrame:
    """Parse one Apron hypergeometric gene table into our scoring schema."""
    a = pd.read_csv(path, sep="\t").rename(columns={
        "Gene Symbol": "gene", "Average LFC": "mean_lfc",
        "Average -log(p-values)": "score", "Number of perturbations": "n_guides"})
    asc = _mean_of_list(a["Individual ascending -log(p-values)"])
    desc = _mean_of_list(a["Individual descending -log(p-values)"])
    a["direction"] = np.where(desc >= asc, "enriched", "depleted")
    a["p_value"] = np.power(10.0, -a["score"])
    a["is_control"] = a["gene"].astype(str).str.startswith(("NO_SITE", "DUMMY"))
    from .core import bh
    a["fdr"] = np.nan
    for d in ("enriched", "depleted"):
        m = a.direction == d
        a.loc[m, "fdr"] = bh(a.loc[m, "p_value"].values)
    return a.set_index("gene")
