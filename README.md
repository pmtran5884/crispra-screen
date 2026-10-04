# crispra-screen

QC and hypergeometric gene scoring for a **px52-amplified CRISPR-activation
screen**, from raw sgRNA FASTQ to a QC report or a scored volcano plot in one
command.

Two commands, no setup required:

```bash
crispra-screen qc    --fastq my_library_R1.fastq.gz --outdir qc_out
crispra-screen score --treatment treated_R1.fastq.gz --mock control_R1.fastq.gz --outdir score_out
```

The gene-level statistic is the rank-based hypergeometric used by the Broad
Institute GPP pooled-screen tool ("Apron"), and this implementation reproduces
the average `-log10(P)` in GPP's published reference output for 98.9 % of genes
to within 1e-6, with identical top-200 gene ranking. `score --method apron` can
also submit to the portal directly instead of scoring locally.

---

## Quick start — nothing installed? Pick one

### Option A — your browser only (easiest)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pmtran5884/crispra-screen/blob/main/notebooks/run_in_colab.ipynb)

Click the badge. It installs everything on Google's free servers, runs both
commands on the bundled example data, and shows the QC report, hit list and
volcano. Then upload your own FASTQ and change the paths. No Python, no
downloads, no account beyond a Google login.

### Option B — Docker (needs nothing but Docker, not even Python)

```bash
docker build -t crispra-screen https://github.com/pmtran5884/crispra-screen.git

# mount your working folder at /data; outputs land back there
docker run --rm -v "$PWD":/data crispra-screen \
    qc --fastq /data/my_library_R1.fastq.gz --outdir /data/qc_out
```

### Option C — you have Python, nothing else

```bash
git clone https://github.com/pmtran5884/crispra-screen.git
cd crispra-screen
./run.sh qc --fastq my_library_R1.fastq.gz --outdir qc_out
```

`run.sh` builds its own local environment the first time and reuses it after.

<details>
<summary>Installing it properly instead</summary>

```bash
pip install .                      # then use the `crispra-screen` command
# conda users:
conda env create -f environment.yml && conda activate crispra-screen && pip install .
```
</details>

---

## `qc` — judge a library from one FASTQ

```bash
crispra-screen qc --fastq my_library_R1.fastq.gz --name my_library --outdir qc_out
```

Uses the sgRNA key bundled with the package by default, so you supply only the
FASTQ. Pass `--key my_key.csv` for a different library (a CSV with `Guide_Seq`
and `Gene` columns).

| output | contents |
| --- | --- |
| `qc_report.pdf` | 3-page report: read fate, primer usage, reads-per-sgRNA histogram, library skew, within-gene guide evenness, depth-subsampling saturation, and a table of up to 10 genes with no or only one detected guide |
| `qc_read_summary.csv` | every percentage in the report as one tidy row |
| `qc_guide_counts.csv` | reads per sgRNA |
| `qc_missing_guide_genes.csv` | every gene with ≤ 1 detected guide |
| `qc_depth_subsampling.csv` | guides detected and count correlation vs sequencing depth |

**What the read-fate numbers mean.** Reads are anchored on the constant U6
sequence `GGAAAGGACGAAACACCG` within the first 40 nt (not at a fixed offset,
because the four px52 forward primers carry 0–3 nt staggers); the next 20 nt is
the barcode, and the tracrRNA scaffold immediately 3' of it is checked as an
independent confirmation that the read is the intended amplicon. A read with
neither the anchor nor the scaffold is contamination.

## `score` — fold change and hypergeometric for a treated arm

```bash
# against your own control library
crispra-screen score --treatment treated_R1.fastq.gz --mock control_R1.fastq.gz \
    --name MyTreatment --outdir score_out

# or against the pre-counted reference library bundled with the package
crispra-screen score --treatment treated_R1.fastq.gz --use-bundled-mock \
    --name MyTreatment --outdir score_out

# or let the Broad GPP portal do the scoring
crispra-screen score --treatment treated_R1.fastq.gz --mock control_R1.fastq.gz \
    --method apron --name MyTreatment --outdir score_out
```

Counts are normalised as `log2(reads per million + 1)` and a fold change is
taken per sgRNA against the control. sgRNAs with no reads in the control are
dropped before ranking: their fold change is undefined, and in a bottlenecked
arm a 0-vs-0 guide otherwise ranks above genuinely depleted ones.

Genes are then ranked by the hypergeometric: for a gene with *n* guides at
ranks r₁ < … < rₙ among all *N* guides, the *i*-th best guide gets the
hypergeometric probability mass `P(X = i)`, and the gene score is the mean
`-log10(p_i)`. Both directions are scored and the more significant kept.
Benjamini–Hochberg FDR is added per direction (the GPP hypergeometric output
itself carries no multiple-testing correction).

Outputs: the full per-gene table (CSV) and the volcano as **PDF + SVG + PNG**,
210 × 130 mm at Arial 10 with no raster elements, so every point and label stays
editable in Illustrator. `--decorated` restores the legend and axis titles,
`--dot-scale`, `--n-label-pos`/`--n-label-neg`, `--width-mm`/`--height-mm` and
`--font-size` adjust the rest.

---

## Example data

`example_data/` holds two **synthetic** FASTQs (40,000 reads each) generated by
`tools/make_example_data.py`. Each read is assembled from the px52 amplicon
model so the QC output resembles a real run (~91 % of reads assign to the key,
~7 % are on-structure amplicons whose barcode is absent from the key, ~2 % are
contamination), but no real sequencing read is included. Guide abundances come
from a real control library's distribution *shape*, shuffled across guides, and
the treated file plants a strong enrichment on four arbitrary genes (`ACIN1`,
`APLP2`, `B4GALT1`, `CNTN1`) so the `score` demo produces a realistic volcano.
Regenerate them with:

```bash
python tools/make_example_data.py
```

## Requirements

Python ≥ 3.9 with numpy, pandas, scipy, matplotlib and requests — all installed
automatically by any of the three quick-start paths. `--method apron` needs
network access to `portals.broadinstitute.org`; everything else runs offline.

## Citation

If this is useful, please cite the Broad Institute GPP pooled-screen analysis
tool for the gene-scoring statistic:
<https://portals.broadinstitute.org/gpp/public/analysis-tools/crispr-gene-scoring>

## License

MIT — see `LICENSE`.
