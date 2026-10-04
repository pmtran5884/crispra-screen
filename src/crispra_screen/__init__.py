"""crispra_screen — read QC, sgRNA counting and hypergeometric gene scoring
for a px52-amplified CRISPR-activation screen (LayV-G / Langya virus G, A375).

Two one-command entry points, both reachable as `crispra-screen qc ...` and
`crispra-screen score ...`:

    qc     full QC of one library FASTQ against the bundled CRISPRa key
    score  fold change + hypergeometric for a treatment FASTQ vs Mock
           (local scorer, or the Broad GPP portal / Apron)

The analysis is the one validated against the Apron reference output.
"""
__version__ = "1.0.0"

from .core import (count_fastq, hypergeom_gene_scores, load_bundled_mock_counts,
                   load_key, score_screen)

__all__ = ["count_fastq", "hypergeom_gene_scores", "load_bundled_mock_counts",
           "load_key", "score_screen", "__version__"]
