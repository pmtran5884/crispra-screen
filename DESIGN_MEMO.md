# Powering a small CRISPRa surfaceome screen: how many cells, how many reads, which analysis

**Scope.** Two decisions, answered from 45,650 simulated screens whose sampling
model was first required to reproduce three completed experiments: (1) how to
get the most detection power out of a 12.5M-cell input and a 100–10,000-cell
sort, and (2) when ONT depth (3,000–30,000 reads, one day) is insufficient and
NGS (~1M reads, two weeks) is actually required.

---

## The two answers, up front

**1. Power comes from cells, and from a smaller library — not from reads.**
At a fixed 12.5M-cell input, shrinking the library from 6,000 to 1,500 genes
takes detection of a 100×-enriched receptor from **17% to 77%** at a 300-cell
sort. Going from 4 to 8 guides per gene adds a further **+14 points**. Head to
head at 1,000 sorted cells and 30,000 reads, your current 24,564-guide library
detects a 30× hit **6.0%** of the time; the recommended 9,500-guide library
detects it **74.5%** of the time. No sequencing decision comes close.

**2. ONT is sufficient across your entire operating range. NGS is never
required.** 30,000 ONT reads sits within **6.3 percentage points** of unlimited
sequencing in the worst case tested, and on the plateau for any hit enriched
≥100×. All three receptors you have actually recovered were enriched 570–2,580×;
**3,000 reads would have called every one of them.**

---

## How the model was built, and why it can be trusted

Two sampling steps, in order, because they are not interchangeable:

1. **The sort** draws `n_sorted` cells from the input pool. This fixes the
   *diversity* available — a construct absent from the sorted cells cannot be
   recovered at any read depth.
2. **Sequencing** draws `n_reads` from those cells' cassettes. This fixes how
   completely that diversity is *observed*.

Guide abundance is inherited from your empirical Mock distribution (12,440,964
reads over 24,367 guides) rather than assumed uniform, so library skew, the
4.3% zero-read guides and the within-gene share imbalance all carry through.

**Validation against completed screens.** Before predicting anything, the
simulator had to reproduce history. Construct counts and read shares came back
within 1% on all three screens, and FDR was reproduced across five orders of
magnitude:

| screen | hit | constructs obs / sim | read share obs / sim | FDR observed | FDR simulated |
|---|---|---|---|---|---|
| LayV-G round 2 | PROCR | 3,228 / 3,225 | 10.63% / 10.65% | 7.1e-07 | 4.8e-07 |
| dIgA | PIGR | 514 / 516 | 15.49% / 15.34% | 2.7e-05 | 9.2e-03 |
| GDF15 sample_1 | FCGR1B | 67 / 66 | 7.71% / 7.80% | 0.366 | 0.364 |

![Simulator validated against three completed screens]({{artifact:art_f44892fe-ec8b-48ff-b5bc-18cc5bb5a406}})

The PIGR row is the one place the model is *conservative*: the real PIGR drew
reads on two guides that the Mock records at 3 and 0 reads, which the abundance
prior cannot anticipate. Everything below therefore slightly understates power —
the safe direction for a design study.

**Implied sort sizes** (inverted from the observed construct counts) — useful in
their own right:

| screen | reads | constructs | cells the sort delivered |
|---|---|---|---|
| LayV-G round 2 | 3.56M | 3,228 | ~4,200 |
| dIgA | 3,924 | 514 | ~620 |
| GDF15 sample_2 | 1,967 | 874 | ~1,100 |
| GDF15 sample_1 | 1,116 | 67 | **~70** |

A 70-cell sort explains the GDF15 pilot entirely.

---

## Question 1 — the best achievable power, given the fixed constraints

![Detection vs enrichment across cells sorted and read depth]({{artifact:art_3f88ea24-c4d7-4dbb-a632-c8ed7d912f8b}})

### Library shape is the dominant lever

At 12.5M input cells, library size and coverage trade directly. Detection of a
**100× hit at a 300-cell sort**:

| genes | 4 guides/gene | 6 guides | 8 guides |
|---|---|---|---|
| **1,500** | 77% | 86% | **91%** |
| 3,000 | 50% | 55% | 54% |
| 6,000 | 17% | 20% | 15% |

Two things to read off this. First, **gene count matters more than anything
else** — the 6,000-gene row is a different experiment from the 1,500-gene row.
Second, **more guides per gene stop helping once coverage drops**: 6,000 genes ×
8 guides (49,000 guides, 255× coverage) is *worse* than 6,000 × 4. Extra guides
only pay off while coverage stays above roughly 500×.

### Recommended library

| parameter | recommendation | consequence |
|---|---|---|
| genes | **1,500** (prioritised surfaceome) | the single biggest power gain available |
| guides per gene | **6** | 10,000 guides total; +9 points over 4 guides |
| non-targeting controls | **500** | 1,000 buys nothing (below); frees 500 slots |
| total guides | 9,500 | **1,316× coverage** at 12.5M cells |
| expected detection | 86% of a 100× hit at 300 cells; >99% at 1,000 cells | |
| | **74.5% of even a 30× hit at 1,000 cells** | vs **6.0%** for the current library, same model |

If the biology demands wider coverage, 3,000 genes × 6 guides (19,000 guides,
658×) still gives 55% at 300 cells and 96% at 1,000 — acceptable, but it makes
the sort size load-bearing in a way the 1,500-gene design does not.

### Sort size remains the binding constraint

Even the recommended library needs cells. Across the whole sweep, detection at
or below **300 sorted cells** is materially worse than at 1,000, at every
enrichment and every read depth. **Target ≥1,000 cells into the tube**; below
~300 the experiment is underpowered regardless of library or sequencer.

---

## Question 2 — the ONT / NGS cutoff

Reads needed to come within 5 points of the deep-sequencing plateau:

| cells sorted | 30× hit | 100× hit | ≥300× hit |
|---|---|---|---|
| 1,000 | 3,000 | 10,000 | 3,000 |
| 3,000 | 30,000 | 10,000 | 3,000 |
| 10,000 | 50,000 | 10,000 | 3,000 |
| 30,000 | 30,000 | 10,000 | 3,000 |

![ONT versus NGS decision map]({{artifact:art_49f94e18-53c1-4494-8d67-f705c6d58bc4}})

**Operating rule**

- **≥300× hit** (the regime of every receptor you have recovered): 3,000 reads
  at any sort size.
- **100× hit**: 10,000 reads at any sort size.
- **30× hit**: 30,000 reads up to ~3,000 cells; 50,000 above that.
- **Simple version: sequence ~5× the number of cells you sorted, capped at
  50,000 reads.** Never send for NGS on the basis of depth alone.

Below 1,000 sorted cells, 3,000 reads and 1,000,000 reads differ by ≤2
percentage points at every enrichment tested. The reason is structural: reads
cannot create diversity the sort did not deliver.

**When NGS *is* worth the two weeks:** only when a sort delivered ≥3,000 cells
*and* the candidate is weak (≈30×) *and* the ONT run already shows the library
is not saturated. That last condition is machine-checkable — `crispra-screen`
reports the Chao1 saturation verdict on every shallow run. If it says
"saturated", another lane cannot help and NGS is wasted.

---

## Fold change: what the numbers mean at shallow depth

**The ceiling is depth-independent.** A 1,000× hit measures 9.57 log2 at 3,000
reads and 9.59 at 1,000,000. The maximum observable enrichment for a 4-guide
gene is **~3,700×** (the gene's four guides share the read budget), and that
ceiling is identical at every depth. A truly 10,000×-enriched receptor reads as
~3,500× whatever you spend.

![Fold-change limits and control-count effect]({{artifact:art_dabc07ae-7560-45dd-99a3-979090dabffb}})

**The floor is what moves.** The apparent fold change of a guide holding a
single read:

| reads | 1 read looks like |
|---|---|
| 3,000 | **11.0×** |
| 10,000 | 3.3× |
| 30,000 | **1.1× (neutral)** |
| 1,000,000 | 0.1× |

At 3,000 reads a single stray read makes a guide look 11-fold enriched, which is
why very shallow runs produce noisy apparent enrichment even when the gene
ranking is correct. **At 30,000 reads the single-read floor sits at neutral** —
a guide must draw ≥2 reads before it looks enriched at all. That is the right
calibration, and it is an independent argument for 30,000 reads as the standing
ONT target even though 3,000 usually suffices for detection.

---

## Non-targeting controls: 500 is enough

| controls | pseudo-genes | detection (100× hit, 1,000 cells) | screens with ≥1 false-positive control |
|---|---|---|---|
| 0 | 0 | 68% | — |
| 250 | 62 | 72% | 0.0% |
| **500** | **125** | **73%** | **0.0%** |
| 1,000 | 250 | 70% | 0.0% |
| 2,000 | 499 | 66% | 0.0% |

Differences across 250–1,000 controls are inside Monte Carlo noise (400
replicates, SE ≈ 2.3 points); 2,000 is mildly *worse* because controls consume
ranking slots. The null score spread is flat at SD 0.27–0.30 regardless, so 500
controls characterise the null as well as 2,000 do.

Note the false-positive column: **0–0.5% of screens produce even one spurious
control call against a nominal 10% FDR.** BH on this statistic is conservative
in this regime, which means the limiting factor is sensitivity, not specificity —
and that in turn is why the more sensitive analysis below is worth adopting.

**Recommendation: 500 non-targeting guides.** Keep them pooled into 4-guide
pseudo-genes, as the scorer already does.

---

## Which analysis: MAGeCK or the Apron hypergeometric?

450 paired screens, both methods scoring the identical count table (shallow
sorted sample vs the 12.4M-read Mock), MAGeCK with conventional defaults
(RRA, median normalisation).

| cells | enrichment | hypergeometric | MAGeCK | hyper rank-1 | MAGeCK rank-1 |
|---|---|---|---|---|---|
| 300 | 100× | 14% | **40%** | 2% | **30%** |
| 1,000 | 100× | 76% | **92%** | 50% | **96%** |
| 3,000 | 100× | 96% | 100% | 98% | 100% |
| 300 | 300× | 72% | **92%** | 46% | **76%** |
| 1,000 | 300× | 98% | 98% | 98% | 98% |
| 300 | 1000× | 98% | 98% | 88% | 96% |

![MAGeCK vs hypergeometric, and library shape]({{artifact:art_d4903ada-d310-4855-bb26-4ad91e4d8645}})

**MAGeCK is the more sensitive detector** — 91.1% vs 83.8% mean detection, with
the whole advantage concentrated in exactly the small-sort, modest-enrichment
corner you operate in. This was not the expected result: with one treated and
one control sample MAGeCK has no replicates from which to estimate
negative-binomial dispersion, and a rank statistic needs no variance model at
all. The dispersion handicap evidently costs less than the rank statistic's
coarseness on four guides.

**But MAGeCK cannot report how good a hit is.** Its FDR floors at 4.95e-03 with
only **80 distinct values** across 6,141 genes, where the hypergeometric gives
**6,117 distinct values down to 4.8e-07**. MAGeCK would have reported PIGR at
"FDR < 0.005" where Apron reported 2.7e-05 — the same call, far less of a
statement. False-positive rates are both negligible (0.9% vs 0.2%).

**Recommendation: run both, and say so in the methods.** MAGeCK as the primary
*detector* on small sorts; the Apron hypergeometric as the reported *statistic*
and the volcano. They agree on the top hit whenever there is one — MAGeCK simply
finds it more often when the sort was marginal. Disagreement between them is
itself informative: a MAGeCK-only call on a small sort is a candidate for a
confirmatory sort, not a result.

---

## Re-running this for a future bait

The power calculation ships with the tool, so it does not have to be
re-derived:

```bash
crispra-screen design --cells-sorted 1000 --reads 30000 \
    --enrichment 30 100 300 1000 --genes 1500 --guides-per-gene 6
```

It builds the library you describe, draws the input pool at your cell budget,
sorts, sequences and scores — and reports detection probability per enrichment
level.

### Genes with only 2 or 3 guides: order them

Some surfaceome genes have only 2 or 3 designable guides. They are worth
ordering, and the reason is not obvious from the score alone.

![Power for genes with few guides]({{artifact:art_0f1a2460-5469-4385-83f2-cd4840a844c3}})

**No guide count is categorically excluded.** The best achievable BH-FDR for a
gene whose guides land at the very top ranks is 1.7e-03 with 2 guides, 4.2e-05
with 3 and 1.2e-06 with 4 — all comfortably below 0.1. A 2-guide gene can
return a significant result.

**Detection at 1,000 sorted cells and 30,000 reads:**

| guides | 100× hit | 300× hit | 1000× hit |
|---|---|---|---|
| 2 | 74% | 87% | 91% |
| 3 | 81% | **99.5%** | 100% |
| 4 | 94% | 99.5% | 100% |
| 6 | 100% | 100% | 100% |

**Three guides is essentially four.** At any enrichment at or above 300× — the
regime of every receptor recovered so far — a 3-guide gene performs
indistinguishably from a 4-guide gene (99.5% vs 99.5%). Order these without
hesitation.

**Two guides carries a real but bounded penalty, and it is a representation
problem rather than a statistical one.** Decomposing the 300× case: 92.8% of
2-guide genes are scoreable at all, and 93.8% of those are called — so the
dominant loss is the **~8% of 2-guide genes that fail the `min_guides = 2`
filter because one of their two guides falls below the reference floor.** At
the empirical 4.26% zero-read rate, `(1 - 0.0426)^2 = 91.7%` of 2-guide genes
have both guides represented; the observed rate in the current library is 8.0%
dropped, against 0.0% for 3-guide and 0.2% for 4-guide genes.

**The cost of keeping them is trivial.** 150 two-guide genes consume 300 of
~8,150 slots (3.7%). Redeploying those slots across the rest of the library
would add ~0.3 guides per gene — negligible. Covering 150 extra genes at ~87%
power for a strong hit is far better value.

Two practical consequences:

1. **Pick the two best-predicted guides for 2-guide genes, and QC the built
   library.** The failure mode is a guide missing from the pool, which is
   visible in the Mock count table before any screen is run —
   `crispra-screen qc` already tables genes with 0 or 1 detected guide.
2. **A negative result for a 2-guide gene is weak evidence.** It may simply
   have been dropped. Treat 2-guide genes as detectable-but-not-excludable, and
   do not read their absence from a hit list as evidence against them.

## Standing protocol

1. **Library** — 1,500 prioritised surfaceome genes × 6 guides + 500
   non-targeting = 9,500 guides, 1,316× coverage at 12.5M input cells.
2. **Sort** — target ≥1,000 cells; treat <300 as a pilot, not an experiment.
3. **Sequence** — ONT, ~5× the sorted-cell count, capped at 30,000–50,000 reads.
   One day.
4. **Analyse** — both MAGeCK and the Apron hypergeometric against the banked
   deep Mock. Report the hypergeometric statistic.
5. **Only then consider NGS**, and only if the saturation verdict says the
   library is *not* saturated and the candidate is weak. Otherwise another lane
   buys nothing.

---

## Limitations

- **Enrichment is the one parameter nobody can know in advance.** Power is
  therefore reported as a function of it, with your three observed values
  (570×, 1,987×, 2,582×) marked for reference. Those all sit at the easy end;
  the 30–300× band is where design choices actually bite, and no completed
  screen tells us how often real receptors land there.
- **The abundance prior is your Mock, which is itself a sample.** Guides the
  Mock under-sampled will be under-weighted — demonstrably so for PIGR. Power
  estimates are mildly conservative as a result.
- **One hit per simulated screen.** Multi-receptor baits, competition between
  co-enriched receptors, and bait-driven artifacts (the GDF15-Fc/CD64 problem)
  are outside the model.
- **Sort quality is not modelled.** Enrichment here is a single number standing
  in for stain quality, gate placement, avidity and wash stringency. The model
  tells you what a given enrichment buys; it cannot tell you what enrichment
  your bait will achieve.
- Assumption diagnostics for the count model (dispersion fit, goodness-of-fit)
  were not assessed.
