# Endpoint support and compression diagnostic

This is a **conditional diagnostic on the compressed pilot's frozen final graph**.
It does not replace the main result or constitute an independently converged
uncompressed analysis. The graph was already selected by the compressed run's
five-iteration RNA coherence loop; the uncompressed comparison below scores that
same selected edge set once, without another coherence loop. A complete
uncompressed rerun from the initial pruned PKN is a separate sensitivity analysis.
That separate run has now completed under `uncompressed_sensitivity/`; its
comparison is summarized at the end of this document and in `fullcomparison.json`.

Reproduce with `Rscript scripts/moon/audit_score_support.R`. The script reads the
main runner's artifacts and creates only diagnostic files. It uses exact input
values from `inputs/factor1_inputs.rds` with the retained IDs from the main pilot.
CSV-roundtripped weights should not be used for a supposedly exact repeat:
printing precision can change the sign of numerically near-zero inferred scores.

## Why the primary extraction is empty

After reachability pruning, the scoring problem contains 74 measured lipid nodes
and 49 upstream TF inputs. There are 398 distinct enzyme/reaction nodes directly
connected to these observations. **Every one of these 398 endpoint enzymes has
exactly one measured lipid target.** This count concerns observations entering
MOON, not the total number of chemically annotated participants of a reaction.

All first-layer ULM scores were reconstructed independently from the signed
coefficient vector and the complete background of 74 observed weights. They
agreed with the exported scores to better than `1e-10`. The reconstruction
confirms that these scores are single-target-versus-background contrasts, not
independent multi-lipid reaction footprints.

Of the 398 enzyme-to-observation edges, 228 agree with the source and target
score signs. Only four of those have endpoint enzyme `abs(ULM score) > 1.7`:

| Measured lipid | Endpoint reaction score | Immediate catalyst | Catalyst score |
|---|---:|---|---:|
| DG(24:0_18:2) | 2.971882 | Smpd1 | 0.867333 |
| DG(24:1_18:2) | 3.988675 | Smpd1 | 0.867333 |
| DG(20:0_18:1) | 1.891168 | Smpd1 | 0.867333 |
| PS(20:0_22:6) | 5.479610 | Ptdss1 | 0.831947 |

The secondary cutoff `abs(score) > 1` excludes both catalysts. It also excludes
L-serine, the other coherent predecessor of the PS reaction, whose score is
0.831947. Water is the remaining predecessor of the three DG reactions; its
score is -0.188424 and its positive input edge disagrees with the positive
reaction scores. Thus **every strong endpoint reaction loses all eligible
incoming explanatory edges**, and unconditional TF-to-observation pruning
correctly returns an empty network.

This is not explained by a shortage of high first-layer reaction scores. It is
the interaction of sparse endpoint support with the score of the immediately
upstream enzyme/cofactor layer. `strong_endpoint_predecessors.csv` contains every
predecessor edge underlying this conclusion.

## What compression changes

The protected-group compatibility fixes ensure correct IDs and lossless graph
projection. They do **not** make ULM scores invariant to node multiplicity.

The compressed final graph contains 26,776 edges; its projected uncompressed
counterpart contains 34,334. The map records 2,720 original nodes in 378 merged
groups. At level 1, 398 original endpoint enzymes correspond to 103 compressed
scoring nodes. Scores for those enzymes are identical between the two runs,
because both score the same 74 observed lipid values. At the next layer,
however, ULM sees **103 rows versus 398 rows**. This changes the repetition of
downstream values, their effective background, the correlation and residual
degrees of freedom. The differences propagate through later layers.

Unconstrained scoring was used here to expose every TF before the package's
upstream sign filter. The compressed repeat reproduces the main inferred values
for every retained main node to better than `1e-10` and reproduces its relaxed
extraction exactly. Both diagnostic runs return a finite inferred score for all
49 TFs.

Comparison in original IDs:

| Level | Shared original nodes | Pearson correlation of scores | Raw sign changes |
|---|---:|---:|---:|
| 0 | 74 | 1.000000 | 0 |
| 1 | 398 | 1.000000 | 0 |
| 2 | 50 | 0.909931 | 18 |
| 3 | 455 | 0.313098 | 149 |
| 4 | 1,261 | 0.566913 | 509 |
| 5 | 3,478 | -0.077781 | 2,542 |
| 6 | 829 | -0.254431 | 535 |

All levels remain the same for shared nodes. The raw sign counts include exact
zero versus tiny floating-point values, so individual near-zero changes should
not be interpreted. The large changes at later levels are not confined to
roundoff: median absolute score changes at levels 4, 5 and 6 are approximately
1.046, 1.469 and 1.473, respectively. Original-ID counts also repeat scores among
compressed group members and therefore are not independent replications.

Among the 49 TFs, 14 have activity-compatible signs in the compressed run and 15
in the uncompressed fixed-graph diagnostic; 11 are shared. Ahr, Esrrg and Gata6
lose agreement; Foxl2, Gata4, Kdm5b and Stat4 gain agreement. Nr1h3 is effectively
zero in both runs and should be described as unsupported, irrespective of the
machine-level sign. Full scores and input signs are in
`tf_unconstrained_compression_audit.csv`.

The illustrative catalyst scores change from 0.867333 to 2.003324 for Smpd1 and
from 0.831947 to 1.869276 for Ptdss1. This explains why the diagnostic uncompressed
scores can pass the primary/secondary extraction thresholds:

| Frozen-graph scoring choice | Primary / secondary | Nodes | Edges | TFs | Measured lipids |
|---|---|---:|---:|---:|---:|
| Compressed | 1 / 0 | 1,476 | 6,320 | 14 | 10 |
| Compressed | 1.7 / 1 | 0 | 0 | 0 | 0 |
| Uncompressed, no new coherence loop | 1 / 0 | 3,697 | 13,282 | 15 | 11 |
| Uncompressed, no new coherence loop | 1.7 / 1 | 467 | 1,650 | 8 | 4 |

These last two rows are **diagnostic extraction counts**, not a converged
alternative solution. They demonstrate that the empty primary result depends on
the compression/scoring choice. The appropriate next check is a separately
converged uncompressed run, with comparison of recurring routes and TFs.
Uncompressed scoring can itself overweight repeated enzyme representations that
share a measured lipid; neither choice removes the fundamental sparse-footprint
limitation.

## Esr1, Dbp and Nfic are inferred but rejected for opposite signs

They are not lost because the algorithm never reaches them:

| TF | Factor-derived input activity | Compressed inferred score | Uncompressed fixed-graph score |
|---|---:|---:|---:|
| Esr1 | -7.088749 | +0.131796 | +0.052402 |
| Dbp | -7.660885 | +0.345156 | +0.868470 |
| Nfic | -4.349231 | +0.718590 | +1.904342 |

All three have inferred signs opposing their RNA-derived TF activities under
both representations of this frozen graph. The package's upstream-input filter
therefore excludes them from the main results. This is specific to the current
PKN, represented lipid observations, factor orientation and scoring assumptions;
it does not refute their biological relevance.

## Files

- `score_support_audit.csv`: one row per directly measured endpoint enzyme,
  measured-target count, coefficients, ULM reconstruction, score and provenance.
- `strong_endpoint_predecessors.csv`: all incoming dependencies of the four
  strong, sign-consistent endpoint reactions.
- `compressed_unconstrained_scores.csv`: compressed-repeat scores expanded to
  original IDs, before TF sign filtering.
- `uncompressed_sensitivity_scores.csv`: one unconstrained ULM score run on the
  frozen uncompressed final graph.
- `tf_unconstrained_compression_audit.csv`: all 49 upstream TF inputs with both
  inferred scores, score availability and sign-compatibility flags.
- `compression_score_comparison.csv` and `compression_level_comparison.csv`:
  complete original-ID and per-level score comparisons.
- `compression_extraction_comparison.csv`: diagnostic extraction counts only.

All numeric conclusions above are calculated from local run artifacts. No
external literature evidence or additional biological measurements were added.

## Subsequent comparison of independently converged runs

The full uncompressed sensitivity starts from the same initially pruned PKN and
converges under its own TF-target coherence loop in four iterations. This should
be distinguished from the earlier frozen-graph diagnostic above. Its primary
extraction has 473 nodes, 1,676 edges, eight TFs and four measured lipids; its
strict extraction has 187 nodes, 537 edges, three TFs and two lipids.

One unconstrained score call on its final graph verifies that all 49 TFs are
inferred before the sign filter, as in the compressed run. Both fully converged
runs retain 14 activity-compatible TFs, of which 11 are shared. The compressed
run alone retains Ahr, Esrrg and Gata6; the uncompressed run alone retains Foxl2,
Kdm5b and Stat4. The relaxed networks share 3,219 edges (Jaccard 0.196917), 1,032
nodes (Jaccard 0.252261) and 10 measured lipids (Jaccard 0.909091). This confirms
substantial dependence of the detailed inferred topology on compression even
after separately reaching coherence convergence.

Esr1, Dbp and Nfic remain opposite to their input activity signs in the fully
converged uncompressed run: their inferred values are +0.05240172, +0.86686593
and +2.07762909, respectively. Their exclusion is thus a sign disagreement,
not failure to score the TF nodes, in both complete runs.

`fullcomparison.json` records run iteration counts, final-graph MD5s, TF identity
sets, key TF values and all extraction overlaps. `fullcomparison_overlap.csv`
contains the edge/node/TF/lipid counts and Jaccards for relaxed, primary and
strict settings. `fullcomparison_tf.csv` records all 49 TF inputs and both
unconstrained inferred scores, including near-zero flags. The additional raw
unconstrained result is
`uncompressed_sensitivity/moon_scores_unconstrained_audit.csv`.
