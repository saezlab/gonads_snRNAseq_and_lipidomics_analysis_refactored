# Editorial and reproducibility notes

These notes accompany `cosmos_results_methods.md`; they are not manuscript Results text. The supplied manuscript was used as a style and scientific-context reference. Its embedded revision bullets were not treated as instructions to modify other analyses or sections.

## Placement and claims

Place the Results paragraph after the closing sentence of “Gonadal sex-biased regulatory circuitry controlling lipid homeostasis”, before the zonation section. This placement does not imply that the preceding siRNA experiments validate the newly inferred MOON network. Figure and supplementary-table numbers have deliberately not been assigned.

The draft uses the previously selected **baseline C03**, located in `results/causal_modelling/cosmos_factor1_compressed_review/threshold_search/C03/`. It does not substitute the different candidate also called C03 in `sensitivity_threshold_search/`. The primary choice has pruning distance 6, scoring depth 4, thresholds 1.85/0.55, 153 nodes and 388 edges. The separate 468-edge network remains a wider-topology sensitivity; it is not the primary network described here.

“Selected configuration” is used instead of “optimal parameter set”. C03 was chosen to balance coverage, cutoffs and local threshold margins. It is not a unique optimum, and higher scores across changed topology/compression backgrounds do not establish stronger biological evidence. The pruning distance and scoring depth are explicitly separate; the selected candidate is a deliberately truncated propagation analysis of the six-step prepared topology, rather than a run with matched pruning/scoring values. Its shortest retained TF–lipid routes contain 7–9 edges.

The Results and Methods distinguish **251,998 curated-foundation graph edges** from the **254,589-edge PKN including species hypotheses**. The additional 149 rule records are not experimentally established reactions on those exact species. Likewise, 33,953 catalyst/reaction/direction nodes are not 33,953 independent biochemical reactions. The broader resource includes metabolite-to-protein activity interactions; direct lipid-to-transcript interactions were not established. The study did not use cached STITCH/BRENDA edges in the selected PKN.

Esr1, Nfic and Dbp were scored and then removed for sign disagreement. Their sparse informative downstream support is compatible with limitations in prior knowledge and in inference from this lipid signature; the audit does not establish prior incompleteness as the only cause. Lowering extraction cutoffs cannot reverse this exclusion. Neither these TFs' omission nor the retained graph establishes physiological absence or causality.

The MOON analysis considered all 1,254 lipid Factor 1 weights, mapping 155 to PKN nodes and retaining 74 after pruning. It was not restricted to the 64 lipids used for the earlier sparse correlation network. The PKN build report's “Factor 1” coverage denominator of 64 refers to that earlier subset. The fitted source model was initialized with 12 factors and retained 10; the saved Factor 1 is directly matched to the model archive. The draft avoids recharacterizing this as a newly trained model.

## Audit reproduced after the move

The following command **was executed successfully** from the repository root:

```sh
Rscript scripts/moon/audit_manuscript_causal_model.R
```

It reads the relocated PKN, input vectors, converged scoring state, compression map and exports. It verifies the PKN/input hashes against the original run, validates the original-to-compressed final graph projection, re-extracts the selected graph from the saved scores at 1.85/0.55, and reproduces its exact 388 signed edges, 153 nodes, four TFs and three lipids. It also checks the complete score-retention counts, expected TF signs and saved threshold/sensitivity summaries. It does not refit MOFA or recompute MOON scores. Results are in `results/causal_modelling/manuscript_audit/{README.md,audit_summary.json,artifact_manifest.csv}`.

**A full replay of the historical generation scripts is a separate outstanding reproducibility task.** Several still hardcode `results/cosmos_*`, and saved RDS objects contain obsolete relative or absolute directory strings. Two current search scripts also differ from their archived source hashes. The frozen PKN, inputs, selected scoring state and compression/extraction helper retain their recorded hashes. These facts permit exact frozen-state extraction auditing, but do not establish byte-identical end-to-end source replay. Before releasing a fully executable reproduction workflow, migrate paths and reconcile the script versions against the recorded manifests. Historical provenance records should be preserved; current locations should be resolved explicitly rather than overwriting archival path strings without a record. No global path migration or scientific refitting was performed during this drafting task.

## Code and data map

All data paths below use the current `results/causal_modelling/` location. These are audit entry points; the historical generation scripts still require the path/version reconciliation described above.

| Analysis component | Relevant code | Current data and audit records |
|---|---|---|
| PKN source assembly | `scripts/pkn/build_mouse_gem.py`, `build_rhea_lipid_reactions.py`, `build_lipid_regulation.py`, `assemble_mouse_pkn.py` | `results/causal_modelling/cosmos_mouse_pkn/pkn_lipid_rules.rds`, `build_manifest.json`, `build_summary.json`; raw downloads in `support/cosmos_mouse/raw/` and `results/causal_modelling/cosmos_pkn_audit/raw/` |
| Lipid identity and reaction rules | `scripts/pkn/map_lipids.py`, `build_lipid_rules.py` | `results/causal_modelling/cosmos_mouse_pkn/mapping_features.csv`, `mapping_candidates.csv`, `feature_coverage.csv`, `lipid_rules_reactions.jsonl`, `lipid_rules_pairs.csv`, `lipid_rules_manifest.json` |
| Full prior evidence | `scripts/pkn/audit_mouse_pkn.py`, `test_pkn.py`, `validate_export.R`, `write_build_manifest.py` | `results/causal_modelling/cosmos_mouse_pkn/edge_evidence.jsonl.gz`, `nodes.jsonl`, `mouse_gem_reactions.jsonl`, `rhea_reactions.jsonl`, `r_validation.json` |
| Factor 1 inputs | `scripts/moon/prepare_factor1_inputs.R`, `audit_factor1_model.py`, `map_factor1_inputs.R` | `results/mofa/main_model_RNA_lipid_weight.RData`, `results/mofa/RNA/mofa_res_12factor.hdf5`; `results/causal_modelling/cosmos_factor1_pilot/inputs/factor1_inputs.rds`, `factor1_lipid_node_mapping.csv`, `factor1_model_audit.json`, `input_source_manifest.csv` |
| Selected compressed scoring state | `scripts/moon/search_factor1_parameter_space.R`, `moon_safe_helpers.R` | `results/causal_modelling/cosmos_factor1_parameter_search/runs/radius6_depth4_compressed/run_state.rds`, `pkn_final_scoring.rds`, `pkn_final_uncompressed.rds`, `coherence_ledger.csv`, `run_manifest.json`, `sessionInfo.txt`; compression map in `results/causal_modelling/cosmos_factor1_parameter_search/topology_radius_6/compressed/compression.rds` |
| Threshold selection | `scripts/moon/search_compressed_threshold_review.R`, `audit_compressed_candidate_neighborhood.R` | `results/causal_modelling/cosmos_factor1_compressed_review/threshold_search/grid_complete.csv`, `selected_candidates.csv`, `selected_validation.csv`, `C03_local_threshold_neighborhood.csv` |
| Selected network and provenance | `scripts/moon/export_factor1_network.py`, `plot_factor1_network.R` | `results/causal_modelling/cosmos_factor1_compressed_review/threshold_search/C03/network.graphml`, `network_SIF.csv`, `network_ATT.csv`, `extraction_audit.json`, `network_evidence_network.csv`, `network_evidence_network_reactions.csv`, `candidate_paths_network.csv` |
| Topology and preprocessing sensitivities | `scripts/moon/run_factor1_compressed_sensitivities.R`, `run_factor1_inputs_only_compression.R`, `search_compressed_threshold_review.R --sensitivity` | `results/causal_modelling/cosmos_factor1_compressed_review/topology_sensitivity/`, `prescore_sensitivity/`, `inputs_only_sensitivity/`, `sensitivity_threshold_search/matched_threshold_sensitivity.csv` |
| Expected TF retention | `scripts/moon/audit_factor1_focus_tf_retention.R`, `validate_factor1_compressed_tf_sensitivities.R`, `export_compressed_score_retention.R` | `results/causal_modelling/cosmos_factor1_compressed_review/tf_audit/focus_tf_stage_audit.csv`, `focus_tf_score_background.csv`, `additional_sensitivity_tf_audit.csv`; `score_retention/radius6_depth4_compressed_scores_all.csv` under the same review directory |
| Relocated frozen-state audit | `scripts/moon/audit_manuscript_causal_model.R` | `results/causal_modelling/manuscript_audit/audit_summary.json`, `artifact_manifest.csv` |

## Citation notes

Author–year citations in the draft are linked to verified primary publications. Import these through the manuscript's reference manager rather than assigning new figure/table or bibliography numbers manually.

- **COSMOS+/MOON:** Dugourd et al., *Modeling causal signal propagation in multi-omic factor space with COSMOS*. bioRxiv [10.1101/2024.07.15.603538](https://doi.org/10.1101/2024.07.15.603538). The draft follows the package/tutorial's 2024 citation convention. Version 3 was posted on 24 April 2026; it remains a preprint in the source checked for this draft. This is the method-specific citation, distinct from the 2021 COSMOS/CARNIVAL paper.
- **OmniPath:** Türei et al. (2021), *Integrated intra- and intercellular signaling knowledge for multicellular omics analysis*. Molecular Systems Biology 17:e9923. [10.15252/msb.20209923](https://doi.org/10.15252/msb.20209923).
- **Mouse-GEM:** Wang et al. (2021), *Genome-scale metabolic network reconstruction of model animals as a platform for translational research*. PNAS 118:e2102344118. [10.1073/pnas.2102344118](https://doi.org/10.1073/pnas.2102344118). Model counts in the draft come from the frozen SBML snapshot, not the live repository.
- **CollecTRI:** Müller-Dott et al. (2023), *Expanding the coverage of regulons from high-confidence prior knowledge for accurate estimation of transcription factor activities*. Nucleic Acids Research 51:10934–10949. [10.1093/nar/gkad841](https://doi.org/10.1093/nar/gkad841).
- **LINEX2:** Rose et al. (2023), *Lipid network and moiety analysis for revealing enzymatic dysregulation and mechanistic alterations from lipidomics data*. Briefings in Bioinformatics, bbac572. [10.1093/bib/bbac572](https://doi.org/10.1093/bib/bbac572). This matches the existing manuscript's Rose et al. citation; our extension uses selected Rhea-derived templates and LINEX2 lipid definitions, rather than importing every LINEX2 network edge.
- **decoupleR:** Badia-i-Mompel et al. (2022), *decoupleR: ensemble of computational methods to infer biological activities from omics data*. Bioinformatics Advances 2:vbac016. [10.1093/bioadv/vbac016](https://doi.org/10.1093/bioadv/vbac016).

Frozen resource manifests document the exact Rhea, SwissLipids, Guide to Pharmacology, ChEBI and UniProt downloads used. Primary resource documentation and version details are also collected in `cosmos_pkn_fact_audit.md`. The saved COSMOS_basic tutorial commit is `7108c6b849c712b1282440dc122c9d783f563da8`. The Methods describe MOON ULM propagation and guarded extraction, not integer-linear-programming optimization.

## Supporting drafting audits

- `cosmos_pkn_fact_audit.md`: independently counted graph composition, provenance, resource versions and coverage denominators.
- `cosmos_reproducibility_audit.md`: exact input preparation, compression/coherence ledgers, scoring/threshold implementation, TF exclusion and relocation status.
- `cosmos_results_style_proposal.md`: manuscript style/placement review; an earlier proposal, not the final text in `cosmos_results_methods.md`.
