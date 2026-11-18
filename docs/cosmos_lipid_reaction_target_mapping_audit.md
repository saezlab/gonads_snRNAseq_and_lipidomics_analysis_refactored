# Measured lipids as reaction-product targets

This audit addresses a narrower question than lipid identity or lipid-to-protein
regulation coverage: which measured lipids are **product-side endpoints** of
the directed lipid reaction graph? A feature merely appearing in a reaction
does not make it a downstream target. In the PKN projection, a substrate is
connected as `metabolite → reaction instance`, whereas a product is connected
as `reaction instance → metabolite`.

The audit is reproducible with:

```bash
python3 scripts/pkn/audit_measured_lipid_reaction_targets.py \
  --output results/causal_modelling/lipid_reaction_target_mapping_audit
```

It reads the frozen `pkn_lipid_rules`, reaction evidence ledger, lipid mapping,
frozen input mapping and selected TF-to-lipid MOON result. It creates a
feature-level inventory, one validation row per species-rule reaction, summary
tables and dependency hashes. It neither changes the PKN nor refits MOON.

## Finding

All 149 measured-species reaction hypotheses are directionally intact. The
builder substitutes the measured source only on the Rhea substrate side and the
measured target only on the product side. The PKN then projects these as
`substrate → enzyme reaction instance → product`. Each of the 449
catalyst-specific LINEX2-rule product edges survives signed-edge conflict
reconciliation.

The 149 features labelled `species_rule_connected` in the PKN coverage table
should therefore not be read as 149 downstream products. They comprise 100
species-rule substrates and 123 species-rule products: 74 occur on both sides,
26 are substrate-only and 49 product-only. The narrower product-side rule
coverage is 123/1,254 measured features (9.8%). Six further features have an
exact SwissLipids-to-ChEBI *candidate* that is an active Rhea product. One of
these, `FA(16:0)`, remains excluded from the frozen input vector because two
exact ChEBI identities cannot be disambiguated. Thus 129/1,254 (10.3%) have a
product candidate under the strict mapping policy, and 128/1,254 (10.2%) are
usable frozen mapped reaction-product endpoints (123 species-rule products plus
five unambiguous ChEBI products).

The historical 64-lipid Factor 1 correlation subset has four product-side rule
features: `DG(24:0_18:2)`, `DG(24:1_18:2)`, `PE(22:1_22:6)` and
`PS(20:0_22:6)`. Seven features in that subset are connected to a rule in any
role, with three substrate-only features. This limited coverage follows from
the deliberately conservative rules: only seven diacyl classes with two
explicit, unmodified ester acyl chains are paired across eight directional
Rhea templates. Sum-composition species, TGs, ethers, oxidized lipids,
sphingolipids, multi-lipid transformations, elongation/desaturation and paths
through unmeasured intermediates are excluded.

## Implication for the selected TF-to-lipid result

The frozen input mapping contained 155 nodes because it accepted every
unambiguous exact chemical node or measured-species rule node, before testing
directed reachability. Of those, 128 are active product endpoints; the 27
remaining nodes are the 26 rule substrate-only features plus
`LPC(20:0/0:0)`, which has exact identity but is only a Rhea substrate.

This does **not** affect the selected solution. Iterative directed topology
pruning retained 74 downstream lipid observations, and every one is a product
node (69 rule products and five exact ChEBI products). The three final
level-0 observations—`DG(20:0_18:1)`, `DG(24:0_18:2)` and
`DG(24:1_18:2)`—are all product-side reaction endpoints. The audit found no
source/product reversal, missing product projection or conflict-induced loss
that could have created the sparse coverage.

The result remains a prior-guided reaction hypothesis. Exact measured-species
rules transfer directional class chemistry and catalyst annotation to a
chain-compatible measured species; they do not establish enzyme specificity,
flux, joint substrate availability or a change in lipid abundance.
