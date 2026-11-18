# Measured lipid reaction-target audit

This directory distinguishes a feature that is merely connected to a lipid rule from a feature that is an active product-side endpoint in the directed PKN. The audit does not modify the PKN or refit MOON.

`feature_reaction_target_inventory.csv` contains one row per measured feature. `rule_direction_validation.csv` confirms the source/product substitution and active catalyst-specific product edges for every rule. `population_summary.csv` compares the complete panel, historical 64-lipid subset, frozen mapped inputs, selected topology and selected solution. `summary.json` contains definitions, headline counts and checks.
