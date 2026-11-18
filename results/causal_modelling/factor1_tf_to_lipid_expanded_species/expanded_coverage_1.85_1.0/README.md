# Expanded Factor 1 threshold re-extraction

This ancillary analysis reuses the completed compressed fit named in results_summary.json. No MOON scoring is performed.
All valid, invalid and empty threshold rows are retained. If no pair meets the compact coverage criteria, network files are empty and no solution.rds is written.
node_retention_diagnostics.csv uses the saved pre-upstream-filter scores. Thresholds are exploratory numerical settings, not calibrated significance measures.
Within-grid representative: 0<edges<500, >=2TFs and >=2observations; descending secondary then primary then observation/TF coverage then fewer edges. Not a calibrated optimum.
