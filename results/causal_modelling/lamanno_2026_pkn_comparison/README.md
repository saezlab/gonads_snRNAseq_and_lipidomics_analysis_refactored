# Mouse brain lipid atlas PKN comparison

Read-only comparison of the recovered Table 10 and deposited annotation-derived lipid graph with our frozen mouse PKN. See `docs/cosmos_lamanno_pkn_comparison.md` for source discovery, interpretation and limitations.

This folder contains resource inventories and conditional mapping audits, not a new PKN or MOON solution. The complete signed CORNETO graph was not located in the cited deposits. Sum-composition overlaps must not be treated as resolved reaction assignments.

Reproduce from the repository root with Python 3.9 or later:

```bash
python3 scripts/pkn/compare_lamanno_lipid_network.py --self-test
python3 scripts/pkn/compare_lamanno_lipid_network.py --output /tmp/lamanno_pkn_comparison_replay
```

Choose an empty output destination. The script verifies source-file hashes and writes dependency hashes. Inputs remain under `support/cosmos_mouse/comparison/lamanno_2026/` and the frozen PKN/product-endpoint audit folders.

Headline audit counts: 350 Table 10 reaction pairs, 113 lipid labels and 30 printed enzymes; 29 enzymes already present in our active PKN. Table 10 has 137 possible measurement matches (46 identical shorthand, 91 conditional sums). The broader deposited lipid draft has 4,562 pairs, 675 nodes and 392 possible measurement matches (125 identical shorthand, 267 conditional sums). The retained-chain audit flags 172 of 288 naive lyso/diacyl assignments across 51 measured features.
