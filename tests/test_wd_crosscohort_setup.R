# Actual-data matching, frozen selection, direct Pearson and pairing invariants.
source("scripts/support_functions.R")
if (dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"), .libPaths()))
for (f in c("cyp2c38_reference.R", "wdsd_intake.R", "wd_intake.R", "wdsd_reference.R",
            "wdsd_comparison_setup.R", "wdsd_network_statistics.R", "wd_crosscohort_setup.R")) source(file.path("scripts/core", f))
a <- wp2_load_inputs(); b <- wd_load_inputs(); ref <- wp2_load_reference()
x <- wp3_prepare_contexts(a, b, ref)
stopifnot(identical(names(x$contexts), c("male_periportal", "male_pericentral")),
          all(x$coverage$n_CF453 == 7), all(x$coverage$n_CF459 == 4),
          all(x$WP2_replay$max_abs_difference < 1e-10))
for (ct in x$contexts) {
  e <- ct$edges
  for (cohort in c("CF453", "CF459")) {
    rna <- ct[[paste0("rna_", cohort)]]; lip <- ct[[paste0("lipid_", cohort)]]
    fast <- wp2_cor_edges(rna, lip, e)
    selected <- unique(c(1L, which(e$gene == "Mtmr11")[1L], nrow(e)))
    selected <- selected[!is.na(selected)]
    slow <- vapply(selected, function(i) stats::cor(rna[e$gene[i], ], lip[e$feature_id[i], ]), numeric(1L))
    stopifnot(max(abs(slow - fast[selected])) < 1e-12)
  }
  frozen <- ref[[e$context[1]]]$edges
  ix <- match(paste(e$gene, e$reference_lipid), paste(frozen$gene, frozen$reference_lipid))
  stopifnot(!anyNA(ix), identical(e$primary_edge, frozen$primary_edge[ix]), identical(e$r_reference, frozen$r_reference[ix]))
}
expect_error <- function(expr, text) {
  z <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(z, "error"), grepl(text, conditionMessage(z), fixed = TRUE))
}
mis <- b; mis$rna$Periportal <- mis$rna$Periportal[, rev(seq_len(8)), drop = FALSE]
expect_error(wp3_prepare_contexts(a, mis, ref), "RNA-lipid mouse order")
# A change confined to the unanalysed SD arm must not affect WP3 eligibility.
sd <- a; sd$lipid[, sd$metadata$diet == "SD"] <- 0
z <- wp3_prepare_contexts(sd, b, ref)
stopifnot(identical(x$coverage, z$coverage), identical(x$contexts, z$contexts))
# A canonical collision is audited/excluded, never forced into a one-to-one join.
amb <- b; key <- x$lipid_crosswalk$canonical_lipid[which(x$lipid_crosswalk$shared_primary)[1]]
row <- amb$lipid_ledger[amb$lipid_ledger$canonical_lipid == key, , drop = FALSE]
row$feature_id <- "WD_TEST_AMBIGUOUS"; row$keep_primary <- FALSE; row$exclusion_reason <- "test_collision"
amb$lipid_ledger <- rbind(amb$lipid_ledger, row)
z <- wp3_prepare_contexts(a, amb, ref)
stopifnot(z$lipid_crosswalk$status_CF459[z$lipid_crosswalk$canonical_lipid == key] == "ambiguous_exact_annotation",
          all(!vapply(z$contexts, function(ct) key %in% wp1_canonical_lipid(ct$edges$reference_lipid), logical(1))))
print(x$coverage); print(x$WP2_replay)
cat("WP3 setup, exact identity, WD-only eligibility and independent Pearson tests passed.\n")
