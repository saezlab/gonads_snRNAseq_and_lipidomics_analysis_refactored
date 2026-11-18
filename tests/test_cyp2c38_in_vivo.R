source("scripts/support_functions.R")
if (dir.exists(".local/R-library")) .libPaths(c(normalizePath(".local/R-library"), .libPaths()))
source("scripts/core/cyp2c38_models.R")
set.seed(20261006)
n <- 1200L
meta <- data.frame(
  sample = paste0("mouse_", seq_len(32L)),
  group = rep(c("F_WT", "F_KO", "M_WT", "M_KO"), c(8L, 8L, 7L, 9L)),
  stringsAsFactors = FALSE
)
meta$sex <- substr(meta$group, 1L, 1L)
baseline <- seq(5, 15, length.out = n)
female_effect <- stats::rnorm(n, sd = 0.4)
male_effect <- stats::rnorm(n, sd = 0.4)
female_effect[1:6] <- c(2, -2, 0, 1, -1, 0.5)
male_effect[1:6] <- c(1, -3, -2, -1, 2, 0.5)
noise <- matrix(stats::rnorm(n * nrow(meta)), nrow = n)
noise <- noise * (1.5 - baseline / 20)
noise[, meta$sex == "M"] <- noise[, meta$sex == "M"] * 1.6
for (group in unique(meta$group)) {
  selected <- meta$group == group
  noise[, selected] <- noise[, selected] - rowMeans(noise[, selected])
}
y <- noise + baseline
y[, meta$sex == "M"] <- y[, meta$sex == "M"] + 3
y[, meta$group == "F_KO"] <- y[, meta$group == "F_KO"] + female_effect
y[, meta$group == "M_KO"] <- y[, meta$group == "M_KO"] + male_effect
dimnames(y) <- list(paste0("lipid_", seq_len(n)), meta$sample)

warnings <- character()
capture_warnings <- function(code) {
  withCallingHandlers(code, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
}
elapsed <- system.time(weighted <- capture_warnings(wp1_fit(y, meta)))[["elapsed"]]
unweighted <- capture_warnings(wp1_fit(y, meta, weighted = FALSE))
separate <- capture_warnings(wp1_fit_separate(y, meta))
coeff <- weighted$fit$coefficients
stopifnot(qr(weighted$design)$rank == 4L,
          all(weighted$fit$df.residual == 28L),
          max(abs(coeff[, "female_KO_vs_WT"] - female_effect)) < 1e-9,
          max(abs(coeff[, "male_KO_vs_WT"] - male_effect)) < 1e-9,
          max(abs(coeff[, "pooled_KO_vs_WT"] - (female_effect + male_effect)/2)) < 1e-9,
          max(abs(coeff[, "interaction_M_minus_F"] - (male_effect - female_effect))) < 1e-9,
          max(abs(coeff - unweighted$fit$coefficients)) < 1e-9,
          max(abs(separate$female$coefficients[,1] - female_effect)) < 1e-9,
          max(abs(separate$male$coefficients[,1] - male_effect)) < 1e-9,
          all(sign(coeff[1:6,"female_KO_vs_WT"])[female_effect[1:6] != 0] ==
              sign(female_effect[1:6])[female_effect[1:6] != 0]))
for (contrast in colnames(coeff)) {
  table <- wp1_table(weighted$fit, contrast)
  stopifnot(nrow(table) == n,
            identical(table$feature_id, rownames(y)),
            identical(table$t, as.numeric(weighted$fit$t[,contrast])),
            all(table$P.Value >= 0 & table$P.Value <= 1),
            all(table$adj.P.Val >= table$P.Value - 1e-12))
}
for (sex in unique(meta$sex)) {
  stopifnot(diff(range(weighted$weights[meta$sex == sex])) < 1e-12)
}
# A large baseline sex effect and unequal male genotype counts must not change
# the known equal-weighted pooled genotype effect. The naive pooled effect is
# intentionally biased in this simulation and must differ.
naive <- rowMeans(y[,grepl("_KO$",meta$group)]) -
  rowMeans(y[,grepl("_WT$",meta$group)])
stopifnot(abs(naive[1] - coeff[1,"pooled_KO_vs_WT"]) > 0.1)

# Reversing the perturbation contrast reverses both t and signed consistency.
reverse <- -weighted$contrasts[,"pooled_KO_vs_WT",drop=FALSE]
reverse_fit <- limma::eBayes(
  limma::contrasts.fit(
    limma::lmFit(y, weighted$design, weights=weighted$weights), reverse
  ), trend=TRUE, robust=TRUE
)
stopifnot(max(abs(reverse_fit$t[,1] + weighted$fit$t[,"pooled_KO_vs_WT"])) < 1e-9)
reference <- seq(-0.8, 0.8, length.out=n)
consistency <- -reference * weighted$fit$t[,"pooled_KO_vs_WT"]
stopifnot(max(abs((-reference * reverse_fit$t[,1]) + consistency)) < 1e-9)

# Metadata scrambling must fail, not silently re-pair mice.
scrambled <- meta[c(2,1,3:nrow(meta)),]
stopifnot(inherits(try(wp1_fit(y,scrambled),silent=TRUE),"try-error"))
bad <- meta
bad$group[1] <- "unknown"
stopifnot(inherits(try(wp1_fit(y,bad),silent=TRUE),"try-error"))

cat("Synthetic checks passed for",n,"features x",nrow(meta),"mice.\n")
cat("Joint weighted fit elapsed:",elapsed,"seconds.\n")
cat("Precision weights by sex:\n")
print(tapply(weighted$weights,meta$sex,unique))
cat("Warnings:",length(warnings),"\n")
if (length(warnings)) print(unique(warnings))


source("scripts/core/cyp2c38_reference.R")
stopifnot(identical(wp1_canonical_lipid(c("PC 16:0/18:1", "PC 18:1_16:0", "PC 34:1|PC 16:0_18:1", NA_character_)),
                    c("PC(16:0/18:1)", "PC(18:1_16:0)", "PC(16:0_18:1)", NA_character_)))
# Formatting must never equate unresolved/resolved species or chain permutations.
stopifnot(length(unique(wp1_canonical_lipid(c("PC 34:1", "PC 16:0/18:1", "PC 16:0_18:1", "PC 18:1_16:0")))) == 4)
# Reordering whole mouse columns together with metadata must preserve estimates.
p <- sample(seq_len(ncol(y)))
reordered <- wp1_fit(y[,p],meta[p,])
stopifnot(max(abs(reordered$fit$coefficients-weighted$fit$coefficients)) < 1e-9)
# One shared genotype permutation is applied to all lipid rows; counts preserved.
mp <- meta
for (sx in c("F", "M")) {
  ii <- which(mp$sex==sx)
  g <- sub("^._", "",mp$group[ii])
  mp$group[ii] <- paste(sx,sample(g),sep="_")
}
stopifnot(identical(table(mp$group),table(meta$group)), identical(colnames(y),mp$sample))
cat("Identity, reordering and permutation-stratum checks passed.\n")
