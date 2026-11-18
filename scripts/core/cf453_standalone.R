# Standalone CF453 statistics. Read workflow plan 14 before editing.
source("scripts/config/analysis_release.R")
cf6_write <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

cf6_tests <- function(fit, families) {
  rows <- lapply(seq_len(ncol(fit$coefficients)), function(j) {
    se <- fit$stdev.unscaled[, j] * sqrt(fit$s2.post)
    crit <- stats::qt(.975, fit$df.total)
    data.frame(gene = rownames(fit$coefficients), contrast = colnames(fit$coefficients)[j],
      effect = fit$coefficients[, j], se = se, t = fit$t[, j], p = fit$p.value[, j],
      low = fit$coefficients[, j] - crit * se, high = fit$coefficients[, j] + crit * se,
      family = families[j], stringsAsFactors = FALSE, row.names = NULL)
  })
  result <- do.call(rbind, rows)
  result$q <- ave(result$p, result$family, FUN = function(v) stats::p.adjust(v, "BH"))
  result
}

cf6_models <- function(root = ".", seed = 20261030L, bootstraps = 999L, input_zonation = NULL) {
  input <- wp2_load_inputs(root, input_zonation = input_zonation)
  input_provenance <- data.frame(path = unname(input$input_paths),
    sha256 = vapply(input$input_paths, digest::digest, character(1L), file = TRUE, algo = "sha256"),
    stringsAsFactors = FALSE, row.names = NULL)
  md <- input$metadata
  pm <- input$rna_profile_metadata
  pm <- pm[pm$zonation_type != "Hepatocytes", , drop = FALSE]
  pm$zone <- ifelse(pm$zonation_type == "Periportal", "PP", "PC")
  pm$cell <- paste(pm$sex, pm$diet, pm$zone, sep = "_")
  rna <- sapply(seq_len(nrow(pm)), function(i) input$rna_data[[pm$zonation_type[i]]][, pm$mouse_id[i]])
  rownames(rna) <- rownames(input$rna$Periportal); colnames(rna) <- pm$profile_id
  counts <- sapply(seq_len(nrow(pm)), function(i) input$rna_counts[[pm$zonation_type[i]]][, pm$mouse_id[i]])
  detected <- rowSums(counts > 0)
  eligible <- data.frame(gene = rownames(rna), detected_profiles = detected,
    variable = apply(rna, 1L, stats::sd) > 1e-8, stringsAsFactors = FALSE)
  eligible$keep <- eligible$detected_profiles >= 7L & eligible$variable
  rna <- rna[eligible$keep, , drop = FALSE]
  cells <- as.vector(outer(as.vector(outer(c("female", "male"), c("SD", "WD"), paste, sep = "_")), c("PC", "PP"), paste, sep = "_"))
  design <- stats::model.matrix(~0 + factor(pm$cell, levels = cells))
  colnames(design) <- cells; rownames(design) <- pm$profile_id
  contrasts <- list(); families <- character()
  add <- function(name, positive, negative, family) {
    v <- stats::setNames(numeric(length(cells)), cells)
    for (p in positive) v[p] <- v[p] + 1
    for (n in negative) v[n] <- v[n] - 1
    contrasts[[name]] <<- v; families[name] <<- family
  }
  for (d in c("SD", "WD")) {
    add(paste0("sex_zone_", d), paste0(c("male_", "female_"), d, c("_PP", "_PC")),
      paste0(c("male_", "female_"), d, c("_PC", "_PP")), "primary_interactions")
  }
  contrasts$sex_diet_zone <- contrasts$sex_zone_WD - contrasts$sex_zone_SD
  families["sex_diet_zone"] <- "primary_interactions"
  for (s in c("female", "male")) for (d in c("SD", "WD")) {
    add(paste("zone", s, d, sep = "_"), paste(s, d, "PP", sep = "_"),
      paste(s, d, "PC", sep = "_"), "secondary_zonation")
  }
  for (s in c("female", "male")) for (z in c("PC", "PP")) {
    add(paste("diet", s, z, sep = "_"), paste(s, "WD", z, sep = "_"),
      paste(s, "SD", z, sep = "_"), "secondary_diet_by_zone")
  }
  for (s in c("female", "male")) {
    contrasts[[paste0("diet_zone_", s)]] <- contrasts[[paste0("zone_", s, "_WD")]] - contrasts[[paste0("zone_", s, "_SD")]]
    families[paste0("diet_zone_", s)] <- "secondary_diet_zonation"
  }
  for (d in c("SD", "WD")) for (z in c("PC", "PP")) {
    add(paste("sex", d, z, sep = "_"), paste("male", d, z, sep = "_"),
      paste("female", d, z, sep = "_"), "secondary_sex_by_zone")
  }
  contrast_matrix <- do.call(cbind, contrasts)
  rho <- limma::duplicateCorrelation(rna, design, block = pm$mouse_id)$consensus.correlation
  fit <- limma::lmFit(rna, design, block = pm$mouse_id, correlation = rho)
  fit <- limma::eBayes(limma::contrasts.fit(fit, contrast_matrix), trend = TRUE, robust = TRUE)
  tests <- cf6_tests(fit, unname(families[colnames(contrast_matrix)]))
  cat("CF453 paired RNA model: ", nrow(rna), " genes x ", ncol(rna), " profiles; rho = ", rho, "\n", sep = "")
  print(stats::aggregate(I(q <= .05) ~ contrast, tests, sum))
  candidates <- list(); k <- 0L; zonation <- list()
  for (d in c("SD", "WD")) {
    interaction <- tests[tests$contrast == paste0("sex_zone_", d), ]
    female <- tests[tests$contrast == paste0("zone_female_", d), ]
    male <- tests[tests$contrast == paste0("zone_male_", d), ]
    stopifnot(identical(interaction$gene, female$gene), identical(male$gene, female$gene))
    category <- rep("unselected", nrow(female))
    for (s in c("female", "male")) {
      simple <- if (s == "female") female else male
      keep <- interaction$q <= .05 & simple$q <= .05
      if (any(keep)) {
        k <- k + 1L
        candidates[[k]] <- data.frame(gene = simple$gene[keep], sex = s, diet = d,
          zone = ifelse(simple$effect[keep] > 0, "PP", "PC"), effect = simple$effect[keep],
          q_interaction = interaction$q[keep], q_zonation = simple$q[keep], stringsAsFactors = FALSE)
      }
      current <- paste(s, ifelse(simple$effect > 0, "PP", "PC"), sep = "_")
      category[keep] <- ifelse(category[keep] == "unselected", current[keep], "shared")
    }
    zonation[[d]] <- data.frame(gene = female$gene, diet = d, t_female = female$t, t_male = male$t,
      effect_female = female$effect, effect_male = male$effect, q_interaction = interaction$q,
      category = category, stringsAsFactors = FALSE)
  }
  candidates <- if (length(candidates)) do.call(rbind, candidates) else data.frame(gene=character(),sex=character(),diet=character(),zone=character(),effect=numeric(),q_interaction=numeric(),q_zonation=numeric())
  ledger <- input$lipid_ledger
  ledger$keep_broad <- ledger$keep_primary
  ledger$keep_primary <- ledger$keep_qc30
  ledger$fill <- as.numeric(ledger[["Fill %"]])
  ledger$keep_fill <- ledger$keep_primary & !is.na(ledger$fill) & ledger$fill >= .5
  ledger$lipid <- ledger$feature_id
  lipid_broad <- input$lipid
  lipid <- lipid_broad[ledger$feature_id[ledger$keep_primary], , drop = FALSE]
  ld <- stats::model.matrix(~0 + factor(md$group, levels = c("female-sd", "female-wd", "male-sd", "male-wd")))
  colnames(ld) <- c("female_SD", "female_WD", "male_SD", "male_WD")
  lc <- cbind(sex_SD = c(-1,0,1,0), sex_WD = c(0,-1,0,1), diet_female = c(-1,1,0,0),
    diet_male = c(0,0,-1,1), sex_diet = c(1,-1,-1,1))
  lfit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(lipid, ld), lc))
  lt <- cf6_tests(lfit, c("lipid_sex", "lipid_sex", "lipid_diet", "lipid_diet", "lipid_interaction"))
  names(lt)[names(lt) == "gene"] <- "lipid"
  lt$lipid_class <- ledger$lipid_class[match(lt$lipid, ledger$feature_id)]
  lt$lipid_label <- ledger$canonical_lipid[match(lt$lipid, ledger$feature_id)]
  extract_effects <- function(prefix, split_zone = TRUE) {
    r <- tests[startsWith(tests$contrast, prefix), ]
    fields <- strsplit(sub(paste0("^", prefix), "", r$contrast), "_")
    r$sex <- vapply(fields, `[`, character(1L), 1L)
    if (split_zone) r$zone <- vapply(fields, `[`, character(1L), 2L)
    r
  }
  qc_pca <- function(x, metadata, kind) {
    p <- stats::prcomp(t(x), center = TRUE, scale. = FALSE)
    data.frame(kind = kind, profile_id = rownames(p$x), mouse_id = metadata$mouse_id,
      sex = metadata$sex, diet = metadata$diet, zone = if ("zone" %in% names(metadata)) metadata$zone else "whole_liver",
      PC1 = p$x[,1], PC2 = p$x[,2], variance_PC1 = p$sdev[1]^2 / sum(p$sdev^2),
      variance_PC2 = p$sdev[2]^2 / sum(p$sdev^2), row.names = NULL)
  }
  list(metadata = md, profile_metadata = pm, modeled_rna = rna, gene_eligibility = eligible,
    design = design, contrasts = contrast_matrix, contrast_families = families,
    rho = rho, fit = fit, rna_results = tests, candidates = candidates,
    zonation = do.call(rbind, zonation), primary_interactions = tests[tests$family == "primary_interactions", ], lipid_matrix = lipid, lipid_broad = lipid_broad,
    lipid_ledger = ledger, lipid_design = ld, lipid_contrasts = lc, lipid_fit = lfit,
    lipid_results = lt, diet_effects = extract_effects("diet_" )[!startsWith(extract_effects("diet_")$contrast, "diet_zone_"), ],
    zonation_diet_effects = extract_effects("diet_zone_", FALSE),
    qc_pca = rbind(qc_pca(rna, pm, "RNA"), qc_pca(lipid, md, "Lipids")),
    input_paths = input$input_paths, input_provenance = input_provenance,
    input_qc = input[c("qc_summary", "rna_profile_qc", "layer_inventory", "marker_by_profile", "zonation_marker_diagnostic", "recorded_commands")],
    settings = list(seed = seed, bootstraps = bootstraps, fdr = .05, display_threshold = .65,
      RNA_scale = "difference in log-normalized expression: ln(1 + 10000 * supplied_counts / profile_library)",
      RNA_method = "limma with mouse duplicateCorrelation, mean-variance trend and robust empirical Bayes; raw aggregation provenance unverified",
      primary_lipid_policy = "repository eligible endogenous classes and analytical QC CV <=30%",
      selection_inference = "RNA-only intersection; correlation inference/resampling conditional on fixed selection",
      primary_family = "all eligible genes x three interactions", zonation_family = "all eligible genes x four simple zonation tests"))
}

# A models-only run is a complete, independently provenance-checked parent
# cache; it does not require the historical all-edge or 30-panel analysis.
cf6_save_models <- function(a, output) {
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  saveRDS(a, file.path(output, "models.rds"))
  for (name in names(a)[vapply(a, is.data.frame, logical(1L))])
    cf6_write(a[[name]], file.path(output, "tables", paste0(name, ".csv")))
  cf6_write(data.frame(profile_id = rownames(a$design), a$design, check.names = FALSE),
    file.path(output, "tables/RNA_design.csv"))
  cf6_write(data.frame(cell = rownames(a$contrasts), a$contrasts, check.names = FALSE),
    file.path(output, "tables/RNA_contrasts.csv"))
  code <- c("scripts/core/wdsd_intake.R", "scripts/core/cf453_standalone.R",
    "scripts/workflows/14_cf453_standalone_figure6.R", "docs/workflow_plans/14_cf453_standalone_figure6.md")
  paths <- unique(c(unname(a$input_paths), file.path(output, "models.rds"), code[file.exists(code)]))
  cf6_write(analysis_provenance(paths),
    file.path(output, "provenance_sha256.csv"))
  writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
  invisible(a)
}

cf6_cor <- function(x, y) {
  center <- function(m) {
    m <- m - rowMeans(m); denom <- sqrt(rowSums(m^2))
    m / ifelse(denom > 1e-12, denom, NA_real_)
  }
  tcrossprod(center(x), center(y))
}

cf6_boot_edges <- function(x, y, gene_idx, lipid_idx, indices, observed, chunk = 384L) {
  n <- ncol(x); b <- ncol(indices)
  weights <- vapply(seq_len(b), function(j) tabulate(indices[,j], nbins = n), integer(n))
  result <- matrix(NA_real_, length(gene_idx), 4L, dimnames = list(NULL, c("bootstrap_low", "bootstrap_high", "bootstrap_sign_fraction", "bootstrap_valid")))
  for (first in seq.int(1L, length(gene_idx), by = chunk)) {
    ix <- first:min(first + chunk - 1L, length(gene_idx))
    xx <- x[gene_idx[ix], , drop = FALSE]; yy <- y[lipid_idx[ix], , drop = FALSE]
    # Center first to reduce cancellation; covariance is invariant to these shifts.
    xx <- xx - rowMeans(xx); yy <- yy - rowMeans(yy)
    sx <- xx %*% weights; sy <- yy %*% weights
    vx <- xx^2 %*% weights - sx^2/n; vy <- yy^2 %*% weights - sy^2/n
    cov <- (xx*yy) %*% weights - sx*sy/n
    rr <- cov / sqrt(pmax(vx,0) * pmax(vy,0))
    rr[vx <= 1e-12 | vy <= 1e-12] <- NA_real_
    rr[] <- pmax(-1, pmin(1, as.vector(rr)))
    result[ix, "bootstrap_valid"] <- rowSums(is.finite(rr))
    result[ix, "bootstrap_sign_fraction"] <- rowMeans(sign(rr) == sign(observed[ix]), na.rm = TRUE)
    qs <- t(apply(rr, 1L, stats::quantile, probs = c(.025,.975), na.rm = TRUE, names = FALSE))
    result[ix, 1:2] <- qs
  }
  result
}

cf6_membership <- function(edges, ledger, threshold = .65, quality = "keep_primary") {
  lip <- ledger[ledger[[quality]], , drop = FALSE]
  order_contexts <- c("female_PP", "female_PC", "male_PP", "male_PC")
  all <- lapply(c("SD", "WD"), function(d) {
    matrix <- matrix(FALSE, nrow(lip), 4L, dimnames = list(lip$feature_id, order_contexts))
    e <- edges[edges$diet == d & is.finite(edges$r) & abs(edges$r) >= threshold, , drop = FALSE]
    if (nrow(e)) {
      e <- e[e$active_candidate, , drop = FALSE]
      at <- match(e$lipid, rownames(matrix)); context <- match(paste(e$sex,e$zone,sep="_"),colnames(matrix))
      keep <- !is.na(at) & !is.na(context)
      matrix[cbind(at[keep],context[keep])] <- TRUE
    }
    n <- rowSums(matrix)
    category <- ifelse(n==0, "unassigned", ifelse(n>1,"shared", colnames(matrix)[max.col(matrix,ties.method="first")]))
    data.frame(lipid=lip$feature_id,lipid_class=lip$lipid_class,diet=d,category=category,
      membership_bits=apply(matrix,1L,function(v) paste(as.integer(v),collapse="")),n_contexts=n,stringsAsFactors=FALSE)
  })
  do.call(rbind, all)
}

cf6_networks <- function(a) {
  set.seed(a$settings$seed)
  groups <- c("female_SD", "female_WD", "male_SD", "male_WD")
  banks <- array(NA_integer_, c(7L,a$settings$bootstraps,4L), dimnames=list(NULL,NULL,groups))
  for (group in groups) banks[,,group] <- replicate(a$settings$bootstraps, sample.int(7L,7L,TRUE))
  all_edges <- list(); all_loo <- list(); eligible <- list(); z <- 0L
  for (s in c("female","male")) for (zone in c("PP","PC")) {
    genes <- sort(unique(a$candidates$gene[a$candidates$sex==s & a$candidates$zone==zone]))
    if (!length(genes)) next
    lipid_ids <- rownames(a$lipid_broad)
    valid_g <- rep(TRUE,length(genes)); valid_l <- rep(TRUE,length(lipid_ids))
    for (d in c("SD","WD")) {
      ids <- a$metadata$mouse_id[a$metadata$sex==s & a$metadata$diet==d]
      p <- a$profile_metadata$profile_id[match(paste(ids,zone),paste(a$profile_metadata$mouse_id,a$profile_metadata$zone))]
      valid_g <- valid_g & apply(a$modeled_rna[genes,p,drop=FALSE],1L,stats::sd)>1e-8
      valid_l <- valid_l & apply(a$lipid_broad[lipid_ids,ids,drop=FALSE],1L,stats::sd)>1e-8
    }
    eligible[[paste(s,zone)]] <- rbind(data.frame(kind="gene",feature=genes,sex=s,zone=zone,keep=valid_g),
      data.frame(kind="lipid",feature=lipid_ids,sex=s,zone=zone,keep=valid_l))
    genes <- genes[valid_g]; lipid_ids <- lipid_ids[valid_l]
    if (!length(genes) || !length(lipid_ids)) next
    for (d in c("SD","WD")) {
      ids <- a$metadata$mouse_id[a$metadata$sex==s & a$metadata$diet==d]
      p <- a$profile_metadata$profile_id[match(paste(ids,zone),paste(a$profile_metadata$mouse_id,a$profile_metadata$zone))]
      x <- a$modeled_rna[genes,p,drop=FALSE]; y <- a$lipid_broad[lipid_ids,ids,drop=FALSE]
      observed <- cf6_cor(x,y)
      e <- expand.grid(gene=genes,lipid=lipid_ids,stringsAsFactors=FALSE)
      e$sex <- s; e$zone <- zone; e$diet <- d; e$r <- as.vector(observed); e$n <- 7L
      ix <- match(e$lipid,a$lipid_ledger$feature_id)
      e$lipid_class <- a$lipid_ledger$lipid_class[ix]; e$lipid_label <- a$lipid_ledger$canonical_lipid[ix]
      e$primary_lipid <- a$lipid_ledger$keep_primary[ix]; e$fill_lipid <- a$lipid_ledger$keep_fill[ix]
      active <- a$candidates$gene[a$candidates$sex==s & a$candidates$zone==zone & a$candidates$diet==d]
      e$active_candidate <- e$gene %in% active
      e$p <- 2*stats::pt(-abs(e$r)*sqrt(5/pmax(1-e$r^2,.Machine$double.eps)),df=5)
      loo <- vapply(seq_len(7L),function(j) as.vector(cf6_cor(x[,-j,drop=FALSE],y[,-j,drop=FALSE])),numeric(nrow(e)))
      e$loo_valid <- rowSums(is.finite(loo))
      e$loo_min <- apply(loo,1L,min,na.rm=TRUE); e$loo_max <- apply(loo,1L,max,na.rm=TRUE)
      e$loo_sign_fraction <- rowMeans(sign(loo)==sign(e$r),na.rm=TRUE)
      e$display <- e$primary_lipid & e$active_candidate & abs(e$r)>=a$settings$display_threshold
      e$stable_display <- e$display & e$loo_valid==7L & e$loo_sign_fraction==1
      # Bootstrap all primary-universe pairs; broad-only pairs retain full/LOO coefficients.
      boot <- matrix(NA_real_,nrow(e),4L,dimnames=list(NULL,c("bootstrap_low","bootstrap_high","bootstrap_sign_fraction","bootstrap_valid")))
      keep <- which(e$primary_lipid)
      boot[keep,] <- cf6_boot_edges(x,y,match(e$gene[keep],genes),match(e$lipid[keep],lipid_ids),banks[,,paste(s,d,sep="_")],e$r[keep])
      e <- cbind(e,boot)
      z <- z+1L; all_edges[[z]] <- e; all_loo[[z]] <- loo
      cat("CF453 network ",s," ",zone," ",d,": ",length(genes)," genes, ",nrow(e)," broad pairs, ",sum(e$stable_display)," stable displayed\n",sep="")
    }
  }
  if (!length(all_edges)) {
    a$edges <- data.frame(gene=character(),lipid=character(),sex=character(),zone=character(),diet=character(),r=numeric(),n=integer(),
      lipid_class=character(),lipid_label=character(),primary_lipid=logical(),fill_lipid=logical(),active_candidate=logical(),p=numeric(),
      loo_valid=integer(),loo_min=numeric(),loo_max=numeric(),loo_sign_fraction=numeric(),display=logical(),stable_display=logical(),
      bootstrap_low=numeric(),bootstrap_high=numeric(),bootstrap_sign_fraction=numeric(),bootstrap_valid=integer())
  } else a$edges <- do.call(rbind,all_edges)
  rownames(a$edges)<-NULL
  a$edges$q <- rep(NA_real_,nrow(a$edges)); k <- which(a$edges$primary_lipid)
  a$edges$q[k] <- stats::p.adjust(a$edges$p[k],"BH")
  a$edge_loo <- if(length(all_loo)) do.call(rbind,all_loo) else matrix(numeric(),0L,7L)
  a$bootstrap_indices <- banks
  a$network_eligibility <- if(length(eligible)) do.call(rbind,eligible) else data.frame(kind=character(),feature=character(),sex=character(),zone=character(),keep=logical())
  a$membership <- cf6_membership(a$edges,a$lipid_ledger)
  cats <- c("female_PP","female_PC","male_PP","male_PC","shared","unassigned")
  counts <- stats::aggregate(lipid~diet+lipid_class+category,a$membership,length)
  names(counts)[4] <- "n"
  totals <- stats::aggregate(lipid~diet+lipid_class,a$membership,length); names(totals)[3] <- "total"
  a$class_summary <- merge(merge(expand.grid(diet=c("SD","WD"),lipid_class=unique(a$membership$lipid_class),category=cats,stringsAsFactors=FALSE),counts,all.x=TRUE),totals)
  a$class_summary$n[is.na(a$class_summary$n)]<-0L
  a$class_summary$proportion<-a$class_summary$n/a$class_summary$total
  primary <- a$edges[a$edges$primary_lipid,,drop=FALSE]
  ranking <- unique(primary[,c("sex","zone","gene")])
  ranking$stable_degree <- vapply(seq_len(nrow(ranking)),function(j) {
    ee<-primary[primary$sex==ranking$sex[j] & primary$zone==ranking$zone[j] & primary$gene==ranking$gene[j],]
    length(unique(ee$lipid[ee$stable_display]))
  },integer(1L))
  ranking$q_interaction <- vapply(seq_len(nrow(ranking)),function(j) min(a$candidates$q_interaction[a$candidates$gene==ranking$gene[j] & a$candidates$sex==ranking$sex[j] & a$candidates$zone==ranking$zone[j]]),numeric(1L))
  ranking<-ranking[order(ranking$sex,ranking$zone,-ranking$stable_degree,ranking$q_interaction,ranking$gene),]
  ranking$rank<-ave(ranking$stable_degree,paste(ranking$sex,ranking$zone),FUN=seq_along)
  a$selection_ledger<-ranking
  a$representatives<-ranking[ranking$rank==1L & ranking$stable_degree>0,,drop=FALSE]
  rep_edges<-list(); examples<-list()
  for (j in seq_len(nrow(a$representatives))) {
    rr<-a$representatives[j,]
    ee<-primary[primary$sex==rr$sex & primary$zone==rr$zone & primary$gene==rr$gene,,drop=FALSE]
    stable_ids<-unique(ee$lipid[ee$stable_display])
    strength<-vapply(stable_ids,function(l) max(abs(ee$r[ee$lipid==l])),numeric(1L))
    partners<-head(stable_ids[order(-strength,stable_ids)],6L)
    ee<-ee[ee$lipid %in% partners,,drop=FALSE]; ee$partner_rank<-match(ee$lipid,partners)
    rep_edges[[j]]<-ee
    best<-partners[1]
    ids<-a$metadata$mouse_id[a$metadata$sex==rr$sex]
    pp<-a$profile_metadata$profile_id[match(paste(ids,rr$zone),paste(a$profile_metadata$mouse_id,a$profile_metadata$zone))]
    examples[[j]]<-data.frame(mouse_id=ids,sex=rr$sex,zone=rr$zone,diet=a$metadata$diet[match(ids,a$metadata$mouse_id)],gene=rr$gene,lipid=best,
      lipid_label=a$lipid_ledger$canonical_lipid[match(best,a$lipid_ledger$feature_id)],rna=as.numeric(a$modeled_rna[rr$gene,pp]),lipid_log2=as.numeric(a$lipid_matrix[best,ids]))
  }
  a$representative_edges<-if(length(rep_edges)) do.call(rbind,rep_edges) else primary[FALSE,]
  a$examples<-if(length(examples)) do.call(rbind,examples) else data.frame()
  sd<-primary[primary$diet=="SD",]; wd<-primary[primary$diet=="WD",]
  key<-function(e) paste(e$sex,e$zone,e$gene,e$lipid,sep="|")
  wd<-wd[match(key(sd),key(wd)),]; stopifnot(identical(key(sd),key(wd)))
  a$diet_edge_changes<-data.frame(sd[,c("sex","zone","gene","lipid","lipid_label","lipid_class")],r_SD=sd$r,r_WD=wd$r,delta_r=wd$r-sd$r,
    active_SD=sd$active_candidate,active_WD=wd$active_candidate,stable_SD=sd$stable_display,stable_WD=wd$stable_display)
  a$thresholds<-do.call(rbind,lapply(c(.5,.6,.65,.7,.8),function(th) {
    m<-cf6_membership(a$edges,a$lipid_ledger,th)
    total<-stats::aggregate(lipid~diet+lipid_class,m,length);names(total)[3]<-"total"
    m$assigned<-as.integer(m$n_contexts>0)
    v<-merge(stats::aggregate(assigned~diet+lipid_class,m,sum),total);v$threshold<-th;v$proportion<-v$assigned/v$total;v
  }))
  a$quality_sensitivity<-do.call(rbind,lapply(c("keep_primary","keep_broad","keep_fill"),function(quality) {
    m<-cf6_membership(a$edges,a$lipid_ledger,.65,quality);m$assigned<-as.integer(m$n_contexts>0)
    x<-stats::aggregate(cbind(assigned,total=rep(1,nrow(m)))~diet+lipid_class,m,sum)
    x$quality<-quality;x$proportion<-x$assigned/x$total;x
  }))
  a
}

cf6_save <- function(a,output) {
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  saveRDS(a,file.path(output,"analysis.rds"),compress="gzip")
  tables<-names(a)[vapply(a,is.data.frame,logical(1L))]
  for (n in tables) cf6_write(a[[n]],file.path(output,"tables",paste0(n,".csv")))
  cf6_write(data.frame(profile_id=rownames(a$design),a$design,check.names=FALSE),file.path(output,"tables","RNA_design.csv"))
  cf6_write(data.frame(cell=rownames(a$contrasts),a$contrasts,check.names=FALSE),file.path(output,"tables","RNA_contrasts.csv"))
  for (n in names(a$input_qc)) if (is.data.frame(a$input_qc[[n]])) cf6_write(a$input_qc[[n]],file.path(output,"QC",paste0(n,".csv")))
  writeLines(capture.output(sessionInfo()),file.path(output,"sessionInfo.txt"))
  files<-c(unname(a$input_paths),"scripts/core/cf453_standalone.R","scripts/core/cf453_figure6_plots.R","scripts/workflows/14_cf453_standalone_figure6.R","scripts/core/wdsd_intake.R","scripts/core/cyp2c38_reference.R")
  cf6_write(analysis_provenance(files),file.path(output,"provenance_sha256.csv"))
  saveRDS(a$bootstrap_indices,file.path(output,"bootstrap_indices.rds"))
  invisible(a)
}

# Exact, reviewable paired-animal draws for the displayed examples; diet groups
# are independent and both zones reuse each group's common index bank.
cf6_representative_draws <- function(a) {
  e <- a$representative_edges
  if (!nrow(e)) return(a)
  draws <- lapply(seq_len(nrow(e)),function(j) {
    ee<-e[j,]; ids<-a$metadata$mouse_id[a$metadata$sex==ee$sex & a$metadata$diet==ee$diet]
    pp<-a$profile_metadata$profile_id[match(paste(ids,ee$zone),paste(a$profile_metadata$mouse_id,a$profile_metadata$zone))]
    x<-as.numeric(a$modeled_rna[ee$gene,pp]);y<-as.numeric(a$lipid_matrix[ee$lipid,ids])
    bank<-a$bootstrap_indices[,,paste(ee$sex,ee$diet,sep="_")]
    values<-vapply(seq_len(ncol(bank)),function(k) {
      ii<-bank[,k];if(length(unique(x[ii]))<2L||length(unique(y[ii]))<2L)NA_real_ else stats::cor(x[ii],y[ii])
    },numeric(1L))
    data.frame(sex=ee$sex,zone=ee$zone,diet=ee$diet,gene=ee$gene,lipid=ee$lipid,draw=seq_along(values),r=values)
  })
  a$representative_bootstrap<-do.call(rbind,draws)
  sd<-a$representative_bootstrap[a$representative_bootstrap$diet=="SD",]
  wd<-a$representative_bootstrap[a$representative_bootstrap$diet=="WD",]
  key<-function(x)paste(x$sex,x$zone,x$gene,x$lipid,x$draw,sep="|")
  wd<-wd[match(key(sd),key(wd)),];stopifnot(identical(key(sd),key(wd)))
  a$representative_delta_draws<-data.frame(sd[,c("sex","zone","gene","lipid","draw")],delta_r=wd$r-sd$r)
  keys<-unique(a$representative_delta_draws[,c("sex","zone","gene","lipid")])
  a$representative_delta_summary<-do.call(rbind,lapply(seq_len(nrow(keys)),function(j) {
    keyj<-keys[j,];v<-a$representative_delta_draws
    v<-v$delta_r[v$sex==keyj$sex & v$zone==keyj$zone & v$gene==keyj$gene & v$lipid==keyj$lipid]
    observed<-a$diet_edge_changes;observed<-observed$delta_r[observed$sex==keyj$sex & observed$zone==keyj$zone & observed$gene==keyj$gene & observed$lipid==keyj$lipid]
    data.frame(keyj,delta_r=observed,low=stats::quantile(v,.025,na.rm=TRUE,names=FALSE),
      high=stats::quantile(v,.975,na.rm=TRUE,names=FALSE),valid=sum(is.finite(v)),requested=length(v),
      positive_fraction=mean(v>0,na.rm=TRUE),stringsAsFactors=FALSE)
  }))
  a
}

# Same gene-lipid universe in both diets; this complements the saturated
# any-edge lipid membership without selecting a new display threshold.
cf6_density_summaries <- function(a) {
  e<-a$edges[a$edges$primary_lipid,c("diet","sex","zone","lipid_class","r","q","loo_valid","loo_sign_fraction"),drop=FALSE]
  if(!nrow(e)) {
    a$class_pair_density<-data.frame(diet=character(),sex=character(),zone=character(),lipid_class=character(),
      total_pairs=integer(),display_pairs=integer(),stable_pairs=integer(),fdr05_pairs=integer(),
      mean_abs_r=numeric(),median_abs_r=numeric(),display_fraction=numeric(),stable_fraction=numeric(),fdr05_fraction=numeric())
    return(a)
  }
  groups<-split(seq_len(nrow(e)),interaction(e$diet,e$sex,e$zone,e$lipid_class,drop=TRUE))
  a$class_pair_density<-do.call(rbind,lapply(groups,function(ii) {
    d<-e[ii,,drop=FALSE];r<-abs(d$r)
    data.frame(d[1L,c("diet","sex","zone","lipid_class")],total_pairs=length(ii),
      display_pairs=sum(r>=.65),stable_pairs=sum(r>=.65 & d$loo_valid==7 & d$loo_sign_fraction==1),
      fdr05_pairs=sum(d$q<=.05,na.rm=TRUE),mean_abs_r=mean(r),median_abs_r=stats::median(r),
      display_fraction=mean(r>=.65),stable_fraction=mean(r>=.65 & d$loo_valid==7 & d$loo_sign_fraction==1),
      fdr05_fraction=mean(d$q<=.05,na.rm=TRUE),row.names=NULL)
  }))
  rownames(a$class_pair_density)<-NULL
  a
}
