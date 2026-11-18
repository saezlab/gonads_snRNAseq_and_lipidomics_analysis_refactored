source("scripts/config/analysis_release.R")
# Read Workflow 16 plan before editing. Pure analysis helpers; no execution on source.
cfdr_anchors <- function(a,direct) {
  tests<-a$rna_results; eligibility<-anchors<-list()
  for(d in c("SD","WD"))for(s in c("female","male")) {
    ix<-tests[tests$contrast==paste0("sex_zone_",d),]
    own<-tests[tests$contrast==paste0("zone_",s,"_",d),]
    other<-tests[tests$contrast==paste0("zone_",if(s=="female")"male" else "female","_",d),]
    own<-own[match(ix$gene,own$gene),];other<-other[match(ix$gene,other$gene),]
    q<-data.frame(gene=ix$gene,anchor_diet=d,sex=s,zone=ifelse(own$effect>0,"PP","PC"),
      anchor_effect=own$effect,anchor_q=ix$q,own_q=own$q,other_t=other$t,
      eligible=ix$q<=.05 & own$q<=.05 & abs(other$t)<1,stringsAsFactors=FALSE)
    q$context<-paste(q$sex,q$zone,sep="_")
    q$rank<-NA_integer_
    for(z in c("PP","PC")) {
      take<-which(q$eligible & q$zone==z)
      take<-take[order(-abs(q$anchor_effect[take]),q$anchor_q[take],q$own_q[take],q$gene[take])]
      q$rank[take]<-seq_along(take)
      anchors[[length(anchors)+1L]]<-q[head(take,5L),]
    }
    eligibility[[length(eligibility)+1L]]<-q
  }
  e<-do.call(rbind,eligibility);x<-do.call(rbind,anchors)
  for(d in c("SD","WD"))for(s in c("female","male")) {
    r<-tests[tests$contrast==paste0("zone_",s,"_",d),];i<-match(x$gene,r$gene)
    x[[paste0(if(s=="female")"x_" else "y_",d)]]<-r$t[i]
    x[[paste0("effect_",s,"_",d)]]<-r$effect[i]
  }
  dr<-direct$overview[direct$overview$contrast=="diet_sex_zone",];i<-match(x$gene,dr$gene)
  x$threeway_effect<-dr$effect[i];x$threeway_p<-dr$p[i];x$threeway_q<-dr$q[i]
  x$conditional_q<-vapply(seq_len(nrow(x)),function(i) {
    r<-direct$conditional; r$q[r$gene==x$gene[i]&r$contrast==paste0("diet_zone_",x$sex[i])]
  },numeric(1L))
  programmes<-unique(x[,c("context","sex","zone","gene")])
  programmes$anchor_diets<-vapply(seq_len(nrow(programmes)),function(i)paste(sort(unique(x$anchor_diet[x$context==programmes$context[i]&x$gene==programmes$gene[i]])),collapse=";"),character(1L))
  programmes$best_anchor_q<-vapply(seq_len(nrow(programmes)),function(i)min(x$anchor_q[x$context==programmes$context[i]&x$gene==programmes$gene[i]]),numeric(1L))
  bins<-expand.grid(anchor_diet=c("SD","WD"),sex=c("female","male"),zone=c("PP","PC"),stringsAsFactors=FALSE)
  bins$eligible<-vapply(seq_len(nrow(bins)),function(i)sum(e$eligible&e$anchor_diet==bins$anchor_diet[i]&e$sex==bins$sex[i]&e$zone==bins$zone[i]),integer(1L))
  bins$selected<-vapply(seq_len(nrow(bins)),function(i)sum(x$anchor_diet==bins$anchor_diet[i]&x$sex==bins$sex[i]&x$zone==bins$zone[i]),integer(1L))
  bins$requested<-5L;bins$status<-ifelse(bins$selected==5L,"Full bin","Underfilled; gates retained")
  stopifnot(all(bins$selected==pmin(bins$requested,bins$eligible)))
  list(anchors=x,eligibility=e,programmes=programmes,anchor_bins=bins)
}

cfdr_class_metrics <- function(r,classes,threshold=.65) {
  ids<-split(seq_along(classes),classes);G<-nrow(r);L<-ncol(r)
  hit<-colSums(abs(r)>=threshold); mag<-colSums(abs(r))
  n<-lengths(ids);h<-vapply(ids,function(j)sum(hit[j]),numeric(1L));m<-vapply(ids,function(j)sum(mag[j]),numeric(1L))
  rate<-h/(G*n);rest<-(sum(hit)-h)/(G*(L-n))
  data.frame(lipid_class=names(ids),n_genes=G,n_lipids=as.integer(n),n_pairs_class=G*n,n_pairs_rest=G*(L-n),
    hits=h,hits_rest=sum(hit)-h,rate=rate,rest=rest,E=rate-rest,
    mean_abs_excess=m/(G*n)-(sum(mag)-m)/(G*(L-n)),row.names=NULL)
}

cfdr_inputs <- function(a,genes,s,d,z) {
  mice<-a$metadata$mouse_id[a$metadata$sex==s&a$metadata$diet==d]
  pm<-a$profile_metadata
  profiles<-pm$profile_id[match(paste(mice,z),paste(pm$mouse_id,pm$zone))]
  stopifnot(length(mice)==7L,!anyNA(profiles),all(mice %in% colnames(a$lipid_matrix)))
  list(x=a$modeled_rna[genes,profiles,drop=FALSE],y=a$lipid_matrix[,mice,drop=FALSE],mice=mice,profiles=profiles)
}

cfdr_analyse <- function(parent=analysis_current("cf453_models/models.rds"),direct_path=analysis_current("cf453_direct_diet/diet_analysis.rds"),B=999L,seed=20261102L,provenance_path=NULL) {
  a<-readRDS(parent);direct<-readRDS(direct_path)
  stopifnot(nrow(a$metadata)==28L,nrow(a$profile_metadata)==56L,nrow(a$lipid_matrix)>0L,
    identical(a$metadata,direct$metadata),identical(a$modeled_rna,direct$modeled_rna))
  if(is.null(provenance_path))provenance_path<-file.path(dirname(parent),"provenance_sha256.csv")
  prior<-utils::read.csv(provenance_path,stringsAsFactors=FALSE)
  paths<-unname(a$input_paths);expected<-prior$sha256[match(paths,prior$path)]
  input_verification<-analysis_verify_inputs(paths,expected)
  out<-cfdr_anchors(a,direct)
  out$input_verification<-input_verification
  if(any(out$anchor_bins$selected<out$anchor_bins$requested))
    message("Some anchor bins are underfilled; all original eligibility gates remain unchanged.")
  print(stats::aggregate(eligible~anchor_diet+sex+zone,out$eligibility,sum))
  cat("Selected programme memberships:",nrow(out$programmes),"unique genes:",length(unique(out$programmes$gene)),"\n")
  ledger<-a$lipid_ledger[match(rownames(a$lipid_matrix),a$lipid_ledger$feature_id),]
  stopifnot(!anyNA(ledger$feature_id),all(ledger$keep_primary))
  classes<-ledger$lipid_class;class_names<-sort(unique(classes))
  contexts<-c("female_PP","female_PC","male_PP","male_PC")
  groups<-c("female_SD","female_WD","male_SD","male_WD")
  set.seed(seed);bank<-array(sample.int(7L,7L*B*4L,replace=TRUE),c(7L,B,4L),dimnames=list(NULL,NULL,groups))
  all_edges<-summaries<-sensitivity<-broad<-boots<-loos<-sample_order<-list()
  for(context in contexts) {
    s<-sub("_[^_]+$","",context);z<-sub("^.*_","",context)
    genes<-sort(out$programmes$gene[out$programmes$context==context])
    G<-length(genes);L<-length(classes)
    edge<-data.frame(context=context,sex=s,zone=z,gene=rep(genes,times=L),lipid=rep(rownames(a$lipid_matrix),each=G),
      lipid_label=rep(ledger$canonical_lipid,each=G),lipid_class=rep(classes,each=G),stringsAsFactors=FALSE)
    inp<-obs<-met<-loo<-list()
    for(d in c("SD","WD")) {
      inp[[d]]<-cfdr_inputs(a,genes,s,d,z);v<-inp[[d]]
      sample_order[[length(sample_order)+1L]]<-data.frame(context=context,diet=d,mouse_id=v$mice,profile_id=v$profiles,index=seq_len(7L))
      obs[[d]]<-cf6_cor(v$x,v$y)
      if(any(!is.finite(obs[[d]])))stop("Constant full-data feature: fixed universe cannot be evaluated")
      edge[[paste0("r_",d)]]<-as.vector(obs[[d]])
      met[[d]]<-cfdr_class_metrics(obs[[d]],classes)
      loo[[d]]<-vapply(seq_len(7L),function(j)as.vector(cf6_cor(v$x[,-j,drop=FALSE],v$y[,-j,drop=FALSE])),numeric(nrow(edge)))
      for(th in c(.60,.65,.70,.80)) {
        mm<-cfdr_class_metrics(obs[[d]],classes,th);mm$context<-context;mm$diet<-d;mm$threshold<-th
        sensitivity[[length(sensitivity)+1L]]<-mm
      }
    }
    edge$delta_r<-edge$r_WD-edge$r_SD;edge$delta_abs<-abs(edge$r_WD)-abs(edge$r_SD)
    gain<-abs(edge$r_WD)>=.8 & abs(edge$r_SD)<=.4 & edge$delta_abs>=.5
    loss<-abs(edge$r_SD)>=.8 & abs(edge$r_WD)<=.4 & edge$delta_abs<=-.5
    edge$direction<-ifelse(gain,"gained",ifelse(loss,"lost","other"))
    robust<-rep(FALSE,nrow(edge));sep<-rep(NA_real_,nrow(edge))
    for(dir in c("gained","lost")) {
      ix<-which(edge$direction==dir);high<-if(dir=="gained")"WD" else "SD";low<-if(dir=="gained")"SD" else "WD"
      for(i in ix) {
        h<-loo[[high]][i,];l<-loo[[low]][i,]
        if(all(is.finite(c(h,l)))) {
          sep[i]<-min(abs(h))-max(abs(l))
          robust[i]<-all(abs(h)>=.65 & sign(h)==sign(edge[[paste0("r_",high)]][i])) && all(abs(l)<=.65)
        }
      }
    }
    edge$loo_robust<-robust;edge$loo_separation_min<-sep
    candidate<-which(robust);R<-length(candidate)
    boot_E<-boot_r<-list();valid<-list()
    for(d in c("SD","WD")) {
      boot_E[[d]]<-matrix(NA_real_,length(class_names),B,dimnames=list(class_names,NULL))
      boot_r[[d]]<-matrix(NA_real_,R,B);valid[[d]]<-rep(FALSE,B)
      v<-inp[[d]];index<-bank[,,paste(s,d,sep="_")]
      for(b in seq_len(B)) {
        rr<-cf6_cor(v$x[,index[,b],drop=FALSE],v$y[,index[,b],drop=FALSE])
        if(R)boot_r[[d]][,b]<-as.vector(rr)[candidate]
        if(all(is.finite(rr))) {
          valid[[d]][b]<-TRUE
          boot_E[[d]][,b]<-cfdr_class_metrics(rr,classes)$E
        }
      }
    }
    delta<-boot_E$WD-boot_E$SD
    sm<-met$SD;names(sm)[names(sm)%in%c("hits","hits_rest","rate","rest","E","mean_abs_excess")]<-paste0(c("hits","hits_rest","rate","rest","E","mean_abs_excess"),"_SD")
    for(nm in c("hits","hits_rest","rate","rest","E","mean_abs_excess"))sm[[paste0(nm,"_WD")]]<-met$WD[[nm]]
    sm$delta_E<-sm$E_WD-sm$E_SD;sm$delta_rate<-sm$rate_WD-sm$rate_SD
    sm$delta_mean_abs_excess<-sm$mean_abs_excess_WD-sm$mean_abs_excess_SD
    sm$context<-context;sm$sex<-s;sm$zone<-z
    qfun<-function(x,p)if(any(is.finite(x)))as.numeric(stats::quantile(x,p,na.rm=TRUE)) else NA_real_
    sm$bootstrap_low<-apply(delta,1,qfun,p=.025);sm$bootstrap_high<-apply(delta,1,qfun,p=.975)
    sm$bootstrap_valid<-rowSums(is.finite(delta));sm$bootstrap_requested<-B
    sm$bootstrap_direction<-vapply(seq_len(nrow(sm)),function(i) {
      x<-delta[i,is.finite(delta[i,])];if(!length(x)||sm$delta_E[i]==0)NA_real_ else mean(sign(x)==sign(sm$delta_E[i]))
    },numeric(1L))
    edge$bootstrap_valid<-NA_integer_;edge$bootstrap_requested<-B;edge$bootstrap_direction<-NA_real_
    edge$bootstrap_low<-edge$bootstrap_high<-NA_real_
    if(R) {
      dd<-abs(boot_r$WD)-abs(boot_r$SD)
      edge$bootstrap_valid[candidate]<-rowSums(is.finite(dd))
      edge$bootstrap_direction[candidate]<-vapply(seq_len(R),function(j) {
        x<-dd[j,is.finite(dd[j,])];if(!length(x))NA_real_ else mean(sign(x)==sign(edge$delta_abs[candidate[j]]))
      },numeric(1L))
      edge$bootstrap_low[candidate]<-apply(dd,1,qfun,p=.025);edge$bootstrap_high[candidate]<-apply(dd,1,qfun,p=.975)
    }
    edge$network_eligible<-edge$loo_robust & !is.na(edge$bootstrap_direction) & edge$bootstrap_direction>=.90 & edge$bootstrap_valid>=ceiling(.90*B)
    for(d in c("SD","WD")) {
      rr<-edge[[paste0("r_",d)]];den<-pmax(1-rr^2,.Machine$double.eps)
      edge[[paste0("p_",d)]]<-2*stats::pt(-abs(rr)*sqrt(5/den),df=5)
    }
    bg<-sort(unique(out$eligibility$gene[out$eligibility$eligible&out$eligibility$context==context]))
    for(d in c("SD","WD")) {
      v<-cfdr_inputs(a,bg,s,d,z);mm<-cfdr_class_metrics(cf6_cor(v$x,v$y),classes);mm$context<-context;mm$diet<-d
      broad[[length(broad)+1L]]<-mm
    }
    all_edges[[context]]<-edge;summaries[[context]]<-sm;loos[[context]]<-loo
    boots[[context]]<-list(E=boot_E,valid=valid,candidate_rows=candidate,candidate_r=boot_r)
    cat(context,":",G,"genes;",sum(gain),"gains",sum(loss),"losses;",R,"LOO candidates;",sum(edge$network_eligible),"bootstrap-supported display edges; class draws",sum(is.finite(delta[1,])),"/",B,"\n")
  }
  edges<-do.call(rbind,all_edges);row.names(edges)<-NULL
  np<-nrow(edges);qq<-stats::p.adjust(c(edges$p_SD,edges$p_WD),"BH")
  edges$q_SD<-qq[seq_len(np)];edges$q_WD<-qq[np+seq_len(np)]
  networks<-ngenes<-list()
  for(context in contexts) {
    e<-edges[edges$context==context&edges$network_eligible,];pr<-out$programmes[out$programmes$context==context,]
    s<-sub("_[^_]+$","",context);z<-sub("^.*_","",context)
    if(!nrow(e)) {
      ngenes[[context]]<-data.frame(context=context,sex=s,zone=z,gene=NA_character_,n_edges=0L,status="No example meets all display gates")
      next
    }
    gs<-unique(e$gene);ranked<-data.frame(gene=gs,n=vapply(gs,function(g)sum(e$gene==g),integer(1L)),
      max_change=vapply(gs,function(g)max(abs(e$delta_abs[e$gene==g])),numeric(1L)),q=pr$best_anchor_q[match(gs,pr$gene)])
    ranked<-ranked[order(-ranked$n,-ranked$max_change,ranked$q,ranked$gene),];g<-ranked$gene[1L]
    v<-e[e$gene==g,];v<-v[order(-abs(v$delta_abs),v$lipid),]
    first<-c(head(which(v$direction=="gained"),3L),head(which(v$direction=="lost"),3L))
    use<-head(c(first,setdiff(seq_len(nrow(v)),first)),6L);v<-v[use,]
    v$display_rank<-seq_len(nrow(v));networks[[context]]<-v
    ngenes[[context]]<-data.frame(context=context,sex=s,zone=z,gene=g,n_edges=nrow(v),status="Exploratory selected example")
  }
  out$edges<-edges;out$class_summary<-do.call(rbind,summaries);out$threshold_sensitivity<-do.call(rbind,sensitivity)
  out$broad_sensitivity<-do.call(rbind,broad);out$network_genes<-do.call(rbind,ngenes)
  out$network_edges<-if(length(networks))do.call(rbind,networks) else edges[FALSE,]
  out$sample_order<-do.call(rbind,sample_order);out$metadata<-a$metadata
  out$bootstrap<-boots;out$bootstrap_indices<-bank;out$loo<-loos
  out$source_paths<-c(parent=normalizePath(parent),direct=normalizePath(direct_path),parent_provenance=normalizePath(provenance_path),a$input_paths)
  out$settings<-list(B=B,seed=seed,n_lipids=nrow(a$lipid_matrix),threshold=.65,anchor_rank="decreasing absolute own zonation effect; interaction q; own q; gene",
    correlation_FDR="BH across all selected programme-context x diet x gene x lipid coefficients",threeway_FDR="Workflow15 genes x three direct diet contrasts")
  out
}

cfdr_save <- function(a,output) {
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  saveRDS(a,file.path(output,"rewiring_analysis.rds"));saveRDS(a$bootstrap_indices,file.path(output,"bootstrap_indices.rds"))
  for(nm in names(a)[vapply(a,is.data.frame,logical(1L))])cf6_write(a[[nm]],file.path(output,"tables",paste0(nm,".csv")))
  writeLines(capture.output(sessionInfo()),file.path(output,"sessionInfo.txt"))
  cfdr_provenance(a,output)
  invisible(a)
}

cfdr_provenance <- function(a,output) {
  code<-c("scripts/core/cf453_diet_rewiring.R","scripts/core/cf453_diet_rewiring_plots.R","scripts/core/cf453_diet_rewiring_report.R",
    "scripts/core/cf453_standalone.R","scripts/workflows/16_cf453_diet_rewiring_figure6.R","tests/qa_cf453_diet_rewiring.R",
    "scripts/reports/package_cf453_diet_rewiring.py","docs/workflow_plans/16_cf453_diet_rewiring_figure6.md")
  files<-c(unname(a$source_paths),code[file.exists(code)])
  cf6_write(analysis_provenance(files),file.path(output,"provenance_sha256.csv"))
}
