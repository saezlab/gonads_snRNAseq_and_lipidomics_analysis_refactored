# Resolve current analyses and immutable baselines through one release registry.
# Paths in saved provenance remain identifiers; resolve them only when reading.
analysis_release <- function(root=".") {
  root<-normalizePath(root,mustWork=TRUE)
  path<-file.path(root,"scripts/config/analysis_release.json")
  if(!file.exists(path))stop("Missing analysis release registry: ",path)
  if(!requireNamespace("jsonlite",quietly=TRUE))stop("jsonlite is required for the analysis release registry")
  x<-jsonlite::fromJSON(path,simplifyVector=FALSE)
  if(is.null(x$schema_version)||as.integer(x$schema_version)!=1L||!is.list(x$paths))
    stop("Unsupported analysis release registry schema")
  attr(x,"root")<-root;x
}

analysis_path <- function(key,root=".",mustWork=FALSE) {
  r<-analysis_release(root);p<-r$paths[[key]]
  if(is.null(p)&&key=="engine_fcg05")p<-"scripts/core/engines/05_zonation_analysis_legacy.R"
  if(is.null(p)&&key=="engine_validation07")p<-"scripts/core/engines/07_validation_experiments_legacy.R"
  if(length(p)!=1L||is.na(p)||!nzchar(p))stop("Missing release path key: ",key)
  if(!grepl("^(/|[A-Za-z]:)",p))p<-file.path(attr(r,"root"),p)
  if(mustWork&&!file.exists(p))stop("Missing registered analysis path ",key,": ",p)
  p
}

analysis_current <- function(component="",root=".",mustWork=FALSE) {
  p<-analysis_path("current_root",root)
  if(nzchar(component)){
    layout<-analysis_release(root)$current_components
    branch<-sub("/.*$","",component)
    if(!is.null(layout[[branch]]))component<-paste0(layout[[branch]],substring(component,nchar(branch)+1L))
    p<-file.path(p,component)
  }
  if(mustWork&&!file.exists(p))stop("Missing current-release artifact: ",p)
  p
}

# Protect the approved release independently of ignored cache/checkpoint files.
# Resolve existing ancestors too, so a new child through a symlink is protected.
# Reject dangling symlinks whose eventual write destination cannot be resolved.
analysis_assert_not_approved_output <- function(output,root=".") {
  stopifnot(length(output)==1L,!is.na(output),nzchar(output))
  canonical<-function(path){
    path<-path.expand(path)
    if(!grepl("^(/|[A-Za-z]:)",path))path<-file.path(normalizePath(root,mustWork=TRUE),path)
    suffix<-character();ancestor<-path
    while(!file.exists(ancestor)&&!dir.exists(ancestor)){
      suffix<-c(basename(ancestor),suffix);parent<-dirname(ancestor)
      if(identical(parent,ancestor))stop("Cannot resolve output path: ",path)
      ancestor<-parent
    }
    resolved<-normalizePath(ancestor,winslash="/",mustWork=TRUE)
    for(part in suffix){
      if(part%in%c("","."))next
      resolved<-if(part=="..")dirname(resolved) else file.path(resolved,part)
      link<-Sys.readlink(resolved)
      if(!is.na(link)&&nzchar(link)&&!file.exists(resolved)&&!dir.exists(resolved))
        stop("Cannot safely resolve dangling output symlink: ",resolved)
      if(file.exists(resolved)||dir.exists(resolved))
        resolved<-normalizePath(resolved,winslash="/",mustWork=TRUE)
    }
    resolved
  }
  target<-canonical(output)
  release<-analysis_release(root)
  protected<-c(release$paths$current_root,release$paths$legacy_results_root)
  for(path in vapply(protected,canonical,character(1L),USE.NAMES=FALSE)){
    if(identical(target,path)||startsWith(target,paste0(path,"/")))
      stop("Preserve the approved release and legacy archive: scientific output must use a fresh directory outside ",path)
  }
  invisible(target)
}

analysis_resolve_recorded_path <- function(path,root=".",mustWork=FALSE) {
  r<-analysis_release(root);repo<-attr(r,"root")
  resolve_one<-function(p){
    p<-gsub("\\\\","/",p)
    rel<-if(startsWith(p,paste0(repo,"/")))substring(p,nchar(repo)+2L) else p
    # Relocate absolute records from a prior checkout without rewriting them.
    if(grepl("^(/|[A-Za-z]:)",rel)){
      candidate<-sub("^.*?/(data|results|scripts|tests|docs|manuscript|miscellaneous)/","\\1/",rel,perl=TRUE)
      if(!identical(candidate,rel))rel<-candidate
    }
    moves<-r$relocations
    if(length(moves)){
      moves<-moves[order(vapply(moves,function(z)nchar(z$from),integer(1L)),decreasing=TRUE)]
      for(move in moves){
        from<-sub("/$","",move$from)
        if(identical(rel,from)||startsWith(rel,paste0(from,"/"))){
          suffix<-substring(rel,nchar(from)+1L)
          rel<-paste0(sub("/$","",move$to),suffix)
          if(!is.null(move$fallback)&&!file.exists(file.path(repo,rel))){
            archived<-paste0(sub("/$","",move$fallback),suffix)
            if(file.exists(file.path(repo,archived)))rel<-archived
          }
          break
        }
      }
    }
    actual<-if(grepl("^(/|[A-Za-z]:)",rel))rel else file.path(repo,rel)
    if(mustWork&&!file.exists(actual))stop("Recorded provenance path cannot be resolved: ",p," -> ",actual)
    actual
  }
  vapply(path,resolve_one,character(1L),USE.NAMES=FALSE)
}

# CSV manifests retain the source RDS identity; that hash is never substituted
# for a CSV file hash. Verify the complete lossless export before using a bridge.
analysis_pseudobulk_bridge <- function(path,expected_sha256=NULL,root=".") {
  if(!grepl("[.]rds$",path,ignore.case=TRUE))return(NULL)
  # Unbound provenance must describe an existing explicit input's actual bytes.
  # Frozen-cache verification supplies its expected identity and remains CSV-first.
  if(is.null(expected_sha256)){
    local_file<-analysis_resolve_recorded_path(path,root)
    if(file.exists(local_file)){
      if(dir.exists(local_file))stop("Recorded RDS input is a directory: ",path)
      expected_sha256<-digest::digest(local_file,file=TRUE,algo="sha256")
    }
  }
  r<-analysis_release(root);keys<-grep("_transport$",names(r$paths),value=TRUE)
  candidates<-lapply(keys,function(k){
    f<-analysis_path(k,root)
    if(!file.exists(f))return(NULL)
    m<-jsonlite::read_json(f,simplifyVector=FALSE)
    if(!identical(m$format,"pseudobulk_csv"))return(NULL)
    list(path=f,manifest=m)
  })
  candidates<-Filter(Negate(is.null),candidates)
  matches<-vapply(candidates,function(z){
    if(!is.null(expected_sha256))return(identical(z$manifest$source$sha256,as.character(expected_sha256)))
    actual<-analysis_resolve_recorded_path(path,root)
    original<-analysis_resolve_recorded_path(z$manifest$source$path,root)
    identical(actual,original)||identical(basename(actual),basename(original))
  },logical(1L))
  candidates<-candidates[matches]
  if(!length(candidates))return(NULL)
  identities<-unique(vapply(candidates,function(z)z$manifest$source$sha256,character(1L)))
  if(length(identities)!=1L)stop("Ambiguous pseudobulk transport identity: ",path)
  bridge<-candidates[[1L]]
  helpers<-new.env(parent=baseenv())
  sys.source(file.path(attr(r,"root"),"scripts/core/pseudobulk_csv.R"),envir=helpers)
  x<-helpers$pb_components(bridge$path,verify=TRUE)
  stopifnot(identical(x$source$sha256,identities))
  bridge$files<-x$source_paths
  bridge
}

analysis_verify_inputs <- function(paths,expected_sha256,root=".") {
  stopifnot(length(paths)==length(expected_sha256))
  rows<-lapply(seq_along(paths),function(i){
    p<-unname(paths[i]);expected<-unname(as.character(expected_sha256[i]))
    if(is.na(expected)||!nzchar(expected))stop("Missing recorded input hash: ",p)
    bridge<-analysis_pseudobulk_bridge(p,expected,root)
    if(!is.null(bridge))return(data.frame(path=p,expected_sha256=expected,verified=TRUE,
      verification="verified_lossless_pseudobulk_CSV_export",opened_path=normalizePath(bridge$path),
      actual_file_sha256=digest::digest(bridge$path,file=TRUE,algo="sha256"),
      source_RDS_sha256=bridge$manifest$source$sha256,stringsAsFactors=FALSE))
    actual<-analysis_resolve_recorded_path(p,root,mustWork=TRUE)
    hash<-digest::digest(actual,file=TRUE,algo="sha256")
    if(!identical(hash,expected))stop("Input hash changed: ",p)
    data.frame(path=p,expected_sha256=expected,verified=TRUE,verification="direct_file_bytes",
      opened_path=normalizePath(actual),actual_file_sha256=hash,source_RDS_sha256=NA_character_,stringsAsFactors=FALSE)
  })
  do.call(rbind,rows)
}

# Keep recorded identities while separately recording the bytes actually opened.
analysis_provenance <- function(paths,root=".") {
  rows<-lapply(unname(as.character(paths)),function(p){
    bridge<-analysis_pseudobulk_bridge(p,root=root)
    if(!is.null(bridge)){
      identity<-data.frame(path=p,resolved_path=normalizePath(bridge$path),sha256=bridge$manifest$source$sha256,
        sha256_scope="recorded_source_RDS_identity",verification="verified_lossless_pseudobulk_CSV_export",stringsAsFactors=FALSE)
      files<-normalizePath(bridge$files,mustWork=TRUE)
      return(rbind(identity,data.frame(path=files,resolved_path=files,
        sha256=unname(vapply(files,digest::digest,character(1L),file=TRUE,algo="sha256")),
        sha256_scope="file_bytes",verification="direct_file_bytes",stringsAsFactors=FALSE)))
    }
    actual<-analysis_resolve_recorded_path(p,root,mustWork=TRUE)
    data.frame(path=p,resolved_path=normalizePath(actual),sha256=digest::digest(actual,file=TRUE,algo="sha256"),
      sha256_scope="file_bytes",verification="direct_file_bytes",stringsAsFactors=FALSE)
  })
  unique(do.call(rbind,rows))
}
