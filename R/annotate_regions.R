# ============================================================================
# annotate_genomic_regions — label ORFs by their mitogenomic context
# ============================================================================

#' Annotate ORFs with Their Mitogenomic Region
#'
#' Adds a `genomic_region` column to an ORF table (output of [find_orfs()])
#' by mapping each ORF's start-codon position against a GenBank feature table
#' (output of [fetch_mito_genbank()]).
#'
#' Region labels follow this priority (highest wins when features overlap):
#' 1. `CDS:<gene>` — protein-coding gene
#' 2. `tRNA:<product>` — transfer RNA
#' 3. `<gene>` for rRNA features (e.g. `12S_rRNA`, `16S_rRNA`)
#' 4. `D-loop` — control / displacement loop
#' 5. Name of any other annotated feature
#' 6. `intergenic` — falls outside all annotated features
#'
#' @param orfs A `data.frame` returned by [find_orfs()].  Must contain
#'   `start` and `strand` columns.
#' @param features A `data.frame` of GenBank features as returned in the
#'   `$features` element of [fetch_mito_genbank()].  Must contain `type`,
#'   `start`, `end`, `strand`, `gene`, and `product` columns.
#' @param seq_len Integer.  Total length of the genome in nucleotides.
#'   Required only to handle ORFs that wrap around the origin
#'   (`wraps_around == TRUE`); ignored otherwise.  Defaults to `NULL`
#'   (wrap-around ORFs are annotated by their `start` coordinate only).
#'
#' @return The input `orfs` data.frame with a `genomic_region` character
#'   column appended.
#'
#' @examples
#' \dontrun{
#' genomes <- fetch_mito_genbank(c("NC_002333.2", "NC_028087.1"))
#'
#' orfs_danio <- find_orfs(genomes[["NC_002333.2"]]$sequence)
#' orfs_danio$seq_id <- "NC_002333.2"
#' orfs_danio <- annotate_genomic_regions(
#'   orfs_danio,
#'   genomes[["NC_002333.2"]]$features,
#'   seq_len = genomes[["NC_002333.2"]]$length
#' )
#' table(orfs_danio$genomic_region)
#' }
#'
#' @export
annotate_genomic_regions <- function(orfs, features, seq_len = NULL) {

  if (nrow(orfs) == 0L) {
    orfs$genomic_region <- character(0L)
    return(orfs)
  }

  if (nrow(features) == 0L) {
    orfs$genomic_region <- "intergenic"
    return(orfs)
  }

  ## ---- build region label for each feature ---------------------------------
  region_label <- function(f) {
    switch(f$type,
      CDS    = paste0("CDS:", if (!is.na(f$gene)) f$gene else
                              if (!is.na(f$product)) f$product else "unknown"),
      tRNA   = paste0("tRNA:", if (!is.na(f$product)) f$product else
                                if (!is.na(f$gene))    f$gene   else "unknown"),
      rRNA   = {
        gene_val <- if (!is.na(f$gene))    f$gene    else ""
        prod_val <- if (!is.na(f$product)) f$product else ""
        combined <- paste(gene_val, prod_val)
        if (grepl("12S|rnr1|RNR1|s-rRNA|small", combined, ignore.case = TRUE))
          "12S_rRNA"
        else if (grepl("16S|rnr2|RNR2|l-rRNA|large", combined, ignore.case = TRUE))
          "16S_rRNA"
        else if (nzchar(gene_val)) gene_val
        else if (nzchar(prod_val)) prod_val
        else "rRNA"
      },
      `D-loop`     = "D-loop",
      `rep_origin` = "D-loop",
      misc_feature = if (!is.na(f$note))    f$note    else "misc_feature",
      gene         = NULL,   # skip bare 'gene' features; CDS/rRNA/tRNA cover them
      source       = NULL,
      if (!is.na(f$gene)) f$gene else f$type
    )
  }

  ## ---- priority score: lower = higher priority -----------------------------
  priority <- function(type) {
    switch(type,
      CDS    = 1L,
      tRNA   = 2L,
      rRNA   = 3L,
      `D-loop`     = 4L,
      `rep_origin` = 4L,
      5L
    )
  }

  ## ---- build lookup: position → (label, priority) -------------------------
  feat_df <- features[features$type != "source" & features$type != "gene", ]

  # Annotate each ORF by its start codon position
  orfs$genomic_region <- vapply(seq_len(nrow(orfs)), function(i) {

    pos <- orfs$start[i]

    # Find all features that contain this position
    hits <- which(feat_df$start <= pos & feat_df$end >= pos)

    if (length(hits) == 0L) {
      # For wrap-around ORFs, also check using end coordinate
      if (isTRUE(orfs$wraps_around[i])) {
        pos2 <- orfs$end[i]
        hits <- which(feat_df$start <= pos2 & feat_df$end >= pos2)
      }
      if (length(hits) == 0L) return("intergenic")
    }

    # Pick highest-priority feature
    prios  <- vapply(hits, function(h) priority(feat_df$type[h]), integer(1L))
    best_h <- hits[which.min(prios)]
    lbl    <- region_label(as.list(feat_df[best_h, ]))
    if (is.null(lbl)) "intergenic" else lbl

  }, character(1L))

  orfs
}
