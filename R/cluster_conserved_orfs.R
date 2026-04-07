# ============================================================================
# cluster_conserved_orfs — group ORFs conserved across multiple mitogenomes
# ============================================================================

#' Cluster ORFs Conserved Across Species
#'
#' Identifies ORFs that are shared across two or more mitochondrial genomes
#' by performing all-vs-all pairwise local alignment of protein sequences
#' from different species.  ORFs connected by sufficient similarity are
#' grouped into clusters; only clusters represented in at least two distinct
#' genomes (identified by `seq_id`) are retained as *conserved*.
#'
#' @param orfs A `data.frame` of ORFs — typically the combined output of
#'   [find_orfs()] run on two or more genomes, with a `seq_id` column added
#'   to identify the source genome.  Must contain `seq_id` and
#'   `protein_sequence` columns.
#' @param min_pct_identity Numeric.  Minimum percent identity of the local
#'   alignment (as returned by [pwalign::pid()]) required to connect two
#'   ORFs.  Default `40`.
#' @param min_coverage Numeric in \[0, 1\].  Minimum fraction of the
#'   *shorter* protein that must be covered by the alignment (gap characters
#'   excluded).  Default `0.5` (50 %).
#' @param substitution_matrix Character.  Substitution matrix name passed to
#'   [pwalign::pairwiseAlignment()].  Default `"BLOSUM62"`.
#' @param gap_opening Numeric.  Gap-opening penalty.  Default `10`.
#' @param gap_extension Numeric.  Gap-extension penalty.  Default `4`.
#' @param verbose Logical.  Print progress?  Default `TRUE`.
#'
#' @return The input `orfs` data.frame, filtered to conserved ORFs only,
#'   with three additional columns:
#' \describe{
#'   \item{`cluster_id`}{Integer cluster identifier.  Clusters are numbered
#'     in order of their first appearance (sorted by `start`, then `seq_id`).}
#'   \item{`n_sequences`}{Number of distinct `seq_id` values in the cluster.}
#'   \item{`cluster_size`}{Total number of ORF rows in the cluster.}
#' }
#'
#' @details
#' ## Algorithm
#' 1. Strip trailing stop-codon asterisks from `protein_sequence`.
#' 2. For every pair of ORFs from *different* `seq_id` values, compute a
#'    Smith–Waterman local alignment.
#' 3. Retain pairs meeting both `min_pct_identity` and `min_coverage`
#'    thresholds.
#' 4. Build an undirected graph; nodes = ORFs, edges = retained pairs.
#' 5. Extract connected components; discard singleton or single-species
#'    components.
#' 6. Return annotated surviving ORFs sorted by `cluster_id`, `seq_id`,
#'    `start`.
#'
#' @examples
#' \dontrun{
#' library(Biostrings)
#'
#' genomes <- fetch_mito_genbank(c("NC_002333.2", "NC_028087.1"))
#'
#' # Find ORFs in each genome and tag with seq_id
#' orfs_list <- lapply(names(genomes), function(acc) {
#'   df        <- find_orfs(genomes[[acc]]$sequence, min_orf_length = 75L)
#'   df$seq_id <- acc
#'   df
#' })
#' orfs_all <- do.call(rbind, orfs_list)
#'
#' conserved <- cluster_conserved_orfs(orfs_all)
#' }
#'
#' @importFrom pwalign pairwiseAlignment pid
#' @importFrom Biostrings AAString
#' @importFrom igraph graph_from_edgelist components
#' @export
cluster_conserved_orfs <- function(orfs,
                                    min_pct_identity    = 40,
                                    min_coverage        = 0.5,
                                    substitution_matrix = "BLOSUM62",
                                    gap_opening         = 10,
                                    gap_extension       = 4,
                                    verbose             = TRUE) {

  ## ---- input checks --------------------------------------------------------
  needed <- c("seq_id", "protein_sequence")
  miss   <- setdiff(needed, names(orfs))
  if (length(miss))
    stop("'orfs' is missing required columns: ",
         paste(miss, collapse = ", "), call. = FALSE)

  seq_ids <- unique(orfs$seq_id)
  if (length(seq_ids) < 2L)
    stop("'orfs' must contain ORFs from at least two distinct seq_id values.",
         call. = FALSE)

  ## ---- strip stop-codon asterisks ------------------------------------------
  proteins <- gsub("\\*$", "", orfs$protein_sequence)
  n        <- nrow(orfs)

  ## ---- index ORFs by seq_id ------------------------------------------------
  id_groups <- split(seq_len(n), orfs$seq_id)

  ## ---- all cross-species pairs --------------------------------------------
  pair_combos <- utils::combn(names(id_groups), 2L, simplify = FALSE)

  edges <- vector("list", length(pair_combos))

  for (pc in seq_along(pair_combos)) {
    s1   <- pair_combos[[pc]][1L]
    s2   <- pair_combos[[pc]][2L]
    idx1 <- id_groups[[s1]]
    idx2 <- id_groups[[s2]]

    if (verbose)
      message(sprintf("  Aligning %s (%d ORFs) vs %s (%d ORFs) ...",
                      s1, length(idx1), s2, length(idx2)))

    edge_rows <- vector("list", length(idx1) * length(idx2))
    k <- 0L

    for (i in idx1) {
      prot_i <- proteins[i]
      len_i  <- nchar(prot_i)
      if (len_i < 2L) next
      aa_i   <- Biostrings::AAString(prot_i)

      for (j in idx2) {
        prot_j <- proteins[j]
        len_j  <- nchar(prot_j)
        if (len_j < 2L) next

        aln <- pwalign::pairwiseAlignment(
          pattern             = Biostrings::AAString(prot_j),
          subject             = aa_i,
          type                = "local",
          substitutionMatrix  = substitution_matrix,
          gapOpening          = gap_opening,
          gapExtension        = gap_extension
        )

        pct_id      <- pwalign::pid(aln, type = "PID1")
        # Non-gap aligned residues for each sequence
        n_sub  <- nchar(gsub("-", "", as.character(pwalign::alignedSubject(aln))))
        n_pat  <- nchar(gsub("-", "", as.character(pwalign::alignedPattern(aln))))
        min_len  <- min(len_i, len_j)
        # Coverage = aligned residues of the shorter sequence / its length
        coverage <- max(n_sub, n_pat) / min_len

        if (pct_id >= min_pct_identity && coverage >= min_coverage) {
          k <- k + 1L
          edge_rows[[k]] <- c(i, j)
        }
      }
    }

    if (k > 0L)
      edges[[pc]] <- do.call(rbind, edge_rows[seq_len(k)])
  }

  edge_mat <- do.call(rbind, Filter(Negate(is.null), edges))

  if (is.null(edge_mat) || nrow(edge_mat) == 0L) {
    message("No conserved ORF pairs found at the specified thresholds.")
    return(invisible(
      cbind(orfs[integer(0L), ],
            cluster_id   = integer(0L),
            n_sequences  = integer(0L),
            cluster_size = integer(0L))
    ))
  }

  ## ---- connected components ------------------------------------------------
  g   <- igraph::graph_from_edgelist(edge_mat, directed = FALSE)
  # Add isolated nodes for all ORFs so membership vector is aligned
  g   <- igraph::add_vertices(g, max(0L, n - igraph::vcount(g)))
  cmp <- igraph::components(g)

  membership <- cmp$membership   # length = number of vertices = n

  ## ---- filter to multi-species clusters ------------------------------------
  # Build a table: component id → set of seq_ids
  comp_seqids <- tapply(orfs$seq_id, membership, function(x) length(unique(x)))
  keep_comps  <- as.integer(names(comp_seqids)[comp_seqids >= 2L])

  keep_rows <- which(membership %in% keep_comps)

  if (length(keep_rows) == 0L) {
    message("No conserved ORF clusters found spanning ≥ 2 genomes.")
    return(invisible(
      cbind(orfs[integer(0L), ],
            cluster_id   = integer(0L),
            n_sequences  = integer(0L),
            cluster_size = integer(0L))
    ))
  }

  out <- orfs[keep_rows, ]
  raw_comp <- membership[keep_rows]

  ## ---- renumber cluster IDs ------------------------------------------------
  # Order clusters by earliest-starting ORF (by start position, then seq_id)
  order_key <- order(out$start, out$seq_id)
  comp_order <- raw_comp[order_key]
  unique_comps <- unique(comp_order)
  new_id <- match(raw_comp, unique_comps)

  out$cluster_id <- new_id

  ## ---- add summary columns -------------------------------------------------
  cluster_stats <- tapply(out$seq_id, out$cluster_id, function(x) {
    c(n_sequences = length(unique(x)), cluster_size = length(x))
  })

  stats_df <- do.call(rbind, cluster_stats)
  out$n_sequences  <- stats_df[as.character(out$cluster_id), "n_sequences"]
  out$cluster_size <- stats_df[as.character(out$cluster_id), "cluster_size"]

  ## ---- sort output ---------------------------------------------------------
  out <- out[order(out$cluster_id, out$seq_id, out$start), ]
  rownames(out) <- NULL

  if (verbose)
    message(sprintf("Found %d conserved cluster(s) across %d genome(s), %d ORF(s) total.",
                    length(unique(out$cluster_id)), length(seq_ids), nrow(out)))
  out
}
