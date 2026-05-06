# ============================================================================
# find_conserved_orfs — amino-acid-identity-based conservation clustering
# ============================================================================

#' Identify Conserved ORFs Across Multiple Sequences
#'
#' Groups ORFs from a [scan_orfs()] result into clusters based on pairwise
#' amino acid sequence identity, then retains clusters that are represented in
#' at least `min_sequences` distinct input sequences.  Amino acid comparison
#' is used rather than nucleotide comparison because it is more sensitive to
#' functional conservation across divergent mitogenomes.
#'
#' @param orfs `data.frame` returned by [scan_orfs()].  Must contain a
#'   `seq_id` column and a `protein_sequence` column.
#' @param min_identity Numeric scalar (0–100).  Minimum percent amino acid
#'   identity for two ORFs to be considered homologous.  Computed as matches /
#'   alignment length (equivalent to `PID3` in Biostrings).  Default `50`.
#' @param min_sequences Numeric scalar (0–1).  Minimum **fraction** of the
#'   input sequences in which a cluster must be represented.  The actual count
#'   threshold is `ceiling(min_sequences * n_sequences)`, so `0.5` requires
#'   presence in at least half of the supplied sequences (rounded up).
#'   Default `0.5`.  Must be in (0, 1].
#' @param alignment_type Character scalar.  Alignment type passed to
#'   [pwalign::pairwiseAlignment()]: `"local"` (Smith-Waterman, default)
#'   or `"global"` (Needleman-Wunsch).  Local alignment is recommended for
#'   ORF sets with variable lengths, as it correctly handles ORFs that are
#'   nested or partially overlapping in sequence space.
#'
#' @return A subset of `orfs` containing only ORFs that belong to a conserved
#'   cluster, with three additional columns appended:
#' \describe{
#'   \item{`cluster_id`}{Integer cluster identifier (arbitrary but consistent
#'     within a single call).}
#'   \item{`n_sequences`}{Number of distinct input sequences represented in
#'     the cluster.}
#'   \item{`cluster_size`}{Total number of ORFs in the cluster (across all
#'     sequences).}
#' }
#' Rows are sorted by `cluster_id`, then `seq_id`, then `start`.
#' Returns an empty `data.frame` invisibly (with a message) when no clusters
#' meet the specified thresholds.
#'
#' @details
#' ## Algorithm
#' Pairwise amino acid alignments are computed only between ORFs from
#' **different** sequences, avoiding self-comparisons.  Clusters are formed
#' using union-find (disjoint-set union): two ORFs are merged into the same
#' cluster whenever their pairwise identity meets `min_identity`.  This means
#' transitivity applies — if A ~ B and B ~ C, all three are clustered together
#' even if A and C do not directly meet the threshold.
#'
#' ## Computational cost
#' Runtime scales quadratically with the number of ORFs.  For large result
#' sets, narrow the search first with tight `min_orf_length` / `max_orf_length`
#' filters in [scan_orfs()].
#'
#' ## Identity metric
#' `PID3` (matches / alignment length, including gap columns) is used as it is
#' the most conservative metric for local alignments and avoids artificially
#' high identities from very short high-similarity regions.
#'
#' @seealso [scan_orfs()], [fetch_sequences()]
#'
#' @examples
#' \dontrun{
#' seqs <- fetch_sequences(c("NC_012920.1", "NC_001643.1", "NC_011137.1"))
#' orfs <- scan_orfs(seqs, min_orf_length = 60, max_orf_length = 900)
#'
#' # Present in >= 50% of input sequences (default)
#' conserved <- find_conserved_orfs(orfs, min_identity = 60)
#'
#' # Require presence in all input sequences
#' conserved_strict <- find_conserved_orfs(orfs, min_identity = 60,
#'                                          min_sequences = 1.0)
#'
#' # Inspect the clusters
#' split(conserved[, c("seq_id", "start", "end", "length_aa")],
#'       conserved$cluster_id)
#' }
#'
#' @importFrom Biostrings AAStringSet
#' @importFrom pwalign pairwiseAlignment pid
#' @export
find_conserved_orfs <- function(orfs,
                                 min_identity   = 50,
                                 min_sequences  = 0.5,
                                 alignment_type = "local") {

  ## ---- validate --------------------------------------------------------------
  if (!is.data.frame(orfs))
    stop("'orfs' must be a data.frame (output of scan_orfs()).", call. = FALSE)
  if (!"seq_id" %in% names(orfs))
    stop("'orfs' must contain a 'seq_id' column (output of scan_orfs()).",
         call. = FALSE)
  if (!"protein_sequence" %in% names(orfs))
    stop("'orfs' must contain a 'protein_sequence' column.", call. = FALSE)
  if (nrow(orfs) == 0L)
    stop("'orfs' is empty.", call. = FALSE)

  min_identity  <- as.numeric(min_identity)
  min_sequences <- as.numeric(min_sequences)

  if (is.na(min_identity) || min_identity < 0 || min_identity > 100)
    stop("'min_identity' must be a number between 0 and 100.", call. = FALSE)
  if (is.na(min_sequences) || min_sequences <= 0 || min_sequences > 1)
    stop("'min_sequences' must be a fraction in (0, 1].", call. = FALSE)
  if (!alignment_type %in% c("local", "global"))
    stop("'alignment_type' must be \"local\" or \"global\".", call. = FALSE)

  seq_ids_unique <- unique(orfs$seq_id)
  n_seqs         <- length(seq_ids_unique)
  if (n_seqs < 2L)
    stop("'orfs' must contain ORFs from at least two sequences; only '",
         seq_ids_unique[1L], "' is present.", call. = FALSE)

  # Resolve fraction to an absolute count (round up to be conservative)
  min_seq_count <- ceiling(min_sequences * n_seqs)
  message("Requiring clusters present in >= ", min_seq_count, " of ",
          n_seqs, " sequence(s) (min_sequences = ", min_sequences, ").")

  ## ---- prepare amino acid sequences ------------------------------------------
  n       <- nrow(orfs)
  aa_raw  <- gsub("\\*", "", orfs$protein_sequence)   # strip stop-codon asterisks
  aa_set  <- Biostrings::AAStringSet(aa_raw)
  seq_id_vec <- orfs$seq_id

  ## ---- union-find ------------------------------------------------------------
  parent <- seq_len(n)

  .find <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]   # path compression
      x <- parent[x]
    }
    x
  }

  .union <- function(x, y) {
    px <- .find(x)
    py <- .find(y)
    if (px != py) parent[px] <<- py
  }

  ## ---- pairwise alignments (cross-sequence pairs only) -----------------------
  message("Aligning ", n, " ORF(s) pairwise across ",
          n_seqs, " sequence(s)...")

  for (i in seq_len(n - 1L)) {
    # Restrict to ORFs that come from a DIFFERENT source sequence
    j_candidates <- which(seq_id_vec[(i + 1L):n] != seq_id_vec[i]) + i
    if (length(j_candidates) == 0L) next

    alns <- pwalign::pairwiseAlignment(
      pattern            = aa_set[j_candidates],
      subject            = aa_set[[i]],
      type               = alignment_type,
      substitutionMatrix = "BLOSUM62",
      gapOpening         = 10,
      gapExtension       = 0.5
    )

    identities <- pwalign::pid(alns, type = "PID3")   # matches / alignment length

    hits <- j_candidates[identities >= min_identity]
    for (j in hits) .union(i, j)
  }

  ## ---- extract clusters meeting min_sequences --------------------------------
  roots <- vapply(seq_len(n), .find, integer(1L))
  cluster_members <- split(seq_len(n), roots)

  keep <- Filter(
    function(idx) length(unique(seq_id_vec[idx])) >= min_seq_count,
    cluster_members
  )

  if (length(keep) == 0L) {
    message("No conserved ORF clusters found at the specified thresholds ",
            "(min_identity = ", min_identity, "%, ",
            "min_sequences = ", min_sequences,
            " [>= ", min_seq_count, " of ", n_seqs, " sequences]).")
    out <- orfs[0L, , drop = FALSE]
    out$cluster_id   <- integer(0L)
    out$n_sequences  <- integer(0L)
    out$cluster_size <- integer(0L)
    return(invisible(out))
  }

  ## ---- annotate and return ---------------------------------------------------
  # Stable sequential cluster IDs (sorted by smallest start position)
  first_start <- vapply(keep, function(idx) min(orfs$start[idx]), numeric(1L))
  keep <- keep[order(first_start)]

  cluster_label <- integer(n)
  n_seq_vec     <- integer(n)
  size_vec      <- integer(n)

  for (k in seq_along(keep)) {
    idx <- keep[[k]]
    cluster_label[idx] <- k
    n_seq_k            <- length(unique(seq_id_vec[idx]))
    n_seq_vec[idx]     <- n_seq_k
    size_vec[idx]      <- length(idx)
  }

  conserved_rows <- which(cluster_label > 0L)
  out <- orfs[conserved_rows, , drop = FALSE]
  out$cluster_id   <- cluster_label[conserved_rows]
  out$n_sequences  <- n_seq_vec[conserved_rows]
  out$cluster_size <- size_vec[conserved_rows]

  rownames(out) <- NULL
  out[order(out$cluster_id, out$seq_id, out$start), ]
}
