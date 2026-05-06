# ============================================================================
# calc_codon_position_rates — per-codon-position substitution rates
# ============================================================================

#' Calculate Per-Codon-Position Substitution Rates
#'
#' For an in-frame codon alignment, computes the proportion of nucleotide
#' differences at first, second, and third codon positions across all
#' pairwise sequence comparisons.  Tests the null hypothesis of equal
#' substitution rates across the three codon positions with a 3 x 2
#' chi-square contingency test on substitutions versus invariant sites.
#'
#' Canonical protein-coding genes show a strongly skewed signature with
#' the third position evolving fastest (most synonymous), the first
#' position intermediate, and the second position most constrained.
#' Departures from this pattern may indicate weakened protein-level
#' constraint, overlapping selection from the antisense strand (as in
#' rRNA-embedded MDP loci), or non-coding origin of the reading frame.
#'
#' @param alignment One of:
#'   - A `character` vector of equal-length, gap-free, in-frame DNA
#'     strings (length must be a multiple of 3).
#'   - A path to a FASTA file containing a codon-aligned set of sequences.
#'   - A [Biostrings::DNAStringSet] of equal-length sequences.
#'   - A `data.frame` from [cluster_conserved_orfs()] with an
#'     `orf_sequence` column; if a `cluster_id` column is present and
#'     more than one cluster is represented, `cluster_id` must also be
#'     supplied.
#' @param cluster_id Optional cluster identifier when `alignment` is a
#'   multi-cluster [cluster_conserved_orfs()] result.
#' @param return_pairs Logical.  If `TRUE`, the per-pair matrix of
#'   position-wise proportions is returned in `$pairs`.  Default `FALSE`.
#'
#' @return A list with:
#' \describe{
#'   \item{`by_position`}{`data.frame` with one row per codon position
#'     (1, 2, 3) and columns `n_sites`, `n_subs`, `prop_diff`, and
#'     `jc_dist` (Jukes-Cantor corrected; `NA` if `prop_diff` >= 0.75).}
#'   \item{`ratios`}{Named numeric vector of position-rate ratios
#'     (`p3_over_p1`, `p3_over_p2`, `p3_over_mean_p1p2`).}
#'   \item{`chisq`}{`htest` object from [chisq.test()] on the 3 x 2
#'     contingency table of substitutions versus invariant sites.}
#'   \item{`n_sequences`, `n_codons`, `n_pairs`}{Integer scalars
#'     describing the alignment.}
#'   \item{`pairs`}{(optional) `n_pairs` x 3 matrix of pairwise
#'     proportions of differences at each codon position; columns
#'     `p1`, `p2`, `p3`.  Returned only when `return_pairs = TRUE`.}
#' }
#'
#' @details
#' Gaps and ambiguous bases (anything outside `{A, C, G, T}`) are treated
#' as missing and excluded from both numerator and denominator at the
#' affected site within each pair.  The Jukes-Cantor correction is
#' \deqn{d = -\tfrac{3}{4}\ln\!\left(1 - \tfrac{4}{3} p\right)}
#' and returns `NA` for `p >= 0.75` (saturation).
#'
#' @examples
#' \dontrun{
#' # From a codon-aligned FASTA
#' res <- calc_codon_position_rates("shlp3_codon_aln.fas")
#' res$by_position
#' res$chisq
#'
#' # From a cluster_conserved_orfs() result
#' clusters <- cluster_conserved_orfs(orfs)
#' calc_codon_position_rates(clusters, cluster_id = "cl9")
#' }
#'
#' @seealso [calc_dnds()], [cluster_conserved_orfs()]
#' @importFrom Biostrings readDNAStringSet
#' @importFrom stats chisq.test
#' @export
calc_codon_position_rates <- function(alignment,
                                      cluster_id    = NULL,
                                      return_pairs  = FALSE) {

  ## ---- coerce input to a character matrix of single-base columns -------------
  seqs <- .cpr_resolve_input(alignment, cluster_id)

  if (length(seqs) < 2L)
    stop("At least two sequences are required.", call. = FALSE)

  widths <- nchar(seqs)
  if (length(unique(widths)) != 1L)
    stop("All sequences must be the same length (gap-free, in-frame).",
         call. = FALSE)

  n_nt     <- widths[[1L]]
  if (n_nt %% 3L != 0L)
    stop("Alignment length (", n_nt, " nt) is not a multiple of 3.",
         call. = FALSE)
  n_codons <- n_nt %/% 3L

  m <- toupper(do.call(rbind, strsplit(seqs, "", fixed = TRUE)))
  rownames(m) <- names(seqs)

  ## ---- per-pair substitution counts at each codon position -------------------
  cp      <- rep(1:3, times = n_codons)
  n_seq   <- nrow(m)
  n_pairs <- choose(n_seq, 2L)

  sub_counts  <- c(p1 = 0L, p2 = 0L, p3 = 0L)
  site_counts <- c(p1 = 0L, p2 = 0L, p3 = 0L)
  per_pair    <- matrix(NA_real_, nrow = n_pairs, ncol = 3L,
                        dimnames = list(NULL, c("p1", "p2", "p3")))

  idx <- 0L
  for (i in seq_len(n_seq - 1L)) {
    for (j in (i + 1L):n_seq) {
      idx     <- idx + 1L
      a       <- m[i, ]
      b       <- m[j, ]
      valid   <- a %in% c("A", "C", "G", "T") &
                 b %in% c("A", "C", "G", "T")
      diff_m  <- (a != b) & valid

      for (p in 1:3) {
        pos_mask <- (cp == p) & valid
        n_sites  <- sum(pos_mask)
        n_diff   <- sum(diff_m & cp == p)
        sub_counts[p]  <- sub_counts[p]  + n_diff
        site_counts[p] <- site_counts[p] + n_sites
        per_pair[idx, p] <- if (n_sites > 0L) n_diff / n_sites else NA_real_
      }
    }
  }

  ## ---- summarize -------------------------------------------------------------
  prop <- sub_counts / site_counts
  jc   <- vapply(prop, .jc_p_dist, numeric(1L))

  by_position <- data.frame(
    position  = 1:3,
    n_sites   = unname(site_counts),
    n_subs    = unname(sub_counts),
    prop_diff = unname(round(prop, 4L)),
    jc_dist   = unname(round(jc,   4L)),
    stringsAsFactors = FALSE
  )

  ratios <- c(
    p3_over_p1        = unname(prop[3L] / prop[1L]),
    p3_over_p2        = unname(prop[3L] / prop[2L]),
    p3_over_mean_p1p2 = unname(prop[3L] / mean(prop[1:2]))
  )

  ct <- rbind(observed  = sub_counts,
              unchanged = site_counts - sub_counts)
  chi <- suppressWarnings(stats::chisq.test(ct))

  out <- list(
    by_position = by_position,
    ratios      = round(ratios, 4L),
    chisq       = chi,
    n_sequences = n_seq,
    n_codons    = n_codons,
    n_pairs     = n_pairs
  )

  if (isTRUE(return_pairs)) out$pairs <- per_pair
  out
}


# --------------------------------------------------------------------------- #
# Internal helpers                                                             #
# --------------------------------------------------------------------------- #

# Resolve an alignment input into a named character vector of equal-length
# DNA strings (gaps stripped column-wise are NOT performed here; caller is
# responsible for supplying a gap-free in-frame alignment).
.cpr_resolve_input <- function(alignment, cluster_id) {

  if (is.character(alignment) && length(alignment) == 1L &&
      file.exists(alignment)) {
    dss <- Biostrings::readDNAStringSet(alignment)
    seqs <- as.character(dss)
    return(seqs)
  }

  if (is.character(alignment)) {
    if (is.null(names(alignment)))
      names(alignment) <- paste0("seq", seq_along(alignment))
    return(alignment)
  }

  if (inherits(alignment, "DNAStringSet")) {
    return(as.character(alignment))
  }

  if (is.data.frame(alignment)) {
    if (!"orf_sequence" %in% names(alignment))
      stop("data.frame input must contain an 'orf_sequence' column.",
           call. = FALSE)
    if ("cluster_id" %in% names(alignment) &&
        length(unique(alignment$cluster_id)) > 1L) {
      if (is.null(cluster_id))
        stop("Multiple clusters present; specify 'cluster_id'.",
             call. = FALSE)
      alignment <- alignment[alignment$cluster_id == cluster_id, ,
                             drop = FALSE]
      if (nrow(alignment) == 0L)
        stop("No rows match cluster_id = '", cluster_id, "'.",
             call. = FALSE)
    }
    seqs <- alignment$orf_sequence
    names(seqs) <- if ("seq_id" %in% names(alignment))
      alignment$seq_id else paste0("seq", seq_len(nrow(alignment)))
    return(seqs)
  }

  stop("Unsupported alignment input class: ",
       paste(class(alignment), collapse = "/"),
       call. = FALSE)
}

# Jukes-Cantor distance correction; returns NA for p >= 0.75 (saturation).
.jc_p_dist <- function(p) {
  if (is.na(p) || p >= 0.75) return(NA_real_)
  if (p == 0)                 return(0)
  -0.75 * log(1 - (4 / 3) * p)
}
