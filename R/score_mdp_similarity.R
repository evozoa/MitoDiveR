# ============================================================================
# score_mdp_similarity — align ORF proteins against known MDPs
# ============================================================================

# Built-in MDP reference sequences.
#
# Every natural peptide below is the UniProt/Swiss-Prot sequence and has been
# verified to be encoded in the human mitochondrial reference (rCRS, NC_012920.1)
# and translated with the STANDARD genetic code (MDPs are made by cytoplasmic
# ribosomes).  Coordinates are rCRS positions of the coding sequence.
# tests/testthat/test-mdp_references.R re-derives each peptide from rCRS.
#
# Only canonical Swiss-Prot peptides are included; synthetic analogues (e.g. the
# S14G humanin variant HNG) are not, since they are not encoded in any genome.
.MDP_REFERENCE <- list(
  Humanin   = "MAPRGFSCLLLLTSEIDLPVKRRA",   # UniProt Q8IVG9; MT-RNR2 rCRS:2633-2704 (+)
  `MOTS-c`  = "MRWQEMGYIFYPRKLR",   # UniProt A0A0C5B5G6; MT-RNR1 rCRS:1343-1390 (+)
  SHLP1     = "MCHWAGGASNTGDARGDVFGKQAG",   # UniProt A0A3G1DJL7; MT-RNR2 rCRS:2488-2559 (-)
  SHLP2     = "MGVKFFTLSTRFFPSVQRAVPLWTNS",   # UniProt A0A3G1DIU6; MT-RNR2 rCRS:2091-2168 (-)
  SHLP3     = "MLGYNFSSFPCGTISIAPGFNFYRLYFIWVNGLAKVVW",   # UniProt A0A3G1DJQ2; MT-RNR2 rCRS:1706-1819 (-)
  SHLP4     = "MLEVMFLVNRRGKICRVPFTFFNLSL",   # UniProt A0A3G1DJK2; MT-RNR2 rCRS:2445-2522 (-)
  SHLP5     = "MYCSEVGFCSEVAPTEIFNAGLVV",   # UniProt A0A3G1DJL1; MT-RNR2 rCRS:2783-2854 (-)
  SHLP6     = "MLDQDIPMVQPLLKVRLFND"   # UniProt A0A3G1DJN1; MT-RNR2 rCRS:2990-3049 (+)
)

#' Score ORF Protein Sequences Against Known Mitochondrial-Derived Peptides
#'
#' Aligns protein sequences from a [find_conserved_orfs()] (or [scan_orfs()])
#' result against a curated reference set of known mitochondrial-derived
#' peptides (MDPs) using local Smith-Waterman alignment.  Returns the best
#' MDP match per query sequence, enabling rapid prioritisation of candidate
#' novel MDPs by similarity to known ones.
#'
#' @param orfs `data.frame` from [find_conserved_orfs()] or [scan_orfs()].
#'   Must contain a `protein_sequence` column.  If a `cluster_id` column is
#'   present, one representative per cluster (the longest sequence) is scored.
#' @param min_identity Numeric (0–100).  Minimum percent identity to report a
#'   hit.  Default `20` (permissive, to capture distant homologues).
#' @param custom_mdps Named character vector of additional reference sequences
#'   to include alongside the built-in MDP set.  Names are used as MDP
#'   identifiers.  `NULL` to use only the built-in references.
#'
#' @return A `data.frame` with one row per query–MDP pair that meets
#'   `min_identity`, with columns:
#' \describe{
#'   \item{`query_id`}{Identifier for the query ORF (cluster_id + seq_id, or
#'     row index).}
#'   \item{`cluster_id`}{Cluster identifier, if present in input.}
#'   \item{`query_length_aa`}{Length of the query protein (aa).}
#'   \item{`query_protein_sequence`}{Full amino acid sequence of the query ORF,
#'     translated using the genetic code selected in [scan_orfs()].}
#'   \item{`mdp_name`}{Name of the matching reference MDP.}
#'   \item{`mdp_length_aa`}{Length of the reference MDP (aa).}
#'   \item{`mdp_sequence`}{Amino acid sequence of the reference MDP.}
#'   \item{`pct_identity`}{Percent identity over the alignment
#'     (matches / alignment length).}
#'   \item{`alignment_score`}{Raw Smith-Waterman alignment score.}
#'   \item{`query_coverage`}{Fraction of query covered by the alignment (%).}
#'   \item{`mdp_coverage`}{Fraction of the reference MDP covered by the
#'     alignment (%).}
#' }
#' Sorted by `cluster_id` then `pct_identity` (descending).  Returns an empty
#' `data.frame` invisibly (with a message) if no matches meet `min_identity`.
#'
#' @details
#' ## Built-in MDP reference set
#' The following human MDPs are included: Humanin, MOTS-c and SHLP1–SHLP6.  Each
#' is the UniProt/Swiss-Prot sequence and has been verified to be encoded in the
#' human mitochondrial reference (rCRS) under the standard genetic code.  Synthetic
#' analogues (e.g. the S14G humanin variant HNG) are not included, since they are
#' not encoded in any genome.  Supply additional sequences via `custom_mdps` to extend the
#' reference, for example with species-specific variants or newly published
#' MDPs.
#'
#' ## No internet required
#' All alignments are performed locally using [pwalign::pairwiseAlignment()]
#' with BLOSUM62 and Smith-Waterman local alignment.
#'
#' @seealso [find_conserved_orfs()], [blast_orfs()]
#'
#' @examples
#' \dontrun{
#' seqs      <- fetch_sequences(c("NC_012920.1", "NC_001643.1"))
#' orfs      <- scan_orfs(seqs, min_orf_length = 60, max_orf_length = 900)
#' conserved <- find_conserved_orfs(orfs, min_identity = 50)
#' hits      <- score_mdp_similarity(conserved)
#' hits[hits$pct_identity >= 40, ]
#'
#' # Add a custom reference peptide
#' hits <- score_mdp_similarity(conserved,
#'           custom_mdps = c(MyPeptide = "MKLLFVIIAVLSACAGAVS"))
#' }
#'
#' @importFrom Biostrings AAStringSet AAString
#' @importFrom pwalign pairwiseAlignment pid
#' @export
score_mdp_similarity <- function(orfs,
                                  min_identity = 20,
                                  custom_mdps  = NULL) {

  ## ---- validate --------------------------------------------------------------
  if (!is.data.frame(orfs) || !"protein_sequence" %in% names(orfs))
    stop("'orfs' must be a data.frame with a 'protein_sequence' column.",
         call. = FALSE)
  if (nrow(orfs) == 0L)
    stop("'orfs' is empty.", call. = FALSE)

  min_identity <- as.numeric(min_identity)
  if (is.na(min_identity) || min_identity < 0 || min_identity > 100)
    stop("'min_identity' must be between 0 and 100.", call. = FALSE)

  has_clusters <- "cluster_id" %in% names(orfs)

  ## ---- build reference set ---------------------------------------------------
  ref_seqs <- .MDP_REFERENCE
  if (!is.null(custom_mdps)) {
    if (is.null(names(custom_mdps)) || any(nchar(names(custom_mdps)) == 0L))
      stop("'custom_mdps' must be a named character vector.", call. = FALSE)
    ref_seqs <- c(ref_seqs, as.list(custom_mdps))
  }

  ref_aa  <- Biostrings::AAStringSet(toupper(unlist(ref_seqs)))
  ref_names <- names(ref_seqs)

  ## ---- select representatives ------------------------------------------------
  if (has_clusters) {
    orfs$..aa_len <- nchar(gsub("\\*", "", orfs$protein_sequence))
    reps <- do.call(rbind, lapply(split(orfs, orfs$cluster_id), function(cl) {
      cl[which.max(cl$..aa_len), , drop = FALSE]
    }))
    reps$..aa_len <- NULL
    orfs$..aa_len <- NULL
  } else {
    reps <- orfs
  }

  seqs_clean <- gsub("\\*", "", reps$protein_sequence)
  if (has_clusters) {
    query_ids <- paste0("cluster_", reps$cluster_id, "_", reps$seq_id)
  } else {
    query_ids <- paste0("seq_", seq_len(nrow(reps)))
  }

  ## ---- align each query against all reference MDPs ---------------------------
  rows <- list()

  for (i in seq_len(nrow(reps))) {
    query_aa  <- Biostrings::AAString(seqs_clean[i])
    query_len <- nchar(seqs_clean[i])

    alns <- pwalign::pairwiseAlignment(
      pattern            = ref_aa,
      subject            = query_aa,
      type               = "local",
      substitutionMatrix = "BLOSUM62",
      gapOpening         = 10,
      gapExtension       = 0.5
    )

    identities <- pwalign::pid(alns, type = "PID3")
    scores     <- pwalign::score(alns)

    for (j in seq_along(ref_names)) {
      if (identities[j] < min_identity) next

      aln      <- alns[j]
      q_start  <- pwalign::start(pwalign::pattern(aln))
      q_end    <- pwalign::end(pwalign::pattern(aln))
      r_start  <- pwalign::start(pwalign::subject(aln))
      r_end    <- pwalign::end(pwalign::subject(aln))
      ref_len  <- nchar(ref_seqs[[j]])

      q_cov <- round((q_end - q_start + 1L) / query_len * 100, 1)
      r_cov <- round((r_end - r_start + 1L) / ref_len  * 100, 1)

      rows[[length(rows) + 1L]] <- data.frame(
        query_id              = query_ids[i],
        cluster_id            = if (has_clusters) reps$cluster_id[i] else NA_integer_,
        query_length_aa       = query_len,
        query_protein_sequence = seqs_clean[i],
        mdp_name              = ref_names[j],
        mdp_length_aa         = ref_len,
        mdp_sequence          = as.character(ref_seqs[[j]]),
        pct_identity          = round(identities[j], 1),
        alignment_score       = scores[j],
        query_coverage        = q_cov,
        mdp_coverage          = r_cov,
        stringsAsFactors      = FALSE
      )
    }
  }

  ## ---- assemble output -------------------------------------------------------
  if (length(rows) == 0L) {
    message("No MDP similarity hits found above min_identity = ", min_identity, "%.")
    return(invisible(data.frame(stringsAsFactors = FALSE)))
  }

  out <- do.call(rbind, rows)
  if (has_clusters) {
    out <- out[order(out$cluster_id, -out$pct_identity), ]
  } else {
    out <- out[order(out$query_id, -out$pct_identity), ]
  }
  rownames(out) <- NULL
  out
}
