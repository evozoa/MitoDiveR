# ============================================================================
# search_mdp_homologs — screen conserved ORFs against known MDPs
# ============================================================================

# --------------------------------------------------------------------------- #
# Built-in MDP reference sequences (human, canonical isoforms)                #
# --------------------------------------------------------------------------- #
.mdp_references <- list(
  "Humanin"   = "MAPRGFSCLLLLTSEIDLPVKRRA",
  "Humanin-G" = "MAPRGFSCLLLLTSEIDLPVKRRAG",
  "MOTS-c"    = "MRWQEMGYIFYPRKLR",
  "SHLP1"     = "MAWRQVGSRFFNLLQALQPIFNTA",
  "SHLP2"     = "MKLFGRGCSIDAFVLQEHLNLKPK",
  "SHLP3"     = "MLKKFKKHKLLEYLSRK",
  "SHLP4"     = "MKLLSSIEQACDICRLKK",
  "SHLP5"     = "MSTPAGSTDFITLGF",
  "SHLP6"     = "MLVFCPFPQTEQALAALQAIGGR"
)


#' Screen Conserved ORFs for Homology to Mitochondria-Derived Peptides
#'
#' Compares a representative protein sequence from each conserved ORF cluster
#' against a built-in library of known mitochondria-derived peptides (MDPs)
#' using Smith–Waterman local alignment.  Returns all pairwise comparisons
#' meeting the minimum identity threshold.
#'
#' @details
#' ## Built-in MDP library
#' The following human MDPs are included by default (canonical sequences):
#'
#' | Name | Length (aa) | Reference |
#' |------|-------------|-----------|
#' | Humanin | 24 | Hashimoto *et al.* 2001 |
#' | Humanin-G | 25 | Humanin G14 variant |
#' | MOTS-c | 16 | Lee *et al.* 2015 |
#' | SHLP1–6 | 15–24 | Cobb *et al.* 2016 |
#'
#' Supply `mdps` to add custom sequences or replace defaults.
#'
#' ## Cluster representative
#' One sequence per cluster is selected: the longest protein (stop-codon
#' asterisk stripped), with ties broken by `seq_id` alphabetically.  The
#' `query_id` column reports `cluster_{id}_{seq_id}` of the representative.
#'
#' ## Coverage metrics
#' - `query_coverage`: percentage of the *query ORF* covered by the
#'   alignment window (alignment width / query length × 100).
#' - `mdp_coverage`: percentage of the *MDP* covered by the alignment
#'   window (alignment width / MDP length × 100).  Values > 100 indicate
#'   that the alignment window (including gap columns) is wider than the MDP.
#'
#' @param conserved_orfs A `data.frame` — the output of
#'   [cluster_conserved_orfs()].  Must contain `cluster_id`, `seq_id`,
#'   `length_aa`, and `protein_sequence` columns.
#' @param min_pct_identity Numeric.  Minimum percent identity to report a hit.
#'   Default `20`.
#' @param min_coverage Numeric.  Minimum `query_coverage` (0–100) to report
#'   a hit.  Default `0` (no filter).
#' @param min_aligned_aa Integer.  Minimum number of aligned amino acid
#'   residues (non-gap positions in the MDP alignment) required to report a
#'   hit.  Prevents spurious 100 %% identity matches over trivially short
#'   windows.  Default `8L`.
#' @param mdps Named character vector or list of additional MDP sequences to
#'   include.  These are appended to (not replacing) the built-in library.
#'   Set to `character(0)` to use built-in MDPs only (default).
#' @param substitution_matrix Character.  Default `"BLOSUM62"`.
#' @param gap_opening Numeric.  Default `10`.
#' @param gap_extension Numeric.  Default `4`.
#'
#' @return A `data.frame` with one row per (cluster representative, MDP) pair
#'   that meets the thresholds, with columns:
#' \describe{
#'   \item{`query_id`}{`cluster_{cluster_id}_{seq_id}` of the representative.}
#'   \item{`cluster_id`}{Integer cluster identifier.}
#'   \item{`query_length_aa`}{Length of the representative protein (without stop codon).}
#'   \item{`query_protein_sequence`}{Representative protein sequence (without `*`).}
#'   \item{`mdp_name`}{MDP name.}
#'   \item{`mdp_length_aa`}{MDP length in amino acids.}
#'   \item{`mdp_sequence`}{MDP sequence.}
#'   \item{`pct_identity`}{Percent identity of the local alignment (PID1: matches / alignment length × 100), rounded to 1 decimal.}
#'   \item{`alignment_score`}{Raw Smith–Waterman alignment score, rounded to 1 decimal.}
#'   \item{`query_coverage`}{Percentage of query ORF covered by the alignment, rounded to 1 decimal.}
#'   \item{`mdp_coverage`}{Percentage of MDP covered by the alignment, rounded to 1 decimal.}
#' }
#'
#' Returns an empty `data.frame` (invisibly) when no hits are found.
#'
#' @examples
#' \dontrun{
#' conserved <- cluster_conserved_orfs(orfs_all)
#' hits <- search_mdp_homologs(conserved)
#'
#' # Best hit per cluster
#' hits[order(hits$cluster_id, -hits$pct_identity), ]
#'
#' # Add a custom peptide
#' hits2 <- search_mdp_homologs(conserved,
#'                               mdps = c(MyPeptide = "MLLSALVGAAVAFVLA"))
#' }
#'
#' @importFrom pwalign pairwiseAlignment pid score
#' @importFrom Biostrings AAString
#' @export
search_mdp_homologs <- function(conserved_orfs,
                                 min_pct_identity    = 20,
                                 min_coverage        = 0,
                                 min_aligned_aa      = 8L,
                                 mdps                = character(0L),
                                 substitution_matrix = "BLOSUM62",
                                 gap_opening         = 10,
                                 gap_extension       = 4) {

  ## ---- input checks --------------------------------------------------------
  needed <- c("cluster_id", "seq_id", "length_aa", "protein_sequence")
  miss   <- setdiff(needed, names(conserved_orfs))
  if (length(miss))
    stop("'conserved_orfs' is missing required columns: ",
         paste(miss, collapse = ", "), call. = FALSE)

  if (nrow(conserved_orfs) == 0L) {
    message("No conserved ORFs supplied.")
    return(invisible(.empty_mdp_hits()))
  }

  ## ---- merge MDP library ---------------------------------------------------
  mdp_lib <- .mdp_references
  if (length(mdps) > 0L) {
    mdps <- as.list(mdps)
    mdp_lib <- c(mdp_lib, mdps[!names(mdps) %in% names(mdp_lib)])
  }

  ## ---- prepare all unique ORFs for alignment --------------------------------
  conserved_orfs$.prot_clean <- gsub("[*]$", "", conserved_orfs$protein_sequence)
  conserved_orfs$.prot_len   <- nchar(conserved_orfs$.prot_clean)

  # De-duplicate by protein sequence — same sequence in multiple rows adds no info
  orfs_unique <- conserved_orfs[!duplicated(conserved_orfs$.prot_clean), ]

  ## ---- pairwise alignment: every unique ORF vs each MDP -------------------
  result_rows <- vector("list", nrow(orfs_unique) * length(mdp_lib))
  k <- 0L

  for (ri in seq_len(nrow(orfs_unique))) {
    rep_row  <- orfs_unique[ri, ]
    prot_seq <- rep_row$.prot_clean
    prot_len <- rep_row$.prot_len
    if (prot_len < 2L) next

    aa_query <- Biostrings::AAString(prot_seq)

    for (mdp_name in names(mdp_lib)) {
      mdp_seq <- as.character(mdp_lib[[mdp_name]])
      mdp_len <- nchar(mdp_seq)

      aln <- pwalign::pairwiseAlignment(
        pattern            = Biostrings::AAString(mdp_seq),
        subject            = aa_query,
        type               = "local",
        substitutionMatrix = substitution_matrix,
        gapOpening         = gap_opening,
        gapExtension       = gap_extension
      )

      pct_id   <- round(pwalign::pid(aln, type = "PID1"), 1L)
      aln_scr  <- round(pwalign::score(aln), 1L)
      aln_wid  <- nchar(as.character(pwalign::alignedPattern(aln)))
      q_cov    <- round(aln_wid / prot_len * 100, 1L)
      mdp_cov  <- round(aln_wid / mdp_len  * 100, 1L)

      # Number of non-gap residues aligned in the MDP
      n_mdp_aligned <- nchar(gsub("-", "", as.character(pwalign::alignedPattern(aln))))

      if (pct_id < min_pct_identity || q_cov < min_coverage ||
          n_mdp_aligned < min_aligned_aa) next

      k <- k + 1L
      result_rows[[k]] <- data.frame(
        query_id               = paste0("cluster_", rep_row$cluster_id,
                                        "_", rep_row$seq_id,
                                        "_", rep_row$start),
        cluster_id             = rep_row$cluster_id,
        seq_id                 = rep_row$seq_id,
        orf_start              = rep_row$start,
        orf_end                = rep_row$end,
        strand                 = rep_row$strand,
        genomic_region         = if ("genomic_region" %in% names(rep_row))
                                   rep_row$genomic_region else NA_character_,
        query_length_aa        = prot_len,
        query_protein_sequence = prot_seq,
        mdp_name               = mdp_name,
        mdp_length_aa          = mdp_len,
        mdp_sequence           = mdp_seq,
        pct_identity           = pct_id,
        alignment_score        = aln_scr,
        query_coverage         = q_cov,
        mdp_coverage           = mdp_cov,
        stringsAsFactors       = FALSE
      )
    }
  }

  ## ---- assemble output -----------------------------------------------------
  if (k == 0L) {
    message("No MDP homologs found at the specified thresholds.")
    return(invisible(.empty_mdp_hits()))
  }

  out <- do.call(rbind, result_rows[seq_len(k)])
  out <- out[order(out$cluster_id, out$orf_start, -out$pct_identity), ]
  rownames(out) <- NULL
  out
}


# --------------------------------------------------------------------------- #
# Internal: empty MDP hits data.frame                                         #
# --------------------------------------------------------------------------- #
.empty_mdp_hits <- function() {
  data.frame(
    query_id               = character(),
    cluster_id             = integer(),
    query_length_aa        = integer(),
    query_protein_sequence = character(),
    mdp_name               = character(),
    mdp_length_aa          = integer(),
    mdp_sequence           = character(),
    pct_identity           = numeric(),
    alignment_score        = numeric(),
    query_coverage         = numeric(),
    mdp_coverage           = numeric(),
    stringsAsFactors       = FALSE
  )
}


#' List Built-in MDP Reference Sequences
#'
#' Returns the named character vector of human MDP sequences that
#' [search_mdp_homologs()] uses by default.
#'
#' @return A named character vector; names are MDP names, values are amino
#'   acid sequences (single-letter code, no stop codon).
#'
#' @examples
#' mdp_sequences()
#'
#' @export
mdp_sequences <- function() {
  vapply(.mdp_references, identity, character(1L))
}
