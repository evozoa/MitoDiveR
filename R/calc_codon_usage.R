# ============================================================================
# calc_codon_usage — RSCU-based codon usage bias for ORF sets
# ============================================================================

#' Calculate Codon Usage Bias for ORF Sets
#'
#' Computes Relative Synonymous Codon Usage (RSCU) for the coding sequences
#' returned by [scan_orfs()] or [find_conserved_orfs()].  RSCU normalises
#' observed codon counts within each synonymous family so that uniform usage
#' gives 1.0; values above 1.0 indicate over-representation and below 1.0
#' indicate under-representation.  Stop codons are excluded from all
#' calculations.
#'
#' @param orfs `data.frame` from [scan_orfs()] or [find_conserved_orfs()].
#'   Must contain an `orf_sequence` column.
#' @param genetic_code Character scalar.  Biostrings genetic-code identifier
#'   passed to [Biostrings::getGeneticCode()].  Default `"SGC1"` (vertebrate
#'   mitochondrial code).  Must match the code used when scanning ORFs.
#' @param per_sequence Logical.  If `TRUE` (default) and a `seq_id` column is
#'   present, RSCU is computed separately for each source sequence.  If
#'   `FALSE`, all ORFs are pooled before counting.
#'
#' @return A `data.frame` with one row per codon (per sequence when
#'   `per_sequence = TRUE`), sorted by `seq_id` then `amino_acid` then
#'   `codon`:
#' \describe{
#'   \item{`seq_id`}{Source sequence identifier (only when
#'     `per_sequence = TRUE` and `seq_id` is present in input).}
#'   \item{`codon`}{Three-letter nucleotide codon.}
#'   \item{`amino_acid`}{Single-letter amino acid encoded by this codon under
#'     the selected genetic code.}
#'   \item{`n_synonymous_codons`}{Number of codons in the synonymous family
#'     (including this codon).}
#'   \item{`count`}{Observed count of this codon across all ORFs in the
#'     group.}
#'   \item{`rscu`}{Relative synonymous codon usage.  Computed as observed
#'     count divided by expected count under uniform synonymous usage:
#'     `RSCU = count / (family_total / n_synonymous_codons)`.  `NA` when
#'     the entire synonymous family has zero counts.}
#' }
#'
#' @details
#' ## RSCU interpretation
#' An RSCU of 1.0 means a codon is used exactly as often as expected if all
#' synonymous codons were equally preferred.  Values markedly different from
#' 1.0 indicate translational selection or mutational bias.  Because RSCU is
#' computed within synonymous families, it is not affected by amino-acid
#' composition.
#'
#' ## Genetic code
#' The genetic code affects which codons are synonymous.  In the vertebrate
#' mitochondrial code (`"SGC1"`), AGA and AGG are stop codons (not Arg), and
#' ATA encodes Met (not Ile), collapsing some synonymous families relative to
#' the standard code.  Always use the code that matches the analysed genome.
#'
#' @seealso [scan_orfs()], [find_conserved_orfs()], [calc_dnds()]
#'
#' @examples
#' \dontrun{
#' seqs <- fetch_sequences(c("NC_012920.1", "NC_011137.1"))
#' orfs <- scan_orfs(seqs, min_orf_length = 60)
#'
#' # Per-sequence RSCU (mito code)
#' rscu <- calc_codon_usage(orfs)
#'
#' # Pooled across all sequences, standard code
#' rscu_std <- calc_codon_usage(orfs, genetic_code = "SGC0",
#'                               per_sequence = FALSE)
#'
#' # Highly biased codons (RSCU > 1.5 or < 0.5)
#' rscu[!is.na(rscu$rscu) & (rscu$rscu > 1.5 | rscu$rscu < 0.5), ]
#' }
#'
#' @importFrom Biostrings getGeneticCode
#' @export
calc_codon_usage <- function(orfs,
                              genetic_code  = "SGC1",
                              per_sequence  = TRUE) {

  ## ---- validate ---------------------------------------------------------------
  if (!is.data.frame(orfs))
    stop("'orfs' must be a data.frame (output of scan_orfs()).", call. = FALSE)
  if (!"orf_sequence" %in% names(orfs))
    stop("'orfs' must contain an 'orf_sequence' column.", call. = FALSE)
  if (nrow(orfs) == 0L)
    stop("'orfs' is empty.", call. = FALSE)

  code     <- Biostrings::getGeneticCode(genetic_code)
  stop_cds <- names(code)[code == "*"]

  ## ---- synonymous-family metadata ---------------------------------------------
  non_stop   <- code[code != "*"]
  aa_groups  <- split(names(non_stop), as.character(non_stop))
  n_syn_lookup <- vapply(names(non_stop), function(cdn) {
    length(aa_groups[[as.character(non_stop[[cdn]])]])
  }, integer(1L))

  ## ---- RSCU computation for one group of ORFs ---------------------------------
  .rscu_for_group <- function(df) {
    all_codons <- unlist(lapply(df$orf_sequence, function(orf) {
      if (nchar(orf) < 3L || nchar(orf) %% 3L != 0L) return(character(0L))
      cdns <- .split_codons(orf)
      cdns[!cdns %in% stop_cds & cdns %in% names(non_stop)]
    }), use.names = FALSE)

    if (length(all_codons) == 0L) return(NULL)

    counts <- table(factor(all_codons, levels = names(non_stop)))

    rscu_vals <- vapply(names(non_stop), function(cdn) {
      aa        <- as.character(non_stop[[cdn]])
      fam       <- aa_groups[[aa]]
      fam_total <- sum(counts[fam])
      if (fam_total == 0L) return(NA_real_)
      n_i       <- length(fam)
      as.numeric(counts[[cdn]]) / (fam_total / n_i)
    }, numeric(1L))

    data.frame(
      codon               = names(non_stop),
      amino_acid          = as.character(non_stop),
      n_synonymous_codons = n_syn_lookup,
      count               = as.integer(counts),
      rscu                = rscu_vals,
      stringsAsFactors    = FALSE,
      row.names           = NULL
    )
  }

  ## ---- dispatch: per-sequence or pooled ---------------------------------------
  if (per_sequence && "seq_id" %in% names(orfs)) {
    parts <- lapply(unique(orfs$seq_id), function(sid) {
      r <- .rscu_for_group(orfs[orfs$seq_id == sid, , drop = FALSE])
      if (is.null(r)) return(NULL)
      cbind(data.frame(seq_id = sid, stringsAsFactors = FALSE), r)
    })
    out <- do.call(rbind, Filter(Negate(is.null), parts))
    if (is.null(out) || nrow(out) == 0L) {
      message("No valid codons found.")
      return(invisible(data.frame()))
    }
    out <- out[order(out$seq_id, out$amino_acid, out$codon), ]
  } else {
    out <- .rscu_for_group(orfs)
    if (is.null(out)) {
      message("No valid codons found.")
      return(invisible(data.frame()))
    }
    out <- out[order(out$amino_acid, out$codon), ]
  }

  rownames(out) <- NULL
  out
}


# --------------------------------------------------------------------------- #
# Internal helper                                                              #
# --------------------------------------------------------------------------- #

# Split a nucleotide string into a character vector of codons
.split_codons <- function(nt) {
  n      <- nchar(nt)
  starts <- seq(1L, n - 2L, by = 3L)
  substring(nt, starts, starts + 2L)
}
