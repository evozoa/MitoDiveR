# ============================================================================
# scan_orfeome — multi-parameter ORF scanning
# ============================================================================

#' Scan Sequences for ORFs Across Multiple Parameter Combinations
#'
#' Wraps [scan_orfs()] to enumerate ORFs under combinations of genetic
#' codes, minimum ORF lengths, and start-codon sets.  Useful for
#' surveying the full ORFeome of a mitogenome (or partial genome)
#' under different translation and detection assumptions.
#'
#' @param sequences NCBI accession(s), path to a FASTA file, raw
#'   sequence string, or [Biostrings::DNAStringSet].  Accessions are
#'   passed to [fetch_sequences()].
#' @param genetic_codes Character vector of Biostrings genetic-code
#'   identifiers.  Required (no default) — e.g. `c("SGC1", "SGC0")` to test
#'   both the vertebrate mitochondrial and standard codes.
#' @param min_orf_lengths Integer vector of minimum ORF lengths in
#'   nucleotides.  Default `60L` (single value).  Provide multiple
#'   values (e.g. `c(60L, 90L, 150L)`) to scan at several thresholds.
#' @param start_codon_sets List of character vectors, each a set of
#'   recognised start codons.  Default `list(canonical = "ATG",
#'   permissive = c("ATG","GTG","ATA","ATT","ATC"))`.
#' @param both_strands Logical. Search the reverse-complement strand?
#'   Default `TRUE`.
#' @param circular Logical. Treat sequences as circular molecules?
#'   Default `TRUE`.
#' @param include_stop_codon Logical. Include the stop codon in
#'   reported sequences and lengths?  Default `TRUE`.
#'
#' @return A `data.frame` combining all [scan_orfs()] outputs with
#'   additional columns identifying each parameter combination:
#'   `genetic_code` (the code identifier), `min_orf_length`,
#'   `start_codon_set` (name of the start-codon-set list element).
#'   The same ORF detected under multiple parameter combinations
#'   appears once per combination.
#'
#' @details
#' This function is intended for exploratory ORFeome characterization
#' — for example, "show me every minus-strand ORF $\geq$ 60 nt in
#' MT-RNR2 under both mitochondrial and standard codes."  When the
#' goal is a single canonical ORF set, prefer [scan_orfs()] directly.
#'
#' @seealso [scan_orfs()], [fetch_sequences()]
#'
#' @examples
#' \dontrun{
#' # Pull rCRS and survey its ORFeome under multiple settings
#' orfs <- scan_orfeome(
#'   "NC_012920.1",
#'   genetic_codes    = c("SGC1", "SGC0"),
#'   min_orf_lengths  = c(60L, 90L)
#' )
#' table(orfs$genetic_code, orfs$min_orf_length)
#' }
#'
#' @importFrom Biostrings DNAStringSet readDNAStringSet
#' @importFrom methods is
#' @export
scan_orfeome <- function(sequences,
                         genetic_codes      = .stop_no_code("genetic_codes"),
                         min_orf_lengths    = 60L,
                         start_codon_sets   = list(
                           canonical  = "ATG",
                           permissive = c("ATG","GTG","ATA","ATT","ATC")
                         ),
                         both_strands       = TRUE,
                         circular           = TRUE,
                         include_stop_codon = TRUE) {

  ## ---- coerce input to DNAStringSet ----------------------------------------
  seqs <- .coerce_to_dnastringset(sequences)
  if (length(seqs) == 0L) stop("'sequences' resolves to zero sequences.",
                                call. = FALSE)

  ## ---- validate parameter vectors ------------------------------------------
  if (length(genetic_codes) == 0L)
    stop("'genetic_codes' must have at least one value.", call. = FALSE)
  min_orf_lengths <- as.integer(min_orf_lengths)
  if (any(min_orf_lengths < 3L))
    stop("All 'min_orf_lengths' must be >= 3.", call. = FALSE)
  if (!is.list(start_codon_sets) || length(start_codon_sets) == 0L)
    stop("'start_codon_sets' must be a non-empty named list.", call. = FALSE)
  if (is.null(names(start_codon_sets)))
    names(start_codon_sets) <- paste0("set_", seq_along(start_codon_sets))

  ## ---- iterate over parameter combinations ---------------------------------
  out_list <- list()
  for (code in genetic_codes) {
    for (mol in min_orf_lengths) {
      for (sname in names(start_codon_sets)) {
        sstart <- start_codon_sets[[sname]]
        message(sprintf("scan_orfs: code=%s, min_orf_length=%d, start_set=%s",
                        code, mol, sname))
        orfs <- scan_orfs(seqs,
                          min_orf_length     = mol,
                          genetic_code       = code,
                          start_codons       = sstart,
                          both_strands       = both_strands,
                          circular           = circular,
                          include_stop_codon = include_stop_codon)
        if (nrow(orfs) > 0L) {
          orfs$genetic_code    <- code
          orfs$min_orf_length  <- mol
          orfs$start_codon_set <- sname
          out_list[[paste(code, mol, sname, sep = "_")]] <- orfs
        }
      }
    }
  }

  if (length(out_list) == 0L) {
    message("No ORFs found under any parameter combination.")
    return(invisible(data.frame()))
  }
  do.call(rbind, out_list)
}

# ---- helper: coerce sequences input to DNAStringSet ------------------------
.coerce_to_dnastringset <- function(x) {
  if (methods::is(x, "DNAStringSet")) return(x)
  if (methods::is(x, "DNAString")) {
    ss <- Biostrings::DNAStringSet(x); names(ss) <- "seq_1"
    return(ss)
  }
  if (is.character(x)) {
    # Treat as FASTA file path if file exists, else as accessions
    if (length(x) == 1L && file.exists(x)) {
      return(Biostrings::readDNAStringSet(x))
    }
    # Treat as raw sequence if predominantly ACGTN
    if (length(x) == 1L && nchar(x) > 100L &&
        nchar(gsub("[^ACGTNacgtn]", "", x)) >= 0.9 * nchar(x)) {
      ss <- Biostrings::DNAStringSet(toupper(x)); names(ss) <- "seq_1"
      return(ss)
    }
    # Otherwise treat as NCBI accessions
    return(fetch_sequences(accessions = x))
  }
  stop("'sequences' must be a DNAStringSet, accession(s), file path, ",
       "or raw sequence string.", call. = FALSE)
}
