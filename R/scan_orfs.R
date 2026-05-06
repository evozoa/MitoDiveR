# ============================================================================
# scan_orfs — ORF scanning across multiple sequences
# ============================================================================

#' Scan Multiple Sequences for Open Reading Frames
#'
#' Applies [find_orfs()] to every sequence in a [Biostrings::DNAStringSet] and
#' returns a single combined `data.frame` with a leading `seq_id` column.
#' Designed to accept output directly from [fetch_sequences()].
#'
#' @param sequences A named [Biostrings::DNAStringSet] or named character
#'   vector of nucleotide sequences.  Typically the output of
#'   [fetch_sequences()].
#' @param min_orf_length Integer. Minimum ORF length in nucleotides (start
#'   codon to stop codon, inclusive).  Default `100L`, matching the NCBI ORF
#'   Finder default.
#' @param max_orf_length Integer or `Inf`. Maximum ORF length in nucleotides.
#'   Default `Inf` (no upper limit).  Set a value to focus on the small ORFs
#'   characteristic of mitochondrial-derived peptides (e.g. `max_orf_length =
#'   900L` for peptides up to 300 aa).
#' @param genetic_code Character scalar. Biostrings genetic-code identifier
#'   passed to [find_orfs()].  Default `"SGC1"` (vertebrate mitochondrial
#'   code).
#' @param start_codons Character vector of recognised start codons passed to
#'   [find_orfs()].  Default `c("ATG", "GTG", "ATT", "ATC", "ATA")`.
#' @param both_strands Logical.  Search the reverse-complement strand?
#'   Default `TRUE`.
#' @param circular Logical.  Treat each sequence as a circular molecule?
#'   Default `TRUE`.
#' @param include_stop_codon Logical.  Include the stop codon in reported
#'   sequences and lengths?  Default `TRUE`.
#'
#' @return A `data.frame` with the same columns as [find_orfs()] plus a
#'   leading `seq_id` column identifying the source sequence.  Rows are
#'   ordered by `seq_id` then genomic `start` position.  Returns an empty
#'   `data.frame` invisibly (with a message) if no ORFs are found.
#'
#' @seealso [find_orfs()], [fetch_sequences()]
#'
#' @examples
#' \dontrun{
#' # Fetch two human and chimp mitogenomes, then scan for small ORFs
#' seqs <- fetch_sequences(c("NC_012920.1", "NC_001643.1"))
#' orfs <- scan_orfs(seqs, min_orf_length = 60, max_orf_length = 900)
#'
#' # Summarise hits per sequence and strand
#' table(orfs$seq_id, orfs$strand)
#'
#' # Restrict to ATG-only, standard code
#' orfs_std <- scan_orfs(seqs,
#'                        genetic_code = "SGC0",
#'                        start_codons = "ATG")
#' }
#'
#' @importFrom Biostrings DNAStringSet
#' @importFrom methods is
#' @export
scan_orfs <- function(sequences,
                      min_orf_length     = 100L,
                      max_orf_length     = Inf,
                      genetic_code       = "SGC1",
                      start_codons       = c("ATG", "GTG", "ATT", "ATC", "ATA"),
                      both_strands       = TRUE,
                      circular           = TRUE,
                      include_stop_codon = TRUE) {

  ## ---- coerce & validate input -----------------------------------------------
  if (is.character(sequences)) {
    sequences <- Biostrings::DNAStringSet(toupper(sequences))
  }
  if (!methods::is(sequences, "DNAStringSet"))
    stop("'sequences' must be a DNAStringSet or named character vector.",
         call. = FALSE)
  if (length(sequences) == 0L)
    stop("'sequences' is empty.", call. = FALSE)

  seq_ids <- names(sequences)
  if (is.null(seq_ids) || any(nchar(seq_ids) == 0L))
    seq_ids <- paste0("seq_", seq_len(length(sequences)))

  ## ---- scan each sequence ----------------------------------------------------
  n <- length(sequences)
  results <- vector("list", n)

  for (i in seq_len(n)) {
    orfs <- find_orfs(
      sequence           = sequences[[i]],
      genetic_code       = genetic_code,
      start_codons       = start_codons,
      min_orf_length     = min_orf_length,
      max_orf_length     = max_orf_length,
      both_strands       = both_strands,
      circular           = circular,
      include_stop_codon = include_stop_codon
    )

    if (nrow(orfs) == 0L) next

    results[[i]] <- cbind(data.frame(seq_id = seq_ids[i], stringsAsFactors = FALSE), orfs)
  }

  ## ---- assemble output -------------------------------------------------------
  out <- do.call(rbind, Filter(Negate(is.null), results))

  if (is.null(out) || nrow(out) == 0L) {
    message("No ORFs found matching the specified criteria.")
    return(invisible(
      data.frame(seq_id = character(), stringsAsFactors = FALSE)
    ))
  }

  rownames(out) <- NULL
  out
}
