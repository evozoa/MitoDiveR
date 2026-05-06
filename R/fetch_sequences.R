# ============================================================================
# fetch_sequences — retrieve nucleotide sequences from NCBI or FASTA input
# ============================================================================

#' Fetch Nucleotide Sequences from NCBI or FASTA Input
#'
#' Retrieves nucleotide sequences by NCBI accession number and/or reads them
#' from a FASTA file or in-memory FASTA text.  All sources are combined into a
#' single named [Biostrings::DNAStringSet] suitable for direct input to
#' [scan_orfs()].
#'
#' @param accessions Character vector of NCBI nucleotide accession numbers
#'   (e.g. `"NC_012920.1"`).  Fetched via the Entrez API using
#'   [rentrez::entrez_fetch()].  `NULL` to skip.
#' @param fasta_file Character scalar.  Path to a FASTA file on disk.
#'   `NULL` to skip.
#' @param fasta_text Character scalar or vector.  Raw FASTA content — either a
#'   single string with embedded newlines or a character vector of lines.
#'   `NULL` to skip.
#'
#' @return A named [Biostrings::DNAStringSet].  Sequence names are set to the
#'   accession (first whitespace-delimited token of the FASTA header), so they
#'   remain stable regardless of the retrieval route and serve as `seq_id`
#'   values in downstream [scan_orfs()] output.
#'
#' @details
#' Requires internet access when `accessions` is supplied.  NCBI rate-limits
#' unauthenticated requests to three per second; [rentrez::entrez_fetch()]
#' handles retry automatically.  For batches of more than ~200 accessions,
#' set an NCBI API key to raise the limit to ten requests per second:
#' ```r
#' Sys.setenv(ENTREZ_KEY = "your_key_here")
#' ```
#'
#' @examples
#' \dontrun{
#' # Human mitochondrial reference genome
#' seqs <- fetch_sequences("NC_012920.1")
#'
#' # Multiple accessions at once
#' seqs <- fetch_sequences(c("NC_012920.1", "NC_011137.1"))
#'
#' # From a local FASTA file
#' seqs <- fetch_sequences(fasta_file = "my_mito.fasta")
#'
#' # Combine accessions and a local file
#' seqs <- fetch_sequences(accessions = "NC_012920.1",
#'                         fasta_file  = "additional.fasta")
#' }
#'
#' @importFrom Biostrings readDNAStringSet
#' @importFrom rentrez entrez_fetch
#' @export
fetch_sequences <- function(accessions = NULL,
                             fasta_file  = NULL,
                             fasta_text  = NULL) {

  if (is.null(accessions) && is.null(fasta_file) && is.null(fasta_text))
    stop("Supply at least one of: 'accessions', 'fasta_file', or 'fasta_text'.",
         call. = FALSE)

  seqs_list <- list()

  ## ---- NCBI accessions -------------------------------------------------------
  if (!is.null(accessions)) {
    accessions <- as.character(accessions)
    message("Fetching ", length(accessions), " sequence(s) from NCBI...")
    raw_fasta <- rentrez::entrez_fetch(
      db      = "nucleotide",
      id      = accessions,
      rettype = "fasta",
      retmode = "text"
    )
    tmp_ncbi <- tempfile(fileext = ".fasta")
    on.exit(unlink(tmp_ncbi), add = TRUE)
    cat(raw_fasta, file = tmp_ncbi)
    seqs_list[["ncbi"]] <- Biostrings::readDNAStringSet(tmp_ncbi)
  }

  ## ---- FASTA file ------------------------------------------------------------
  if (!is.null(fasta_file)) {
    if (!file.exists(fasta_file))
      stop("'fasta_file' not found: ", fasta_file, call. = FALSE)
    seqs_list[["file"]] <- Biostrings::readDNAStringSet(fasta_file)
  }

  ## ---- Raw FASTA text --------------------------------------------------------
  if (!is.null(fasta_text)) {
    tmp_txt <- tempfile(fileext = ".fasta")
    on.exit(unlink(tmp_txt), add = TRUE)
    cat(fasta_text, file = tmp_txt, sep = "\n")
    seqs_list[["text"]] <- Biostrings::readDNAStringSet(tmp_txt)
  }

  ## ---- Combine & normalise names ---------------------------------------------
  out <- do.call(c, unname(seqs_list))

  # Trim to first whitespace-delimited token (the accession / sequence ID)
  names(out) <- sub("\\s.*$", "", names(out))

  out
}
