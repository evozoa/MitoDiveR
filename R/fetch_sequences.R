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
#'   (e.g. `c("NC_012920.1", "NC_028285")`).  Comma- or space-separated
#'   strings are split, and case and surrounding whitespace are ignored.
#'   Fetched via the Entrez API using [rentrez::entrez_fetch()].  `NULL` to
#'   skip.  To search by species and gene name instead, use [fetch_gene()].
#' @param fasta_file Character vector of paths to FASTA files on disk.
#'   `NULL` to skip.
#' @param fasta_text Character scalar or vector.  Raw FASTA content — either a
#'   single string with embedded newlines or a character vector of lines.
#'   `NULL` to skip.
#'
#' @return A named [Biostrings::DNAStringSet].  NCBI sequences are named by
#'   the accession as requested (so `fetch_sequences("NC_028285")` returns a
#'   sequence named `"NC_028285"`, not `"NC_028285.1"`) and returned in the
#'   requested order; FASTA sequences are named by the first
#'   whitespace-delimited token of the header.  These names serve as `seq_id`
#'   values in downstream [scan_orfs()] output.
#'
#' @details
#' Requires internet access when `accessions` is supplied.  Accessions are
#' requested in batches of 200; transient failures (rate limiting, NCBI server
#' errors, dropped connections) are retried with back-off.  Accessions NCBI
#' does not return are skipped with a warning, and an error is raised if none
#' are returned.  Inputs that are not accession numbers (e.g. species names)
#' are rejected before any request is made.
#'
#' NCBI rate-limits unauthenticated requests to three per second.  For large
#' jobs, set an NCBI API key to raise the limit to ten requests per second:
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
    accessions <- .parse_accession_arg(accessions)
    if (length(accessions) == 0L)
      stop("'accessions' contains no accession numbers.", call. = FALSE)

    message("Fetching ", length(accessions), " sequence(s) from NCBI...")
    seqs_list[["ncbi"]] <- .fetch_fasta_batch(accessions)
  }

  ## ---- FASTA file ------------------------------------------------------------
  if (!is.null(fasta_file)) {
    fasta_file <- as.character(fasta_file)
    for (f in fasta_file) {
      if (file.exists(f)) next
      if (.is_accession(toupper(trimws(f))))
        stop("'", f, "' was passed as 'fasta_file' but looks like an accession.\n",
             "To fetch several accessions, wrap them in c(): ",
             "fetch_sequences(c(\"",
             if (length(accessions)) accessions[[1L]] else "NC_012920.1",
             "\", \"", f, "\"))", call. = FALSE)
      stop("'fasta_file' not found: ", f,
           " (working directory: ", getwd(), ")", call. = FALSE)
    }
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


# --------------------------------------------------------------------------- #
# Internal: fetch FASTA for a vector of accessions in batches, name results   #
# by the accession as requested, and warn about any NCBI did not return.      #
# --------------------------------------------------------------------------- #
.fetch_fasta_batch <- function(accessions, batch_size = 200L) {

  fetch_one_batch <- function(ids) {
    tryCatch(
      .entrez_retry(.efetch, db = "nuccore", id = ids,
                    rettype = "fasta", retmode = "text"),
      error = function(e) {
        # NCBI rejects a batch outright when it holds only bad IDs; fall back
        # to one-by-one so the good ones still come through.
        if (length(ids) == 1L) return("")
        paste(vapply(ids, fetch_one_batch, character(1)), collapse = "\n")
      })
  }

  batches <- split(accessions, ceiling(seq_along(accessions) / batch_size))
  raw <- character(length(batches))
  for (b in seq_along(batches)) {
    raw[[b]] <- fetch_one_batch(batches[[b]])
    if (b < length(batches)) .entrez_pause()
  }
  raw <- paste(raw, collapse = "\n")

  # Anything before the first record is an NCBI error message, not sequence.
  first <- regexpr("(^|\n)>", raw)
  got <- if (first > 0L) {
    tmp <- tempfile(fileext = ".fasta")
    on.exit(unlink(tmp), add = TRUE)
    cat(substring(raw, first), file = tmp)
    Biostrings::readDNAStringSet(tmp)
  } else Biostrings::DNAStringSet()

  # Match requested accessions to returned ones, with or without version.
  ret_full <- sub("\\s.*$", "", names(got))
  ret_base <- sub("\\.\\d+$", "", ret_full)
  idx <- vapply(accessions, function(a) {
    hit <- if (grepl("\\.\\d+$", a)) which(ret_full == a) else which(ret_base == a)
    if (length(hit)) hit[[1L]] else NA_integer_
  }, integer(1))

  missing <- accessions[is.na(idx)]
  if (length(missing) == length(accessions))
    stop("NCBI returned no sequences for: ", paste(missing, collapse = ", "),
         ". Check that these are valid nucleotide accessions.", call. = FALSE)
  if (length(missing))
    warning(length(missing), " of ", length(accessions),
            " accession(s) not returned by NCBI and skipped: ",
            paste(missing, collapse = ", "), call. = FALSE)

  out <- got[idx[!is.na(idx)]]
  names(out) <- accessions[!is.na(idx)]
  out
}
