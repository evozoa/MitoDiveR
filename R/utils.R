# ============================================================================
# Internal utilities shared across MitoDiveR functions
# ============================================================================

# --------------------------------------------------------------------------- #
# .resolve_dna_input                                                           #
#                                                                              #
# Accepts any of:                                                              #
#   - NCBI accession string  (e.g. "NC_002333.2")                             #
#   - Path to a FASTA file   (single- or multi-sequence)                      #
#   - Raw nucleotide string                                                    #
#   - DNAString / DNAStringSet                                                 #
#                                                                              #
# Returns a named list of DNAStrings, one entry per sequence.                 #
# Names are: accession ID, FASTA header (first word), or "seq_1".            #
# --------------------------------------------------------------------------- #
.resolve_dna_input <- function(x, api_key = NULL) {

  ## --- already Bioconductor objects -----------------------------------------
  if (methods::is(x, "DNAStringSet")) {
    if (length(x) == 0L) stop("DNAStringSet is empty.", call. = FALSE)
    nms <- if (!is.null(names(x)) && any(nzchar(names(x))))
             sub("\\s.*", "", names(x))   # first word of FASTA header
           else paste0("seq_", seq_along(x))
    seqs <- as.list(x)
    names(seqs) <- nms
    return(seqs)
  }

  if (methods::is(x, "DNAString")) {
    return(list(seq_1 = x))
  }

  ## --- character input ------------------------------------------------------
  if (!is.character(x) || length(x) != 1L)
    stop("'sequence' must be a length-1 character string, DNAString, or DNAStringSet.",
         call. = FALSE)

  x_trim <- trimws(x)

  # (1) FASTA file path -------------------------------------------------------
  if (nchar(x_trim) < 500L && file.exists(x_trim)) {
    seqs <- Biostrings::readDNAStringSet(x_trim)
    if (length(seqs) == 0L)
      stop("No sequences found in file: ", x_trim, call. = FALSE)
    nms  <- sub("\\s.*", "", names(seqs))
    nms  <- ifelse(nzchar(nms), nms, paste0("seq_", seq_along(seqs)))
    seqs_list <- as.list(seqs)
    names(seqs_list) <- nms
    return(seqs_list)
  }

  # (2) NCBI accession --------------------------------------------------------
  # Pattern covers RefSeq (NC_, NM_, NR_, XM_, ...) and GenBank (2-letter
  # prefix + digits, e.g. AY123456) with optional version suffix (.1, .2, ...)
  acc_pattern <- "^[A-Z]{1,2}_?[A-Z0-9]*\\d{4,9}(\\.\\d+)?$"
  if (grepl(acc_pattern, x_trim)) {
    if (!requireNamespace("rentrez", quietly = TRUE))
      stop("Package 'rentrez' is required to fetch accessions. ",
           "Install with: install.packages('rentrez')", call. = FALSE)

    raw <- rentrez::entrez_fetch(
      db      = "nuccore",
      id      = x_trim,
      rettype = "fasta",
      retmode = "text",
      api_key = api_key
    )

    # Parse FASTA text directly (avoids temp file)
    lines    <- strsplit(raw, "\n", fixed = TRUE)[[1L]]
    hdr_idx  <- grep("^>", lines)
    if (length(hdr_idx) == 0L)
      stop("NCBI returned no FASTA sequence for accession: ", x_trim, call. = FALSE)

    seqs_list <- vector("list", length(hdr_idx))
    for (h in seq_along(hdr_idx)) {
      seq_start <- hdr_idx[h] + 1L
      seq_end   <- if (h < length(hdr_idx)) hdr_idx[h + 1L] - 1L else length(lines)
      seq_str   <- paste(lines[seq_start:seq_end], collapse = "")
      seq_str   <- toupper(gsub("[^ACGTNacgtn]", "", seq_str))
      nm        <- sub("^>([^\\s]+).*", "\\1", lines[hdr_idx[h]])
      seqs_list[[h]] <- Biostrings::DNAString(seq_str)
      names(seqs_list)[h] <- nm
    }
    return(seqs_list)
  }

  # (3) Raw nucleotide string -------------------------------------------------
  seq_clean <- toupper(gsub("[^ACGTNacgtn]", "", x_trim))
  if (nchar(seq_clean) < 3L)
    stop("Input does not appear to be a valid sequence, file path, or accession.",
         call. = FALSE)
  list(seq_1 = Biostrings::DNAString(seq_clean))
}
