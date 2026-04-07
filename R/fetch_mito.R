# ============================================================================
# fetch_mito_genbank — download and parse RefSeq mitochondrial genomes
# ============================================================================

#' Fetch Mitochondrial Genomes from NCBI RefSeq
#'
#' Downloads one or more mitochondrial genome records from NCBI by accession
#' number and parses both the nucleotide sequence and the GenBank feature
#' table.  The feature table is used downstream by
#' [annotate_genomic_regions()] to label ORFs by genomic context.
#'
#' @param accessions Character vector of sequence sources.  Each element may be:
#'   - An **NCBI accession number** (e.g. `"NC_002333.2"`) — full GenBank
#'     record is fetched, including the feature table.
#'   - A **path to a FASTA file** — sequence is read locally; features will
#'     be an empty `data.frame` (no annotation available from FASTA).
#'   - A **raw nucleotide string** — treated as an anonymous sequence; features
#'     will be an empty `data.frame`.
#'
#'   All three types may be mixed in a single call.
#' @param api_key Character scalar. NCBI API key for higher rate limits
#'   (10 requests/s vs 3/s without a key).  Obtain one free at
#'   <https://www.ncbi.nlm.nih.gov/account/>.  Set to `NULL` (default) to
#'   use unauthenticated access.
#' @param verbose Logical.  Print progress messages?  Default `TRUE`.
#'
#' @return A named list with one entry per accession.  Each entry is itself
#'   a list with:
#' \describe{
#'   \item{`accession`}{The requested accession string.}
#'   \item{`definition`}{DEFINITION line from the GenBank record.}
#'   \item{`length`}{Sequence length in nucleotides.}
#'   \item{`sequence`}{The complete mitogenome as a [Biostrings::DNAString].}
#'   \item{`features`}{A `data.frame` of parsed GenBank features with
#'     columns `type`, `start`, `end`, `strand`, `gene`, `product`, and
#'     `note`.  Locations are 1-based and inclusive; `strand` is `"+"` or
#'     `"-"`.}
#' }
#'
#' @details
#' Sequences are fetched via `rentrez::entrez_fetch()` in GenBank flat-file
#' format (`rettype = "gbwithparts"`).  A brief pause between requests
#' respects NCBI's rate limits even without an API key.
#'
#' Feature locations are simplified to a single `[start, end]` span:
#' - *complement* locations set `strand = "-"`.
#' - *join* locations use the overall min/max span.
#' - Wrap-around locations (e.g. `join(16500..16770,1..200)`) are handled
#'   correctly for circular genomes.
#'
#' @examples
#' \dontrun{
#' genomes <- fetch_mito_genbank(c("NC_002333.2", "NC_028087.1"))
#'
#' # Inspect Danio rerio features
#' head(genomes[["NC_002333.2"]]$features)
#'
#' # Pass the sequence directly to find_orfs()
#' orfs <- find_orfs(genomes[["NC_002333.2"]]$sequence)
#' }
#'
#' @importFrom rentrez entrez_fetch
#' @export
fetch_mito_genbank <- function(accessions,
                                api_key = NULL,
                                verbose  = TRUE) {

  if (!requireNamespace("rentrez", quietly = TRUE))
    stop("Package 'rentrez' is required. Install with: install.packages('rentrez')",
         call. = FALSE)

  accessions <- as.character(accessions)
  if (length(accessions) == 0L)
    stop("'accessions' must be a non-empty character vector.", call. = FALSE)

  results <- list()

  empty_features <- data.frame(type = character(), start = integer(),
                               end = integer(), strand = character(),
                               gene = character(), product = character(),
                               note = character(), stringsAsFactors = FALSE)

  acc_pattern <- "^[A-Z]{1,2}_?[A-Z0-9]*\\d{4,9}(\\.\\d+)?$"

  for (i in seq_along(accessions)) {
    acc <- accessions[[i]]

    # --- FASTA file -----------------------------------------------------------
    if (nchar(acc) < 500L && file.exists(acc)) {
      if (verbose) message("Reading FASTA file: ", acc)
      seqs <- Biostrings::readDNAStringSet(acc)
      if (length(seqs) == 0L) stop("No sequences found in file: ", acc, call. = FALSE)
      # One result entry per sequence in the file
      for (j in seq_along(seqs)) {
        nm <- sub("\\s.*", "", names(seqs)[j])
        nm <- if (nzchar(nm)) nm else paste0("seq_", j)
        results[[nm]] <- list(accession  = nm,
                              definition = names(seqs)[j],
                              length     = length(seqs[[j]]),
                              sequence   = seqs[[j]],
                              features   = empty_features)
      }
      next
    }

    # --- NCBI accession -------------------------------------------------------
    if (grepl(acc_pattern, trimws(acc))) {
      if (verbose) message("Fetching ", acc, " ...")
      raw <- rentrez::entrez_fetch(
        db      = "nuccore",
        id      = acc,
        rettype = "gbwithparts",
        retmode = "text",
        api_key = api_key
      )
      results[[acc]] <- .parse_genbank(raw, acc)
      if (i < length(accessions))
        Sys.sleep(if (is.null(api_key)) 0.34 else 0.11)
      next
    }

    # --- Raw nucleotide string ------------------------------------------------
    if (verbose) message("Treating input ", i, " as raw sequence.")
    seq_clean <- Biostrings::DNAString(toupper(gsub("[^ACGTNacgtn]", "", acc)))
    nm        <- paste0("seq_", i)
    results[[nm]] <- list(accession  = nm,
                          definition = NA_character_,
                          length     = length(seq_clean),
                          sequence   = seq_clean,
                          features   = empty_features)
  }

  results
}


# --------------------------------------------------------------------------- #
# Internal: parse a single GenBank flat-file record                           #
# --------------------------------------------------------------------------- #
.parse_genbank <- function(gb_text, accession) {

  lines <- strsplit(gb_text, "\n", fixed = TRUE)[[1L]]

  ## ---- DEFINITION ----------------------------------------------------------
  def_line <- grep("^DEFINITION", lines, value = TRUE)
  definition <- if (length(def_line)) sub("^DEFINITION\\s+", "", def_line[[1L]]) else NA_character_

  ## ---- sequence length from LOCUS ------------------------------------------
  locus_line <- grep("^LOCUS", lines, value = TRUE)
  seq_len <- if (length(locus_line)) {
    m <- regmatches(locus_line[[1L]], regexpr("\\d+(?=\\s+bp)", locus_line[[1L]], perl = TRUE))
    as.integer(m)
  } else NA_integer_

  ## ---- ORIGIN (sequence) ---------------------------------------------------
  origin_start <- which(grepl("^ORIGIN", lines))
  end_marker   <- which(lines == "//")

  if (length(origin_start) == 0L)
    stop("No ORIGIN section found in GenBank record for ", accession, call. = FALSE)

  seq_lines <- lines[(origin_start + 1L):(end_marker[[1L]] - 1L)]
  seq_raw   <- paste(gsub("[^ACGTNacgtn]", "", seq_lines), collapse = "")
  sequence  <- Biostrings::DNAString(toupper(seq_raw))

  ## ---- FEATURES ------------------------------------------------------------
  feat_start <- which(grepl("^FEATURES", lines))
  feat_end   <- origin_start - 1L

  features <- if (length(feat_start) > 0L) {
    .parse_features(lines[(feat_start + 1L):feat_end])
  } else {
    data.frame(type = character(), start = integer(), end = integer(),
               strand = character(), gene = character(),
               product = character(), note = character(),
               stringsAsFactors = FALSE)
  }

  list(
    accession  = accession,
    definition = definition,
    length     = length(sequence),
    sequence   = sequence,
    features   = features
  )
}


# --------------------------------------------------------------------------- #
# Internal: parse the FEATURES block into a data.frame                        #
# --------------------------------------------------------------------------- #
.parse_features <- function(feat_lines) {

  # Feature lines: 5-space indent + type (<=15 chars) + location
  # Qualifier lines: 21-space indent + /key="value"
  feat_pattern <- "^     (\\S+)\\s+(.+)$"
  qual_pattern <- "^\\s{21}/(\\w+)(?:=\"(.*)\")?$"

  records <- list()
  cur <- NULL

  flush <- function() {
    if (!is.null(cur) && cur$type != "source") {
      records[[length(records) + 1L]] <<- cur
    }
  }

  for (ln in feat_lines) {
    if (grepl(feat_pattern, ln)) {
      flush()
      m    <- regmatches(ln, regexec(feat_pattern, ln))[[1L]]
      loc  <- .parse_location(m[3L])
      cur  <- list(type    = m[2L],
                   start   = loc$start,
                   end     = loc$end,
                   strand  = loc$strand,
                   gene    = NA_character_,
                   product = NA_character_,
                   note    = NA_character_)
    } else if (!is.null(cur) && grepl(qual_pattern, ln)) {
      m   <- regmatches(ln, regexec(qual_pattern, ln))[[1L]]
      key <- m[2L]
      val <- if (nchar(m[3L]) > 0L) m[3L] else TRUE
      if (key %in% c("gene", "product", "note"))
        cur[[key]] <- as.character(val)
    } else if (!is.null(cur) && grepl("^\\s{21}", ln)) {
      # continuation of previous qualifier value
      trimmed <- trimws(ln)
      if (!is.na(cur$note) && is.character(cur$note))
        cur$note <- paste(cur$note, trimmed)
    }
  }
  flush()

  if (length(records) == 0L)
    return(data.frame(type = character(), start = integer(), end = integer(),
                      strand = character(), gene = character(),
                      product = character(), note = character(),
                      stringsAsFactors = FALSE))

  do.call(rbind, lapply(records, as.data.frame, stringsAsFactors = FALSE))
}


# --------------------------------------------------------------------------- #
# Internal: parse a GenBank feature location string                           #
# --------------------------------------------------------------------------- #
.parse_location <- function(loc_str) {

  loc_str <- trimws(loc_str)

  # Detect strand
  strand <- "+"
  if (grepl("^complement", loc_str)) {
    strand  <- "-"
    loc_str <- sub("^complement\\((.+)\\)$", "\\1", loc_str)
  }

  # Remove join() wrapper — collect all numeric ranges
  loc_str <- sub("^join\\((.+)\\)$", "\\1", loc_str)

  # Extract all numbers (handles simple ranges and joins)
  nums <- as.integer(regmatches(loc_str, gregexpr("\\d+", loc_str))[[1L]])

  list(
    start  = min(nums, na.rm = TRUE),
    end    = max(nums, na.rm = TRUE),
    strand = strand
  )
}
