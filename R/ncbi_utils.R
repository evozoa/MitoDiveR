# ============================================================================
# Internal helpers shared by the NCBI fetch functions
# (fetch_sequences, fetch_mito_genbank, fetch_gene)
# ============================================================================

# RefSeq (NC_, NM_, ...) and GenBank (1-2 letter prefix + digits, WGS-style
# prefixes) accessions, with optional version suffix.
.acc_pattern <- "^[A-Z]{1,2}_?[A-Z0-9]*\\d{4,9}(\\.\\d+)?$"

# Thin wrappers so tests can mock the network layer.
.efetch  <- function(...) rentrez::entrez_fetch(...)
.esearch <- function(...) rentrez::entrez_search(...)


# --------------------------------------------------------------------------- #
# Split, trim and upper-case accession input.                                 #
# "nc_028285, NC_012920.1" -> c("NC_028285", "NC_012920.1")                   #
# --------------------------------------------------------------------------- #
.split_accessions <- function(x) {
  x <- as.character(unlist(x, use.names = FALSE))
  tokens <- unlist(strsplit(trimws(x), "[,;[:space:]]+"), use.names = FALSE)
  toupper(tokens[nzchar(tokens)])
}

.is_accession <- function(x) grepl(.acc_pattern, x)

# Clean an 'accessions' argument, erroring on the first element (as the user
# typed it) that does not consist of accession numbers.
.parse_accession_arg <- function(x) {
  x <- as.character(unlist(x, use.names = FALSE))
  out <- character()
  for (el in x[!is.na(x)]) {
    toks <- .split_accessions(el)
    if (length(toks) && !all(.is_accession(toks))) .stop_bad_input(el)
    out <- c(out, toks)
  }
  unique(out)
}

# Mostly nucleotide characters (after dropping FASTA headers and whitespace)?
.looks_like_dna <- function(x) {
  lines <- strsplit(x, "\n", fixed = TRUE)[[1L]]
  s <- gsub("[[:space:][:digit:]]", "",
            paste(lines[!grepl("^\\s*>", lines)], collapse = ""))
  nchar(s) > 0L &&
    nchar(gsub("[^ACGTUNRYKMSWBDHVacgtunrykmswbdhv-]", "", s)) >= 0.9 * nchar(s)
}

# Helpful error for input that is neither an accession, a file, nor DNA.
.stop_bad_input <- function(x, what = "an NCBI accession") {
  shown <- if (nchar(x) > 60L) paste0(substr(x, 1L, 57L), "...") else x
  hint <- if (grepl("[/\\\\]|\\.(fa|fas|fasta|fna|txt|gb|gbk)$", x,
                    ignore.case = TRUE)) {
    sprintf("\nIt looks like a file path, but no such file exists (working directory: %s).",
            getwd())
  } else if (grepl("^[A-Za-z]+( [A-Za-z.]+)+$", trimws(x))) {
    paste0("\nNCBI fetches need accession numbers, not names. To search by ",
           "species and gene, use fetch_gene(), e.g.\n",
           "  fetch_gene(\"COI\", species = \"", trimws(x), "\")")
  } else ""
  stop(sprintf("'%s' is not %s.%s", shown, what, hint), call. = FALSE)
}


# --------------------------------------------------------------------------- #
# Call an Entrez function, retrying transient failures (rate limits, server   #
# errors, dropped connections) with exponential back-off.                     #
# --------------------------------------------------------------------------- #
.entrez_retry <- function(fun, ..., .tries = 4L, .wait = 1) {
  transient <- paste0("HTTP failure: (429|5\\d\\d)|timed? ?out|Could not resolve|",
                      "Failed to connect|Connection reset|Recv failure|",
                      "SSL|empty reply")
  for (attempt in seq_len(.tries)) {
    res <- tryCatch(fun(...), error = function(e) e)
    if (!inherits(res, "error")) return(res)
    if (attempt == .tries ||
        !grepl(transient, conditionMessage(res), ignore.case = TRUE))
      stop(res)
    Sys.sleep(.wait * 2^(attempt - 1L))
  }
}

# Turn rentrez's URL-encoded HTTP errors into something readable.
.ncbi_error_reason <- function(e) {
  msg  <- conditionMessage(e)
  code <- regmatches(msg, regexpr("(?<=HTTP failure: )\\d{3}", msg, perl = TRUE))
  if (!length(code)) return(msg)
  switch(substr(code, 1L, 1L),
         "4" = if (code == "429")
                 "HTTP 429: rate-limited by NCBI; set an API key or retry later"
               else sprintf("HTTP %s: accession not found or invalid", code),
         "5" = sprintf("HTTP %s: NCBI server error; retry later", code),
         sprintf("HTTP %s", code))
}

# Pause between successive Entrez requests (NCBI allows 3/s, 10/s with a key).
.entrez_pause <- function(api_key = NULL) {
  has_key <- !is.null(api_key) || nzchar(Sys.getenv("ENTREZ_KEY"))
  Sys.sleep(if (has_key) 0.11 else 0.34)
}
