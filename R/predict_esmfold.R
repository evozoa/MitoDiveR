# ============================================================================
# predict_esmfold — structure prediction via ESMFold public API
# ============================================================================

#' Predict Protein Structure with ESMFold
#'
#' Submits one or more protein sequences to the ESMFold public REST API
#' (Meta AI) and returns per-sequence confidence statistics derived from the
#' predicted pLDDT scores.  Accepts the same inputs as [predict_tmhmm()] so
#' the two functions can be applied sequentially in a pipeline.
#'
#' @param x Either a **character vector** of protein sequences or a
#'   **`data.frame`** returned by [find_orfs()], [cluster_conserved_orfs()],
#'   or [find_conserved_windows()].  Stop-codon characters (`*`, `X`) and
#'   gap characters (`-`) are stripped before submission.
#' @param ids Character vector of sequence identifiers, one per sequence.
#'   Ignored when `x` is a `data.frame` (identifiers are constructed
#'   automatically).  Defaults to `seq_1`, `seq_2`, ... when `NULL`.
#' @param sequence_col Character.  When `x` is a `data.frame`, the name of
#'   the column containing protein sequences.  Defaults to `NULL`, which
#'   triggers automatic detection: `protein_sequence` (find_orfs output),
#'   then `window_aa_a` (find_conserved_windows output).
#' @param plddt_threshold Numeric (0–100).  Residues at or above this pLDDT
#'   are counted as "confident" for the `frac_confident` column.
#'   Default `70`.
#' @param save_pdb Logical.  If `TRUE`, write each predicted PDB to
#'   `pdb_dir` rather than storing the full string in the result.
#'   Default `FALSE`.
#' @param pdb_dir Character.  Directory to write PDB files when
#'   `save_pdb = TRUE`.  Created if it does not exist.  Default `NULL`
#'   (uses the current working directory).
#' @param api_url Character.  ESMFold API endpoint.
#'   Default `"https://api.esmatlas.com/foldSequence/v1/pdb/"`.
#' @param verbose Logical.  Print progress messages?  Default `TRUE`.
#'
#' @return A `data.frame` with one row per input sequence:
#' \describe{
#'   \item{`seq_id`}{Sequence identifier.}
#'   \item{`seq_length`}{Number of amino acids submitted.}
#'   \item{`plddt_mean`}{Mean pLDDT across all residues (0–100).}
#'   \item{`plddt_min`}{Minimum per-residue pLDDT.}
#'   \item{`plddt_max`}{Maximum per-residue pLDDT.}
#'   \item{`frac_confident`}{Fraction of residues with pLDDT >=
#'     `plddt_threshold`.}
#'   \item{`structure_quality`}{Ordered factor summarising overall
#'     confidence: `"high"` (mean >= 90), `"medium"` (70–90),
#'     `"low"` (50–70), or `"disordered"` (< 50).}
#'   \item{`pdb`}{Full PDB text returned by ESMFold, or the file path if
#'     `save_pdb = TRUE`.  `NA` on API failure for that sequence.}
#' }
#' When `x` is a `data.frame`, these columns are appended to the right of
#' the original table.
#'
#' The returned object carries a `"plddt"` attribute: a long-format
#' `data.frame` of per-residue pLDDT values.  Access it with
#' [esmfold_plddt()].
#'
#' @details
#' ## API limits
#' The ESMFold public endpoint is rate-limited and does not require
#' authentication.  Sequences longer than ~400 aa may time out or return an
#' error; consider breaking long sequences into shorter segments.  A 1-second
#' pause is inserted between requests.
#'
#' ## pLDDT interpretation
#' Per AlphaFold/ESMFold conventions:
#' \itemize{
#'   \item pLDDT >= 90 — very high confidence, well-folded
#'   \item 70–90 — generally good backbone accuracy
#'   \item 50–70 — low confidence, likely disordered or flexible
#'   \item < 50 — should not be interpreted as a folded structure
#' }
#'
#' @references
#' Lin Z, Akin H, Rao R, *et al.* (2023) Evolutionary-scale prediction of
#' atomic-level protein structure with a language model.
#' *Science* **379**, 1123–1130.
#'
#' @seealso [esmfold_plddt()] to extract per-residue confidence scores;
#'   [predict_tmhmm()] for transmembrane topology prediction.
#'
#' @examples
#' \dontrun{
#' # Humanin sequence
#' sf <- predict_esmfold("MAPRGFSCLLLLTSEIDLPVKRRA")
#' sf$plddt_mean
#'
#' # Pipeline: conserved ORFs -> TMHMM -> ESMFold
#' orfs  <- find_orfs(genome)
#' topo  <- predict_tmhmm(orfs)
#' tm_only <- topo[topo$topology_class == "TM", ]
#' sf    <- predict_esmfold(tm_only)
#'
#' # Per-residue pLDDT
#' esmfold_plddt(sf)
#' }
#'
#' @export
predict_esmfold <- function(x,
                             ids             = NULL,
                             sequence_col    = NULL,
                             plddt_threshold = 70,
                             save_pdb        = FALSE,
                             pdb_dir         = NULL,
                             api_url         = "https://api.esmatlas.com/foldSequence/v1/pdb/",
                             verbose         = TRUE) {

  ## ---- extract sequences ---------------------------------------------------
  is_df <- is.data.frame(x)

  if (is_df) {
    col <- .detect_seq_col(x, sequence_col)
    seqs <- x[[col]]
    ids  <- paste0("seq_", seq_len(nrow(x)))
  } else {
    seqs <- as.character(x)
    if (is.null(ids))
      ids <- paste0("seq_", seq_len(length(seqs)))
  }

  if (length(seqs) == 0L)
    stop("No sequences supplied.", call. = FALSE)

  # Strip stop codons, gaps, and ambiguous X characters
  seqs_clean <- gsub("[*X\\-]", "", seqs)

  zero_len <- nchar(seqs_clean) == 0L
  if (any(zero_len))
    warning(sum(zero_len), " sequence(s) are empty after stripping and will ",
            "be returned as NA.", call. = FALSE)

  ## ---- PDB output directory ------------------------------------------------
  if (save_pdb) {
    if (is.null(pdb_dir)) pdb_dir <- getwd()
    if (!dir.exists(pdb_dir)) dir.create(pdb_dir, recursive = TRUE)
  }

  ## ---- query ESMFold per sequence ------------------------------------------
  n         <- length(seqs_clean)
  result_rows <- vector("list", n)
  plddt_rows  <- list()

  for (i in seq_len(n)) {
    seq_i <- seqs_clean[[i]]
    id_i  <- ids[[i]]

    if (verbose)
      message(sprintf("  [%d/%d] %s (%d aa)...", i, n, id_i, nchar(seq_i)))

    if (nchar(seq_i) == 0L) {
      result_rows[[i]] <- .esmfold_na_row(id_i, nchar(seqs[[i]]))
      if (i < n) Sys.sleep(1)
      next
    }

    pdb_text <- tryCatch(
      .call_esmfold_api(seq_i, api_url),
      error = function(e) {
        warning("ESMFold API error for ", id_i, ": ", conditionMessage(e),
                call. = FALSE)
        NA_character_
      }
    )

    if (is.na(pdb_text)) {
      result_rows[[i]] <- .esmfold_na_row(id_i, nchar(seq_i))
      if (i < n) Sys.sleep(1)
      next
    }

    ## -- extract per-residue pLDDT from B-factor column --------------------
    plddt_vec <- .extract_plddt(pdb_text)

    plddt_mean <- mean(plddt_vec)
    plddt_min  <- min(plddt_vec)
    plddt_max  <- max(plddt_vec)
    frac_conf  <- mean(plddt_vec >= plddt_threshold)

    quality <- .plddt_quality(plddt_mean)

    ## -- save or store PDB ------------------------------------------------
    pdb_out <- if (save_pdb) {
      pdb_path <- file.path(pdb_dir, paste0(id_i, ".pdb"))
      writeLines(pdb_text, pdb_path)
      pdb_path
    } else {
      pdb_text
    }

    result_rows[[i]] <- data.frame(
      seq_id           = id_i,
      seq_length       = nchar(seq_i),
      plddt_mean       = round(plddt_mean, 1L),
      plddt_min        = round(plddt_min,  1L),
      plddt_max        = round(plddt_max,  1L),
      frac_confident   = round(frac_conf,  3L),
      structure_quality = quality,
      pdb              = pdb_out,
      stringsAsFactors = FALSE
    )

    ## -- per-residue pLDDT table ------------------------------------------
    plddt_rows[[length(plddt_rows) + 1L]] <- data.frame(
      seq_id   = id_i,
      position = seq_along(plddt_vec),
      residue  = strsplit(seq_i, "")[[1L]],
      plddt    = round(plddt_vec, 1L),
      stringsAsFactors = FALSE
    )

    if (i < n) Sys.sleep(1)
  }

  ## ---- assemble output -----------------------------------------------------
  out <- do.call(rbind, result_rows)

  # structure_quality as ordered factor
  lvls <- c("disordered", "low", "medium", "high")
  out$structure_quality <- factor(out$structure_quality,
                                   levels = lvls, ordered = TRUE)

  plddt_tbl <- if (length(plddt_rows) > 0L)
    do.call(rbind, plddt_rows) else NULL

  if (!is.null(plddt_tbl))
    attr(out, "plddt") <- plddt_tbl

  if (is_df) {
    out[["seq_id"]] <- NULL   # row order matches; drop redundant key
    out <- cbind(x, out)
    if (!is.null(plddt_tbl)) attr(out, "plddt") <- plddt_tbl
  }

  out
}


# --------------------------------------------------------------------------- #
# Exported accessor: per-residue pLDDT                                         #
# --------------------------------------------------------------------------- #

#' Extract Per-Residue pLDDT Scores from a predict_esmfold() Result
#'
#' Returns the long-format per-residue confidence table attached to the
#' [predict_esmfold()] result as the `"plddt"` attribute.
#'
#' @param x The `data.frame` returned by [predict_esmfold()].
#'
#' @return A `data.frame` with columns:
#' \describe{
#'   \item{`seq_id`}{Sequence identifier.}
#'   \item{`position`}{1-based residue position.}
#'   \item{`residue`}{Single-letter amino acid.}
#'   \item{`plddt`}{Per-residue pLDDT score (0–100).}
#' }
#' Returns `NULL` invisibly (with a message) if no pLDDT data is present.
#'
#' @seealso [predict_esmfold()]
#'
#' @examples
#' \dontrun{
#' sf <- predict_esmfold("MAPRGFSCLLLLTSEIDLPVKRRA")
#' esmfold_plddt(sf)
#' }
#'
#' @export
esmfold_plddt <- function(x) {
  tbl <- attr(x, "plddt")
  if (is.null(tbl)) {
    message("No pLDDT data attached. Was predict_esmfold() run successfully?")
    return(invisible(NULL))
  }
  tbl
}


# --------------------------------------------------------------------------- #
# Internals                                                                    #
# --------------------------------------------------------------------------- #

.call_esmfold_api <- function(sequence, api_url) {
  tmp <- tempfile(fileext = ".pdb")
  on.exit(unlink(tmp), add = TRUE)

  # system2 passes each element of args as a separate argument — no quoting needed
  exit_code <- system2(
    "curl",
    args   = c("-s", "-X", "POST", api_url,
               "-H", "Content-Type: text/plain",
               "--data-raw", sequence,
               "-o", tmp),
    stdout = FALSE,
    stderr = FALSE
  )

  if (!file.exists(tmp) || file.info(tmp)$size == 0L)
    stop("curl produced no output (exit code ", exit_code, ").", call. = FALSE)

  pdb <- paste(readLines(tmp, warn = FALSE), collapse = "\n")

  if (!grepl("ATOM", pdb, fixed = TRUE))
    stop("Unexpected API response (no ATOM records). ",
         "Possible rate limit or invalid sequence.", call. = FALSE)

  pdb
}

.extract_plddt <- function(pdb_text) {
  lines <- strsplit(pdb_text, "\n", fixed = TRUE)[[1L]]
  # Alpha-carbon (CA) atoms only — one per residue
  ca_lines <- lines[startsWith(lines, "ATOM") &
                      trimws(substr(lines, 13L, 16L)) == "CA"]
  if (length(ca_lines) == 0L)
    stop("No CA atoms found in PDB — cannot extract pLDDT.", call. = FALSE)

  # PDB columns 61-66 = B-factor (pLDDT stored as 0-1 by ESMFold)
  bfac <- as.numeric(trimws(substr(ca_lines, 61L, 66L)))
  bfac * 100   # convert to 0-100 scale
}

.plddt_quality <- function(mean_plddt) {
  if (mean_plddt >= 90)       "high"
  else if (mean_plddt >= 70)  "medium"
  else if (mean_plddt >= 50)  "low"
  else                         "disordered"
}

.esmfold_na_row <- function(id, len) {
  data.frame(
    seq_id            = id,
    seq_length        = len,
    plddt_mean        = NA_real_,
    plddt_min         = NA_real_,
    plddt_max         = NA_real_,
    frac_confident    = NA_real_,
    structure_quality = NA_character_,
    pdb               = NA_character_,
    stringsAsFactors  = FALSE
  )
}

.detect_seq_col <- function(df, sequence_col) {
  if (!is.null(sequence_col)) {
    if (!sequence_col %in% names(df))
      stop("Column '", sequence_col, "' not found in data.frame.", call. = FALSE)
    return(sequence_col)
  }
  # Auto-detect: find_orfs output
  if ("protein_sequence" %in% names(df)) return("protein_sequence")
  # find_conserved_windows output
  if ("window_aa_a" %in% names(df)) return("window_aa_a")
  stop(
    "Cannot auto-detect sequence column. Provide 'sequence_col' explicitly.\n",
    "Available columns: ", paste(names(df), collapse = ", "),
    call. = FALSE
  )
}
