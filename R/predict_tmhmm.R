# ============================================================================
# predict_tmhmm — transmembrane topology prediction via DeepTMHMM
# ============================================================================

#' Predict Transmembrane Topology with DeepTMHMM
#'
#' Submits protein sequences to DeepTMHMM (via the local `biolib` CLI) and
#' returns a tidy summary of predicted transmembrane topology.  Accepts either
#' a plain character vector of sequences or the `data.frame` returned by
#' [find_orfs()]; in the latter case, results are merged back onto the ORF
#' table so that every ORF — canonical or novel — carries a topology
#' annotation.
#'
#' @param x Either a **character vector** of protein sequences or the
#'   **`data.frame`** returned by [find_orfs()].  Stop-codon asterisks
#'   (`"*"`) are stripped automatically before submission.
#' @param ids Character vector of sequence identifiers, one per sequence.
#'   Ignored when `x` is a `data.frame` (row indices are used).  If `NULL`
#'   and `x` is a character vector, identifiers default to `seq_1`, `seq_2`,
#'   ...
#' @param biolib_bin Path to the `biolib` executable.  Defaults to
#'   `"biolib"`, assuming it is on `PATH`.  Install with:
#'   `pip install pybiolib`.
#'
#' @return A `data.frame` with one row per input sequence and the following
#'   columns:
#' \describe{
#'   \item{`seq_id`}{Sequence identifier.}
#'   \item{`seq_length`}{Number of amino acids submitted.}
#'   \item{`topology_class`}{Overall predicted class:
#'     `"TM"` (alpha-helical transmembrane),
#'     `"SP"` (signal peptide, no TM),
#'     `"SP+TM"` (signal peptide plus TM helices),
#'     `"BETA"` (beta-barrel membrane protein), or
#'     `"GLOB"` (globular / soluble).}
#'   \item{`n_tm_segments`}{Number of predicted TM segments (helices or
#'     beta-strands).}
#'   \item{`topology_string`}{Per-residue label string returned by
#'     DeepTMHMM: `M` = membrane helix, `i` = inside (cytoplasmic),
#'     `o` = outside (lumenal/IMS), `S` = signal peptide, `B` = beta-barrel
#'     strand.}
#' }
#' When `x` is a [find_orfs()] `data.frame` these columns are appended to
#' the right of the original table.
#'
#' The returned object also carries a `"domains"` attribute: a long-format
#' `data.frame` of per-residue domain boundaries parsed from the DeepTMHMM
#' GFF3 output.  Access it with [tmhmm_domains()].
#'
#' @details
#' ## Requirements
#' Install the BioLib CLI once into any Python ≥ 3.8 environment visible
#' to R:
#' ```
#' pip install pybiolib
#' ```
#' The first run downloads the DeepTMHMM model weights (~100 MB) into a
#' local cache; subsequent runs are faster.  A free BioLib account is
#' required for authenticated runs; create one at <https://biolib.com>.
#' Run `biolib login` in a terminal before the first use.
#'
#' ## Why DeepTMHMM over TMHMM 2.0?
#' DeepTMHMM (Hallgren *et al.*, 2022) jointly predicts alpha-helical TM
#' proteins, beta-barrel outer membrane proteins, and signal peptides in a
#' single model, achieving higher accuracy than TMHMM 2.0, particularly for
#' short and non-canonical sequences — the primary target of this package.
#'
#' @references
#' Hallgren J, Tsirigos KD, Pedersen MD, *et al.* (2022) DeepTMHMM predicts
#' alpha and beta transmembrane proteins using a deep neural network.
#' *bioRxiv*. <https://doi.org/10.1101/2022.04.08.487609>
#'
#' @seealso [tmhmm_domains()] to extract per-residue domain boundaries.
#'
#' @examples
#' \dontrun{
#' library(Biostrings)
#'
#' # Annotate all ORFs from a mitogenome in one pipeline
#' mito <- readDNAStringSet("my_mito.fasta")[[1]]
#' orfs <- find_orfs(mito)
#' topo <- predict_tmhmm(orfs)
#'
#' # Inspect only TM-containing ORFs
#' topo[topo$topology_class %in% c("TM", "SP+TM"), ]
#'
#' # View domain boundaries for every sequence
#' tmhmm_domains(topo)
#' }
#'
#' @export
predict_tmhmm <- function(x,
                           ids        = NULL,
                           biolib_bin = "biolib") {

  ## ---- extract sequences ---------------------------------------------------
  if (is.data.frame(x)) {
    if (!"protein_sequence" %in% names(x))
      stop("'x' must contain a 'protein_sequence' column (output of find_orfs()).",
           call. = FALSE)
    seqs  <- x$protein_sequence
    is_df <- TRUE
    ids   <- paste0("orf_", seq_len(nrow(x)))
  } else {
    seqs  <- as.character(x)
    is_df <- FALSE
    if (is.null(ids))
      ids <- paste0("seq_", seq_len(length(seqs)))
  }

  if (length(seqs) == 0L)
    stop("No sequences supplied.", call. = FALSE)

  # Strip trailing stop-codon asterisk; DeepTMHMM rejects non-AA characters
  seqs_clean <- gsub("\\*$", "", seqs)

  ## ---- write temp FASTA ----------------------------------------------------
  fasta_in <- tempfile(fileext = ".fasta")
  on.exit(unlink(fasta_in), add = TRUE)

  fasta_lines <- unlist(
    mapply(function(id, sq) c(paste0(">", id), sq),
           ids, seqs_clean,
           SIMPLIFY = FALSE, USE.NAMES = FALSE)
  )
  writeLines(fasta_lines, fasta_in)

  ## ---- run DeepTMHMM -------------------------------------------------------
  result <- .run_deeptmhmm(fasta_in, biolib_bin)

  ## ---- merge with ORF table if applicable ----------------------------------
  if (is_df) {
    doms <- attr(result, "domains")
    result[["seq_id"]] <- NULL          # row order matches; drop key column

    out <- cbind(x, result)
    if (!is.null(doms)) attr(out, "domains") <- doms
    return(out)
  }

  result
}


# --------------------------------------------------------------------------- #
# Internal: invoke biolib CLI and collect output                               #
# --------------------------------------------------------------------------- #
.run_deeptmhmm <- function(fasta_in, biolib_bin) {

  run_dir <- tempfile(pattern = "MitoDiveR_deeptmhmm_")
  dir.create(run_dir)
  on.exit(unlink(run_dir, recursive = TRUE), add = TRUE)

  # Run biolib from within run_dir so output lands there
  cmd <- sprintf(
    'cd "%s" && "%s" run DTU/DeepTMHMM --fasta "%s" 2>&1',
    run_dir, biolib_bin, normalizePath(fasta_in, winslash = "/")
  )
  exit_code <- system(cmd, ignore.stdout = FALSE)

  if (exit_code != 0L)
    stop(
      "DeepTMHMM returned a non-zero exit code (", exit_code, ").\n",
      "  - Is pybiolib installed?  pip install pybiolib\n",
      "  - Are you logged in?      biolib login\n",
      "  - Test manually:          biolib run DTU/DeepTMHMM --help",
      call. = FALSE
    )

  # Locate output files — biolib nests them under biolib_results/<app>/<job>/
  all_files <- list.files(run_dir, recursive = TRUE, full.names = TRUE)

  topo_file <- grep("predicted_topologies\\.3line$", all_files, value = TRUE)
  gff_file  <- grep("TMRs\\.gff3$",                 all_files, value = TRUE)

  if (length(topo_file) == 0L)
    stop(
      "Expected output file 'predicted_topologies.3line' not found.\n",
      "Files present: ", paste(basename(all_files), collapse = ", "),
      call. = FALSE
    )

  result <- .parse_deeptmhmm_3line(topo_file[[1L]])

  if (length(gff_file) > 0L)
    attr(result, "domains") <- .parse_tmhmm_gff3(gff_file[[1L]])

  result
}


# --------------------------------------------------------------------------- #
# Internal: parse the 3-line topology file                                     #
# --------------------------------------------------------------------------- #
# Format (three lines per sequence, blank lines ignored):
#   >sequence_id
#   iiiMMMMMMMMMMMMMMooooooo...   <- per-residue topology labels
#   MKLLLLLGAAAFVVVVVV...         <- amino acid sequence
.parse_deeptmhmm_3line <- function(path) {

  lines <- readLines(path, warn = FALSE)
  lines <- lines[nchar(trimws(lines)) > 0L]   # drop blanks

  if (length(lines) %% 3L != 0L)
    stop("Malformed 3-line file: ", length(lines),
         " non-blank lines (expected a multiple of 3).", call. = FALSE)

  n_seqs <- length(lines) %/% 3L
  rows   <- vector("list", n_seqs)

  for (i in seq_len(n_seqs)) {
    base     <- (i - 1L) * 3L
    seq_id   <- sub("^>", "", lines[base + 1L])
    topo_str <- lines[base + 2L]
    topo_chr <- strsplit(topo_str, "")[[1L]]

    has_tm   <- any(topo_chr == "M")
    has_beta <- any(topo_chr == "B")
    has_sp   <- any(topo_chr == "S")

    topo_class <- if (has_beta)            "BETA"
                  else if (has_tm && has_sp) "SP+TM"
                  else if (has_tm)          "TM"
                  else if (has_sp)          "SP"
                  else                      "GLOB"

    rle_res      <- rle(topo_chr)
    n_tm_segs    <- sum(rle_res$values %in% c("M", "B"))

    rows[[i]] <- data.frame(
      seq_id         = seq_id,
      seq_length     = length(topo_chr),
      topology_class = topo_class,
      n_tm_segments  = n_tm_segs,
      topology_string = topo_str,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}


# --------------------------------------------------------------------------- #
# Internal: parse the GFF3 domain-boundary file                               #
# --------------------------------------------------------------------------- #
# Retains only data lines (skips ## comment lines) and extracts
# seq_id, domain type, start, and end.
.parse_tmhmm_gff3 <- function(path) {

  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^#", lines) & nchar(trimws(lines)) > 0L]

  if (length(lines) == 0L)
    return(data.frame(seq_id = character(), type = character(),
                      start  = integer(),   end  = integer(),
                      stringsAsFactors = FALSE))

  fields <- strsplit(lines, "\t")

  data.frame(
    seq_id = vapply(fields, `[[`, character(1L), 1L),
    type   = vapply(fields, `[[`, character(1L), 3L),
    start  = as.integer(vapply(fields, `[[`, character(1L), 4L)),
    end    = as.integer(vapply(fields, `[[`, character(1L), 5L)),
    stringsAsFactors = FALSE
  )
}


# --------------------------------------------------------------------------- #
# Exported accessor: domain boundaries                                         #
# --------------------------------------------------------------------------- #

#' Extract Per-Residue Domain Boundaries from a predict_tmhmm() Result
#'
#' Returns the long-format domain-boundary table that DeepTMHMM emits in its
#' GFF3 output, parsed and attached to the [predict_tmhmm()] result as an
#' attribute.
#'
#' @param x The `data.frame` returned by [predict_tmhmm()].
#'
#' @return A `data.frame` with columns:
#' \describe{
#'   \item{`seq_id`}{Sequence identifier.}
#'   \item{`type`}{Domain type: e.g. `"TMhelix"`, `"inside"`, `"outside"`,
#'     `"signal"`, `"Beta strand"`.}
#'   \item{`start`}{1-based start residue of the domain.}
#'   \item{`end`}{1-based end residue of the domain (inclusive).}
#' }
#' Returns `NULL` invisibly (with a message) if no domain data is present.
#'
#' @seealso [predict_tmhmm()]
#'
#' @examples
#' \dontrun{
#' topo <- predict_tmhmm(orfs)
#' tmhmm_domains(topo)
#' }
#'
#' @export
tmhmm_domains <- function(x) {
  doms <- attr(x, "domains")
  if (is.null(doms)) {
    message("No domain data attached. Was predict_tmhmm() run successfully?")
    return(invisible(NULL))
  }
  doms
}
