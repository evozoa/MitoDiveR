# ============================================================================
# fold_rna — RNA secondary structure prediction via ViennaRNA (RNAfold)
# ============================================================================

# ── internal helpers ──────────────────────────────────────────────────────────

# Parse dot-bracket notation into a base-pair table
.parse_dot_bracket <- function(db) {
  chars <- strsplit(db, "")[[1]]
  stack <- integer(0)
  pos5  <- integer(0)
  pos3  <- integer(0)
  for (i in seq_along(chars)) {
    if      (chars[i] == "(") { stack <- c(stack, i) }
    else if (chars[i] == ")") {
      if (length(stack) == 0L) stop("Unbalanced dot-bracket structure.", call. = FALSE)
      j     <- stack[length(stack)]
      stack <- stack[-length(stack)]
      pos5  <- c(pos5, j)
      pos3  <- c(pos3, i)
    }
  }
  data.frame(pos5 = pos5, pos3 = pos3, stringsAsFactors = FALSE)
}

# Count hairpin loops: innermost stem-loops of the form (.....) with no
# nested parentheses inside
.count_hairpins <- function(db) {
  m <- gregexpr("\\([.]+\\)", db)[[1]]
  if (m[1L] == -1L) 0L else length(m)
}

# Convert a sequence string to RNA (T -> U, lowercase tolerated)
.to_rna <- function(seq) toupper(chartr("Tt", "UU", seq))

# Check for RNAfold binary and return its path, or stop with instructions
.check_rnafold <- function() {
  path <- Sys.which("RNAfold")
  if (nchar(path) == 0L)
    stop(
      "RNAfold (ViennaRNA) not found on PATH.\n",
      "Install ViennaRNA: https://www.tbi.univie.ac.at/RNA/#download\n",
      "  conda: conda install -c bioconda viennarna\n",
      "  apt  : sudo apt-get install vienna-rna",
      call. = FALSE
    )
  path
}

# ── main function ─────────────────────────────────────────────────────────────

#' Predict RNA Secondary Structure
#'
#' Folds a nucleotide sequence into its minimum free-energy (MFE) secondary
#' structure using ViennaRNA's `RNAfold` program, which implements the same
#' Turner-energy-based algorithm as Mfold and UNAFold.  DNA sequences are
#' automatically converted to RNA (T → U) before folding.
#'
#' @param sequence A nucleotide sequence as a character string, or a
#'   [Biostrings::DNAString], [Biostrings::RNAString], or
#'   [Biostrings::DNAStringSet] / [Biostrings::RNAStringSet] object.  If a
#'   `StringSet` is supplied, each element is folded separately and a list of
#'   results is returned.  Ambiguity codes are passed through to `RNAfold`
#'   unchanged.
#' @param temperature Numeric.  Folding temperature in °C.  Default `37`.
#' @param circular Logical.  Treat the sequence as circular (uses
#'   `--circ` flag).  Default `FALSE`.
#' @param constraints Character string of dot-bracket constraints to enforce
#'   (same length as sequence; `.` = unconstrained, `x` = unpaired,
#'   `<`/`>` = paired).  `NULL` for unconstrained folding (default).
#'
#' @return For a single sequence, a list of class `"rna_fold"` with:
#' \describe{
#'   \item{`sequence`}{The RNA sequence (U-alphabet) that was folded.}
#'   \item{`structure`}{Dot-bracket secondary structure string.}
#'   \item{`mfe`}{Minimum free energy (kcal/mol, numeric).}
#'   \item{`base_pairs`}{`data.frame` with columns `pos5` and `pos3` giving
#'     the 1-based positions of each base pair.}
#'   \item{`n_pairs`}{Total number of base pairs.}
#'   \item{`n_hairpins`}{Number of hairpin (stem-loop) motifs.}
#'   \item{`length_nt`}{Sequence length (nt).}
#' }
#' For a `StringSet` input, a named list of such objects.
#'
#' @details
#' ## Dependency
#' `fold_rna()` requires the **ViennaRNA** package (`RNAfold` binary) to be
#' installed and on the system `PATH`.  ViennaRNA is freely available at
#' <https://www.tbi.univie.ac.at/RNA/> and via conda (`conda install -c
#' bioconda viennarna`) or apt (`sudo apt-get install vienna-rna`).
#'
#' ## Algorithm
#' `RNAfold` minimises the Gibbs free energy of the structure using the
#' nearest-neighbour Turner energy parameters, identical in principle to the
#' Mfold / UNAFold / RNAstructure family of programs.  Only the MFE structure
#' is returned; for the full partition function / base-pair probabilities use
#' `RNAfold -p` directly.
#'
#' @seealso [scan_orfs()], [find_conserved_orfs()]
#'
#' @examples
#' \dontrun{
#' # Fold a single ORF nucleotide sequence
#' fold <- fold_rna("ATGAAACCCGGGTTTTAA")
#' fold$structure
#' fold$mfe
#' fold$n_hairpins
#'
#' # Fold at a different temperature
#' fold_rna("ATGAAACCCGGGTTTTAA", temperature = 25)
#'
#' # Fold all ORF sequences from scan_orfs() output
#' orfs  <- scan_orfs(seqs, min_orf_length = 60)
#' folds <- fold_rna(Biostrings::DNAStringSet(
#'            setNames(orfs$orf_sequence, orfs$seq_id)))
#' }
#'
#' @importFrom Biostrings DNAString RNAString DNAStringSet RNAStringSet
#' @export
fold_rna <- function(sequence,
                     temperature = 37,
                     circular    = FALSE,
                     constraints = NULL) {

  rnafold_bin <- .check_rnafold()

  ## ── normalise input to a named character vector ──────────────────────────
  if (inherits(sequence, c("DNAStringSet", "RNAStringSet"))) {
    seqs <- as.character(sequence)
    multi <- TRUE
  } else if (inherits(sequence, c("DNAString", "RNAString"))) {
    seqs  <- c(seq1 = as.character(sequence))
    multi <- FALSE
  } else if (is.character(sequence)) {
    if (length(sequence) == 1L) {
      seqs  <- c(seq1 = sequence)
      multi <- FALSE
    } else {
      seqs  <- sequence
      multi <- TRUE
      if (is.null(names(seqs))) names(seqs) <- paste0("seq", seq_along(seqs))
    }
  } else {
    stop("'sequence' must be a character string or Biostrings object.", call. = FALSE)
  }

  if (any(nchar(seqs) == 0L))
    stop("One or more sequences are empty.", call. = FALSE)

  temperature <- as.numeric(temperature)
  if (is.na(temperature))
    stop("'temperature' must be numeric.", call. = FALSE)

  ## ── build RNAfold argument list ──────────────────────────────────────────
  args <- c("--noPS",                              # suppress PostScript output
            paste0("--temp=", temperature))
  if (isTRUE(circular)) args <- c(args, "--circ")

  ## ── fold each sequence ───────────────────────────────────────────────────
  fold_one <- function(nm, nt_seq) {
    rna_seq <- .to_rna(nt_seq)

    # Build stdin: FASTA-format header + sequence (+ optional constraints)
    stdin_lines <- c(paste0(">", nm), rna_seq)
    if (!is.null(constraints)) {
      if (nchar(constraints) != nchar(rna_seq))
        stop("'constraints' length must equal sequence length.", call. = FALSE)
      stdin_lines <- c(stdin_lines, constraints)
      args_run <- c(args, "-C")
    } else {
      args_run <- args
    }

    tmp_in  <- tempfile(fileext = ".fa")
    tmp_out <- tempfile(fileext = ".txt")
    on.exit({ unlink(tmp_in); unlink(tmp_out) }, add = TRUE)
    writeLines(stdin_lines, tmp_in)

    exit_code <- system2(
      rnafold_bin,
      args   = args_run,
      stdin  = tmp_in,
      stdout = tmp_out,
      stderr = FALSE
    )
    if (exit_code != 0L)
      stop("RNAfold exited with code ", exit_code,
           " for sequence '", nm, "'.", call. = FALSE)

    raw <- readLines(tmp_out)
    # Output lines: ">name", sequence, "structure  (mfe)"
    struct_line <- raw[grepl("^[.()+x<>]+\\s+\\(", raw)][1L]
    if (is.na(struct_line))
      stop("Could not parse RNAfold output for '", nm, "'.", call. = FALSE)

    parts     <- regmatches(struct_line,
                            regexec("^([.()+x<>]+)\\s+\\(\\s*([+-]?[0-9.]+)\\)",
                                    struct_line))[[1L]]
    structure <- parts[2L]
    mfe       <- as.numeric(parts[3L])

    bp <- .parse_dot_bracket(structure)

    structure_obj <- list(
      sequence   = rna_seq,
      structure  = structure,
      mfe        = mfe,
      base_pairs = bp,
      n_pairs    = nrow(bp),
      n_hairpins = .count_hairpins(structure),
      length_nt  = nchar(rna_seq)
    )
    class(structure_obj) <- "rna_fold"
    structure_obj
  }

  results <- mapply(fold_one, names(seqs), seqs, SIMPLIFY = FALSE)

  if (multi) results else results[[1L]]
}

# ── S3 print method ───────────────────────────────────────────────────────────

#' @export
print.rna_fold <- function(x, ...) {
  cat(sprintf(
    "RNA fold  |  %d nt  |  MFE = %.2f kcal/mol  |  %d bp  |  %d hairpin(s)\n",
    x$length_nt, x$mfe, x$n_pairs, x$n_hairpins
  ))
  cat("Sequence : ", x$sequence,  "\n", sep = "")
  cat("Structure: ", x$structure, "\n", sep = "")
  invisible(x)
}
