# ============================================================================
# analyze_snp — single-nucleotide variant impact on noncanonical ORFs
# ============================================================================

#' Analyze a Single Nucleotide Variant for Noncanonical ORF Impact
#'
#' Introduces a single-base substitution into a reference mitogenome
#' (rCRS by default) and compares ORFs spanning the variant position
#' between reference and synthetic mutant under one or more genetic
#' codes.  Classifies each affected ORF as a true disruption, a novel
#' (read-through or new-start) ORF, or a substitution-only change.
#'
#' Useful for assessing whether a longevity- or disease-associated
#' mtDNA SNP carries noncanonical coding consequences invisible to
#' standard variant-effect annotation.
#'
#' @param variant Variant specification. Either a character scalar in
#'   the form `"m.16189T>C"` / `"16189T>C"` / `"T16189C"`, or a list
#'   with elements `pos` (integer), `ref` (character), `alt` (character).
#' @param reference Character scalar.  NCBI accession of the reference
#'   genome (default `"NC_012920.1"` = rCRS) or a [Biostrings::DNAString]
#'   / [Biostrings::DNAStringSet] containing a single sequence.
#' @param genetic_codes Character vector of Biostrings genetic-code
#'   identifiers to test.  Default `c("SGC1", "SGC0")` — vertebrate
#'   mitochondrial code and standard code.
#' @param min_orf_length Integer. Minimum ORF length in nucleotides.
#'   Default `60L`.
#' @param start_codons Character vector of recognised start codons.
#'   Default `c("ATG","GTG","ATA","ATT","ATC")`.
#' @param substitution_tolerance Integer. Maximum number of amino-acid
#'   differences allowed between a lost and gained ORF of equal length
#'   for the pair to be classified as substitution-only (and excluded
#'   from the disrupted/novel tally).  Default `3L`.
#' @param rna_structure Logical.  If `TRUE`, also fold the reference and
#'   mutant RNA around the variant and report the secondary-structure
#'   change.  Requires the ViennaRNA `RNAfold` binary (see [fold_rna()]).
#'   Default `FALSE`.
#' @param rna_flank Integer.  Number of nucleotides included on each side of
#'   the variant when folding (window length `= 2 * rna_flank + 1`, clamped at
#'   genome ends).  Default `100L`.
#' @param rna_strand Either `"sense"` (fold the reference strand as deposited)
#'   or `"antisense"` (fold the reverse complement).  Default `"sense"`.
#' @param rna_temperature Numeric folding temperature in degrees Celsius,
#'   passed to [fold_rna()].  Default `37`.
#' @param rna_max_bulge Integer.  `max_bulge` passed to [longest_stem()] when
#'   measuring the longest-stem change.  Default `2L`.
#'
#' @return A list with components:
#'   \describe{
#'     \item{`variant`}{Parsed variant (pos, ref, alt).}
#'     \item{`reference`}{Reference accession or name used.}
#'     \item{`summary`}{A `data.frame` with one row per genetic code,
#'       columns `genetic_code`, `n_ref_orfs`, `n_alt_orfs`, `n_lost`,
#'       `n_gained`, `n_disrupted`, `n_novel`.}
#'     \item{`detail`}{A named list (keyed by genetic code) of lists
#'       containing `ref_orfs`, `alt_orfs`, `disrupted` (character
#'       vector of protein sequences lost without substitution-only
#'       counterpart), and `novel` (character vector of gained ORFs
#'       without a substitution-only counterpart).}
#'     \item{`rna`}{`NULL` unless `rna_structure = TRUE`, otherwise a list with
#'       the folded window (`window`, `variant_offset`, `strand`), folding free
#'       energies (`mfe_ref`, `mfe_alt`, `delta_mfe`), base-pair topology
#'       (`n_pairs_ref`, `n_pairs_alt`, `shared_pairs`, `bp_distance`,
#'       `bp_distance_frac`), longest-stem change (`longest_stem_ref`,
#'       `longest_stem_alt`, `delta_stem`), and the two dot-bracket structures
#'       (`structure_ref`, `structure_alt`).}
#'   }
#'
#' @details
#' The reference base at the specified position is verified against
#' the reference genome; mismatches raise an error.  For each genetic
#' code, [scan_orfs()] is run on both the reference and a synthetic
#' single-mutant constructed by substituting the alt allele at the
#' variant position.  ORFs spanning the variant position are extracted
#' from each set; proteins unique to the reference are candidate
#' disruptions, proteins unique to the mutant are candidate novel
#' ORFs.  A lost/gained pair of equal protein length differing by no
#' more than `substitution_tolerance` residues is classified as
#' substitution-only and excluded from both tallies.
#'
#' When `rna_structure = TRUE`, a window of `rna_flank` nt either side of the
#' variant is folded on the chosen strand for both reference and mutant via
#' [fold_rna()] (which requires the ViennaRNA `RNAfold` binary on `PATH`), and
#' the change is quantified as the shift in folding free energy (`delta_mfe`),
#' the base-pair distance (pairs gained + lost between the two structures), and
#' the change in the longest stem.  The window is clamped at the linear genome
#' ends — circularity is not wrapped, so variants within `rna_flank` of the
#' deposit boundary get a truncated structural context.
#'
#' @seealso [scan_orfs()], [compare_mitogenomes()], [fold_rna()],
#'   [pair_partners()], [longest_stem()]
#'
#' @examples
#' \dontrun{
#' # The classic D4a longevity variant (Japanese centenarians)
#' res <- analyze_snp("m.3206C>T")
#' res$summary
#'
#' # m.T16189C — extends a minus-strand D-loop peptide under mt code
#' res16189 <- analyze_snp("m.16189T>C")
#' res16189$detail[["SGC1"]]$novel
#'
#' # RNA structural impact of m.3206C>T (Liu et al. 2007 reported minimal effect)
#' res_rna <- analyze_snp("m.3206C>T", rna_structure = TRUE)
#' res_rna$rna$delta_mfe      # small shift in folding free energy
#' res_rna$rna$bp_distance    # few base pairs gained / lost
#' }
#'
#' @importFrom Biostrings DNAString DNAStringSet reverseComplement
#' @export
analyze_snp <- function(variant,
                        reference              = "NC_012920.1",
                        genetic_codes          = c("SGC1", "SGC0"),
                        min_orf_length         = 60L,
                        start_codons           = c("ATG","GTG","ATA","ATT","ATC"),
                        substitution_tolerance = 3L,
                        rna_structure          = FALSE,
                        rna_flank              = 100L,
                        rna_strand             = c("sense", "antisense"),
                        rna_temperature        = 37,
                        rna_max_bulge          = 2L) {

  rna_strand <- match.arg(rna_strand)

  ## ---- 1. Parse the variant ------------------------------------------------
  parsed <- .parse_variant(variant)
  pos <- parsed$pos
  ref <- parsed$ref
  alt <- parsed$alt
  if (!grepl("^[ACGT]$", ref) || !grepl("^[ACGT]$", alt) || ref == alt)
    stop("Variant ref/alt must be single non-identical A/C/G/T bases.",
         call. = FALSE)

  ## ---- 2. Load reference ---------------------------------------------------
  ref_seq <- .resolve_reference(reference)
  ref_str <- as.character(ref_seq)
  ref_name <- if (is.character(reference)) reference else "user-supplied"

  if (pos < 1L || pos > nchar(ref_str))
    stop(sprintf("Variant position %d outside reference length (%d nt).",
                 pos, nchar(ref_str)), call. = FALSE)
  observed <- substr(ref_str, pos, pos)
  if (toupper(observed) != ref)
    stop(sprintf("Reference base at position %d is %s, not %s as specified.",
                 pos, observed, ref), call. = FALSE)

  ## ---- 3. Build synthetic single-mutant ------------------------------------
  mut_str <- ref_str
  substr(mut_str, pos, pos) <- alt
  mut_seq <- Biostrings::DNAString(mut_str)

  ## ---- 4. Per-code analysis ------------------------------------------------
  detail <- list()
  summary_rows <- list()

  for (code in genetic_codes) {
    ref_ss <- Biostrings::DNAStringSet(ref_seq); names(ref_ss) <- "reference"
    mut_ss <- Biostrings::DNAStringSet(mut_seq)
    names(mut_ss) <- sprintf("m%d%s>%s", pos, ref, alt)

    ref_orfs <- scan_orfs(ref_ss,
                          min_orf_length = min_orf_length,
                          genetic_code   = code,
                          start_codons   = start_codons,
                          both_strands   = TRUE,
                          circular       = TRUE)
    mut_orfs <- scan_orfs(mut_ss,
                          min_orf_length = min_orf_length,
                          genetic_code   = code,
                          start_codons   = start_codons,
                          both_strands   = TRUE,
                          circular       = TRUE)

    ref_at <- ref_orfs[ref_orfs$start <= pos & ref_orfs$end >= pos &
                       !ref_orfs$wraps_around, , drop = FALSE]
    mut_at <- mut_orfs[mut_orfs$start <= pos & mut_orfs$end >= pos &
                       !mut_orfs$wraps_around, , drop = FALSE]

    lost   <- setdiff(ref_at$protein_sequence, mut_at$protein_sequence)
    gained <- setdiff(mut_at$protein_sequence, ref_at$protein_sequence)

    disrupted <- .filter_non_substitution(lost,   gained, substitution_tolerance)
    novel     <- .filter_non_substitution(gained, lost,   substitution_tolerance)

    detail[[code]] <- list(
      ref_orfs  = ref_at,
      alt_orfs  = mut_at,
      lost      = lost,
      gained    = gained,
      disrupted = disrupted,
      novel     = novel
    )
    summary_rows[[code]] <- data.frame(
      genetic_code = code,
      n_ref_orfs   = nrow(ref_at),
      n_alt_orfs   = nrow(mut_at),
      n_lost       = length(lost),
      n_gained     = length(gained),
      n_disrupted  = length(disrupted),
      n_novel      = length(novel),
      stringsAsFactors = FALSE
    )
  }

  ## ---- 5. Optional RNA secondary-structure impact --------------------------
  rna <- NULL
  if (isTRUE(rna_structure)) {
    rna <- .rna_structure_impact(ref_str, mut_str, pos,
                                 flank       = as.integer(rna_flank),
                                 strand      = rna_strand,
                                 temperature = rna_temperature,
                                 max_bulge   = as.integer(rna_max_bulge))
  }

  list(
    variant   = list(pos = pos, ref = ref, alt = alt,
                     label = sprintf("m.%d%s>%s", pos, ref, alt)),
    reference = ref_name,
    summary   = do.call(rbind, summary_rows),
    detail    = detail,
    rna       = rna
  )
}

# ---- helpers ----------------------------------------------------------------

# Parse variant input — accept "m.16189T>C", "16189T>C", "T16189C", or list
.parse_variant <- function(x) {
  if (is.list(x) && !is.null(x$pos) && !is.null(x$ref) && !is.null(x$alt)) {
    return(list(pos = as.integer(x$pos), ref = toupper(x$ref), alt = toupper(x$alt)))
  }
  if (!is.character(x) || length(x) != 1L)
    stop("'variant' must be a character scalar or a list(pos, ref, alt).",
         call. = FALSE)
  s <- sub("^M\\.", "", toupper(x))
  # form 1: NNN<pos>NNN  e.g. T16189C
  m1 <- regmatches(s, regexec("^([ACGT])([0-9]+)([ACGT])$", s))[[1L]]
  if (length(m1) == 4L) {
    return(list(pos = as.integer(m1[3L]), ref = m1[2L], alt = m1[4L]))
  }
  # form 2: <pos>NNN>NNN  e.g. 16189T>C
  m2 <- regmatches(s, regexec("^([0-9]+)([ACGT])>([ACGT])$", s))[[1L]]
  if (length(m2) == 4L) {
    return(list(pos = as.integer(m2[2L]), ref = m2[3L], alt = m2[4L]))
  }
  stop("Could not parse variant '", x, "'. Use forms like 'm.16189T>C', ",
       "'16189T>C', or 'T16189C'.", call. = FALSE)
}

# Resolve the reference argument to a single DNAString
.resolve_reference <- function(reference) {
  if (methods::is(reference, "DNAString")) return(reference)
  if (methods::is(reference, "DNAStringSet")) {
    if (length(reference) != 1L)
      stop("'reference' DNAStringSet must contain exactly one sequence.",
           call. = FALSE)
    return(reference[[1L]])
  }
  if (is.character(reference) && length(reference) == 1L) {
    seqs <- fetch_sequences(accessions = reference)
    return(seqs[[1L]])
  }
  stop("'reference' must be an NCBI accession, DNAString, or DNAStringSet ",
       "of length 1.", call. = FALSE)
}

# For each protein in 'losts', remove any with an equal-length counterpart
# in 'gaineds' differing by <= tol residues (= substitution-only equivalent).
.filter_non_substitution <- function(losts, gaineds, tol) {
  if (length(losts) == 0L) return(character(0L))
  keep <- logical(length(losts))
  for (i in seq_along(losts)) {
    lp <- losts[i]
    matched <- FALSE
    for (gp in gaineds) {
      if (nchar(lp) == nchar(gp)) {
        n_diff <- sum(strsplit(lp, "")[[1L]] != strsplit(gp, "")[[1L]])
        if (n_diff <= tol) { matched <- TRUE; break }
      }
    }
    keep[i] <- !matched
  }
  losts[keep]
}

# Quantify the RNA secondary-structure impact of the substitution.
# Folds a window of `flank` nt either side of the variant on the chosen strand
# for both reference and mutant, then reports the change in folding free energy
# and base-pair topology.  Requires the ViennaRNA RNAfold binary (via fold_rna).
.rna_structure_impact <- function(ref_str, mut_str, pos, flank, strand,
                                  temperature, max_bulge) {
  glen    <- nchar(ref_str)
  w_start <- max(1L, pos - flank)
  w_end   <- min(glen, pos + flank)
  ref_win <- substr(ref_str, w_start, w_end)
  mut_win <- substr(mut_str, w_start, w_end)
  var_off <- pos - w_start + 1L            # variant offset within the window

  if (strand == "antisense") {
    ref_win <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(ref_win)))
    mut_win <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(mut_win)))
    var_off <- nchar(ref_win) - var_off + 1L
  }

  f_ref <- fold_rna(ref_win, temperature = temperature)
  f_mut <- fold_rna(mut_win, temperature = temperature)

  # Represent each structure as the set of "i-j" base pairs (i < j)
  pair_set <- function(structure) {
    p   <- pair_partners(structure)
    idx <- which(!is.na(p) & p > seq_along(p))
    paste(idx, p[idx], sep = "-")
  }
  s_ref  <- pair_set(f_ref$structure)
  s_mut  <- pair_set(f_mut$structure)
  shared <- length(intersect(s_ref, s_mut))
  bp_dist <- length(setdiff(s_ref, s_mut)) + length(setdiff(s_mut, s_ref))
  denom   <- length(s_ref) + length(s_mut)

  ls_ref <- longest_stem(f_ref$structure, max_bulge = max_bulge)
  ls_mut <- longest_stem(f_mut$structure, max_bulge = max_bulge)

  list(
    strand           = strand,
    window           = c(start = w_start, end = w_end,
                         length = w_end - w_start + 1L),
    variant_offset   = var_off,
    mfe_ref          = f_ref$mfe,
    mfe_alt          = f_mut$mfe,
    delta_mfe        = f_mut$mfe - f_ref$mfe,        # >0 = destabilising
    n_pairs_ref      = f_ref$n_pairs,
    n_pairs_alt      = f_mut$n_pairs,
    shared_pairs     = shared,
    bp_distance      = bp_dist,                      # base pairs gained + lost
    bp_distance_frac = if (denom > 0L) bp_dist / denom else 0,  # 0 = identical
    longest_stem_ref = ls_ref$stem_length,
    longest_stem_alt = ls_mut$stem_length,
    delta_stem       = ls_mut$stem_length - ls_ref$stem_length,
    structure_ref    = f_ref$structure,
    structure_alt    = f_mut$structure
  )
}
