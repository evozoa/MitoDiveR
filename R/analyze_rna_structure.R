# ============================================================================
# analyze_rna_structure — post-processing of dot-bracket structures from
# fold_rna() to identify stem-loops, megahairpins, and other diagnostic
# features (e.g., asmtRNA-class transcripts)
# ============================================================================

#' Map base-pair partners from a dot-bracket structure
#'
#' Walks a dot-bracket secondary-structure string (as produced by
#' [fold_rna()]) and returns, for each position, the index of its
#' base-pair partner — or `NA` if the position is unpaired.
#'
#' @param structure A character string of dot-bracket notation
#'   (`(`, `)`, `.`).  Must be balanced.  May come from `fold_rna()$structure`
#'   or any other RNA-folding output.
#'
#' @return An integer vector of length `nchar(structure)`. Position `i` holds
#'   the partner position of base `i`, or `NA_integer_` if base `i` is
#'   unpaired (`.`) or pseudoknotted to a partner outside this notation.
#'
#' @examples
#' # Simple stem-loop: ((((....))))
#' partners <- pair_partners("((((....))))")
#' partners       # 12 11 10  9 NA NA NA NA  4  3  2  1
#'
#' \dontrun{
#' # Pair with fold_rna()
#' f <- fold_rna("GGGAAACCC")
#' pair_partners(f$structure)
#' }
#'
#' @seealso [fold_rna()], [longest_stem()]
#' @export
pair_partners <- function(structure) {
  if (!is.character(structure) || length(structure) != 1L)
    stop("'structure' must be a single character string.", call. = FALSE)
  ch <- strsplit(structure, "")[[1L]]
  L  <- length(ch)
  partners <- rep(NA_integer_, L)
  stack <- integer(0)
  for (i in seq_len(L)) {
    c_i <- ch[i]
    if (c_i == "(") {
      stack <- c(stack, i)
    } else if (c_i == ")") {
      if (length(stack) == 0L)
        stop("Unbalanced dot-bracket structure at position ", i, ".",
             call. = FALSE)
      j <- stack[length(stack)]
      stack <- stack[-length(stack)]
      partners[i] <- j
      partners[j] <- i
    }
    # all other characters (".", whitespace, IUPAC pseudoknot codes) → unpaired
  }
  if (length(stack) > 0L)
    stop("Unbalanced dot-bracket structure: ", length(stack),
         " unmatched '('.", call. = FALSE)
  partners
}

#' Find the longest contiguous base-pair stem in a structure
#'
#' Identifies the longest consecutive base-pair "ladder" in a dot-bracket
#' structure — i.e., the longest run of positions where `i, i+1, i+2, ...`
#' on the 5' side pair with `j, j-1, j-2, ...` on the 3' side, defining a
#' true stem rather than aggregate base-pair counts.
#'
#' Optionally tolerates small bulges (single-strand insertions of up to
#' `max_bulge` residues on either strand) so that biologically realistic
#' stems with internal loops are still captured as a single stem.  This is
#' the recommended setting for detecting long inverted-repeat-driven
#' hairpins such as those found in mammalian ASncmtRNAs (Burzio lab; Gao
#' et al. 2017), where the diagnostic feature is a megahairpin of >=40
#' paired residues, typically interrupted by a few small bulges.
#'
#' @param structure Dot-bracket structure string (single character; see
#'   [fold_rna()]).
#' @param max_bulge Integer >= 0.  Maximum number of unpaired residues
#'   allowed on either the 5' or 3' strand between two consecutive base
#'   pairs that are still considered part of the same stem.  Default `0`
#'   (strict — no bulges, only perfect ladders).  Set to `2`–`3` for
#'   biologically realistic stem detection.
#'
#' @return A list with elements:
#' \describe{
#'   \item{`stem_length`}{Integer.  Number of base pairs in the longest
#'     stem (counting through allowed bulges).}
#'   \item{`stem_5p_pos`}{Integer.  1-based 5' position of the first base
#'     pair of the stem, or `NA` if no stem found.}
#'   \item{`stem_3p_pos`}{Integer.  1-based 3' position of the partner of
#'     the first base pair (i.e., the outermost 3' end of the stem).}
#'   \item{`stem_span_nt`}{Integer.  Number of nucleotides between
#'     `stem_5p_pos` and `stem_3p_pos` inclusive (the "footprint" of the
#'     stem on the sequence).}
#' }
#'
#' @details
#' The strict mode (`max_bulge = 0`) requires perfect consecutive pairing.
#' This is overly restrictive for biological RNAs because real stems
#' typically contain small bulges.  For most asmtRNA / lncRNA hairpin
#' detection use `max_bulge = 2L` or `3L`.
#'
#' Note that the longest stem may not be the most thermodynamically stable
#' or biologically functional one — for those questions consult
#' `fold_rna()$mfe` and the full base-pair table.  This function is
#' designed for *topology* inference (presence/absence of long
#' inverted-repeat-driven megahairpins), not energetics.
#'
#' @examples
#' # Strict stem detection — perfect ladder
#' s1 <- longest_stem("(((((....)))))")
#' s1$stem_length    # 5
#'
#' # Two stems separated by an unpaired junction:
#' #  ((((..(((....))))..))))
#' s2 <- longest_stem("((((..(((....))))..))))")
#' s2$stem_length    # 4 (the outer stem)
#'
#' # Permissive mode handles small bulges in a single conceptual stem:
#' #  (((((.((....)).)))))
#' s3 <- longest_stem("(((((.((....)).)))))", max_bulge = 1L)
#' s3$stem_length    # 7 (bulged stem extends through the .((..)).)
#'
#' \dontrun{
#' # Test for asmtRNA-class megahairpin in a folded transcript
#' f <- fold_rna(my_lstrand_rna_seq)
#' s <- longest_stem(f$structure, max_bulge = 3L)
#' is_megahairpin <- s$stem_length >= 40L
#' }
#'
#' @seealso [fold_rna()], [pair_partners()]
#' @export
longest_stem <- function(structure, max_bulge = 0L) {
  if (!is.character(structure) || length(structure) != 1L)
    stop("'structure' must be a single character string.", call. = FALSE)
  if (!is.numeric(max_bulge) || length(max_bulge) != 1L || max_bulge < 0)
    stop("'max_bulge' must be a non-negative integer.", call. = FALSE)
  max_bulge <- as.integer(max_bulge)

  partners <- pair_partners(structure)
  L <- length(partners)
  if (max_bulge == 0L) {
    # strict perfect ladder
    best_len <- 0L; best_5p <- NA_integer_; best_3p <- NA_integer_
    i <- 1L
    while (i <= L) {
      if (!is.na(partners[i]) && partners[i] > i) {
        run <- 1L
        while (i + run <= L &&
               !is.na(partners[i + run]) &&
               partners[i + run] == partners[i] - run) {
          run <- run + 1L
        }
        if (run > best_len) {
          best_len <- run; best_5p <- i; best_3p <- partners[i]
        }
        i <- i + run
      } else {
        i <- i + 1L
      }
    }
  } else {
    # permissive: walk along the structure, allowing up to max_bulge unpaired
    # residues per strand between consecutive paired residues counted as part
    # of the same stem
    visited <- logical(L)
    best_len <- 0L; best_5p <- NA_integer_; best_3p <- NA_integer_
    for (start in seq_len(L)) {
      if (visited[start]) next
      if (is.na(partners[start]) || partners[start] <= start) next
      count   <- 1L
      last_5p <- start
      last_3p <- partners[start]
      visited[start] <- TRUE
      pos <- start + 1L
      while (pos <= L && pos < last_3p) {
        if (visited[pos]) { pos <- pos + 1L; next }
        if (!is.na(partners[pos]) &&
            partners[pos] < last_3p &&
            partners[pos] > pos &&
            (pos - last_5p) <= (max_bulge + 1L) &&
            (last_3p - partners[pos]) <= (max_bulge + 1L)) {
          count    <- count + 1L
          last_5p  <- pos
          last_3p  <- partners[pos]
          visited[pos] <- TRUE
        }
        pos <- pos + 1L
      }
      if (count > best_len) {
        best_len <- count; best_5p <- start; best_3p <- partners[start]
      }
    }
  }

  list(
    stem_length  = best_len,
    stem_5p_pos  = best_5p,
    stem_3p_pos  = best_3p,
    stem_span_nt = if (!is.na(best_3p)) best_3p - best_5p + 1L else NA_integer_
  )
}

#' Construct a circular L-strand transcript spanning the mt control region
#'
#' Helper for testing antisense mitochondrial RNA (asmtRNA) hypotheses.
#' Mitochondrial genomes are circular, so an L-strand transcript can
#' originate in the control region (D-loop) and extend continuously across
#' the linear-deposit "wrap" boundary into MT-RNR1 / MT-RNR2.  Standard
#' linear extraction misses this — and the inverted-repeat element that
#' drives the mammalian ASncmtRNA-1 megahairpin sits exactly at the
#' D-loop / 12S boundary.
#'
#' This function takes a mitogenome H-strand sequence and constructs the
#' L-strand 5'→3' of the transcript by:
#' \enumerate{
#'   \item Concatenating the last `dloop_upstream_nt` of the linear
#'     deposit (= the D-loop / control region tail) onto the start of the
#'     genome through `rrna_end`.
#'   \item Reverse-complementing the result.
#' }
#'
#' @param sequence Mitogenome H-strand as a [Biostrings::DNAString] or a
#'   character string.  Assumed to be the deposited rRNA-coding strand
#'   (standard mtDNA convention).
#' @param rrna_end Integer.  Last position of the rRNA region on the
#'   linear deposit (typically the end of MT-RNR2 / 16S).
#' @param dloop_upstream_nt Integer.  Number of nucleotides at the
#'   genome's 3' end (= D-loop / control region) to prepend before the
#'   rRNA-region extraction.  Default `600L`.  Use ~600–1000 to capture
#'   the inverted-repeat element responsible for the mammalian
#'   ASncmtRNA-1 megahairpin.
#'
#' @return A [Biostrings::DNAString] containing the L-strand transcript
#'   in 5'→3' orientation, length =
#'   `dloop_upstream_nt + rrna_end`.
#'
#' @examples
#' \dontrun{
#' # Test the human ASncmtRNA-1 megahairpin
#' rec <- fetch_mito_genbank("NC_012920.1")[[1]]
#' lstrand <- circular_lstrand_transcript(
#'   rec$sequence,
#'   rrna_end          = 3229L,   # end of MT-RNR2
#'   dloop_upstream_nt = 600L
#' )
#' f <- fold_rna(lstrand)
#' s <- longest_stem(f$structure, max_bulge = 3L)
#' s$stem_length    # ~42 paired residues for human ASncmtRNA-1
#' }
#'
#' @seealso [fold_rna()], [longest_stem()]
#' @importFrom Biostrings DNAString reverseComplement subseq
#' @export
circular_lstrand_transcript <- function(sequence,
                                         rrna_end,
                                         dloop_upstream_nt = 600L) {
  if (inherits(sequence, "DNAString")) {
    s <- sequence
  } else if (is.character(sequence) && length(sequence) == 1L) {
    s <- Biostrings::DNAString(sequence)
  } else {
    stop("'sequence' must be a DNAString or single character string.",
         call. = FALSE)
  }
  glen <- length(s)
  if (!is.numeric(rrna_end) || length(rrna_end) != 1L ||
      rrna_end < 1 || rrna_end > glen)
    stop("'rrna_end' must be an integer in [1, length(sequence)].",
         call. = FALSE)
  if (!is.numeric(dloop_upstream_nt) || length(dloop_upstream_nt) != 1L ||
      dloop_upstream_nt < 0)
    stop("'dloop_upstream_nt' must be a non-negative integer.",
         call. = FALSE)
  rrna_end          <- as.integer(rrna_end)
  dloop_upstream_nt <- as.integer(dloop_upstream_nt)

  back_start <- max(1L, glen - dloop_upstream_nt + 1L)
  back_part  <- as.character(Biostrings::subseq(s, back_start, glen))
  front_part <- as.character(Biostrings::subseq(s, 1L, rrna_end))
  full_h     <- Biostrings::DNAString(paste0(back_part, front_part))
  Biostrings::reverseComplement(full_h)
}
