# ============================================================================
# annotate_reference_frame — anchor ORF reading frames to a shared coordinate
# system, so "same frame" is comparable across genomes.
# ============================================================================

# Build a query-position -> anchor-position lookup (integer vector, one entry
# per base of the query, NA where the query base aligns before the anchor's
# first base). baseOnly = FALSE because real deposits carry IUPAC ambiguity
# codes (R, Y, ...) that a baseOnly matrix cannot score.
.build_anchor_pos_map <- function(query_str, anchor_str, gap_opening, gap_extension) {
  aln <- pwalign::pairwiseAlignment(
    query_str, anchor_str, type = "global",
    substitutionMatrix = pwalign::nucleotideSubstitutionMatrix(
      match = 1, mismatch = -1, baseOnly = FALSE
    ),
    gapOpening = gap_opening, gapExtension = gap_extension
  )
  patt <- strsplit(as.character(pwalign::alignedPattern(aln)), "")[[1L]]
  subj <- strsplit(as.character(pwalign::alignedSubject(aln)), "")[[1L]]
  subj_pos <- cumsum(subj != "-")
  map <- subj_pos[patt != "-"]        # one entry per query position, in order
  map[map == 0L] <- NA_integer_       # query base precedes the anchor's start
  map
}

#' Anchor ORF Reading Frames to a Shared Reference Coordinate System
#'
#' Reading `frame` in [find_orfs()] / [scan_orfs()] output is `start position
#' mod 3` in each genome's own coordinates, so it drifts across genomes with
#' every indel accumulated since their common ancestor — two truly orthologous
#' ORFs routinely end up with different frame numbers, and two unrelated
#' overlapping ORFs can just as easily share one. Frame numbers are therefore
#' not comparable across sequences on their own.
#'
#' This function fixes that by aligning every sequence to a single shared
#' anchor sequence once, mapping each ORF's start-codon (5') position through
#' that alignment, and labelling it with `anchor_position mod 3` (combined
#' with strand). Two ORFs from different genomes with the same
#' `<prefix>_frame` value sit in the same phase of the same anchor-relative
#' reading frame — a cheap, alignment-verified fact that a raw coordinate
#' overlap cannot confirm on its own. It is a filter for candidate positional
#' orthologs, not a homology call: overlapping, same-frame ORFs still need an
#' amino-acid comparison (e.g. [find_conserved_orfs()], or a direct
#' [pwalign::pairwiseAlignment()] on `protein_sequence`) to confirm they are
#' actually the same peptide, since two unrelated ORFs can coincidentally
#' share a frame.
#'
#' @param orfs `data.frame` from [scan_orfs()] / [find_orfs()] /
#'   [collapse_nested_orfs()]. Must contain `seq_id`, `strand`, `start`, and
#'   `end` columns. ORFs from more than one genetic-code scan of the same
#'   genomes can be combined (e.g. `rbind()`ed) before calling — the alignment
#'   step is independent of genetic code and translation, so it only needs to
#'   run once either way.
#' @param genomes Named [Biostrings::DNAStringSet] or named character vector
#'   of the genome sequences `orfs$seq_id` refers to — the same object passed
#'   to [scan_orfs()]. Must include every value in `orfs$seq_id`, including
#'   `anchor`.
#' @param anchor Character scalar naming the sequence in `genomes` to use as
#'   the shared coordinate system. Its choice is arbitrary — the anchor is a
#'   ruler, not a biological reference — but it must be the same anchor for
#'   every call being compared. Defaults to `names(genomes)[1]`, with a
#'   message identifying the sequence used.
#' @param gap_opening,gap_extension Numeric. Gap penalties passed to
#'   [pwalign::pairwiseAlignment()] for the nucleotide alignment of each
#'   non-anchor genome against `anchor`. Defaults `5` and `2`.
#' @param prefix Character scalar prepended to the names of the added columns.
#'   Default `"ref"`.
#' @param verbose Logical. Report progress as each genome is aligned to the
#'   anchor? Default `TRUE`.
#'
#' @return `orfs` with four columns added (named using `prefix`):
#' \describe{
#'   \item{`ref_5p`}{Anchor-coordinate position of the ORF's start codon (the
#'     `start` column on the `+` strand, `end` on the `-` strand, per
#'     [collapse_nested_orfs()]'s stop/start convention).}
#'   \item{`ref_start`, `ref_end`}{Anchor-coordinate span of the ORF
#'     (`start`/`end` mapped through the alignment; always `ref_start <=
#'     ref_end`).}
#'   \item{`ref_frame`}{`"<strand>_<(ref_5p - 1) %% 3>"`, e.g. `"-_1"`. Two
#'     ORFs sharing this value are in the same anchor-relative frame.}
#' }
#' All four are `NA` for an ORF whose start codon aligns before the anchor's
#' first base (no anchor coordinate exists there). Row order and all existing
#' columns are preserved.
#'
#' @details
#' One global nucleotide alignment is computed per non-anchor sequence in
#' `genomes` (the anchor maps to itself with no alignment needed), not one per
#' ORF or per pair of ORFs — cost scales with the number of genomes, not the
#' number of ORFs or ORF pairs. For `n` genomes that is `n - 1` alignments
#' regardless of how many times "same frame?" is subsequently asked, which is
#' the point: it replaces a pairwise, quadratic-in-genomes alignment burden
#' with a linear one.
#'
#' @examples
#' \dontrun{
#' orfs <- scan_orfs(seqs, genetic_code = "SGC0", min_orf_length = 60,
#'                    max_orf_length = 900, start_codons = "ATG")
#' orfs <- collapse_nested_orfs(orfs)
#' orfs <- annotate_reference_frame(orfs, genomes = seqs,
#'                                   anchor = names(seqs)[1])
#'
#' # candidate positional orthologs: same strand, same anchor frame, overlapping
#' same_frame <- split(orfs, orfs$ref_frame)
#' }
#'
#' @seealso [scan_orfs()], [find_orfs()], [collapse_nested_orfs()],
#'   [find_conserved_orfs()]
#' @export
annotate_reference_frame <- function(orfs,
                                      genomes,
                                      anchor         = NULL,
                                      gap_opening    = 5,
                                      gap_extension  = 2,
                                      prefix         = "ref",
                                      verbose        = TRUE) {

  required_cols <- c("seq_id", "strand", "start", "end")
  missing_cols  <- setdiff(required_cols, names(orfs))
  if (length(missing_cols) > 0L)
    stop("'orfs' is missing required column(s): ", paste(missing_cols, collapse = ", "),
         call. = FALSE)

  if (!methods::is(genomes, "DNAStringSet"))
    genomes <- Biostrings::DNAStringSet(genomes)
  if (is.null(names(genomes)) || any(!nzchar(names(genomes))))
    stop("'genomes' must be named.", call. = FALSE)

  missing_seqs <- setdiff(unique(orfs$seq_id), names(genomes))
  if (length(missing_seqs) > 0L)
    stop("'genomes' is missing sequence(s) referenced in orfs$seq_id: ",
         paste(missing_seqs, collapse = ", "), call. = FALSE)

  if (is.null(anchor)) {
    anchor <- names(genomes)[1L]
    if (verbose) message("No anchor specified -- using '", anchor, "' as the reference frame anchor.")
  } else if (!anchor %in% names(genomes)) {
    stop("'anchor' ('", anchor, "') is not a name in 'genomes'.", call. = FALSE)
  }

  genome_len <- stats::setNames(BiocGenerics::width(genomes), names(genomes))
  anchor_str <- as.character(genomes[[anchor]])

  pos_maps <- list()
  pos_maps[[anchor]] <- seq_len(genome_len[[anchor]])
  for (nm in setdiff(unique(orfs$seq_id), anchor)) {
    if (verbose) message("Aligning '", nm, "' to anchor '", anchor, "'...")
    pos_maps[[nm]] <- .build_anchor_pos_map(as.character(genomes[[nm]]), anchor_str,
                                             gap_opening, gap_extension)
  }

  to_ref <- function(seq_id, pos) {
    pos <- pmin(pmax(pos, 1L), genome_len[[seq_id]])
    pos_maps[[seq_id]][pos]
  }

  five_prime <- ifelse(orfs$strand == "+", orfs$start, orfs$end)
  ref_5p    <- mapply(to_ref, orfs$seq_id, five_prime)
  ref_start <- mapply(to_ref, orfs$seq_id, orfs$start)
  ref_end   <- mapply(to_ref, orfs$seq_id, orfs$end)

  orfs[[paste0(prefix, "_5p")]]    <- unname(ref_5p)
  orfs[[paste0(prefix, "_start")]] <- unname(pmin(ref_start, ref_end))
  orfs[[paste0(prefix, "_end")]]   <- unname(pmax(ref_start, ref_end))
  orfs[[paste0(prefix, "_frame")]] <- ifelse(
    is.na(ref_5p), NA_character_,
    paste0(orfs$strand, "_", (unname(ref_5p) - 1L) %% 3L)
  )

  orfs
}
