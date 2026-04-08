# ============================================================================
# find_conserved_windows — interspecific conserved AA windows via frame translation
# ============================================================================

#' Find Conserved Amino Acid Windows Across Mitochondrial Genomes
#'
#' Compares full six-frame translations of two or more mitochondrial genomes
#' to identify windows of conserved amino acid sequence.  Unlike
#' [cluster_conserved_orfs()], this approach does not require ORF calling:
#' the entire genome is translated as a continuous string under each genetic
#' code.  Each inter-stop peptide segment is then compared pairwise between
#' species, revealing all regions of inter-species homology — canonical ETC
#' genes, rRNA-embedded ORFs, D-loop peptides, and novel micropeptide
#' candidates alike.
#'
#' Canonical protein-coding gene windows will typically appear at the top of
#' the output (highest identity); novel MDP candidates appear in non-CDS
#' regions.  Use the `genomic_region` column to stratify results.
#'
#' @param genomes Named list of genome records as returned by
#'   [fetch_mito_genbank()].  Must contain at least two elements.
#' @param genetic_codes Character vector of genetic-code identifiers passed to
#'   [translate_all_frames()].  Default `c("SGC0", "SGC1")`.
#' @param min_pct_identity Numeric.  Minimum percent identity of the local
#'   alignment to report a window.  Default `40`.
#' @param min_window_aa Integer.  Minimum number of amino acid residues
#'   (non-gap) that must be present in the reported aligned window.
#'   Default `10L`.
#' @param min_segment_aa Integer.  Minimum length of a peptide segment
#'   (inter-stop run) to be included in alignment.  Larger values reduce
#'   runtime by excluding very short segments unlikely to produce informative
#'   alignments.  Must be >= `min_window_aa`.  Default `30L`.
#'   Set to `min_window_aa` for exhaustive search (slower).
#' @param same_strand_only Logical.  If `TRUE` (default), only compare
#'   forward-strand segments of genome A with forward-strand segments of genome
#'   B, and reverse-strand with reverse-strand.  Cross-strand comparisons are
#'   biologically uncommon for conserved ORFs and roughly double runtime.  Set
#'   to `FALSE` to enable exhaustive cross-strand search.
#' @param substitution_matrix Character.  Default `"BLOSUM62"`.
#' @param gap_opening Numeric.  Default `10`.
#' @param gap_extension Numeric.  Default `4`.
#' @param verbose Logical.  Print progress?  Default `TRUE`.
#'
#' @return A `data.frame` with one row per conserved window, sorted by
#'   descending `pct_identity`, with columns:
#' \describe{
#'   \item{`acc_a`, `acc_b`}{Accessions of the two genomes compared.}
#'   \item{`code`}{Genetic code used for translation (e.g. `"SGC0"`).}
#'   \item{`strand_a`, `frame_a`}{Strand (`"+"` / `"-"`) and frame (1–3) of
#'     the window in genome A.}
#'   \item{`strand_b`, `frame_b`}{Strand and frame in genome B.}
#'   \item{`genome_start_a`, `genome_end_a`}{1-based genomic coordinates of
#'     the window in genome A.}
#'   \item{`genome_start_b`, `genome_end_b`}{1-based genomic coordinates in
#'     genome B.}
#'   \item{`window_aa_a`}{Amino acid sequence of the window in genome A
#'     (no gaps).}
#'   \item{`window_aa_b`}{Amino acid sequence of the window in genome B.}
#'   \item{`window_length_aa`}{Number of aligned amino acid columns
#'     (including gaps).}
#'   \item{`pct_identity`}{Percent identity over the aligned window.}
#'   \item{`alignment_score`}{Raw Smith-Waterman score.}
#'   \item{`genomic_region_a`, `genomic_region_b`}{Genomic context from the
#'     GenBank feature table (if features are present in the genome records).}
#' }
#'
#' @examples
#' \dontrun{
#' genomes <- fetch_mito_genbank(c("NC_002333.2", "NC_011379.1"))
#' windows <- find_conserved_windows(genomes, min_pct_identity = 50)
#'
#' # Novel candidates: not in canonical CDS
#' novel <- windows[!grepl("^CDS:", windows$genomic_region_a), ]
#'
#' # Exhaustive search including very short segments (slower)
#' windows_all <- find_conserved_windows(genomes, min_segment_aa = 10L)
#' }
#'
#' @importFrom pwalign pairwiseAlignment pid score alignedPattern alignedSubject
#'   subject pattern
#' @importFrom Biostrings AAString
#' @importFrom BiocGenerics start end
#' @importFrom methods is
#' @export
find_conserved_windows <- function(genomes,
                                    genetic_codes       = c("SGC0", "SGC1"),
                                    min_pct_identity    = 40,
                                    min_window_aa       = 10L,
                                    min_segment_aa      = 30L,
                                    same_strand_only    = TRUE,
                                    substitution_matrix = "BLOSUM62",
                                    gap_opening         = 10,
                                    gap_extension       = 4,
                                    verbose             = TRUE) {

  ## ---- input checks --------------------------------------------------------
  if (!is.list(genomes) || length(genomes) < 2L)
    stop("'genomes' must be a named list of at least two genome records.",
         call. = FALSE)

  accessions <- names(genomes)
  if (is.null(accessions) || any(!nzchar(accessions)))
    stop("'genomes' must be a *named* list (names = accession strings).",
         call. = FALSE)

  min_segment_aa <- max(min_segment_aa, min_window_aa)

  ## ---- translate all genomes -----------------------------------------------
  if (verbose) message("Translating genomes in all frames...")

  translations <- lapply(genomes, function(rec) {
    translate_all_frames(rec, genetic_codes = genetic_codes)
  })
  names(translations) <- accessions

  ## ---- split each frame string into peptide segments -----------------------
  if (verbose) message("Splitting translations into peptide segments...")

  segments <- lapply(translations, function(tr) {
    lapply(tr, function(aa_str) {
      .split_by_stops(as.character(aa_str), min_segment_aa)
    })
  })

  ## ---- all pairwise genome comparisons ------------------------------------
  acc_pairs   <- utils::combn(accessions, 2L, simplify = FALSE)
  result_rows <- list()
  row_k       <- 0L

  for (pair in acc_pairs) {
    acc_a <- pair[1L];  acc_b <- pair[2L]
    seg_a <- segments[[acc_a]]
    seg_b <- segments[[acc_b]]

    ## -- features for region annotation --------------------------------------
    feat_a <- genomes[[acc_a]]$features
    feat_b <- genomes[[acc_b]]$features
    len_a  <- genomes[[acc_a]]$length
    len_b  <- genomes[[acc_b]]$length

    n_seg_a <- sum(vapply(seg_a, length, integer(1L)))
    n_seg_b <- sum(vapply(seg_b, length, integer(1L)))

    if (verbose)
      message(sprintf("  Comparing %s vs %s (%d x %d segments)...",
                      acc_a, acc_b, n_seg_a, n_seg_b))

    for (key_a in names(seg_a)) {
      segs_a <- seg_a[[key_a]]
      if (length(segs_a) == 0L) next

      meta_a <- .parse_frame_key(key_a, len_a)
      code_a <- .extract_code(key_a)

      for (key_b in names(seg_b)) {
        # Only compare matching genetic codes
        if (.extract_code(key_b) != code_a) next
        # Strand filter: skip cross-strand comparisons when same_strand_only=TRUE
        if (same_strand_only && .extract_strand(key_b) != meta_a$strand) next

        segs_b <- seg_b[[key_b]]
        if (length(segs_b) == 0L) next

        meta_b <- .parse_frame_key(key_b, len_b)

        # Precompute k-mers for all segments in B (reused across seg_a loop)
        km_b_list <- lapply(segs_b, function(sg) .kmers4(sg$seq))

        for (sg_a in segs_a) {
          # Precompute k-mers for this segment (k=4)
          km_a <- .kmers4(sg_a$seq)

          for (ib in seq_along(segs_b)) {
            sg_b  <- segs_b[[ib]]
            km_b  <- km_b_list[[ib]]
            # k-mer pre-filter: skip if no shared 4-mers (fast rejection)
            if (length(km_a) > 0L && length(km_b) > 0L &&
                !any(km_a %in% km_b)) next

            aln <- pwalign::pairwiseAlignment(
              pattern            = sg_b$seq,
              subject            = sg_a$seq,
              type               = "local",
              substitutionMatrix = substitution_matrix,
              gapOpening         = gap_opening,
              gapExtension       = gap_extension
            )

            pct_id <- pwalign::pid(aln, type = "PID1")
            if (pct_id < min_pct_identity) next

            pat_str <- as.character(pwalign::alignedPattern(aln))
            sub_str <- as.character(pwalign::alignedSubject(aln))
            n_aligned_a <- nchar(gsub("-", "", sub_str))
            n_aligned_b <- nchar(gsub("-", "", pat_str))
            if (min(n_aligned_a, n_aligned_b) < min_window_aa) next

            aln_score  <- pwalign::score(aln)
            window_len <- nchar(sub_str)   # alignment columns incl. gaps

            ## -- map AA positions back to genomic coordinates ----------------
            sub_ranges <- pwalign::subject(aln)
            pat_ranges <- pwalign::pattern(aln)

            aa_start_a <- sg_a$aa_start + as.integer(BiocGenerics::start(sub_ranges)) - 1L
            aa_end_a   <- sg_a$aa_start + as.integer(BiocGenerics::end(sub_ranges))   - 1L
            aa_start_b <- sg_b$aa_start + as.integer(BiocGenerics::start(pat_ranges)) - 1L
            aa_end_b   <- sg_b$aa_start + as.integer(BiocGenerics::end(pat_ranges))   - 1L

            gstart_a <- .aa_to_genome(aa_start_a, meta_a, len_a)
            gend_a   <- .aa_to_genome(aa_end_a,   meta_a, len_a)
            gstart_b <- .aa_to_genome(aa_start_b, meta_b, len_b)
            gend_b   <- .aa_to_genome(aa_end_b,   meta_b, len_b)

            ## -- genomic region annotation -----------------------------------
            reg_a <- .region_for_pos(gstart_a, feat_a)
            reg_b <- .region_for_pos(gstart_b, feat_b)

            row_k <- row_k + 1L
            result_rows[[row_k]] <- data.frame(
              acc_a            = acc_a,
              acc_b            = acc_b,
              code             = code_a,
              strand_a         = meta_a$strand,
              frame_a          = meta_a$frame,
              strand_b         = meta_b$strand,
              frame_b          = meta_b$frame,
              genome_start_a   = gstart_a,
              genome_end_a     = gend_a,
              genome_start_b   = gstart_b,
              genome_end_b     = gend_b,
              window_aa_a      = gsub("-", "", sub_str),
              window_aa_b      = gsub("-", "", pat_str),
              window_length_aa = window_len,
              pct_identity     = round(pct_id, 1L),
              alignment_score  = round(aln_score, 1L),
              genomic_region_a = reg_a,
              genomic_region_b = reg_b,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }

  ## ---- assemble output -----------------------------------------------------
  if (row_k == 0L) {
    message("No conserved windows found at the specified thresholds.")
    return(invisible(.empty_windows()))
  }

  out <- do.call(rbind, result_rows)
  out <- out[order(-out$pct_identity, -out$alignment_score), ]
  rownames(out) <- NULL

  if (verbose)
    message(sprintf("Found %d conserved window(s).", nrow(out)))

  out
}


# --------------------------------------------------------------------------- #
# Internals                                                                    #
# --------------------------------------------------------------------------- #

#' Split a translated string into inter-stop peptide segments
#'
#' @param aa_char Character string (translated genome, X = stop).
#' @param min_len Integer minimum segment length.
#' @return List of lists with: aa_start (1-based), aa_end, seq (character).
.split_by_stops <- function(aa_char, min_len) {
  if (!nzchar(aa_char)) return(list())

  x_pos  <- c(0L, which(strsplit(aa_char, "", fixed = TRUE)[[1L]] == "X"),
               nchar(aa_char) + 1L)

  out <- list()
  for (i in seq_len(length(x_pos) - 1L)) {
    seg_start <- x_pos[i] + 1L
    seg_end   <- x_pos[i + 1L] - 1L
    seg_len   <- seg_end - seg_start + 1L
    if (seg_len < min_len) next
    out[[length(out) + 1L]] <- list(
      aa_start = seg_start,
      aa_end   = seg_end,
      seq      = substr(aa_char, seg_start, seg_end)
    )
  }
  out
}

#' Compute all 4-mers of an amino acid string
.kmers4 <- function(seq, k = 4L) {
  n <- nchar(seq)
  if (n < k) return(character(0L))
  vapply(seq_len(n - k + 1L),
         function(i) substr(seq, i, i + k - 1L),
         character(1L))
}

.extract_code <- function(key) {
  parts <- strsplit(key, "_", fixed = TRUE)[[1L]]
  parts[grepl("^SGC", parts)]
}

.extract_strand <- function(key) {
  # Last token is strand+frame e.g. "+1", "-3"
  substr(regmatches(key, regexpr("[+-][123]$", key)), 1L, 1L)
}

.parse_frame_key <- function(key, genome_len) {
  sf    <- regmatches(key, regexpr("[+-][123]$", key))
  strand <- substr(sf, 1L, 1L)
  frame  <- as.integer(substr(sf, 2L, 2L))
  list(strand = strand, frame = frame, genome_len = genome_len)
}

.aa_to_genome <- function(aa_pos, meta, genome_len) {
  nt_in_chunk <- (aa_pos - 1L) * 3L + 1L

  if (meta$strand == "+") {
    genomic <- meta$frame + nt_in_chunk - 1L
    min(genomic, genome_len)
  } else {
    rc_pos  <- meta$frame + nt_in_chunk - 1L
    genomic <- genome_len - rc_pos + 1L
    max(1L, genomic)
  }
}

.region_for_pos <- function(pos, features) {
  if (is.null(features) || nrow(features) == 0L) return("unknown")

  priority <- function(type) {
    switch(type, CDS = 1L, tRNA = 2L, rRNA = 3L,
           `D-loop` = 4L, `rep_origin` = 4L, 5L)
  }

  hits <- which(features$start <= pos & features$end >= pos &
                  features$type != "source")
  if (length(hits) == 0L) return("intergenic")

  prios  <- vapply(hits, function(h) priority(features$type[h]), integer(1L))
  best   <- hits[which.min(prios)]
  f      <- features[best, ]

  switch(f$type,
    CDS  = paste0("CDS:", if (!is.na(f$gene)) f$gene else "unknown"),
    rRNA = {
      lbl <- paste(if (!is.na(f$gene)) f$gene else "",
                   if (!is.na(f$product)) f$product else "")
      if (grepl("12S|rnr1|small", lbl, ignore.case = TRUE)) "12S_rRNA"
      else if (grepl("16S|rnr2|large", lbl, ignore.case = TRUE)) "16S_rRNA"
      else trimws(lbl)
    },
    tRNA        = paste0("tRNA:", if (!is.na(f$product)) f$product else f$gene),
    `D-loop`    = "D-loop",
    `rep_origin`= "D-loop",
    if (!is.na(f$gene)) f$gene else f$type
  )
}

.empty_windows <- function() {
  data.frame(
    acc_a = character(), acc_b = character(), code = character(),
    strand_a = character(), frame_a = integer(),
    strand_b = character(), frame_b = integer(),
    genome_start_a = integer(), genome_end_a = integer(),
    genome_start_b = integer(), genome_end_b = integer(),
    window_aa_a = character(), window_aa_b = character(),
    window_length_aa = integer(),
    pct_identity = numeric(), alignment_score = numeric(),
    genomic_region_a = character(), genomic_region_b = character(),
    stringsAsFactors = FALSE
  )
}
