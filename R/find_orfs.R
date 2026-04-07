# ============================================================================
# find_orfs — ORF detection across all six reading frames
# ============================================================================

# (Input normalisation is handled by .resolve_dna_input() in utils.R)


# --------------------------------------------------------------------------- #
# Main exported function                                                        #
# --------------------------------------------------------------------------- #

#' Find Open Reading Frames in a (Mitochondrial) Genome
#'
#' Searches all six reading frames (three forward, three reverse-complement)
#' of a nucleotide sequence and returns every ORF — defined as a run from a
#' recognised start codon to the first in-frame stop codon — that meets the
#' minimum length threshold.
#'
#' @param sequence Nucleotide sequence to search. Accepts any of:
#'   - A raw nucleotide `character` string (case-insensitive)
#'   - A path to a FASTA file (single- or multi-sequence)
#'   - An NCBI accession number (e.g. `"NC_002333.2"`) — sequence is
#'     fetched automatically via the NCBI API
#'   - A [Biostrings::DNAString] or [Biostrings::DNAStringSet]
#'
#'   When the input resolves to multiple sequences (multi-record FASTA or
#'   multi-sequence DNAStringSet), ORFs are searched in each sequence and
#'   results are combined into a single `data.frame` with a `seq_id` column
#'   identifying the source sequence.
#' @param genetic_code Character scalar. Biostrings genetic-code identifier
#'   passed to [Biostrings::getGeneticCode()]. Defaults to `"SGC0"`
#'   (standard vertebrate nuclear code, NCBI table 1). Use `"SGC1"` for the
#'   vertebrate mitochondrial code (TGA = Trp; AGA/AGG = stop).
#' @param start_codons Character vector of codons that initiate an ORF.
#'   Defaults to `"ATG"`. For vertebrate mitochondrial genomes you may wish
#'   to add alternative starts: `c("ATG", "GTG", "ATT", "ATC", "ATA")`.
#' @param min_orf_length Integer. Minimum ORF length in nucleotides,
#'   measured from the first base of the start codon to the last base of
#'   the stop codon (inclusive when `include_stop_codon = TRUE`).
#'   Default `100L`.
#' @param max_orf_length Integer. Maximum ORF length in nucleotides
#'   (same measurement as `min_orf_length`). ORFs longer than this value
#'   are silently excluded. Default `Inf` (no upper limit).
#' @param both_strands Logical. Search the reverse-complement strand as well
#'   as the forward strand? Default `TRUE`.
#' @param circular Logical. Treat the sequence as a circular molecule and
#'   detect ORFs that span the linearisation boundary? Default `TRUE`.
#' @param include_stop_codon Logical. Include the stop codon in the reported
#'   `orf_sequence` and `length_nt`? Default `TRUE`. When `FALSE`, the stop
#'   codon column is still reported but the sequence/length exclude it.
#'
#' @return A `data.frame` sorted by `start` then `strand`, with columns:
#' \describe{
#'   \item{`strand`}{`"+"` or `"-"`.}
#'   \item{`frame`}{Reading frame: `+1`, `+2`, `+3`, `-1`, `-2`, or `-3`.
#'     On the minus strand, frame `-1` starts from the last nucleotide of
#'     the genome, `-2` from the second-to-last, etc.}
#'   \item{`start`}{1-based genomic start coordinate (always `<=` `end`
#'     unless `wraps_around` is `TRUE`).}
#'   \item{`end`}{1-based genomic end coordinate (inclusive).}
#'   \item{`length_nt`}{ORF length in nucleotides.}
#'   \item{`length_aa`}{Number of residues in the translated sequence
#'     (including `"*"` for the stop codon when `include_stop_codon =
#'     TRUE`).}
#'   \item{`start_codon`}{Initiating codon.}
#'   \item{`stop_codon`}{Terminating codon (`NA` when
#'     `include_stop_codon = FALSE`).}
#'   \item{`wraps_around`}{`TRUE` when the ORF spans the origin of a
#'     circular genome (only possible when `circular = TRUE`).}
#'   \item{`orf_sequence`}{Nucleotide sequence of the ORF on the sense
#'     strand.}
#'   \item{`protein_sequence`}{Amino-acid translation. The first residue
#'     reflects the actual start codon (not forced to Met); subsequent
#'     residues follow the selected genetic code.}
#' }
#'
#' Returns an empty `data.frame` (invisibly) when no ORFs are found.
#'
#' @details
#' ## Genetic codes
#' The vertebrate mitochondrial code (`"SGC1"`) differs from the standard
#' code: TGA = Trp; AGA = stop; AGG = stop; ATA = Met.  Always choose the
#' code matching the taxon whose mitogenome you are analysing — results
#' obtained with the wrong code will include spurious truncations (stop
#' codons read through) or missed ORFs.
#'
#' ## Circular genomes
#' When `circular = TRUE` the genome is temporarily duplicated.  An ORF is
#' retained only when its start codon lies in the first copy (positions
#' 1–*L*), eliminating wholly duplicated hits while capturing boundary-
#' spanning ORFs.  For such ORFs `wraps_around = TRUE` and `end < start`.
#'
#' ## Overlapping ORFs
#' All ORFs meeting the length threshold are reported, including nested and
#' overlapping entries within the same reading frame.  To retain only the
#' longest non-overlapping set, filter the returned table downstream.
#'
#' @examples
#' \dontrun{
#' library(Biostrings)
#'
#' # Quick synthetic test (standard code, ATG only, linear)
#' seq <- paste0("ATGAAACCCGGGTTTTAA", strrep("A", 50))
#' orfs <- find_orfs(seq,
#'                   genetic_code   = "SGC0",
#'                   start_codons   = "ATG",
#'                   min_orf_length = 15,
#'                   circular       = FALSE)
#' print(orfs[, c("strand","frame","start","end","length_nt","protein_sequence")])
#'
#' # Real mitochondrial genome from a FASTA file
#' mito <- readDNAStringSet("my_mito.fasta")[[1]]
#' orfs <- find_orfs(mito)
#' }
#'
#' @importFrom Biostrings DNAString DNAStringSet getGeneticCode
#'   reverseComplement translate xscat
#' @importFrom methods is
#' @export
find_orfs <- function(sequence,
                      genetic_code       = "SGC0",
                      start_codons       = "ATG",
                      min_orf_length     = 100L,
                      max_orf_length     = Inf,
                      both_strands       = TRUE,
                      circular           = TRUE,
                      include_stop_codon = TRUE) {

  ## ---- resolve input to one or more DNAStrings ----------------------------
  seq_list <- .resolve_dna_input(sequence)

  # If multiple sequences, run recursively and combine
  if (length(seq_list) > 1L) {
    parts <- lapply(names(seq_list), function(nm) {
      df        <- find_orfs(seq_list[[nm]],
                             genetic_code       = genetic_code,
                             start_codons       = start_codons,
                             min_orf_length     = min_orf_length,
                             max_orf_length     = max_orf_length,
                             both_strands       = both_strands,
                             circular           = circular,
                             include_stop_codon = include_stop_codon)
      df$seq_id <- nm
      df
    })
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    return(out)
  }

  sequence <- seq_list[[1L]]

  ## ---- input validation & normalisation ------------------------------------
  genome_len     <- length(sequence)          # Biostrings length() = nt count
  code           <- Biostrings::getGeneticCode(genetic_code)
  stop_cdns      <- names(code)[code == "*"]
  start_codons   <- toupper(start_codons)
  min_orf_length <- as.integer(min_orf_length)
  max_orf_length <- if (is.infinite(max_orf_length)) Inf else as.integer(max_orf_length)

  if (genome_len < 3L)
    stop("'sequence' is too short to contain any codons.", call. = FALSE)

  ## ---- extend for circular search ------------------------------------------
  # Doubling ensures ORFs spanning the origin are found.
  work_seq <- if (circular) Biostrings::xscat(sequence, sequence) else sequence
  work_len <- length(work_seq)   # genome_len or 2 * genome_len

  ## ---- inner search: one strand --------------------------------------------
  .search_strand <- function(strand_seq, strand_label) {
    sq  <- as.character(strand_seq)
    len <- nchar(sq)

    frame_results <- vector("list", 3L)

    for (fr in 1L:3L) {
      chunk   <- substr(sq, fr, len)
      n_cod   <- nchar(chunk) %/% 3L
      if (n_cod < 2L) next

      # Vectorised codon extraction
      cod_s  <- seq(1L, by = 3L, length.out = n_cod)
      cod_e  <- seq(3L, by = 3L, length.out = n_cod)
      codons <- substring(chunk, cod_s, cod_e)

      is_start <- codons %in% start_codons
      is_stop  <- codons %in% stop_cdns

      start_idxs <- which(is_start)
      if (length(start_idxs) == 0L) next

      orf_rows <- vector("list", length(start_idxs))
      k <- 0L

      for (s_ci in start_idxs) {
        # First stop codon strictly after the start codon
        e_ci_candidates <- which(is_stop & seq_along(codons) > s_ci)
        if (length(e_ci_candidates) == 0L) next
        e_ci <- e_ci_candidates[[1L]]

        # Nucleotide positions within strand_seq (1-based, inclusive)
        nt_start <- fr + (s_ci - 1L) * 3L
        nt_end   <- fr + e_ci * 3L - 1L
        if (!include_stop_codon) nt_end <- nt_end - 3L

        len_nt <- nt_end - nt_start + 1L
        if (len_nt < min_orf_length || len_nt > max_orf_length) next

        # Skip ORFs whose start codon lies entirely in the duplicated half
        if (circular && nt_start > genome_len) next

        orf_str <- substr(sq, nt_start, nt_end)
        prot    <- as.character(
          Biostrings::translate(
            Biostrings::DNAString(orf_str),
            genetic.code   = code,
            no.init.codon  = TRUE,      # translate start codon as-is
            if.fuzzy.codon = "solve"    # handle N-containing codons
          )
        )

        ## -- genomic coordinate mapping ------------------------------------
        wraps <- FALSE

        if (strand_label == "+") {
          g_start <- nt_start   # always in [1, genome_len] (filtered above)
          if (circular && nt_end > genome_len) {
            g_end <- ((nt_end - 1L) %% genome_len) + 1L
            wraps <- TRUE
          } else {
            g_end <- nt_end
          }
          frame_val <- fr

        } else {
          # strand_seq is RC of work_seq; map back to forward coordinates.
          # RC position p in a sequence of length W maps to forward position
          # W - p + 1.  Here W = work_len.
          g_end   <- ((work_len - nt_start) %% genome_len) + 1L
          g_start <- ((work_len - nt_end)   %% genome_len) + 1L
          wraps   <- g_start > g_end
          frame_val <- -fr
        }

        k <- k + 1L
        orf_rows[[k]] <- data.frame(
          strand           = strand_label,
          frame            = frame_val,
          start            = g_start,
          end              = g_end,
          length_nt        = len_nt,
          length_aa        = nchar(prot),
          start_codon      = codons[s_ci],
          stop_codon       = if (include_stop_codon) codons[e_ci] else NA_character_,
          wraps_around     = wraps,
          orf_sequence     = orf_str,
          protein_sequence = prot,
          stringsAsFactors = FALSE
        )
      }

      if (k > 0L)
        frame_results[[fr]] <- do.call(rbind, orf_rows[seq_len(k)])
    }

    do.call(rbind, Filter(Negate(is.null), frame_results))
  }

  ## ---- search both strands -------------------------------------------------
  results <- list(.search_strand(work_seq, "+"))

  if (both_strands) {
    rc_work    <- Biostrings::reverseComplement(work_seq)
    results[[2L]] <- .search_strand(rc_work, "-")
  }

  ## ---- assemble output -----------------------------------------------------
  out <- do.call(rbind, Filter(Negate(is.null), results))

  if (is.null(out) || nrow(out) == 0L) {
    message("No ORFs found matching the specified criteria.")
    return(invisible(
      data.frame(
        strand = character(), frame = integer(),
        start = integer(), end = integer(),
        length_nt = integer(), length_aa = integer(),
        start_codon = character(), stop_codon = character(),
        wraps_around = logical(),
        orf_sequence = character(), protein_sequence = character(),
        stringsAsFactors = FALSE
      )
    ))
  }

  out <- out[order(out$start, out$strand), ]
  rownames(out) <- NULL
  out
}
