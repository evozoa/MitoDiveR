# ============================================================================
# translate_all_frames — translate a mitogenome in all 6 frames × N codes
# ============================================================================

#' Translate a Mitochondrial Genome in All Six Reading Frames
#'
#' Translates a mitogenome sequence in all three forward and three
#' reverse-complement reading frames, optionally under multiple genetic codes.
#' Stop codons are replaced with `X` so that the resulting strings can be used
#' as continuous amino acid sequences for window-based homology searches with
#' [find_conserved_windows()].
#'
#' @param genome_record A single genome record — one element of the list
#'   returned by [fetch_mito_genbank()].  Must contain `$sequence`
#'   (a [Biostrings::DNAString]) and `$accession`.
#' @param genetic_codes Character vector of Biostrings genetic-code
#'   identifiers passed to [Biostrings::getGeneticCode()].
#'   Default `c("SGC0", "SGC1")` translates under both the standard
#'   vertebrate nuclear code and the vertebrate mitochondrial code,
#'   yielding 12 strings total.
#'
#' @return A named list of [Biostrings::AAString] objects.  Names follow the
#'   pattern `{accession}_{code}_{strand}{frame}`, e.g.
#'   `NC_002333.2_SGC0_+1`, `NC_002333.2_SGC1_-3`.  Stop codons are
#'   represented as `X`.
#'
#' @details
#' Each reading frame is translated as a single uninterrupted string — stop
#' codons do **not** terminate translation.  This preserves the full-length
#' context needed for inter-species window comparison and allows
#' [find_conserved_windows()] to detect conserved regions that span or
#' adjoin canonical stop codons in one species but not another.
#'
#' Reverse-complement frames are numbered `-1`, `-2`, `-3` following the
#' same convention as [find_orfs()]: `-1` starts from the last nucleotide
#' of the genome, `-2` from the second-to-last, etc.
#'
#' @examples
#' \dontrun{
#' genomes <- fetch_mito_genbank("NC_002333.2")
#' tr <- translate_all_frames(genomes[["NC_002333.2"]])
#' names(tr)
#' nchar(tr[["NC_002333.2_SGC0_+1"]])
#' }
#'
#' @importFrom Biostrings getGeneticCode reverseComplement DNAString AAString
#'   translate
#' @importFrom methods is
#' @export
translate_all_frames <- function(genome_record,
                                  genetic_codes = c("SGC0", "SGC1")) {

  ## ---- input checks --------------------------------------------------------
  if (!is.list(genome_record) || !"sequence" %in% names(genome_record))
    stop("'genome_record' must be a single element from fetch_mito_genbank().",
         call. = FALSE)

  seq     <- genome_record$sequence
  acc     <- genome_record$accession
  seq_len <- length(seq)

  if (!methods::is(seq, "DNAString"))
    stop("'genome_record$sequence' must be a DNAString.", call. = FALSE)

  rc_seq <- Biostrings::reverseComplement(seq)

  results <- list()

  for (code_id in genetic_codes) {
    code <- Biostrings::getGeneticCode(code_id)

    for (fr in 1L:3L) {
      ## -- forward strand ----------------------------------------------------
      fwd_chunk  <- Biostrings::DNAString(substr(as.character(seq), fr, seq_len))
      fwd_aa     <- .translate_no_stop(fwd_chunk, code)
      fwd_key    <- paste0(acc, "_", code_id, "_+", fr)
      results[[fwd_key]] <- fwd_aa

      ## -- reverse strand ----------------------------------------------------
      rev_chunk  <- Biostrings::DNAString(substr(as.character(rc_seq), fr, seq_len))
      rev_aa     <- .translate_no_stop(rev_chunk, code)
      rev_key    <- paste0(acc, "_", code_id, "_-", fr)
      results[[rev_key]] <- rev_aa
    }
  }

  results
}


# --------------------------------------------------------------------------- #
# Internal: translate a DNA chunk, replacing stop codons with X               #
# --------------------------------------------------------------------------- #
.translate_no_stop <- function(dna_chunk, code) {
  # Trim to codon-complete length
  n_nt    <- length(dna_chunk)
  trim    <- n_nt - (n_nt %% 3L)
  if (trim < 3L) return(Biostrings::AAString(""))

  dna_trim <- Biostrings::DNAString(substr(as.character(dna_chunk), 1L, trim))

  aa <- Biostrings::translate(
    dna_trim,
    genetic.code   = code,
    no.init.codon  = TRUE,
    if.fuzzy.codon = "solve"
  )

  # Replace stop (*) with X
  aa_str <- gsub("[*]", "X", as.character(aa))
  Biostrings::AAString(aa_str)
}
