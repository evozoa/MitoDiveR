# ============================================================================
# compare_mitogenomes — pairwise mitochondrial genome comparison
# ============================================================================

# ── internal helpers ──────────────────────────────────────────────────────────

# Load a sequence from an accession, file path, DNAString, or raw string.
# Returns list(seq = DNAString, name = character).
.mito_load_seq <- function(x, label) {

  ## Already a DNAString
  if (methods::is(x, "DNAString")) {
    nm <- if (!is.null(label) && nchar(label) > 0L) label else "sequence"
    return(list(seq = x, name = nm))
  }

  ## DNAStringSet (single-entry)
  if (methods::is(x, "DNAStringSet")) {
    if (length(x) != 1L)
      stop("DNAStringSet input must contain exactly one sequence.", call. = FALSE)
    nm <- if (!is.null(label) && nchar(label) > 0L) label else
          if (!is.null(names(x)) && nchar(names(x)[1L]) > 0L) names(x)[1L] else "sequence"
    return(list(seq = x[[1L]], name = nm))
  }

  if (!is.character(x) || length(x) != 1L)
    stop("Sequence input must be a character string, DNAString, or DNAStringSet.",
         call. = FALSE)

  x_trimmed <- trimws(x)

  ## NCBI accession: 6-12 non-whitespace characters, dots/underscores allowed
  if (grepl("^[A-Za-z0-9_.]{6,20}$", x_trimmed)) {
    seqset <- fetch_sequences(accessions = x_trimmed)
    nm <- if (!is.null(label) && nchar(label) > 0L) label else names(seqset)[1L]
    return(list(seq = seqset[[1L]], name = nm))
  }

  ## File path
  if (file.exists(x_trimmed)) {
    seqset <- Biostrings::readDNAStringSet(x_trimmed)
    if (length(seqset) == 0L)
      stop("FASTA file '", x_trimmed, "' contains no sequences.", call. = FALSE)
    nm <- if (!is.null(label) && nchar(label) > 0L) label else
          sub("\\s.*$", "", names(seqset)[1L])
    return(list(seq = seqset[[1L]], name = nm))
  }

  ## Long raw nucleotide string (>100 chars, predominantly ACGTN)
  clean <- toupper(gsub("[^ACGTNacgtn]", "", x_trimmed))
  if (nchar(clean) > 100L && nchar(clean) >= 0.9 * nchar(x_trimmed)) {
    nm <- if (!is.null(label) && nchar(label) > 0L) label else "sequence"
    return(list(seq = Biostrings::DNAString(clean), name = nm))
  }

  stop("Cannot interpret sequence input: '", substr(x_trimmed, 1L, 60L), "...'",
       call. = FALSE)
}


# Parse a GenBank flat-file record for ALL feature types (CDS, rRNA, tRNA).
# Returns data.frame(feature, gene, start, end, strand).
.parse_gb_features <- function(acc) {

  gb_text <- tryCatch(
    rentrez::entrez_fetch(db      = "nuccore",
                          id      = acc,
                          rettype = "gb",
                          retmode = "text"),
    error = function(e) {
      warning("Could not fetch GenBank record for '", acc, "': ", conditionMessage(e),
              call. = FALSE)
      NULL
    }
  )
  if (is.null(gb_text)) return(NULL)

  lines <- strsplit(gb_text, "\n")[[1L]]

  ## ---- locate the FEATURES table -------------------------------------------
  feat_start <- which(grepl("^FEATURES", lines))[1L]
  feat_end   <- which(grepl("^ORIGIN", lines))[1L]
  if (is.na(feat_start) || is.na(feat_end)) return(NULL)

  feat_lines <- lines[seq(feat_start + 1L, feat_end - 1L)]

  ## ---- state-machine parser -------------------------------------------------
  rows           <- list()
  cur_feature    <- NULL   # "CDS" / "rRNA" / "tRNA"
  cur_loc        <- ""
  cur_gene       <- ""
  cur_product    <- ""
  in_qualifiers  <- FALSE  # TRUE once the first /qualifier line is seen

  flush_feature <- function() {
    if (is.null(cur_feature)) return()
    if (!cur_feature %in% c("CDS", "rRNA", "tRNA")) return()

    loc   <- cur_loc
    gene  <- cur_gene
    prod  <- cur_product

    ## Parse complement / join notation
    strand <- if (grepl("complement", loc, fixed = TRUE)) "-" else "+"
    ## Strip complement( )
    loc_clean <- gsub("complement\\(|\\)", "", loc)
    ## Strip join( )
    loc_clean <- gsub("join\\(|\\)", "", loc_clean)
    ## Take the min/max of all coordinate ranges
    all_nums <- as.integer(unlist(regmatches(loc_clean,
                                              gregexpr("[0-9]+", loc_clean))))
    if (length(all_nums) < 2L) return()
    feat_start_coord <- min(all_nums)
    feat_end_coord   <- max(all_nums)

    ## Gene / name detection
    if (nchar(gene) == 0L) {
      ## Try to pull a gene name from the product
      gene <- prod
    }

    ## For rRNA: normalise to 12S / 16S labels
    ## Match on product ("s-rRNA"/"l-rRNA"), gene name (RNR1/RNR2), or
    ## explicit "12S"/"16S" strings anywhere in gene or product
    feat_label <- cur_feature
    if (cur_feature == "rRNA") {
      is_12s <- grepl("12S|12 S|small subunit|s-rRNA|RNR1", prod, ignore.case = TRUE) ||
                grepl("12S|12 S|small subunit|s-rRNA|RNR1", gene, ignore.case = TRUE)
      is_16s <- grepl("16S|16 S|large subunit|l-rRNA|RNR2", prod, ignore.case = TRUE) ||
                grepl("16S|16 S|large subunit|l-rRNA|RNR2", gene, ignore.case = TRUE)
      feat_label <- if (is_12s) "12S_rRNA" else if (is_16s) "16S_rRNA" else "rRNA"
    }
    ## For tRNA: extract amino acid from "tRNA-Xxx"
    if (cur_feature == "tRNA") {
      m <- regexpr("tRNA-[A-Za-z]{3}", prod)
      gene <- if (m != -1L) regmatches(prod, m) else if (nchar(gene) > 0L) gene else "tRNA"
    }

    rows[[length(rows) + 1L]] <<- data.frame(
      feature = feat_label,
      gene    = gene,
      start   = feat_start_coord,
      end     = feat_end_coord,
      strand  = strand,
      stringsAsFactors = FALSE
    )
    cur_feature   <<- NULL
    cur_loc       <<- ""
    cur_gene      <<- ""
    cur_product   <<- ""
    in_qualifiers <<- FALSE
  }

  for (line in feat_lines) {
    ## Feature key line: 5 spaces, key (left-justified), then location
    if (grepl("^     [A-Za-z_]", line) && !grepl("^      ", line)) {
      flush_feature()
      parts         <- strsplit(trimws(line), "\\s+")[[1L]]
      cur_feature   <- parts[1L]
      cur_loc       <- if (length(parts) >= 2L) parts[2L] else ""
      cur_gene      <- ""
      cur_product   <- ""
      in_qualifiers <- FALSE
    } else if (grepl("^      ", line)) {
      trimmed <- trimws(line)
      if (grepl("^/", trimmed)) {
        ## Qualifier line — stop treating subsequent lines as location
        in_qualifiers <- TRUE
        if      (grepl("^/gene=",    trimmed)) cur_gene    <- gsub('^/gene="|"$',    "", trimmed)
        else if (grepl("^/product=", trimmed)) cur_product <- gsub('^/product="|"$', "", trimmed)
      } else if (!in_qualifiers && nchar(cur_loc) > 0L) {
        ## True location continuation line (before any qualifier)
        cur_loc <- paste0(cur_loc, trimmed)
      }
      ## Lines after first qualifier are qualifier value continuations — ignore
    }
  }
  flush_feature()   # catch last feature

  if (length(rows) == 0L) return(NULL)

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


# Align two genome strings (character) and extract SNPs and indels.
# Returns data.frame(pos_ref, pos_alt, ref_nt, alt_nt, mutation_type).
.align_and_extract_mutations <- function(seq1_str, seq2_str) {

  ## Replace non-ACGT characters with A — the substitution matrix is 4-base
  ## only; we need position mapping not exact scoring so any base works here
  seq1_clean <- gsub("[^ACGTacgt]", "A", seq1_str)
  seq2_clean <- gsub("[^ACGTacgt]", "A", seq2_str)

  ## Build a simple +1/-1 substitution matrix
  subst_mat <- pwalign::nucleotideSubstitutionMatrix(match    = 1L,
                                                      mismatch = -1L,
                                                      baseOnly = TRUE)

  aln <- pwalign::pairwiseAlignment(
    pattern            = seq2_clean,
    subject            = seq1_clean,
    type               = "global",
    substitutionMatrix = subst_mat,
    gapOpening         = 10,
    gapExtension       = 0.5
  )

  ref_aln <- as.character(pwalign::subject(aln))   # reference (seq1) in alignment
  alt_aln <- as.character(pwalign::pattern(aln))   # query     (seq2) in alignment

  ref_chars <- strsplit(ref_aln, "")[[1L]]
  alt_chars <- strsplit(alt_aln, "")[[1L]]
  n_col     <- length(ref_chars)

  ## Walk alignment columns and collect raw events
  pos_ref <- 0L
  pos_alt <- 0L
  events  <- list()

  for (col in seq_len(n_col)) {
    r <- ref_chars[col]
    a <- alt_chars[col]

    is_ref_gap <- (r == "-")
    is_alt_gap <- (a == "-")

    if (!is_ref_gap) pos_ref <- pos_ref + 1L
    if (!is_alt_gap) pos_alt <- pos_alt + 1L

    if (is_ref_gap) {
      events[[length(events) + 1L]] <- list(
        pos_ref  = NA_integer_,
        pos_alt  = pos_alt,
        ref_nt   = "-",
        alt_nt   = a,
        type     = "insertion"
      )
    } else if (is_alt_gap) {
      events[[length(events) + 1L]] <- list(
        pos_ref  = pos_ref,
        pos_alt  = NA_integer_,
        ref_nt   = r,
        alt_nt   = "-",
        type     = "deletion"
      )
    } else if (toupper(r) != toupper(a)) {
      events[[length(events) + 1L]] <- list(
        pos_ref  = pos_ref,
        pos_alt  = pos_alt,
        ref_nt   = toupper(r),
        alt_nt   = toupper(a),
        type     = "SNP"
      )
    }
  }

  if (length(events) == 0L) {
    return(data.frame(
      pos_ref       = integer(),
      pos_alt       = integer(),
      ref_nt        = character(),
      alt_nt        = character(),
      mutation_type = character(),
      stringsAsFactors = FALSE
    ))
  }

  ## ---- merge consecutive indels of the same type ----------------------------
  merged <- list()
  cur    <- events[[1L]]

  merge_ok <- function(prev, next_ev) {
    if (prev$type != next_ev$type) return(FALSE)
    if (prev$type == "insertion") {
      return(!is.na(next_ev$pos_alt) && next_ev$pos_alt == prev$pos_alt + nchar(prev$alt_nt))
    }
    if (prev$type == "deletion") {
      return(!is.na(next_ev$pos_ref) && next_ev$pos_ref == prev$pos_ref + nchar(prev$ref_nt))
    }
    FALSE
  }

  for (i in seq(2L, length(events))) {
    ev <- events[[i]]
    if (merge_ok(cur, ev)) {
      cur$ref_nt <- paste0(cur$ref_nt, ev$ref_nt)
      cur$alt_nt <- paste0(cur$alt_nt, ev$alt_nt)
    } else {
      merged[[length(merged) + 1L]] <- cur
      cur <- ev
    }
  }
  merged[[length(merged) + 1L]] <- cur

  ## Assemble data.frame
  pos_ref_v  <- vapply(merged, function(e) if (is.na(e$pos_ref)) NA_integer_ else e$pos_ref, integer(1L))
  pos_alt_v  <- vapply(merged, function(e) if (is.na(e$pos_alt)) NA_integer_ else e$pos_alt, integer(1L))
  ref_nt_v   <- vapply(merged, `[[`, character(1L), "ref_nt")
  alt_nt_v   <- vapply(merged, `[[`, character(1L), "alt_nt")
  type_v     <- vapply(merged, `[[`, character(1L), "type")

  data.frame(
    pos_ref       = pos_ref_v,
    pos_alt       = pos_alt_v,
    ref_nt        = ref_nt_v,
    alt_nt        = alt_nt_v,
    mutation_type = type_v,
    stringsAsFactors = FALSE
  )
}


# Annotate a single mutation position (pos_ref) with its genomic region.
# Priority: CDS > rRNA > tRNA > D-loop > intergenic.
# Returns a list(region, canonical_gene, feat_row) where feat_row is the
# matching features row (or NULL) for downstream codon computation.
.annotate_position <- function(pos_ref, features) {

  if (is.null(features) || nrow(features) == 0L || is.na(pos_ref)) {
    return(list(region = "intergenic", canonical_gene = NA_character_, feat_row = NULL))
  }

  hits <- features[features$start <= pos_ref & features$end >= pos_ref, , drop = FALSE]

  if (nrow(hits) == 0L) {
    return(list(region = "intergenic", canonical_gene = NA_character_, feat_row = NULL))
  }

  ## CDS first
  cds_hits <- hits[hits$feature == "CDS", , drop = FALSE]
  if (nrow(cds_hits) > 0L) {
    row_i <- cds_hits[1L, , drop = FALSE]
    gene  <- toupper(row_i$gene)
    return(list(
      region         = paste0("CDS:", gene),
      canonical_gene = gene,
      feat_row       = row_i
    ))
  }

  ## rRNA (covers 12S_rRNA, 16S_rRNA, rRNA)
  rrna_hits <- hits[hits$feature %in% c("12S_rRNA", "16S_rRNA", "rRNA"), , drop = FALSE]
  if (nrow(rrna_hits) > 0L) {
    lbl <- rrna_hits$feature[1L]
    return(list(region = lbl, canonical_gene = NA_character_, feat_row = rrna_hits[1L, , drop = FALSE]))
  }

  ## tRNA
  trna_hits <- hits[hits$feature == "tRNA", , drop = FALSE]
  if (nrow(trna_hits) > 0L) {
    gene <- trna_hits$gene[1L]
    lbl  <- if (nchar(gene) > 0L) paste0("tRNA:", gene) else "tRNA"
    return(list(region = lbl, canonical_gene = NA_character_, feat_row = trna_hits[1L, , drop = FALSE]))
  }

  ## D-loop (feature == "D-loop" etc. — handled as anything whose gene contains D-loop)
  dloop_hits <- hits[grepl("D.loop|control|d-loop", hits$gene, ignore.case = TRUE) |
                     grepl("D.loop|control|d-loop", hits$feature, ignore.case = TRUE), , drop = FALSE]
  if (nrow(dloop_hits) > 0L) {
    return(list(region = "D-loop", canonical_gene = NA_character_, feat_row = NULL))
  }

  list(region = "intergenic", canonical_gene = NA_character_, feat_row = NULL)
}


# Compute codon-level effect of a SNP at pos_ref within a CDS feature.
# Returns list(codon_ref, codon_alt, aa_ref, aa_alt, effect).
.snp_codon_effect <- function(pos_ref, ref_nt, alt_nt, cds_row, seq_str, genetic_code) {

  code  <- Biostrings::getGeneticCode(genetic_code)

  ## Offset within the CDS (0-based)
  if (cds_row$strand == "+") {
    offset <- pos_ref - cds_row$start   # 0-based
  } else {
    offset <- cds_row$end - pos_ref     # 0-based
  }

  if (offset < 0L) {
    return(list(codon_ref = NA_character_, codon_alt = NA_character_,
                aa_ref = NA_character_, aa_alt = NA_character_,
                effect = NA_character_))
  }

  codon_idx <- offset %/% 3L           # 0-based codon number
  codon_pos <- offset %%  3L           # 0-based position within codon (0,1,2)

  ## Genomic coordinates of that codon (always on forward strand of seq_str)
  if (cds_row$strand == "+") {
    cod_g_start <- cds_row$start + codon_idx * 3L       # 1-based
    cod_g_end   <- cod_g_start + 2L
  } else {
    cod_g_end   <- cds_row$end - codon_idx * 3L         # 1-based
    cod_g_start <- cod_g_end - 2L
  }

  if (cod_g_start < 1L || cod_g_end > nchar(seq_str)) {
    return(list(codon_ref = NA_character_, codon_alt = NA_character_,
                aa_ref = NA_character_, aa_alt = NA_character_,
                effect = NA_character_))
  }

  codon_fwd <- substr(seq_str, cod_g_start, cod_g_end)

  if (cds_row$strand == "-") {
    codon_sense <- as.character(
      Biostrings::reverseComplement(Biostrings::DNAString(codon_fwd))
    )
  } else {
    codon_sense <- codon_fwd
  }

  ## Build alt codon
  codon_sense_chars        <- strsplit(codon_sense, "")[[1L]]

  ## On minus strand the mutant base is the complement of alt_nt at the
  ## genomic position; codon_pos in the sense-strand codon is (2 - codon_pos)
  if (cds_row$strand == "-") {
    sense_codon_pos <- 2L - codon_pos
    alt_sense_nt    <- as.character(
      Biostrings::reverseComplement(Biostrings::DNAString(toupper(alt_nt)))
    )
  } else {
    sense_codon_pos <- codon_pos
    alt_sense_nt    <- toupper(alt_nt)
  }

  codon_alt_chars                    <- codon_sense_chars
  codon_alt_chars[sense_codon_pos + 1L] <- alt_sense_nt
  codon_alt_sense <- paste(codon_alt_chars, collapse = "")

  ## Translate
  safe_translate <- function(cdn) {
    tryCatch(
      as.character(Biostrings::translate(
        Biostrings::DNAString(toupper(cdn)),
        genetic.code   = code,
        no.init.codon  = TRUE,
        if.fuzzy.codon = "solve"
      )),
      error = function(e) "?"
    )
  }

  aa_ref <- safe_translate(codon_sense)
  aa_alt <- safe_translate(codon_alt_sense)

  effect <- if (aa_ref == aa_alt) {
    "synonymous"
  } else if (aa_alt == "*") {
    "nonsense"
  } else {
    "missense"
  }

  list(
    codon_ref = codon_sense,
    codon_alt = codon_alt_sense,
    aa_ref    = aa_ref,
    aa_alt    = aa_alt,
    effect    = effect
  )
}


# Fold the rRNA sequence from a genome string using annotated coordinates.
# Returns the rna_fold object or NULL on failure.
.fold_rrna_region <- function(seq_str, feat_row, temperature) {
  if (is.null(feat_row) || nrow(feat_row) == 0L) return(NULL)
  rrna_seq <- substr(seq_str, feat_row$start[1L], feat_row$end[1L])
  if (nchar(rrna_seq) == 0L) return(NULL)
  tryCatch(
    fold_rna(rrna_seq, temperature = temperature),
    error = function(e) NULL
  )
}


# Check whether RNAfold is available without stopping.
.rnafold_available <- function() {
  nchar(Sys.which("RNAfold")) > 0L
}


# Extract the mRNA sequence for a CDS feature from a genome string.
# For minus-strand genes the reverse complement is returned.
.extract_cds_seq <- function(seq_str, feat_row) {
  raw <- substr(seq_str, as.integer(feat_row$start[1L]), as.integer(feat_row$end[1L]))
  if (feat_row$strand[1L] == "-")
    raw <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(raw)))
  raw
}


# Introduce a single mutation (SNP or indel) into a CDS string.
# pos_ref    : 1-based genomic position of the mutation in the reference
# ref_nt/alt_nt : genomic-strand alleles ("-" or "--" = gap for indels)
# feat_row   : one-row data.frame with $start, $end, $strand
# Returns the mutated CDS string, or NULL if the position is out of range.
.introduce_mutation_in_cds <- function(cds_seq, pos_ref, ref_nt, alt_nt, feat_row) {
  s <- as.integer(feat_row$start[1L])
  e <- as.integer(feat_row$end[1L])

  if (feat_row$strand[1L] == "+") {
    cds_pos <- pos_ref - s + 1L
    ins_nt  <- gsub("-", "", toupper(alt_nt))   # deletion → ""
    del_len <- nchar(gsub("-", "", ref_nt))     # insertion ref side → 0
  } else {
    ## Minus-strand: CDS string is the RC; position flips and base complements
    cds_pos <- e - pos_ref + 1L
    ins_raw <- gsub("-", "", toupper(alt_nt))
    ins_nt  <- if (nchar(ins_raw) > 0L)
                 as.character(Biostrings::reverseComplement(Biostrings::DNAString(ins_raw)))
               else ""
    del_len <- nchar(gsub("-", "", ref_nt))
  }

  if (cds_pos < 1L || cds_pos > nchar(cds_seq)) return(NULL)
  paste0(substr(cds_seq, 1L, cds_pos - 1L),
         ins_nt,
         substr(cds_seq, cds_pos + del_len, nchar(cds_seq)))
}


# Determine if an ORF (start, end, strand) substantially overlaps a CDS feature
# (>50% of ORF length).
.orf_overlaps_cds <- function(orf_start, orf_end, orf_strand, orf_len, cds_df) {
  if (is.null(cds_df) || nrow(cds_df) == 0L) return(FALSE)
  cds <- cds_df[cds_df$feature == "CDS" & cds_df$strand == orf_strand, , drop = FALSE]
  if (nrow(cds) == 0L) return(FALSE)
  any(vapply(seq_len(nrow(cds)), function(i) {
    overlap <- min(orf_end, cds$end[i]) - max(orf_start, cds$start[i]) + 1L
    overlap > 0L && (overlap / orf_len) > 0.5
  }, logical(1L)))
}


# Build a human-readable mdp_effect string from lists of lost/gained/changed ORFs.
.mdp_effect_string <- function(lost, gained, changed) {
  parts <- character(0L)

  if (!is.null(lost) && nrow(lost) > 0L) {
    for (i in seq_len(nrow(lost))) {
      parts <- c(parts, sprintf("ORF lost (frame%+d, pos %d-%d, %d nt)",
                                lost$frame[i], lost$start[i], lost$end[i], lost$length_nt[i]))
    }
  }
  if (!is.null(gained) && nrow(gained) > 0L) {
    for (i in seq_len(nrow(gained))) {
      parts <- c(parts, sprintf("ORF gained (frame%+d, pos %d-%d, %d nt)",
                                gained$frame[i], gained$start[i], gained$end[i], gained$length_nt[i]))
    }
  }
  if (!is.null(changed) && nrow(changed) > 0L) {
    for (i in seq_len(nrow(changed))) {
      parts <- c(parts, sprintf("ORF changed (frame%+d, pos %d-%d, %d nt)",
                                changed$frame[i], changed$start[i], changed$end[i], changed$length_nt[i]))
    }
  }

  if (length(parts) == 0L) return(NA_character_)
  paste(parts, collapse = "; ")
}


# ── main exported function ────────────────────────────────────────────────────

#' Compare Two Mitochondrial Genome Sequences
#'
#' @title Compare Two Mitochondrial Genome Sequences
#'
#' @description
#' Given two mitochondrial genome sequences, performs a global pairwise
#' alignment and produces a mutation table with one row per SNP or indel,
#' annotated with the overlapping genomic region, canonical protein effect (for
#' CDS mutations), RNA secondary-structure impact (for rRNA mutations), and
#' MDP ORF impact (gained / lost / changed small ORFs).
#'
#' @param seq1 Reference sequence.  Accepts an NCBI accession string (e.g.
#'   `"NC_012920.1"`), a path to a FASTA file, a
#'   [Biostrings::DNAString], or a raw nucleotide character string (>100 nt).
#' @param seq2 Query sequence.  Same accepted types as `seq1`.
#' @param name1 Character scalar.  Display label for `seq1`.  Auto-detected
#'   from the accession / FASTA header / object name when `NULL` (default).
#' @param name2 Character scalar.  Display label for `seq2`.  Same rules.
#' @param genetic_code Character scalar.  Biostrings genetic-code identifier.
#'   Default `"SGC1"` (vertebrate mitochondrial code, NCBI table 2).
#' @param min_mdp_nt Integer.  Minimum ORF length in nucleotides for MDP
#'   scanning.  Default `60L`.
#' @param fold_rrna Logical.  Compute RNA secondary-structure delta metrics
#'   (ΔMfe, Δbase-pairs, Δhairpins) for mutations that fall in annotated rRNA
#'   features?  Requires ViennaRNA `RNAfold` on the system `PATH`.  Default
#'   `TRUE`.
#' @param scan_mdps Logical.  Scan both genomes for small ORFs (MDP candidates)
#'   and report gained / lost / changed ORFs at each mutation position?
#'   Default `TRUE`.
#' @param temperature Numeric.  Folding temperature in °C passed to
#'   [fold_rna()].  Default `37`.
#'
#' @return A `data.frame` of class `c("mito_comparison", "data.frame")` with
#'   one row per mutation and the following columns:
#' \describe{
#'   \item{`pos_ref`}{1-based position in the reference genome.  `NA` for pure
#'     insertions.}
#'   \item{`pos_alt`}{1-based position in the query genome.  `NA` for pure
#'     deletions.}
#'   \item{`ref_nt`}{Reference nucleotide(s); `"-"` for insertions.}
#'   \item{`alt_nt`}{Query nucleotide(s); `"-"` for deletions.}
#'   \item{`mutation_type`}{`"SNP"`, `"insertion"`, or `"deletion"`.}
#'   \item{`region`}{Genomic region, e.g. `"CDS:ND1"`, `"16S_rRNA"`,
#'     `"tRNA:Phe"`, `"intergenic"`, `"D-loop"`.}
#'   \item{`canonical_gene`}{Gene name when in an annotated CDS; otherwise
#'     `NA`.}
#'   \item{`codon_ref`}{Three-letter reference codon (sense strand).  `NA` for
#'     non-CDS mutations or frameshift-causing indels.}
#'   \item{`codon_alt`}{Alternate codon.}
#'   \item{`aa_ref`}{Reference amino acid (single-letter IUPAC).}
#'   \item{`aa_alt`}{Alternate amino acid.}
#'   \item{`aa_effect`}{`"synonymous"`, `"missense"`, `"nonsense"`,
#'     `"frameshift"`, `"inframe_indel"`, or `NA`.}
#'   \item{`rrna_mfe_ref`}{Minimum free energy of the reference rRNA sequence
#'     (kcal/mol).  `NA` if the mutation is not in an rRNA region, or if
#'     `fold_rrna = FALSE`, or if ViennaRNA is unavailable.}
#'   \item{`rrna_mfe_alt`}{MFE of the query rRNA sequence.}
#'   \item{`rrna_delta_mfe`}{`rrna_mfe_alt − rrna_mfe_ref`.}
#'   \item{`rrna_delta_pairs`}{Change in base-pair count (alt − ref).}
#'   \item{`rrna_delta_hairpins`}{Change in hairpin count (alt − ref).}
#'   \item{`mdp_n_ref`}{Number of MDP ORF candidates overlapping this position
#'     in the reference.  `NA` if `scan_mdps = FALSE`.}
#'   \item{`mdp_n_alt`}{Number in the query.}
#'   \item{`mdp_effect`}{Human-readable summary of gained / lost / changed MDP
#'     ORFs, or `NA`.}
#' }
#' Returns an empty `data.frame` with the correct column names (invisibly, with
#' a message) when no mutations are found.
#'
#' @details
#' ## Alignment
#' The two genome sequences are aligned with [pwalign::pairwiseAlignment()]
#' using a global Needleman-Wunsch alignment with a +1/−1 nucleotide
#' substitution matrix, gap-open penalty 10, and gap-extension penalty 0.5.
#' Non-ACGT characters are replaced with N before alignment to avoid IUPAC
#' ambiguity-code lookup errors.  Consecutive gap columns of the same type are
#' merged into a single indel event.
#'
#' ## Annotation
#' Genomic features (CDS, rRNA, tRNA) are obtained by downloading the GenBank
#' flat-file record for `seq1` via [rentrez::entrez_fetch()].  If the download
#' fails (e.g. `seq1` is a raw sequence without an accession), the function
#' falls back to [fetch_annotations()] for CDS-only information; rRNA columns
#' are then set to `NA`.
#'
#' ## Protein effect
#' For SNPs inside annotated CDS features the reference and alternate codons are
#' extracted from the reference sequence and classified as `"synonymous"`,
#' `"missense"`, or `"nonsense"`.  For indels inside CDS features, the effect
#' is `"frameshift"` when `length(indel) %% 3 != 0`, else `"inframe_indel"`.
#'
#' ## rRNA structure
#' When `fold_rrna = TRUE` and ViennaRNA `RNAfold` is on the PATH, the entire
#' rRNA sequence (12S or 16S, separately) is extracted from each genome using
#' the annotated coordinates and folded with [fold_rna()].  All mutations
#' falling in the same rRNA are annotated with the same delta values
#' (alt − ref).  If `fold_rrna = FALSE` or `RNAfold` is not found, a warning
#' is issued and rRNA columns are set to `NA`.
#'
#' ## MDP ORFs
#' When `scan_mdps = TRUE`, [scan_orfs()] is called on both genomes
#' (`both_strands = TRUE`, `circular = TRUE`) with `min_orf_length = min_mdp_nt`.
#' ORFs that overlap an annotated CDS by more than 50 % of the ORF length are
#' excluded.  For each mutation, the ORFs overlapping that position in the
#' reference and query are compared: ORFs present only in the reference are
#' "lost"; those present only in the query are "gained"; those present in both
#' but with different protein sequences are "changed".
#'
#' @examples
#' \dontrun{
#' # Compare human reference mitogenome to chimpanzee
#' result <- compare_mitogenomes("NC_012920.1", "NC_001643.1")
#' head(result)
#'
#' # Only structural annotation, no folding or MDP scanning
#' result <- compare_mitogenomes("NC_012920.1", "NC_001643.1",
#'                               fold_rrna = FALSE, scan_mdps = FALSE)
#'
#' # Supply sequences as DNAString objects
#' library(Biostrings)
#' s1 <- readDNAStringSet("human_mito.fasta")[[1]]
#' s2 <- readDNAStringSet("chimp_mito.fasta")[[1]]
#' result <- compare_mitogenomes(s1, s2, name1 = "human", name2 = "chimp")
#' }
#'
#' @importFrom Biostrings DNAString DNAStringSet reverseComplement translate
#'   getGeneticCode subseq readDNAStringSet
#' @importFrom pwalign pairwiseAlignment nucleotideSubstitutionMatrix
#'   pattern subject
#' @importFrom rentrez entrez_fetch
#' @importFrom methods is
#' @export
compare_mitogenomes <- function(
    seq1,
    seq2,
    name1        = NULL,
    name2        = NULL,
    genetic_code = "SGC1",
    min_mdp_nt   = 60L,
    fold_rrna    = TRUE,
    scan_mdps    = TRUE,
    temperature  = 37,
    mrna_fold    = c("per_snp", "joint")
) {
  mrna_fold <- match.arg(mrna_fold)

  ## =========================================================================
  ## 0. Empty result template
  ## =========================================================================
  .empty_result <- function() {
    df <- data.frame(
      pos_ref            = integer(),
      pos_alt            = integer(),
      ref_nt             = character(),
      alt_nt             = character(),
      mutation_type      = character(),
      region             = character(),
      canonical_gene     = character(),
      codon_ref          = character(),
      codon_alt          = character(),
      aa_ref             = character(),
      aa_alt             = character(),
      aa_effect          = character(),
      rrna_mfe_ref       = numeric(),
      rrna_mfe_alt       = numeric(),
      rrna_delta_mfe     = numeric(),
      rrna_delta_pairs   = integer(),
      rrna_delta_hairpins = integer(),
      mdp_n_ref          = integer(),
      mdp_n_alt          = integer(),
      mdp_effect         = character(),
      stringsAsFactors   = FALSE
    )
    class(df) <- c("mito_comparison", "data.frame")
    df
  }

  ## =========================================================================
  ## 1. Load sequences
  ## =========================================================================
  message("Loading sequences...")
  s1 <- .mito_load_seq(seq1, name1)
  s2 <- .mito_load_seq(seq2, name2)

  seq1_str <- as.character(s1$seq)
  seq2_str <- as.character(s2$seq)
  label1   <- s1$name
  label2   <- s2$name

  ## Detect accession for annotation fetching (only attempt if looks like one)
  acc1 <- if (is.character(seq1) && grepl("^[A-Za-z0-9_.]{6,20}$", trimws(seq1)))
    trimws(seq1) else NULL

  ## =========================================================================
  ## 2. Fetch GenBank features
  ## =========================================================================
  message("Fetching GenBank annotations for '", label1, "'...")

  features <- NULL
  if (!is.null(acc1)) {
    features <- .parse_gb_features(acc1)
  }

  if (is.null(features)) {
    ## Fall back to fetch_annotations (CDS only, no rRNA coords)
    if (!is.null(acc1)) {
      ann_fallback <- tryCatch(
        fetch_annotations(acc1),
        error = function(e) {
          warning("fetch_annotations() also failed for '", acc1, "': ",
                  conditionMessage(e), call. = FALSE)
          NULL
        }
      )
      if (!is.null(ann_fallback) && nrow(ann_fallback) > 0L) {
        ## Convert fetch_annotations output to .parse_gb_features schema
        feature_col <- sub(":.*$", "", ann_fallback$region_type)
        gene_col    <- sub("^[^:]+:", "", ann_fallback$region_type)
        gene_col[gene_col == feature_col] <- ""
        features <- data.frame(
          feature  = feature_col,
          gene     = gene_col,
          start    = ann_fallback$start,
          end      = ann_fallback$end,
          strand   = ann_fallback$strand,
          stringsAsFactors = FALSE
        )
      }
    }
    if (is.null(features)) {
      warning("No genomic annotations available; mutations will be labelled 'intergenic'.",
              call. = FALSE)
    }
  }

  ## =========================================================================
  ## 3. Pairwise alignment → mutation table
  ## =========================================================================
  message("Aligning genomes (", nchar(seq1_str), " nt vs ", nchar(seq2_str), " nt)...")

  muts <- .align_and_extract_mutations(seq1_str, seq2_str)

  if (nrow(muts) == 0L) {
    message("No mutations found between '", label1, "' and '", label2, "'.")
    return(invisible(.empty_result()))
  }

  message(nrow(muts), " mutation(s) identified.")

  ## =========================================================================
  ## 4. rRNA folding (pre-compute once per rRNA type)
  ## =========================================================================
  rrna_folds <- list()   # named by feature label (e.g. "16S_rRNA")

  do_fold <- isTRUE(fold_rrna)
  if (do_fold && !.rnafold_available()) {
    warning("RNAfold not found on PATH; rRNA folding disabled. ",
            "Install ViennaRNA: https://www.tbi.univie.ac.at/RNA/",
            call. = FALSE)
    do_fold <- FALSE
  }

  if (do_fold && !is.null(features)) {
    rrna_types <- unique(features$feature[features$feature %in%
                                            c("12S_rRNA", "16S_rRNA", "rRNA")])
    for (rtype in rrna_types) {
      feat_row_r <- features[features$feature == rtype, , drop = FALSE][1L, , drop = FALSE]
      fold1 <- .fold_rrna_region(seq1_str, feat_row_r, temperature)
      fold2 <- .fold_rrna_region(seq2_str, feat_row_r, temperature)
      rrna_folds[[rtype]] <- list(ref = fold1, alt = fold2)
    }
  }

  ## =========================================================================
  ## 4b. mRNA folding for CDS genes (pre-compute once per gene)
  ## =========================================================================
  mrna_folds <- list()   # named by gene (e.g. "ND1")

  if (do_fold && !is.null(features)) {
    cds_feats <- features[features$feature == "CDS", , drop = FALSE]
    for (j in seq_len(nrow(cds_feats))) {
      gene_nm <- cds_feats$gene[j]
      if (nchar(gene_nm) == 0L || gene_nm %in% names(mrna_folds)) next
      feat_j <- cds_feats[j, , drop = FALSE]
      f1 <- tryCatch(.fold_rrna_region(seq1_str, feat_j, temperature),
                     error = function(e) NULL)
      ## In per_snp mode the alt fold is computed per mutation, not per gene
      f2 <- if (mrna_fold == "joint")
               tryCatch(.fold_rrna_region(seq2_str, feat_j, temperature),
                        error = function(e) NULL)
            else NULL
      if (!is.null(f1))
        mrna_folds[[gene_nm]] <- list(ref = f1, alt = f2, feat_row = feat_j)
    }
    if (length(mrna_folds) > 0L)
      message("  Folded mRNA for ", length(mrna_folds), " CDS gene(s) [",
              mrna_fold, " mode].")
  }

  ## =========================================================================
  ## 5. MDP ORF scanning (pre-compute once)
  ## =========================================================================
  mdp_ref <- NULL
  mdp_alt <- NULL

  if (isTRUE(scan_mdps)) {
    message("Scanning MDP ORFs in reference...")
    seqset1 <- Biostrings::DNAStringSet(setNames(seq1_str, label1))
    mdp_ref_raw <- tryCatch(
      scan_orfs(seqset1, min_orf_length = min_mdp_nt, both_strands = TRUE,
                circular = TRUE, genetic_code = genetic_code),
      error = function(e) {
        warning("scan_orfs() failed for reference: ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    message("Scanning MDP ORFs in query...")
    seqset2 <- Biostrings::DNAStringSet(setNames(seq2_str, label2))
    mdp_alt_raw <- tryCatch(
      scan_orfs(seqset2, min_orf_length = min_mdp_nt, both_strands = TRUE,
                circular = TRUE, genetic_code = genetic_code),
      error = function(e) {
        warning("scan_orfs() failed for query: ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )

    ## Filter out ORFs that are essentially canonical CDS
    filter_canonical <- function(orfs_df) {
      if (is.null(orfs_df) || nrow(orfs_df) == 0L) return(orfs_df)
      keep <- vapply(seq_len(nrow(orfs_df)), function(i) {
        !.orf_overlaps_cds(orfs_df$start[i], orfs_df$end[i],
                           orfs_df$strand[i], orfs_df$length_nt[i],
                           features)
      }, logical(1L))
      orfs_df[keep, , drop = FALSE]
    }

    mdp_ref <- filter_canonical(mdp_ref_raw)
    mdp_alt <- filter_canonical(mdp_alt_raw)
  }

  ## =========================================================================
  ## 6. Per-mutation annotation
  ## =========================================================================
  n_muts <- nrow(muts)

  out_region          <- character(n_muts)
  out_canonical_gene  <- character(n_muts)
  out_codon_ref       <- character(n_muts)
  out_codon_alt       <- character(n_muts)
  out_aa_ref          <- character(n_muts)
  out_aa_alt          <- character(n_muts)
  out_aa_effect       <- character(n_muts)
  out_rrna_mfe_ref    <- rep(NA_real_,    n_muts)
  out_rrna_mfe_alt    <- rep(NA_real_,    n_muts)
  out_rrna_delta_mfe  <- rep(NA_real_,    n_muts)
  out_rrna_delta_pairs    <- rep(NA_integer_, n_muts)
  out_rrna_delta_hairpins <- rep(NA_integer_, n_muts)
  out_mrna_mfe_ref    <- rep(NA_real_,    n_muts)
  out_mrna_mfe_alt    <- rep(NA_real_,    n_muts)
  out_mrna_delta_mfe  <- rep(NA_real_,    n_muts)
  out_mrna_delta_pairs    <- rep(NA_integer_, n_muts)
  out_mrna_delta_hairpins <- rep(NA_integer_, n_muts)
  out_mdp_n_ref       <- rep(NA_integer_, n_muts)
  out_mdp_n_alt       <- rep(NA_integer_, n_muts)
  out_mdp_effect      <- rep(NA_character_, n_muts)

  for (i in seq_len(n_muts)) {

    pos_r <- muts$pos_ref[i]
    pos_a <- muts$pos_alt[i]
    rnt   <- muts$ref_nt[i]
    ant   <- muts$alt_nt[i]
    mtype <- muts$mutation_type[i]

    ## ---- 6a. Region annotation ----
    ann_i <- .annotate_position(pos_r, features)
    out_region[i]         <- ann_i$region
    out_canonical_gene[i] <- if (is.na(ann_i$canonical_gene)) NA_character_
                              else ann_i$canonical_gene

    ## ---- 6b. Protein effect ----
    if (grepl("^CDS:", ann_i$region) && !is.null(ann_i$feat_row)) {
      if (mtype == "SNP") {
        ce <- .snp_codon_effect(pos_r, rnt, ant, ann_i$feat_row,
                                seq1_str, genetic_code)
        out_codon_ref[i] <- if (is.na(ce$codon_ref)) NA_character_ else ce$codon_ref
        out_codon_alt[i] <- if (is.na(ce$codon_alt)) NA_character_ else ce$codon_alt
        out_aa_ref[i]    <- if (is.na(ce$aa_ref))    NA_character_ else ce$aa_ref
        out_aa_alt[i]    <- if (is.na(ce$aa_alt))    NA_character_ else ce$aa_alt
        out_aa_effect[i] <- if (is.na(ce$effect))    NA_character_ else ce$effect
      } else {
        ## indel in CDS
        indel_len        <- if (mtype == "deletion") nchar(rnt) else nchar(ant)
        out_codon_ref[i] <- NA_character_
        out_codon_alt[i] <- NA_character_
        out_aa_ref[i]    <- NA_character_
        out_aa_alt[i]    <- NA_character_
        out_aa_effect[i] <- if (indel_len %% 3L != 0L) "frameshift" else "inframe_indel"
      }
    } else {
      out_codon_ref[i] <- NA_character_
      out_codon_alt[i] <- NA_character_
      out_aa_ref[i]    <- NA_character_
      out_aa_alt[i]    <- NA_character_
      out_aa_effect[i] <- NA_character_
    }

    ## ---- 6c. RNA structure (rRNA for rRNA-region mutations; mRNA for CDS) ----
    region_label <- ann_i$region

    if (region_label %in% names(rrna_folds)) {
      fl <- rrna_folds[[region_label]]
      f1 <- fl$ref;  f2 <- fl$alt
      if (!is.null(f1) && !is.null(f2)) {
        out_rrna_mfe_ref[i]        <- f1$mfe
        out_rrna_mfe_alt[i]        <- f2$mfe
        out_rrna_delta_mfe[i]      <- f2$mfe - f1$mfe
        out_rrna_delta_pairs[i]    <- as.integer(f2$n_pairs    - f1$n_pairs)
        out_rrna_delta_hairpins[i] <- as.integer(f2$n_hairpins - f1$n_hairpins)
      }
    }

    gene_label <- out_canonical_gene[i]
    if (!is.na(gene_label) && gene_label %in% names(mrna_folds)) {
      fl <- mrna_folds[[gene_label]]
      f_ref <- fl$ref

      if (mrna_fold == "per_snp") {
        ## Build a single-mutant CDS and fold it
        cds_ref_seq <- .extract_cds_seq(seq1_str, fl$feat_row)
        mut_seq     <- .introduce_mutation_in_cds(cds_ref_seq, pos_r,
                                                   rnt, ant, fl$feat_row)
        f_alt <- if (!is.null(mut_seq))
                   tryCatch(fold_rna(mut_seq, temperature = temperature),
                            error = function(e) NULL)
                 else NULL
      } else {
        f_alt <- fl$alt   # pre-computed joint alt fold
      }

      if (!is.null(f_ref) && !is.null(f_alt)) {
        out_mrna_mfe_ref[i]        <- f_ref$mfe
        out_mrna_mfe_alt[i]        <- f_alt$mfe
        out_mrna_delta_mfe[i]      <- f_alt$mfe - f_ref$mfe
        out_mrna_delta_pairs[i]    <- as.integer(f_alt$n_pairs    - f_ref$n_pairs)
        out_mrna_delta_hairpins[i] <- as.integer(f_alt$n_hairpins - f_ref$n_hairpins)
      }
    }

    ## ---- 6d. MDP ORF impact ----
    if (isTRUE(scan_mdps)) {

      ## ORFs overlapping this ref position
      if (!is.null(mdp_ref) && nrow(mdp_ref) > 0L && !is.na(pos_r)) {
        ref_over <- mdp_ref[mdp_ref$start <= pos_r & mdp_ref$end >= pos_r, , drop = FALSE]
      } else {
        ref_over <- NULL
      }

      ## ORFs overlapping this alt position
      if (!is.null(mdp_alt) && nrow(mdp_alt) > 0L && !is.na(pos_a)) {
        alt_over <- mdp_alt[mdp_alt$start <= pos_a & mdp_alt$end >= pos_a, , drop = FALSE]
      } else {
        alt_over <- NULL
      }

      n_ref_over <- if (is.null(ref_over)) 0L else nrow(ref_over)
      n_alt_over <- if (is.null(alt_over)) 0L else nrow(alt_over)

      out_mdp_n_ref[i] <- as.integer(n_ref_over)
      out_mdp_n_alt[i] <- as.integer(n_alt_over)

      ## Classify: lost (in ref, not alt), gained (in alt, not ref), changed
      ## Use (strand, frame, length_nt) as an ORF "identity" key; compare proteins
      lost    <- NULL
      gained  <- NULL
      changed <- NULL

      if (n_ref_over > 0L || n_alt_over > 0L) {
        make_key <- function(df) {
          if (is.null(df) || nrow(df) == 0L) return(character(0L))
          paste0(df$strand, "_", df$frame, "_", df$start, "_", df$end)
        }
        ref_keys <- make_key(ref_over)
        alt_keys <- make_key(alt_over)

        ## Lost: in ref but protein absent from alt
        lost_keys <- setdiff(ref_keys, alt_keys)
        if (length(lost_keys) > 0L)
          lost <- ref_over[ref_keys %in% lost_keys, , drop = FALSE]

        ## Gained: in alt but not ref
        gained_keys <- setdiff(alt_keys, ref_keys)
        if (length(gained_keys) > 0L)
          gained <- alt_over[alt_keys %in% gained_keys, , drop = FALSE]

        ## Changed: same key but different protein
        shared_keys <- intersect(ref_keys, alt_keys)
        if (length(shared_keys) > 0L) {
          ref_shared <- ref_over[ref_keys %in% shared_keys, , drop = FALSE]
          alt_shared <- alt_over[alt_keys %in% shared_keys, , drop = FALSE]
          changed_mask <- vapply(seq_len(nrow(ref_shared)), function(j) {
            ref_shared$protein_sequence[j] != alt_shared$protein_sequence[j]
          }, logical(1L))
          if (any(changed_mask))
            changed <- ref_shared[changed_mask, , drop = FALSE]
        }
      }

      out_mdp_effect[i] <- .mdp_effect_string(lost, gained, changed)
    }
  }

  ## =========================================================================
  ## 7. Assemble output data.frame
  ## =========================================================================
  result <- data.frame(
    pos_ref             = muts$pos_ref,
    pos_alt             = muts$pos_alt,
    ref_nt              = muts$ref_nt,
    alt_nt              = muts$alt_nt,
    mutation_type       = muts$mutation_type,
    region              = out_region,
    canonical_gene      = out_canonical_gene,
    codon_ref           = out_codon_ref,
    codon_alt           = out_codon_alt,
    aa_ref              = out_aa_ref,
    aa_alt              = out_aa_alt,
    aa_effect           = out_aa_effect,
    rrna_mfe_ref        = out_rrna_mfe_ref,
    rrna_mfe_alt        = out_rrna_mfe_alt,
    rrna_delta_mfe      = out_rrna_delta_mfe,
    rrna_delta_pairs    = out_rrna_delta_pairs,
    rrna_delta_hairpins = out_rrna_delta_hairpins,
    mrna_mfe_ref        = out_mrna_mfe_ref,
    mrna_mfe_alt        = out_mrna_mfe_alt,
    mrna_delta_mfe      = out_mrna_delta_mfe,
    mrna_delta_pairs    = out_mrna_delta_pairs,
    mrna_delta_hairpins = out_mrna_delta_hairpins,
    mdp_n_ref           = out_mdp_n_ref,
    mdp_n_alt           = out_mdp_n_alt,
    mdp_effect          = out_mdp_effect,
    stringsAsFactors    = FALSE
  )

  ## Replace empty strings with NA for consistency
  chr_cols <- c("canonical_gene", "codon_ref", "codon_alt", "aa_ref", "aa_alt",
                "aa_effect", "mdp_effect")
  for (col in chr_cols) {
    result[[col]][result[[col]] == ""] <- NA_character_
  }

  ## Sort by pos_ref (NA insertions go to end)
  result <- result[order(result$pos_ref, na.last = TRUE), ]
  rownames(result) <- NULL

  ## Attach metadata and raw data as attributes for use by analyze_protein_impacts()
  attr(result, "label1")   <- label1
  attr(result, "label2")   <- label2
  attr(result, "seq1_str") <- seq1_str
  attr(result, "seq2_str") <- seq2_str
  attr(result, "features") <- features
  attr(result, "mdp_ref")  <- mdp_ref
  attr(result, "mdp_alt")  <- mdp_alt

  class(result) <- c("mito_comparison", "data.frame")
  result
}


# ── S3 print method ───────────────────────────────────────────────────────────

#' Print a mito_comparison Object
#'
#' Displays a brief summary header followed by the full mutation table.
#'
#' @param x A `mito_comparison` object returned by [compare_mitogenomes()].
#' @param ... Additional arguments passed to [print.data.frame()].
#'
#' @return `x`, invisibly.
#'
#' @export
print.mito_comparison <- function(x, ...) {
  label1 <- attr(x, "label1")
  label2 <- attr(x, "label2")
  if (is.null(label1)) label1 <- "ref"
  if (is.null(label2)) label2 <- "alt"

  n_total  <- nrow(x)
  n_cds    <- sum(grepl("^CDS:", x$region), na.rm = TRUE)
  n_rrna   <- sum(grepl("rRNA$", x$region),  na.rm = TRUE)
  n_mdp    <- if ("mdp_effect" %in% names(x))
                sum(!is.na(x$mdp_effect), na.rm = TRUE) else NA_integer_

  cat(sprintf(
    "── mito_comparison: %s vs %s ────────────────────────────────\n",
    label1, label2
  ))
  cat(sprintf("   Mutations total : %d\n", n_total))
  cat(sprintf("   In CDS          : %d\n", n_cds))
  cat(sprintf("   In rRNA         : %d\n", n_rrna))
  if (!is.na(n_mdp)) {
    cat(sprintf("   Affecting MDPs  : %d\n", n_mdp))
  }
  cat("────────────────────────────────────────────────────────────\n")

  NextMethod()
  invisible(x)
}
