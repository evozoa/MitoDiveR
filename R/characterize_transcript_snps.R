# ============================================================================
# characterize_transcript_snps — transcript-first variant survey.
# Given a noncanonical mtRNA transcript region, enumerate the SNPs observed
# within it for a taxon (human via MITOMAP, other taxa via NCBI) and
# characterize each one's peptide and RNA-structure consequences with
# analyze_snp().
# ============================================================================

#' Characterize SNPs within a noncanonical mtRNA transcript (transcript-first)
#'
#' The transcript-first complement to [analyze_snp()].  Rather than starting
#' from a single variant, you supply a transcript region; the function
#' enumerates the single-nucleotide variants observed within it for a taxon and
#' runs each through [analyze_snp()] to report peptide and (optionally) RNA
#' secondary-structure consequences alongside the variant's population
#' frequency.
#'
#' Variant sources by taxon:
#' \itemize{
#'   \item **Human** (`taxon = "Homo sapiens"`): the MITOMAP GenBank
#'     polymorphisms table (pre-computed allele frequencies from ~60k+ GenBank
#'     mitogenomes), in rCRS coordinates.  No genome download.
#'   \item **Other taxa**: a sample of up to `max_genomes` mitogenomes is
#'     fetched from NCBI, the region is extracted, and substitutions relative to
#'     `reference` are tabulated with their sample frequency.
#' }
#'
#' @param region Either a preset transcript name (`"MOTS-c"`, `"humanin"`) or a
#'   list with `start`, `end`, and optional `strand` (`"+"`/`"-"`) giving the
#'   region's coordinates on `reference`.
#' @param taxon NCBI organism name.  Default `"Homo sapiens"`.
#' @param reference Reference genome: NCBI accession (default `"NC_012920.1"` =
#'   rCRS) or a [Biostrings::DNAString].
#' @param source One of `"auto"` (MITOMAP for *Homo sapiens*, NCBI otherwise),
#'   `"mitomap"`, or `"ncbi"`.
#' @param min_freq Numeric.  Only characterize variants at or above this
#'   population/sample frequency, expressed as a fraction in `[0, 1]` (MITOMAP
#'   percentages are converted to fractions internally).  Default `0` (all).
#' @param max_genomes Integer.  Cap on NCBI mitogenomes sampled (non-human
#'   path).  Default `1000L`.  Set `Inf` for the full set (warns if very large).
#' @param genetic_codes,min_orf_length,rna_structure Passed to [analyze_snp()].
#'   `min_orf_length` defaults to `45L` here so short MDP ORFs (MOTS-c, SHLPs)
#'   are detected.
#' @param analysis_flank Integer.  Each variant is characterized against a
#'   window of the reference spanning the transcript +/- this many nt rather
#'   than the whole genome (the key speedup — full-genome rescans per variant
#'   are prohibitively slow).  Keep it larger than the longest ORF of interest
#'   and >= `analyze_snp()`'s `rna_flank`.  Default `300L`.
#' @param mitomap_file Optional path to a previously downloaded MITOMAP
#'   polymorphisms TSV; if `NULL` it is downloaded once and cached in the
#'   session temp directory.
#' @param ... Further arguments passed to [analyze_snp()].
#'
#' @return A list with `transcript` (resolved region), `source`, `n_genomes`
#'   (NCBI path only), and `variants` — a `data.frame` with one row per
#'   characterized SNP: `pos`, `ref`, `alt`, `freq`, `count`, `consequence`
#'   (`"disrupts_orf"` / `"novel_orf"` / `"substitution"` / `"no_orf"`),
#'   `n_disrupted`, `n_novel`, and (if `rna_structure`) `delta_mfe`,
#'   `bp_distance`, `delta_stem`.  Variants that error in [analyze_snp()] (e.g.
#'   reference-base mismatch) are reported in the `failed` attribute.
#'
#' @details
#' MITOMAP coordinates and `ref` alleles are rCRS; if you supply a non-rCRS
#' `reference` the human path will mismatch.  The NCBI path fetches each record
#' by coordinate (`seq_start`/`seq_stop`) and currently skips records whose
#' region length differs from the reference (indels) — adequate for
#' within-species surveys where indels are rare in coding/rRNA regions.
#'
#' Preset coordinates are rCRS-based and **provisional** — verify against your
#' annotation before publication.  ASncmtRNA-1/-2 are intentionally not
#' presets (their transcript boundaries vary in the literature); supply them as
#' explicit `list(start, end, strand = "-")`.
#'
#' @seealso [analyze_snp()], [fold_rna()], [fetch_sequences()]
#'
#' @examples
#' \dontrun{
#' # Human MOTS-c region — should surface the K14Q variant (m.1382A>C)
#' res <- characterize_transcript_snps("MOTS-c", min_freq = 0.001)
#' res$variants
#'
#' # An explicit antisense transcript in another species
#' characterize_transcript_snps(
#'   list(start = 1671, end = 3229, strand = "-"),
#'   taxon = "Sebastes", reference = "NC_005450.1", max_genomes = 50)
#' }
#'
#' @importFrom Biostrings DNAString subseq
#' @export
characterize_transcript_snps <- function(region,
                                          taxon          = "Homo sapiens",
                                          reference      = "NC_012920.1",
                                          source         = c("auto", "mitomap", "ncbi"),
                                          min_freq       = 0,
                                          max_genomes    = 1000L,
                                          genetic_codes  = c("SGC1", "SGC0"),
                                          min_orf_length = 45L,
                                          analysis_flank = 300L,
                                          rna_structure  = TRUE,
                                          mitomap_file   = NULL,
                                          ...) {
  source <- match.arg(source)
  reg    <- .resolve_region(region)
  rna_strand <- if (identical(reg$strand, "-")) "antisense" else "sense"

  use_mitomap <- source == "mitomap" ||
    (source == "auto" && grepl("homo sapiens", tolower(taxon)))

  ## ---- 1. Fetch the reference once (reused for every variant) ---------------
  ref_seq <- .resolve_reference_once(reference)

  # Window the reference around the transcript so each per-variant analyze_snp()
  # scans a few hundred nt rather than the whole genome (the dominant cost).
  # Variant genome positions are translated to window-local coordinates for the
  # call; reported positions remain genomic.  ORFs spanning the variant farther
  # than `analysis_flank` are clipped at the window edge, so keep it comfortably
  # larger than the longest ORF of interest and >= analyze_snp()'s rna_flank.
  glen       <- length(ref_seq)
  w_start    <- max(1L, reg$start - as.integer(analysis_flank))
  w_end      <- min(glen, reg$end + as.integer(analysis_flank))
  window_seq <- Biostrings::subseq(ref_seq, w_start, w_end)

  ## ---- 2. Enumerate variants in the region ---------------------------------
  if (use_mitomap) {
    vars      <- .mitomap_region_variants(reg, min_freq, mitomap_file)
    n_genomes <- NA_integer_
  } else {
    vars      <- .ncbi_region_variants(reg, taxon, max_genomes, min_freq, ref_seq)
    n_genomes <- attr(vars, "n_genomes")
  }
  if (nrow(vars) == 0L)
    return(list(transcript = reg, source = if (use_mitomap) "mitomap" else "ncbi",
                n_genomes = n_genomes, variants = vars))

  ## ---- 3. Characterize each variant via analyze_snp ------------------------
  rows   <- vector("list", nrow(vars))
  failed <- list()
  for (i in seq_len(nrow(vars))) {
    v   <- vars[i, ]
    lbl <- sprintf("m.%d%s>%s", v$pos - w_start + 1L, v$ref, v$alt)  # window-local pos
    res <- tryCatch(
      analyze_snp(lbl, reference = window_seq, genetic_codes = genetic_codes,
                  min_orf_length = min_orf_length, rna_structure = rna_structure,
                  rna_strand = rna_strand, ...),
      error = function(e) e)

    if (inherits(res, "error")) {
      failed[[lbl]] <- conditionMessage(res)
      next
    }

    s   <- res$summary
    nd  <- sum(s$n_disrupted); nv <- sum(s$n_novel)
    nl  <- sum(s$n_lost);      ng <- sum(s$n_gained)
    cons <- if (nd > 0L) "disrupts_orf"
            else if (nv > 0L) "novel_orf"
            else if (nl > 0L || ng > 0L) "substitution"
            else "no_orf"

    row <- data.frame(
      pos = v$pos, ref = v$ref, alt = v$alt,
      freq = v$freq, count = v$count,
      consequence = cons, n_disrupted = nd, n_novel = nv,
      mitomap_aachange = if (!is.null(v$mitomap_aachange)) v$mitomap_aachange else NA_character_,
      stringsAsFactors = FALSE)
    if (rna_structure && !is.null(res$rna)) {
      row$delta_mfe   <- res$rna$delta_mfe
      row$bp_distance <- res$rna$bp_distance
      row$delta_stem  <- res$rna$delta_stem
    }
    rows[[i]] <- row
  }

  variants <- do.call(rbind, rows)
  if (!is.null(variants))
    variants <- variants[order(-variants$freq), , drop = FALSE]

  out <- list(transcript = reg,
              source      = if (use_mitomap) "mitomap" else "ncbi",
              n_genomes   = n_genomes,
              variants    = variants)
  if (length(failed)) attr(out, "failed") <- failed
  out
}

# ---- helpers ---------------------------------------------------------------

# rCRS (NC_012920.1) transcript presets. PROVISIONAL coordinates — verify.
.transcript_presets <- function() {
  list(
    `MOTS-c`  = list(start = 1343L, end = 1390L, strand = "+", note = "12S rRNA (MT-RNR1)"),
    humanin   = list(start = 2634L, end = 2707L, strand = "+", note = "16S rRNA (MT-RNR2)")
  )
}

.resolve_region <- function(region) {
  if (is.character(region) && length(region) == 1L) {
    p   <- .transcript_presets()
    key <- names(p)[tolower(names(p)) == tolower(region)]
    if (length(key) == 0L)
      stop("Unknown transcript preset '", region, "'. Known presets: ",
           paste(names(p), collapse = ", "),
           ". For anything else (incl. ASncmtRNA-1/-2) pass ",
           "list(start=, end=, strand=).", call. = FALSE)
    r <- p[[key]]; r$name <- key
    return(r)
  }
  if (is.list(region) && !is.null(region$start) && !is.null(region$end)) {
    strand <- if (is.null(region$strand)) "+" else region$strand
    return(list(start  = as.integer(region$start),
                end    = as.integer(region$end),
                strand = strand,
                name   = if (!is.null(region$name)) region$name else "custom",
                note   = NA_character_))
  }
  stop("'region' must be a preset name or list(start, end, strand).", call. = FALSE)
}

.resolve_reference_once <- function(reference) {
  if (methods::is(reference, "DNAString"))    return(reference)
  if (methods::is(reference, "DNAStringSet")) return(reference[[1L]])
  fetch_sequences(accessions = reference)[[1L]]
}

# Download (and cache) the MITOMAP GenBank polymorphisms table.
.get_mitomap_polymorphisms <- function(cache_path = NULL) {
  if (is.null(cache_path))
    cache_path <- file.path(tempdir(), "mitomap_polymorphisms.tsv")
  if (!file.exists(cache_path)) {
    url <- "https://www.mitomap.org/cgi-bin/polymorphisms.cgi?format=text"
    old <- getOption("HTTPUserAgent")
    options(HTTPUserAgent =
      "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 Chrome/120 Safari/537.36")
    on.exit(options(HTTPUserAgent = old), add = TRUE)
    message("Downloading MITOMAP polymorphisms table to ", cache_path, " ...")
    utils::download.file(url, destfile = cache_path, quiet = TRUE)
  }
  utils::read.delim(cache_path, sep = "\t", header = TRUE,
                    stringsAsFactors = FALSE, quote = "")
}

.mitomap_region_variants <- function(reg, min_freq, mitomap_file) {
  df <- .get_mitomap_polymorphisms(mitomap_file)
  df <- df[df$pos >= reg$start & df$pos <= reg$end, , drop = FALSE]
  snv <- nchar(df$ref) == 1L & nchar(df$alt) == 1L &
         grepl("^[ACGT]$", df$ref) & grepl("^[ACGT]$", df$alt) & df$ref != df$alt
  df <- df[snv, , drop = FALSE]
  # MITOMAP gbfreq is a PERCENTAGE of the GenBank set; convert to a fraction so
  # it matches the NCBI path's count / n_sampled.
  freq  <- suppressWarnings(as.numeric(df$gbfreq)) / 100; freq[is.na(freq)] <- 0
  count <- suppressWarnings(as.integer(df$gbcnt))
  keep  <- freq >= min_freq
  data.frame(pos = df$pos[keep], ref = df$ref[keep], alt = df$alt[keep],
             freq = freq[keep], count = count[keep],
             mitomap_aachange = df$aachange[keep], stringsAsFactors = FALSE)
}

# Non-human: sample mitogenomes from NCBI and tabulate substitutions in region.
.ncbi_region_variants <- function(reg, taxon, max_genomes, min_freq, ref_seq) {
  term <- sprintf("%s[Organism] AND mitochondrion[Filter] AND 14000:20000[SLEN]", taxon)
  total <- rentrez::entrez_search(db = "nuccore", term = term, retmax = 0)$count
  if (total == 0L) stop("No mitogenomes found for taxon '", taxon, "'.", call. = FALSE)
  n   <- if (is.infinite(max_genomes)) total else min(as.integer(max_genomes), total)
  ids <- rentrez::entrez_search(db = "nuccore", term = term, retmax = n)$ids

  ref_region <- as.character(Biostrings::subseq(ref_seq, reg$start, reg$end))
  L          <- nchar(ref_region)
  ref_chars  <- strsplit(ref_region, "")[[1L]]
  counts <- new.env(parent = emptyenv())
  n_used <- 0L

  for (id in ids) {
    fa <- tryCatch(
      rentrez::entrez_fetch(db = "nuccore", id = id, rettype = "fasta",
                            retmode = "text",
                            seq_start = reg$start, seq_stop = reg$end),
      error = function(e) NA_character_)
    if (is.na(fa)) next
    lines <- strsplit(fa, "\n")[[1L]]
    s     <- toupper(paste0(lines[-1L], collapse = ""))
    if (nchar(s) != L) next                       # skip indel/length mismatch (v1)
    sc   <- strsplit(s, "")[[1L]]
    diff <- which(sc != ref_chars & sc %in% c("A", "C", "G", "T"))
    for (k in diff) {
      key <- paste(reg$start + k - 1L, ref_chars[k], sc[k], sep = "|")
      assign(key, get0(key, envir = counts, ifnotfound = 0L) + 1L, envir = counts)
    }
    n_used <- n_used + 1L
  }

  keys <- ls(counts)
  if (length(keys) == 0L || n_used == 0L) {
    out <- data.frame(); attr(out, "n_genomes") <- n_used; return(out)
  }
  parts <- do.call(rbind, strsplit(keys, "\\|"))
  cnt   <- vapply(keys, function(k) get(k, envir = counts), integer(1L))
  df <- data.frame(pos = as.integer(parts[, 1L]), ref = parts[, 2L], alt = parts[, 3L],
                   count = cnt, freq = cnt / n_used,
                   mitomap_aachange = NA_character_, stringsAsFactors = FALSE)
  df <- df[df$freq >= min_freq, , drop = FALSE]
  df <- df[order(-df$freq), , drop = FALSE]
  attr(df, "n_genomes") <- n_used
  df
}
