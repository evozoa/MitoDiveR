# ============================================================================
# compare_clade_orfeomes — taxon-first ORFeome comparison.
# For a clade, fetch mitogenomes, scan ORFs, optionally restrict to
# noncanonical ORFs, cluster across genomes into families, align each family
# (MAFFT), concatenate to a supermatrix, and return a pairwise identity matrix
# (+ Jaccard presence/absence and an hclust tree).
# ============================================================================

#' Compare ORFeomes across a clade (taxon-first)
#'
#' Builds a pairwise genome × genome **identity matrix** from a concatenated
#' protein alignment of ORF families shared across a clade's mitogenomes, and
#' (optionally) a complementary **Jaccard presence/absence** matrix.  The
#' identity matrix measures within-shared-family sequence divergence; the
#' Jaccard matrix measures lineage-specific ORF family gain/loss.  Defaults to
#' restricting to **noncanonical ORFs** (i.e. ORFs whose 5' position does not
#' coincide with a CDS feature in the same strand), so the comparison reflects
#' the non-canonical coding repertoire rather than the standard 13 mt CDSs.
#'
#' Pipeline: NCBI search for `taxon` mitogenomes → [fetch_mito_genbank()] for
#' sequence + features → [scan_orfs()] per genome under the (required)
#' `genetic_code` → noncanonical filter → [find_conserved_orfs()] clustering
#' into ORF families → per-family MAFFT alignment → concatenated supermatrix →
#' pairwise identity + Jaccard.
#'
#' @param taxon NCBI organism / taxon name (e.g. `"Hominidae"`).
#' @param genetic_code Required (no default) — Biostrings genetic-code id, e.g.
#'   `"SGC1"` (vertebrate mitochondrial) or `"SGC0"` (standard).
#' @param max_genomes Integer cap on mitogenomes fetched.  Default `20L`.
#' @param refseq_only Restrict the NCBI search to RefSeq?  Default `TRUE`
#'   (one representative per species; balances species representation).
#' @param min_length Integer.  Minimum genome length (nt) to retain.
#'   Default `13000L`.
#' @param min_orf_length Integer.  Passed to [scan_orfs()].  Default `60L`
#'   (catches MDP-scale ORFs).
#' @param start_codons Character vector of recognised start codons passed to
#'   [scan_orfs()].
#' @param noncanonical_only Logical.  If `TRUE` (default), drop ORFs whose 5'
#'   position coincides with the 5' end of a CDS feature on the same strand
#'   (i.e. the canonical 13 mt CDSs).  All other ORFs — minus-strand,
#'   frame-shifted same-strand, rRNA/D-loop-internal — are retained.
#' @param collapse_nested Logical.  If `TRUE` (default), call
#'   [collapse_nested_orfs()] before cross-genome clustering so each (frame,
#'   stop) per genome contributes a single longest ORF rather than every
#'   nested internal start codon.  Cuts the clustering input from ~10³ to ~10²
#'   ORFs per genome and is essential for tractable runtime; turn off only if
#'   you need every nested start.
#' @param min_identity Numeric percent aa identity passed to
#'   [find_conserved_orfs()] for cross-genome family clustering.  Default `70`.
#' @param gap_handling How to treat gap-vs-residue positions in pairwise
#'   identity: `"pairwise_complete"` (default; ignore positions where either
#'   sequence has a gap — measures sequence divergence within shared ORF
#'   regions) or `"mismatch"` (count gaps as mismatches — folds in
#'   lineage-specific gain/loss).
#' @param also_jaccard Logical.  Also compute a genome × genome Jaccard
#'   distance matrix on family presence/absence?  Default `TRUE`.
#' @param mafft_bin Path to the MAFFT binary.  Default `"mafft"` (must be on
#'   `PATH`).
#' @param mafft_opts Character options passed to MAFFT.  Default `"--auto"`.
#' @param api_key Optional NCBI API key.
#'
#' @return A list with components:
#'   \describe{
#'     \item{`taxon`, `accessions`, `metadata`}{Inputs / record metadata.}
#'     \item{`all_orfs`}{Combined [scan_orfs()] table across genomes.}
#'     \item{`noncanonical_orfs`}{Filtered subset used downstream.}
#'     \item{`orf_families`}{Per-cluster summary (`cluster_id`, `n_genomes`,
#'       `n_orfs`, `longest_aa`, `longest_len`).}
#'     \item{`presence_absence`}{Integer matrix `[family x genome]` of 0/1.}
#'     \item{`msa_concat`}{[Biostrings::AAStringSet] of length `n_genomes` —
#'       the concatenated per-family alignments.}
#'     \item{`identity_matrix`}{Genome × genome `% identity` (symmetric;
#'       diagonal = `100`).}
#'     \item{`jaccard_matrix`}{Genome × genome Jaccard distance on
#'       presence/absence (only when `also_jaccard = TRUE`; `0` =
#'       identical family sets, `1` = disjoint).}
#'     \item{`gap_fraction_matrix`}{Fraction of supermatrix positions where
#'       at least one of the pair carries a gap — the lineage-specific
#'       gain/loss signal not captured by `pairwise_complete` identity.}
#'     \item{`tree`}{`hclust` on `as.dist(1 - identity_matrix/100)` (average
#'       linkage).}
#'   }
#'
#' @details
#' Requires the **MAFFT** binary on `PATH` for the per-family MSA step.
#' Families present in only one genome (singletons) are produced by the
#' clustering step *only* when `min_identity` would group them — by default
#' [find_conserved_orfs()] keeps families with `>= 2` genomes, so singletons
#' are not included; the comparison thus measures shared-ORF divergence.
#'
#' @seealso [find_clade_conserved_orfs()], [find_conserved_orfs()],
#'   [scan_orfs()], [fetch_mito_genbank()]
#'
#' @examples
#' \dontrun{
#' res <- compare_clade_orfeomes(
#'   taxon          = "Hominidae",
#'   genetic_code   = "SGC1",
#'   refseq_only    = TRUE,
#'   min_orf_length = 60
#' )
#' res$identity_matrix
#' plot(res$tree, main = "Hominidae ORFeome (noncanonical) — identity tree")
#' }
#'
#' @importFrom Biostrings AAStringSet readAAStringSet
#' @importFrom stats as.dist hclust
#' @export
compare_clade_orfeomes <- function(
    taxon,
    genetic_code      = .stop_no_code("genetic_code"),
    max_genomes       = 20L,
    refseq_only       = TRUE,
    min_length        = 13000L,
    min_orf_length    = 60L,
    start_codons      = c("ATG","GTG","ATA","ATT","ATC"),
    noncanonical_only = TRUE,
    collapse_nested   = TRUE,
    min_identity      = 70,
    gap_handling      = c("pairwise_complete","mismatch"),
    also_jaccard      = TRUE,
    mafft_bin         = "mafft",
    mafft_opts        = "--auto",
    api_key           = NULL) {

  gap_handling <- match.arg(gap_handling)
  if (Sys.which(mafft_bin) == "")
    stop("MAFFT binary '", mafft_bin, "' not found on PATH. ",
         "Install MAFFT (https://mafft.cbrc.jp/) or pass mafft_bin = '/path/to/mafft'.",
         call. = FALSE)
  if (!is.null(api_key)) rentrez::set_entrez_key(api_key)

  ## ---- 1. NCBI search + metadata + length filter ---------------------------
  term <- paste(c(sprintf("%s[ORGN]", taxon), "mitochondrion[FILT]",
                  if (refseq_only) "srcdb_refseq[PROP]"), collapse = " AND ")
  message(sprintf("Searching NCBI nucleotide: %s", term))
  ids <- rentrez::entrez_search(db = "nucleotide", term = term,
                                retmax = max_genomes)$ids
  if (length(ids) == 0L)
    stop(sprintf("No NCBI nucleotide hits for taxon '%s'.", taxon),
         call. = FALSE)
  s <- rentrez::entrez_summary(db = "nucleotide", id = ids)
  if (length(ids) == 1L) s <- list(s)
  meta <- data.frame(
    acc      = vapply(s, function(x) x$accessionversion, character(1L)),
    organism = vapply(s, function(x) if (!is.null(x$organism)) x$organism else NA_character_, character(1L)),
    title    = vapply(s, function(x) x$title, character(1L)),
    length   = vapply(s, function(x) as.integer(x$slen), integer(1L)),
    stringsAsFactors = FALSE)
  meta <- meta[meta$length >= min_length, , drop = FALSE]
  if (nrow(meta) == 0L)
    stop(sprintf("No mitogenomes >= %d nt for taxon '%s'.", min_length, taxon),
         call. = FALSE)
  message(sprintf("Retained %d mitogenome(s).", nrow(meta)))

  ## ---- 2. Fetch GenBank per accession (sequence + features) ----------------
  message("Fetching GenBank records (sequence + features) for each genome...")
  gb_list <- fetch_mito_genbank(meta$acc)
  meta    <- meta[meta$acc %in% names(gb_list), , drop = FALSE]
  features_by_acc <- lapply(gb_list[meta$acc], `[[`, "features")
  seqs <- Biostrings::DNAStringSet(vapply(gb_list[meta$acc],
                                          function(g) as.character(g$sequence),
                                          character(1L)))
  names(seqs) <- meta$acc

  ## ---- 3. Scan ORFs per genome ---------------------------------------------
  message(sprintf("Scanning ORFs under code %s (min %d nt)...",
                  genetic_code, min_orf_length))
  all_orfs <- scan_orfs(seqs, min_orf_length = min_orf_length,
                        genetic_code = genetic_code, start_codons = start_codons,
                        both_strands = TRUE, circular = TRUE)
  message(sprintf("  %d ORFs total across %d genomes.",
                  nrow(all_orfs), nrow(meta)))

  ## ---- 4. Restrict to noncanonical (drop ORFs matching a CDS 5' end) -------
  if (isTRUE(noncanonical_only)) {
    keep <- !.orf_is_canonical(all_orfs, features_by_acc)
    nc_orfs <- all_orfs[keep, , drop = FALSE]
    message(sprintf("  %d noncanonical ORFs after filtering.", nrow(nc_orfs)))
  } else {
    nc_orfs <- all_orfs
  }
  if (nrow(nc_orfs) == 0L)
    stop("No ORFs left after filtering.", call. = FALSE)

  ## ---- 4b. Collapse nested ORFs per (frame, stop) per genome ---------------
  if (isTRUE(collapse_nested)) {
    nc_split <- split(nc_orfs, nc_orfs$seq_id)
    nc_orfs  <- do.call(rbind, lapply(nc_split, collapse_nested_orfs))
    rownames(nc_orfs) <- NULL
    message(sprintf("  %d noncanonical ORFs after nested-locus collapse.",
                    nrow(nc_orfs)))
  }

  ## ---- 5. Cluster into families across genomes (find_conserved_orfs) -------
  message(sprintf("Clustering into ORF families at >= %g%% aa identity...",
                  min_identity))
  conserved <- find_conserved_orfs(nc_orfs, min_identity = min_identity,
                                   min_sequences = 2 / nrow(meta))
  if (nrow(conserved) == 0L)
    stop("No multi-genome ORF families found at min_identity=", min_identity,
         ". Lower min_identity or relax filters.", call. = FALSE)
  fam_ids <- sort(unique(conserved$cluster_id))
  message(sprintf("  %d ORF families across %d genomes.",
                  length(fam_ids), nrow(meta)))

  ## ---- 6. Per-family MAFFT alignment + concatenation -----------------------
  message(sprintf("Aligning %d families with %s...",
                  length(fam_ids), mafft_bin))
  fam_summary <- vector("list", length(fam_ids))
  per_genome_concat <- setNames(rep("", nrow(meta)), meta$acc)
  pa <- matrix(0L, nrow = length(fam_ids), ncol = nrow(meta),
               dimnames = list(paste0("fam_", fam_ids), meta$acc))
  for (k in seq_along(fam_ids)) {
    fid <- fam_ids[k]
    rows <- conserved[conserved$cluster_id == fid, , drop = FALSE]
    # one representative per genome (longest ORF)
    rep <- do.call(rbind, lapply(split(rows, rows$seq_id), function(d)
      d[which.max(nchar(d$protein_sequence)), , drop = FALSE]))
    seqs_aa <- sub("\\*+$", "", rep$protein_sequence)   # strip trailing stop
    names(seqs_aa) <- rep$seq_id
    pa[k, names(seqs_aa)] <- 1L
    aln <- .run_mafft(seqs_aa, mafft_bin, mafft_opts)
    aln_len <- nchar(aln[1])
    for (g in meta$acc) {
      per_genome_concat[g] <- paste0(per_genome_concat[g],
        if (g %in% names(aln)) aln[g] else strrep("-", aln_len))
    }
    fam_summary[[k]] <- data.frame(
      cluster_id  = fid,
      n_genomes   = length(seqs_aa),
      n_orfs      = nrow(rows),
      aln_len     = aln_len,
      longest_len = max(nchar(seqs_aa)),
      longest_aa  = seqs_aa[[which.max(nchar(seqs_aa))]],
      stringsAsFactors = FALSE)
  }
  fam_summary <- do.call(rbind, fam_summary)
  msa_concat  <- Biostrings::AAStringSet(per_genome_concat)
  names(msa_concat) <- names(per_genome_concat)

  ## ---- 7. Pairwise distances + tree ----------------------------------------
  message("Computing pairwise identity / Jaccard / gap-fraction matrices...")
  ident <- .pairwise_identity_matrix(per_genome_concat, gap_handling)
  gapf  <- .pairwise_gap_fraction(per_genome_concat)
  jacc  <- if (isTRUE(also_jaccard)) .jaccard_distance(pa) else NULL
  tree  <- stats::hclust(stats::as.dist(1 - ident / 100), method = "average")

  list(
    taxon               = taxon,
    accessions          = meta$acc,
    metadata            = meta,
    all_orfs            = all_orfs,
    noncanonical_orfs   = nc_orfs,
    orf_families        = fam_summary,
    presence_absence    = pa,
    msa_concat          = msa_concat,
    identity_matrix     = ident,
    jaccard_matrix      = jacc,
    gap_fraction_matrix = gapf,
    tree                = tree
  )
}

# ---- helpers ---------------------------------------------------------------

# ORF is "canonical" if a CDS feature with the same strand shares its 5'
# coordinate (start for +, end for -).
.orf_is_canonical <- function(orfs, features_by_acc) {
  out <- logical(nrow(orfs))
  for (i in seq_len(nrow(orfs))) {
    feats <- features_by_acc[[orfs$seq_id[i]]]
    if (is.null(feats) || nrow(feats) == 0L) next
    cds <- feats[!is.na(feats$type) & feats$type == "CDS" &
                 !is.na(feats$strand) & feats$strand == orfs$strand[i], ,
                 drop = FALSE]
    if (nrow(cds) == 0L) next
    out[i] <- if (orfs$strand[i] == "+")
                any(cds$start == orfs$start[i])
              else
                any(cds$end == orfs$end[i])
  }
  out
}

# Run MAFFT on a named character vector of AA sequences; return aligned
# sequences as a named character vector (input order).
.run_mafft <- function(seqs, mafft_bin = "mafft", opts = "--auto") {
  n <- length(seqs)
  if (n == 0L) return(setNames(character(0L), character(0L)))
  if (n == 1L || length(unique(seqs)) == 1L) {
    L <- max(nchar(seqs))
    pad <- vapply(seqs, function(s) paste0(s, strrep("-", L - nchar(s))),
                  character(1L))
    return(setNames(pad, names(seqs)))
  }
  tf_in  <- tempfile(fileext = ".fa")
  tf_out <- tempfile(fileext = ".aln.fa")
  on.exit(unlink(c(tf_in, tf_out)), add = TRUE)
  writeLines(paste(paste0(">", names(seqs)), seqs, sep = "\n"), tf_in)
  args <- c(strsplit(opts, " ")[[1L]], tf_in)
  rc <- system2(mafft_bin, args, stdout = tf_out, stderr = FALSE)
  if (rc != 0L) stop("MAFFT failed (exit ", rc, ").", call. = FALSE)
  aln <- Biostrings::readAAStringSet(tf_out)
  as.character(aln)[names(seqs)]
}

# Pairwise % identity from a named character vector of equal-length aligned
# sequences (the concatenated supermatrix, one string per genome).
.pairwise_identity_matrix <- function(concat, gap_handling) {
  g <- names(concat)
  L <- nchar(concat[1L])
  mats <- do.call(rbind, lapply(concat, function(s) strsplit(s, "")[[1L]]))
  rownames(mats) <- g
  out <- matrix(NA_real_, length(g), length(g), dimnames = list(g, g))
  for (i in seq_along(g)) {
    out[i, i] <- 100
    for (j in seq_len(i - 1L)) {
      a <- mats[i, ]; b <- mats[j, ]
      both_aa <- a != "-" & b != "-"
      if (gap_handling == "pairwise_complete") {
        denom <- sum(both_aa)
        pid   <- if (denom > 0) 100 * sum(a == b & both_aa) / denom else NA_real_
      } else {  # mismatch: gaps count as mismatches over all positions
        pid <- 100 * sum(a == b & both_aa) / L
      }
      out[i, j] <- out[j, i] <- pid
    }
  }
  out
}

.pairwise_gap_fraction <- function(concat) {
  g <- names(concat); L <- nchar(concat[1L])
  mats <- do.call(rbind, lapply(concat, function(s) strsplit(s, "")[[1L]]))
  rownames(mats) <- g
  out <- matrix(0, length(g), length(g), dimnames = list(g, g))
  for (i in seq_along(g)) for (j in seq_len(i - 1L)) {
    out[i, j] <- out[j, i] <-
      sum(mats[i, ] == "-" | mats[j, ] == "-") / L
  }
  out
}

# Jaccard DISTANCE on a presence/absence matrix [family x genome].
.jaccard_distance <- function(pa) {
  g <- colnames(pa)
  out <- matrix(0, length(g), length(g), dimnames = list(g, g))
  for (i in seq_along(g)) for (j in seq_len(i - 1L)) {
    a <- pa[, i] == 1L; b <- pa[, j] == 1L
    u <- sum(a | b); inter <- sum(a & b)
    out[i, j] <- out[j, i] <- if (u > 0L) 1 - inter / u else 0
  }
  out
}
