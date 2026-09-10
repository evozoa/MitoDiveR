# ---------------------------------------------------------------------------
# plot_orf_depth() — taxon-first per-reading-frame ORF "depth" across species.
#
# Given a taxon name, fetch every mitogenome, scan ORFs in each, align the
# genomes to a reference, map ORF coordinates onto reference bp, and draw seven
# stacked depth tracks: one per signed reading frame (+1,+2,+3,-1,-2,-3) plus a
# 7th summing the six.  A peak in a frame track marks a reference locus where
# many species carry an ORF in that frame — a candidate conserved noncanonical
# ORF.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Internal: fetch all mitogenomes for a taxon as GenBank records.
# Mirrors the search/filter/fetch block in compare_clade_orfeomes() but keeps
# the organism column (needed to deduplicate species).  Returns a list with
# `meta` (acc, organism, length), `gb_list`, `features_by_acc`, and `seqs`.
# ---------------------------------------------------------------------------
.fetch_clade_genbank <- function(taxon, max_genomes, refseq_only, min_length,
                                 api_key = NULL) {
  if (!is.null(api_key)) rentrez::set_entrez_key(api_key)
  term <- paste(c(sprintf("%s[ORGN]", taxon), "mitochondrion[FILT]",
                  if (refseq_only) "srcdb_refseq[PROP]"), collapse = " AND ")
  message(sprintf("Searching NCBI nucleotide: %s", term))
  ids <- rentrez::entrez_search(db = "nucleotide", term = term,
                                retmax = max_genomes)$ids
  if (length(ids) == 0L)
    stop(sprintf("No NCBI nucleotide hits for taxon '%s'.", taxon), call. = FALSE)
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

  message("Fetching GenBank records (sequence + features) for each genome...")
  gb_list <- fetch_mito_genbank(meta$acc)
  meta    <- meta[meta$acc %in% names(gb_list), , drop = FALSE]
  features_by_acc <- lapply(gb_list[meta$acc], `[[`, "features")
  seqs <- Biostrings::DNAStringSet(vapply(gb_list[meta$acc],
                                          function(g) as.character(g$sequence),
                                          character(1L)))
  names(seqs) <- meta$acc
  list(meta = meta, gb_list = gb_list,
       features_by_acc = features_by_acc, seqs = seqs)
}

# ---------------------------------------------------------------------------
# Internal (pure, offline-testable): build the reference-anchored depth matrix.
#
# aligned    : named character vector of equal-length aligned sequences
#              (gaps as "-"); names identify each genome (seq_id).
# orfs       : data.frame with seq_id, frame (signed 1,2,3,-1,-2,-3), start,
#              end (1-based, in each genome's OWN ungapped coordinates), and
#              optionally wraps_around (logical).
# reference  : name in `aligned` whose ungapped coordinates define the bp axis.
# depth_unit : "species" -> each group contributes 0/1 per position per frame;
#              "orf"     -> every ORF adds 1 (raw pileup).
# group_of   : named vector mapping seq_id -> group label (species).  Only used
#              for "species"; defaults to one group per seq_id.
#
# Returns a numeric matrix [ref_len x 7] with columns
# "+1","+2","+3","-1","-2","-3","sum"; row i = reference bp i.
# ---------------------------------------------------------------------------
.orf_depth_profile <- function(aligned, orfs,
                               reference  = names(aligned)[1L],
                               depth_unit = c("species", "orf"),
                               group_of   = NULL) {
  depth_unit <- match.arg(depth_unit)
  frames     <- c(1L, 2L, 3L, -1L, -2L, -3L)
  frame_lab  <- c("+1", "+2", "+3", "-1", "-2", "-3")
  if (length(aligned) == 0L) stop("'aligned' is empty.", call. = FALSE)
  if (!reference %in% names(aligned))
    stop("'reference' is not among the aligned sequences.", call. = FALSE)

  L <- nchar(aligned[[1L]])
  # genome bp -> alignment column, per genome
  col_of_pos <- lapply(aligned, function(a) which(strsplit(a, "", fixed = TRUE)[[1L]] != "-"))
  names(col_of_pos) <- names(aligned)

  # alignment column -> reference bp (insertion columns fold onto preceding bp)
  ref_nongap   <- strsplit(aligned[[reference]], "", fixed = TRUE)[[1L]] != "-"
  ref_bp_at_col <- cumsum(ref_nongap)               # length L; 0 before ref starts
  ref_len      <- sum(ref_nongap)
  if (ref_len == 0L) stop("Reference has zero ungapped length.", call. = FALSE)

  if (is.null(group_of))
    group_of <- stats::setNames(names(aligned), names(aligned))

  # alignment columns covered by one ORF (handles circular wrap-around)
  orf_cols <- function(seq_id, start, end, wraps) {
    cp   <- col_of_pos[[seq_id]]
    if (is.null(cp) || length(cp) == 0L) return(integer(0))
    glen <- length(cp)
    start <- max(1L, min(glen, start)); end <- max(1L, min(glen, end))
    pos <- if (isTRUE(wraps) && end < start) c(start:glen, 1:end) else start:end
    cp[pos]
  }

  wraps_vec <- if ("wraps_around" %in% names(orfs)) orfs$wraps_around else rep(FALSE, nrow(orfs))

  M <- matrix(0, nrow = ref_len, ncol = 7L,
              dimnames = list(NULL, c(frame_lab, "sum")))

  for (fi in seq_along(frames)) {
    fr   <- frames[fi]
    rows <- which(orfs$frame == fr)
    col_depth <- numeric(L)
    if (length(rows)) {
      if (depth_unit == "orf") {
        for (r in rows) {
          cols <- orf_cols(orfs$seq_id[r], orfs$start[r], orfs$end[r], wraps_vec[r])
          if (length(cols)) col_depth[cols] <- col_depth[cols] + 1
        }
      } else {                                   # "species": 0/1 per group
        grp <- group_of[orfs$seq_id[rows]]
        for (g in unique(grp)) {
          mask <- logical(L)
          for (r in rows[grp == g]) {
            cols <- orf_cols(orfs$seq_id[r], orfs$start[r], orfs$end[r], wraps_vec[r])
            if (length(cols)) mask[cols] <- TRUE
          }
          col_depth <- col_depth + as.integer(mask)
        }
      }
    }
    # collapse alignment columns onto reference bp (max within each bp)
    keep <- ref_bp_at_col >= 1L
    if (any(keep)) {
      rb  <- ref_bp_at_col[keep]
      dv  <- col_depth[keep]
      agg <- tapply(dv, rb, max)
      M[as.integer(names(agg)), fi] <- agg
    }
  }
  M[, "sum"] <- rowSums(M[, seq_along(frames), drop = FALSE])
  M
}

#' Per-reading-frame ORF depth across a taxon
#'
#' Taxon-first companion to [plot_orfeome_circos()].  Fetches every mitogenome
#' in a taxon, scans each for ORFs, aligns the genomes to a reference, maps ORF
#' coordinates onto reference base pairs, and draws **seven stacked depth
#' tracks**: one per signed reading frame (`+1, +2, +3, -1, -2, -3`) and a
#' seventh summing the six.  Each track shows, at every reference position, how
#' many species (or ORFs) carry an ORF in that frame there.  A peak in a frame
#' track flags a reference locus where the same noncanonical ORF recurs across
#' species — a conserved-ORF candidate.
#'
#' Cross-species positions are made comparable by a reference-anchored
#' alignment: all genomes are aligned with MAFFT and each ORF's coordinates are
#' projected onto the chosen reference genome's base pairs (columns inserted
#' relative to the reference fold onto the preceding reference position).  The
#' depth profile is therefore in homologous coordinates, not raw bp.
#'
#' @param taxon Character scalar.  Taxon name searched as `<taxon>[ORGN]`
#'   (e.g. `"Sebastes"`, `"Cyprinidae"`).
#' @param genetic_code Genetic-code identifier passed to [scan_orfs()].
#'   Required (no default).  Use `"SGC0"` for a noncanonical/MDP scan or
#'   `"SGC1"` for canonical vertebrate-mitochondrial genes.
#' @param reference Optional accession (must be one of the fetched genomes)
#'   whose coordinates define the depth axis.  `NULL` (default) uses the
#'   longest genome.
#' @param depth_unit `"species"` (default) counts each species at most once per
#'   position per frame (genomes are deduplicated by organism); `"orf"` is a
#'   raw ORF pileup where every ORF adds one and no deduplication is done.
#' @param noncanonical_only Logical.  If `TRUE` (default), drop ORFs whose 5'
#'   end matches an annotated CDS (via the same test as
#'   [compare_clade_orfeomes()]), keeping the maps focused on noncanonical ORFs.
#' @param max_genomes Integer.  Maximum genomes to fetch.  Default `20L`.
#' @param refseq_only Logical.  Restrict the NCBI search to RefSeq
#'   (`srcdb_refseq[PROP]`).  Default `TRUE`.
#' @param min_length Integer.  Minimum genome length (nt) to retain.
#'   Default `13000L`.
#' @param min_orf_length Integer.  Minimum ORF length (nt) for [scan_orfs()].
#'   Default `60L`.
#' @param start_codons Character vector of recognised start codons.
#' @param annotate_genes Logical.  Draw a reference gene ruler (CDS/rRNA/tRNA)
#'   beneath the bottom track.  Default `TRUE`.
#' @param frame_colors Optional named character vector of six colours keyed by
#'   `"1","2","3","-1","-2","-3"`.  `NULL` uses the [plot_orfeome_circos()]
#'   palette.
#' @param mafft_bin Path to the MAFFT binary.  Default `"mafft"` (must be on
#'   `PATH`).
#' @param mafft_opts Character options passed to MAFFT.  Default `"--auto"`.
#' @param title Optional plot title.  `NULL` derives one from the inputs.
#' @param file Optional output path ending in `.png` or `.pdf`.  `NULL`
#'   (default) draws to the active device.
#' @param api_key Optional NCBI API key.
#'
#' @return Invisibly, a list with components `depth` (a numeric matrix
#'   `reference_length x 7`; columns `+1,+2,+3,-1,-2,-3,sum`), `reference`
#'   (accession used), `n_species`, `meta` (the genome table), and `params`.
#'
#' @section External tools:
#' Requires MAFFT on `PATH` (as for [compare_clade_orfeomes()]) and network
#' access to NCBI.
#'
#' @seealso [plot_orfeome_circos()], [compare_clade_orfeomes()],
#'   [find_clade_conserved_orfs()], [scan_orfs()], [fetch_mito_genbank()].
#'
#' @examples
#' \dontrun{
#' # Seven depth tracks of noncanonical ORFs across rockfishes, on screen
#' plot_orf_depth("Sebastes", genetic_code = "SGC0")
#'
#' # Raw ORF pileup, saved to a file, anchored on a chosen reference
#' res <- plot_orf_depth("Sebastes", genetic_code = "SGC0",
#'                       depth_unit = "orf", reference = "NC_005450.1",
#'                       file = "sebastes_orf_depth.png")
#' str(res$depth)
#' }
#'
#' @importFrom Biostrings DNAStringSet
#' @importFrom rentrez entrez_search entrez_summary set_entrez_key
#' @importFrom stats setNames
#' @importFrom graphics par layout plot rect axis mtext text abline box
#' @importFrom grDevices png pdf dev.off
#' @export
plot_orf_depth <- function(taxon,
                           genetic_code      = .stop_no_code("genetic_code"),
                           reference         = NULL,
                           depth_unit        = c("species", "orf"),
                           noncanonical_only = TRUE,
                           max_genomes       = 20L,
                           refseq_only       = TRUE,
                           min_length        = 13000L,
                           min_orf_length    = 60L,
                           start_codons      = c("ATG","GTG","ATA","ATT","ATC"),
                           annotate_genes    = TRUE,
                           frame_colors      = NULL,
                           mafft_bin         = "mafft",
                           mafft_opts        = "--auto",
                           title             = NULL,
                           file              = NULL,
                           api_key           = NULL) {
  depth_unit <- match.arg(depth_unit)
  if (!is.character(taxon) || length(taxon) != 1L || nchar(taxon) == 0L)
    stop("'taxon' must be a single non-empty character string.", call. = FALSE)
  if (Sys.which(mafft_bin) == "")
    stop("MAFFT binary '", mafft_bin, "' not found on PATH. ",
         "Install MAFFT (https://mafft.cbrc.jp/) or pass mafft_bin = '/path/to/mafft'.",
         call. = FALSE)

  ## ---- 1. fetch genomes ----------------------------------------------------
  fc   <- .fetch_clade_genbank(taxon, max_genomes, refseq_only, min_length, api_key)
  meta <- fc$meta; seqs <- fc$seqs

  ## ---- 2. scan ORFs, optional noncanonical filter --------------------------
  message(sprintf("Scanning ORFs under code %s (min %d nt)...",
                  genetic_code, min_orf_length))
  orfs <- scan_orfs(seqs, min_orf_length = min_orf_length,
                    genetic_code = genetic_code, start_codons = start_codons,
                    both_strands = TRUE, circular = TRUE)
  if (is.null(orfs) || nrow(orfs) == 0L) {
    message("No ORFs found; nothing to plot."); return(invisible(NULL))
  }
  if (isTRUE(noncanonical_only)) {
    orfs <- orfs[!.orf_is_canonical(orfs, fc$features_by_acc), , drop = FALSE]
    message(sprintf("  %d noncanonical ORFs retained.", nrow(orfs)))
    if (nrow(orfs) == 0L) {
      message("No noncanonical ORFs; nothing to plot."); return(invisible(NULL))
    }
  }

  ## ---- 3. reference + whole-genome alignment -------------------------------
  if (is.null(reference)) {
    reference <- meta$acc[which.max(meta$length)]
  } else if (!reference %in% names(seqs)) {
    stop("'reference' (", reference, ") is not among the fetched genomes: ",
         paste(names(seqs), collapse = ", "), call. = FALSE)
  }
  message(sprintf("Aligning %d genome(s) with MAFFT (reference %s)...",
                  length(seqs), reference))
  aligned <- .run_mafft(stats::setNames(as.character(seqs), names(seqs)),
                        mafft_bin, mafft_opts)

  ## ---- 4. depth profile in reference coordinates ---------------------------
  group_of  <- stats::setNames(meta$organism, meta$acc)
  group_of[is.na(group_of)] <- names(group_of)[is.na(group_of)]   # NA organism -> own group
  n_species <- length(unique(group_of))
  depth <- .orf_depth_profile(aligned, orfs, reference = reference,
                              depth_unit = depth_unit, group_of = group_of)

  ## ---- 5. draw -------------------------------------------------------------
  if (is.null(frame_colors))
    frame_colors <- c("1" = "#08519c", "2" = "#3182bd", "3" = "#6baed6",
                      "-1" = "#a50f15", "-2" = "#de2d26", "-3" = "#fb6a4a")
  if (is.null(title))
    title <- sprintf("%s ORF depth — %d species, ref %s, code %s, by %s",
                     taxon, n_species, reference, genetic_code,
                     if (depth_unit == "species") "species" else "ORF")
  genes <- if (isTRUE(annotate_genes)) fc$features_by_acc[[reference]] else NULL

  if (!is.null(file)) {
    ext <- tolower(tools::file_ext(file))
    if (ext == "png")      grDevices::png(file, width = 2200, height = 1600, res = 200)
    else if (ext == "pdf") grDevices::pdf(file, width = 11, height = 8)
    else stop("'file' must end in .png or .pdf", call. = FALSE)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  .draw_orf_depth(depth, frame_colors, title, depth_unit, genes,
                  has_ruler = !is.null(genes) && nrow(genes) > 0L)

  invisible(list(depth = depth, reference = reference, n_species = n_species,
                 meta = meta,
                 params = list(taxon = taxon, genetic_code = genetic_code,
                               depth_unit = depth_unit,
                               noncanonical_only = noncanonical_only)))
}

# ---------------------------------------------------------------------------
# Internal: render the seven stacked tracks (+ optional reference gene ruler).
# ---------------------------------------------------------------------------
.draw_orf_depth <- function(depth, frame_colors, title, depth_unit,
                            genes = NULL, has_ruler = FALSE) {
  frames    <- c(1L, 2L, 3L, -1L, -2L, -3L)
  frame_lab <- c("+1", "+2", "+3", "-1", "-2", "-3")
  ref_len   <- nrow(depth)
  x         <- seq_len(ref_len)
  ymax_fr   <- max(1, max(depth[, frame_lab]))
  ymax_sum  <- max(1, max(depth[, "sum"]))
  ylab      <- if (depth_unit == "species") "species" else "ORFs"

  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op), add = TRUE)
  n_rows  <- 7L + as.integer(has_ruler)
  heights <- c(rep(1, 6L), 1.3, if (has_ruler) 0.9 else NULL)
  graphics::layout(matrix(seq_len(n_rows), ncol = 1L), heights = heights)
  graphics::par(oma = c(3.5, 1, 3, 1), mar = c(0.4, 5, 0.4, 1))

  draw_track <- function(y, col, lab, ymax, bottom = FALSE) {
    graphics::plot(x, y, type = "h", col = col, lwd = 1, xaxs = "i", yaxs = "i",
                   xlim = c(1, ref_len), ylim = c(0, ymax * 1.08),
                   axes = FALSE, xlab = "", ylab = "")
    graphics::axis(2, at = c(0, ymax), las = 1, cex.axis = 0.7, tcl = -0.3,
                   mgp = c(2, 0.5, 0))
    graphics::mtext(lab, side = 2, line = 2.6, las = 1, cex = 0.8, font = 2)
    if (bottom)
      graphics::axis(1, cex.axis = 0.8, mgp = c(2, 0.6, 0))
    graphics::box(col = "grey85")
  }

  for (i in seq_along(frames))
    draw_track(depth[, frame_lab[i]], frame_colors[as.character(frames[i])],
               frame_lab[i], ymax_fr)
  draw_track(depth[, "sum"], "grey25", "Σ", ymax_sum,
             bottom = !has_ruler)

  if (has_ruler) {
    type_col <- c(CDS = "#4575b4", rRNA = "#fd8d3c", tRNA = "grey75",
                  `D-loop` = "#ffd700", rep_origin = "#ffd700")
    keep <- genes$type %in% names(type_col)
    genes <- genes[keep, , drop = FALSE]
    graphics::plot(x, rep(0, ref_len), type = "n", xaxs = "i",
                   xlim = c(1, ref_len), ylim = c(0, 1), axes = FALSE,
                   xlab = "", ylab = "")
    graphics::mtext("genes", side = 2, line = 2.6, las = 1, cex = 0.8, font = 2)
    if (nrow(genes)) {
      st <- pmax(1L, pmin(ref_len, genes$start))
      en <- pmax(1L, pmin(ref_len, genes$end))
      gc <- type_col[genes$type]
      graphics::rect(st, 0.25, en, 0.75, col = gc, border = "grey30", lwd = 0.5)
      lab_fld <- if ("gene" %in% names(genes)) "gene" else
                 if ("product" %in% names(genes)) "product" else NULL
      if (!is.null(lab_fld)) {
        show <- genes$type %in% c("CDS", "rRNA") & !is.na(genes[[lab_fld]])
        if (any(show))
          graphics::text((st[show] + en[show]) / 2, 0.5, genes[[lab_fld]][show],
                         cex = 0.45, srt = 90, col = "grey20")
      }
    }
    graphics::axis(1, cex.axis = 0.8, mgp = c(2, 0.6, 0))
  }

  graphics::mtext(sprintf("reference position (bp)  —  depth = number of %s with an ORF",
                          ylab), side = 1, outer = TRUE, line = 2, cex = 0.8)
  graphics::mtext(title, side = 3, outer = TRUE, line = 0.8, cex = 1.0, font = 2)
  invisible(NULL)
}
