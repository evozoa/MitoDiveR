# ============================================================================
# plot_orfeome_circos — Circos-style plot of a six-frame ORFeome around a
# central ring of nucleotide coordinates and canonical genes.
# ============================================================================

#' Circos plot of a six-frame ORFeome
#'
#' Renders the ORFs returned by [scan_orfs()] / [scan_orfeome()] as six
#' concentric rings — one per reading frame (`+1, +2, +3` then `-1, -2, -3`,
#' outermost to innermost) — arranged *outside* a central ring that shows the
#' nucleotide coordinate axis and the canonical gene annotation.  This is the
#' genome-first companion view: every noncanonical ORF is placed in its frame
#' and strand context against the standard mitochondrial gene map.
#'
#' Pure plotting: it does no network access.  Fetch the genome/annotation
#' yourself (e.g. with [fetch_mito_genbank()]) and pass `genome_length` and
#' `genes`.
#'
#' @param orfs A `data.frame` of ORFs with integer `frame` (`1,2,3,-1,-2,-3`),
#'   `start`, `end`, and (optionally) `wraps_around` columns — as returned by
#'   [scan_orfs()].  If it has a `seq_id` column with more than one value, pass
#'   `seq_id` to select one.
#' @param genome_length Integer length of the (circular) genome in nt.
#' @param genes Optional `data.frame` of canonical features for the central
#'   ring, with `start`, `end`, and a label column (`gene` or `product`); a
#'   `type` column (e.g. `CDS`, `rRNA`, `tRNA`, `D-loop`) is used for colouring
#'   when present.  Pass `fetch_mito_genbank(acc)[[acc]]$features`.
#' @param seq_id Optional value to subset `orfs$seq_id` when multiple sequences
#'   are present.
#' @param file Optional output path ending in `.png` or `.pdf`.  If `NULL` the
#'   plot is drawn on the current graphics device.
#' @param title Optional plot title.
#' @param frame_colors Named character vector of six colours keyed by
#'   `"1","2","3","-1","-2","-3"`.  Defaults to blues (forward) / reds (reverse).
#' @param label_genes Character vector of feature `type`s to label (others are
#'   drawn but unlabelled, to reduce clutter).  Default `c("CDS","rRNA")`.
#' @param track_height Height of each ORF ring (fraction of radius).  Default
#'   `0.07`.
#'
#' @return Invisibly `NULL`; called for its side effect (the plot).
#'
#' @details
#' Requires the **circlize** package.  ORFs flagged `wraps_around` (spanning the
#' circular origin) are drawn as two arcs.  Rings are labelled at the top gap.
#'
#' @seealso [scan_orfs()], [scan_orfeome()], [fetch_mito_genbank()]
#'
#' @examples
#' \dontrun{
#' g    <- fetch_mito_genbank("NC_012920.1")[["NC_012920.1"]]
#' orfs <- scan_orfs(Biostrings::DNAStringSet(g$sequence),
#'                   min_orf_length = 150, circular = TRUE)
#' plot_orfeome_circos(orfs, genome_length = g$length, genes = g$features,
#'                     title = "Human rCRS ORFeome", file = "orfeome.png")
#' }
#'
#' @export
plot_orfeome_circos <- function(orfs,
                                genome_length,
                                genes        = NULL,
                                seq_id       = NULL,
                                file         = NULL,
                                title        = NULL,
                                frame_colors = NULL,
                                label_genes  = c("CDS", "rRNA"),
                                track_height = 0.07) {
  if (!requireNamespace("circlize", quietly = TRUE))
    stop("Package 'circlize' is required for plot_orfeome_circos(). ",
         "Install it with install.packages('circlize').", call. = FALSE)
  if (!all(c("frame", "start", "end") %in% names(orfs)))
    stop("'orfs' must have columns 'frame', 'start', 'end'.", call. = FALSE)
  if (!is.null(seq_id) && "seq_id" %in% names(orfs))
    orfs <- orfs[orfs$seq_id == seq_id, , drop = FALSE]
  genome_length <- as.numeric(genome_length)

  if (is.null(frame_colors))
    frame_colors <- c("1" = "#08519c", "2" = "#3182bd", "3" = "#6baed6",
                      "-1" = "#a50f15", "-2" = "#de2d26", "-3" = "#fb6a4a")
  ring_order   <- c(1L, 2L, 3L, -1L, -2L, -3L)   # outermost -> innermost
  frame_label  <- function(fr) if (fr > 0) paste0("+", fr) else as.character(fr)
  wrap_col     <- if ("wraps_around" %in% names(orfs)) orfs$wraps_around else rep(FALSE, nrow(orfs))

  if (!is.null(file)) {
    ext <- tolower(tools::file_ext(file))
    if (ext == "png")      grDevices::png(file, width = 2200, height = 2200, res = 220)
    else if (ext == "pdf") grDevices::pdf(file, width = 10, height = 10)
    else stop("'file' must end in .png or .pdf", call. = FALSE)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  circlize::circos.clear()
  circlize::circos.par(start.degree   = 90,
                       gap.degree      = 10,
                       cell.padding    = c(0, 0, 0, 0),
                       track.margin    = c(0.004, 0.004),
                       points.overflow.warning = FALSE)
  circlize::circos.initialize(sectors = "MT", xlim = c(0, genome_length))

  ## ---- six ORF rings (outermost -> innermost) ------------------------------
  for (fr in ring_order) {
    circlize::circos.track(ylim = c(0, 1), track.height = track_height,
                           bg.border = "grey92")
    sub <- orfs[orfs$frame == fr, , drop = FALSE]
    sub_wrap <- wrap_col[orfs$frame == fr]
    col <- grDevices::adjustcolor(frame_colors[as.character(fr)], alpha.f = 0.7)
    if (nrow(sub)) for (i in seq_len(nrow(sub))) {
      s <- sub$start[i]; e <- sub$end[i]
      if (isTRUE(sub_wrap[i]) && e < s) {
        circlize::circos.rect(s, 0.15, genome_length, 0.85, sector.index = "MT",
                              col = col, border = NA)
        circlize::circos.rect(0, 0.15, e, 0.85, sector.index = "MT",
                              col = col, border = NA)
      } else {
        circlize::circos.rect(s, 0.15, e, 0.85, sector.index = "MT",
                              col = col, border = NA)
      }
    }
    # ring label, just inside the top gap
    circlize::circos.text(0, 0.5, labels = frame_label(fr), sector.index = "MT",
                          adj = c(1.4, 0.5), cex = 0.6, col = col, facing = "downward")
  }

  ## ---- central ring: coordinate axis + canonical genes ---------------------
  circlize::circos.track(ylim = c(0, 1), track.height = 0.14, bg.border = "grey60")
  circlize::circos.axis(h = "top",
                        major.at   = seq(0, genome_length, by = 2000),
                        labels     = paste0(seq(0, genome_length, by = 2000) / 1000, "k"),
                        labels.cex = 0.45, major.tick.length = 0.4,
                        sector.index = "MT")
  # Drop redundant 'gene' wrappers and non-feature rows so the central ring is
  # a clean product map (GenBank lists both a 'gene' and a CDS/tRNA/rRNA row per
  # locus); keep only the canonical product types when a 'type' column exists.
  if (!is.null(genes) && "type" %in% names(genes))
    genes <- genes[genes$type %in%
                   c("CDS", "rRNA", "tRNA", "D-loop", "rep_origin", "ncRNA"),
                   , drop = FALSE]
  if (!is.null(genes) && nrow(genes)) {
    type_col <- c(CDS = "#4575b4", rRNA = "#fd8d3c", tRNA = "grey75",
                  `D-loop` = "#ffd700", rep_origin = "#ffd700")
    label_fld <- if ("gene" %in% names(genes)) "gene" else
                 if ("product" %in% names(genes)) "product" else NULL
    for (i in seq_len(nrow(genes))) {
      ty  <- if ("type" %in% names(genes)) genes$type[i] else "CDS"
      gc  <- if (!is.na(type_col[ty])) type_col[ty] else "grey60"
      circlize::circos.rect(genes$start[i], 0.25, genes$end[i], 0.75,
                            sector.index = "MT", col = gc, border = "grey30")
      if (!is.null(label_fld) && ty %in% label_genes && !is.na(genes[[label_fld]][i])) {
        mid <- (genes$start[i] + genes$end[i]) / 2
        circlize::circos.text(mid, 1.25, labels = genes[[label_fld]][i],
                              sector.index = "MT", cex = 0.4,
                              facing = "clockwise", niceFacing = TRUE)
      }
    }
  }

  if (!is.null(title)) graphics::title(main = title, cex.main = 1.1)
  # legend (base graphics over the circlize plot)
  graphics::legend("bottomleft",
                   legend = vapply(ring_order, frame_label, character(1)),
                   fill   = frame_colors[as.character(ring_order)],
                   title  = "frame", bty = "n", cex = 0.7)
  circlize::circos.clear()
  invisible(NULL)
}
