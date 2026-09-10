# ============================================================================
# plot_rna_structure — 2D RNA secondary-structure drawing with ORF highlights
# ============================================================================

# ── internal helpers ──────────────────────────────────────────────────────────

# Check for the RNAplot binary and return its path, or stop with instructions
.check_rnaplot <- function() {
  path <- Sys.which("RNAplot")
  if (nchar(path) == 0L)
    stop(
      "RNAplot (ViennaRNA) not found on PATH.\n",
      "Install ViennaRNA: https://www.tbi.univie.ac.at/RNA/#download\n",
      "  conda: conda install -c bioconda viennarna\n",
      "  apt  : sudo apt-get install vienna-rna",
      call. = FALSE
    )
  path
}

# Run RNAplot and parse the /coor array out of its PostScript output.
# Returns a two-column matrix (x, y) with one row per nucleotide.
.rnaplot_layout <- function(sequence, structure, layout_type = 1L) {
  .check_rnaplot()
  if (nchar(sequence) != nchar(structure))
    stop("sequence and structure must be the same length (", nchar(sequence),
         " vs ", nchar(structure), ").", call. = FALSE)

  tag <- "mitodiver"
  wd  <- file.path(tempdir(), paste0("rnaplot_", as.integer(Sys.getpid()), "_",
                                     sample.int(1e6, 1L)))
  dir.create(wd, showWarnings = FALSE, recursive = TRUE)
  old <- setwd(wd)
  on.exit({ setwd(old); unlink(wd, recursive = TRUE) }, add = TRUE)

  input <- sprintf(">%s\n%s\n%s\n", tag, .to_rna(sequence), structure)
  out   <- suppressWarnings(
    system2("RNAplot", args = c("-t", layout_type),
            stdout = TRUE, stderr = TRUE, input = input)
  )

  ps <- file.path(wd, paste0(tag, "_ss.ps"))
  if (!file.exists(ps))
    stop("RNAplot did not produce a layout: ", paste(out, collapse = " "),
         call. = FALSE)

  ln <- readLines(ps, warn = FALSE)
  i0 <- grep("^/coor\\s*\\[", ln)
  if (length(i0) == 0L)
    stop("Could not find coordinates in RNAplot output.", call. = FALSE)
  i0 <- i0[1L]
  i1 <- which(grepl("^\\]\\s*def", ln) & seq_along(ln) > i0)[1L]
  rows <- ln[(i0 + 1L):(i1 - 1L)]
  rows <- rows[grepl("^\\s*\\[", rows)]

  xy <- do.call(rbind, lapply(strsplit(trimws(gsub("[][]", "", rows)), "\\s+"),
                              function(v) as.numeric(v[nzchar(v)][1:2])))
  colnames(xy) <- c("x", "y")
  if (nrow(xy) != nchar(structure))
    stop("RNAplot returned ", nrow(xy), " coordinates for a ",
         nchar(structure), " nt structure.", call. = FALSE)
  xy
}

# Add an alpha channel to a colour
.alpha_col <- function(col, alpha) {
  m <- grDevices::col2rgb(col)
  grDevices::rgb(m[1L, ], m[2L, ], m[3L, ],
                 alpha = alpha * 255, maxColorValue = 255)
}

# Resolve colours for the named MDPs present in a plot.
.mdp_palette <- function(names_present, mdp_colors) {
  u <- unique(stats::na.omit(names_present))
  if (length(u) == 0L) return(character(0))

  if (is.function(mdp_colors)) {
    cols <- mdp_colors(length(u))
    if (length(cols) < length(u))
      stop("`mdp_colors` returned ", length(cols), " colours for ", length(u),
           " MDPs.", call. = FALSE)
    return(stats::setNames(cols[seq_along(u)], u))
  }

  if (!is.null(mdp_colors)) {
    if (is.null(names(mdp_colors)))
      stop("`mdp_colors` must be named (keyed on MDP name) or a function.",
           call. = FALSE)
    miss <- setdiff(u, names(mdp_colors))
    if (length(miss))
      stop("`mdp_colors` has no entry for: ", paste(miss, collapse = ", "),
           call. = FALSE)
    return(stats::setNames(unname(mdp_colors[u]), u))
  }

  # Okabe-Ito minus the greys, recycled
  pal <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
           "#E69F00", "#56B4E9", "#F0E442", "#7B3294")
  stats::setNames(pal[(seq_along(u) - 1L) %% length(pal) + 1L], u)
}

# Place MDP labels around the structure: anchor on the span midpoint, push
# radially outward from the layout centroid, then separate any that collide.
.draw_mdp_labels <- function(xy, idx, labels, col, cex, offset) {
  if (length(idx) == 0L) return(invisible(NULL))

  cx <- mean(range(xy[, "x"])); cy <- mean(range(xy[, "y"]))
  scale <- max(diff(range(xy[, "x"])), diff(range(xy[, "y"])))
  ax <- xy[idx, "x"]; ay <- xy[idx, "y"]

  dx <- ax - cx; dy <- ay - cy
  r  <- sqrt(dx^2 + dy^2)
  r[r == 0] <- 1e-9
  ux <- dx / r; uy <- dy / r

  # An ORF sitting near the centroid has a small radius, so a fixed outward
  # nudge would leave its label lying on top of the drawing.  Lift such labels
  # most of the way to the structure's extent along their own direction;
  # peripheral labels, where r already ~ the extent, are left where they are.
  px <- xy[, "x"] - cx; py <- xy[, "y"] - cy
  reach <- vapply(seq_along(idx), function(i) max(px * ux[i] + py * uy[i]),
                  numeric(1))
  target <- pmax(r, 0.75 * reach) + offset * scale

  push <- target - r
  lx <- ax + ux * push; ly <- ay + uy * push

  # greedy separation: nudge outward until labels are not on top of one another
  minsep <- 0.052 * scale
  if (length(idx) > 1L) {
    for (pass in 1:60) {
      moved <- FALSE
      o <- order(r, decreasing = TRUE)
      for (a in seq_along(o)) for (b in seq_along(o)) {
        if (b <= a) next
        i <- o[a]; j <- o[b]
        d <- sqrt((lx[i] - lx[j])^2 + (ly[i] - ly[j])^2)
        if (d < minsep) {
          push[j] <- push[j] + 0.22 * minsep
          lx[j] <- ax[j] + ux[j] * push[j]
          ly[j] <- ay[j] + uy[j] * push[j]
          moved <- TRUE
        }
      }
      if (!moved) break
    }
  }

  segments(ax, ay, lx, ly, col = .alpha_col(col, 0.75), lwd = 0.8)
  points(ax, ay, pch = 16, col = col, cex = 0.45)
  text(lx, ly, labels = labels, col = col, cex = cex, font = 2,
       adj = ifelse(ux >= 0, -0.06, 1.06))
  invisible(NULL)
}

# ── main function ─────────────────────────────────────────────────────────────

#' Draw an RNA Secondary Structure, Optionally Highlighting ORFs
#'
#' Renders the classic 2D "squiggle plot" of an RNA secondary structure --
#' backbone, base pairs, loops -- and optionally paints open reading frames onto
#' the backbone as coloured highlights, with a separate colour per strand.  This
#' makes it possible to see directly where a noncanonical ORF (for example a
#' mitochondrial-derived peptide) sits relative to the stems and loops of the
#' RNA that encodes it, and whether it reads in the sense or antisense
#' direction.
#'
#' @param x An `"rna_fold"` object as returned by [fold_rna()], or any list with
#'   `sequence` and `structure` elements.  A dot-bracket structure string may
#'   also be passed directly, in which case `sequence` must be supplied too.
#' @param orfs Optional `data.frame` of ORFs to highlight, normally the output
#'   of [scan_orfs()] or [find_orfs()].  Must contain `start` and `end` columns;
#'   a `strand` column (`"+"` / `"-"`) selects the highlight colour.  `NULL`
#'   (default) draws the structure with no highlights.
#' @param sequence Nucleotide sequence, required only when `x` is a bare
#'   dot-bracket string.
#' @param seq_id Optional character scalar.  When `orfs` has a `seq_id` column
#'   spanning several sequences, keep only the rows for this one.
#' @param offset Integer.  Coordinate of the first folded base in the coordinate
#'   system used by `orfs`.  Use this when the fold covers a sub-region of a
#'   larger genome and `orfs` carries genomic coordinates -- e.g. `offset =
#'   1671` for the human 16S rRNA within the rCRS.  Default `1` (ORF and fold
#'   coordinates agree).  ORFs are clipped to the folded region and any that
#'   fall entirely outside it are dropped.
#' @param strand_colors Named character vector of colours for the `"+"` and
#'   `"-"` strands.  Defaults to a colourblind-safe blue/vermillion pair.
#' @param strand_labels Named character vector of legend labels for the `"+"`
#'   and `"-"` strands.
#' @param mdp_colors Colours for named MDPs, used when `orfs` carries an
#'   `mdp_name` column (see [annotate_mdp_orfs()]).  Either a named character
#'   vector keyed on MDP name, or a single palette function taking the number of
#'   MDPs and returning that many colours.  Defaults to a colourblind-safe
#'   qualitative palette recycled as needed.
#' @param mdp_label How named ORFs are labelled on the drawing: `"name"`
#'   (default), `"name_pct"` to append the percent identity -- useful when the
#'   names came from `method = "homology"` and are inferences rather than
#'   identities -- or `"none"` to colour them without on-figure text.
#' @param mdp_label_cex Character expansion for MDP labels.  Default `0.72`.
#' @param mdp_label_offset Radial distance of the label from the structure, as a
#'   fraction of the plot's larger dimension.  Default `0.075`.
#' @param unnamed_alpha Opacity multiplier applied to ORFs that did *not* match
#'   an MDP, so that named ones read as foreground.  Default `0.45`.  Set to `1`
#'   to draw all ORFs at `highlight_alpha`.
#' @param mark Optional integer vector of positions to mark with a point (same
#'   coordinate system as `orfs`, i.e. subject to `offset`) -- useful for
#'   flagging variant sites.
#' @param mark_labels Optional character vector of labels drawn next to `mark`.
#' @param highlight_lwd Line width of the ORF highlights.  Default `4.2`.
#' @param highlight_alpha Opacity of the ORF highlights, `0`--`1`.  Default
#'   `0.55`, so that overlapping ORFs remain visible through one another.
#' @param backbone_col Colour of the sequence backbone.
#' @param pair_col Colour of the base pairs.  Either a single colour, one colour
#'   per base pair, or a **function** taking the base-pair `data.frame`
#'   (columns `pos5`, `pos3`, as returned by [fold_rna()]) and returning a
#'   colour vector.  The function form is the safe way to colour pairs by
#'   identity -- e.g. shared vs genome-specific helices -- without depending on
#'   the order in which pairs are enumerated.
#' @param layout Layout algorithm passed to `RNAplot`: `"naview"` (default),
#'   `"simple"`, `"circular"`, `"turtle"` or `"puzzler"`.
#' @param main,sub Optional title and subtitle.
#' @param legend Logical.  Draw a strand legend when ORFs are highlighted.
#'   Default `TRUE`.
#' @param legend_pos Legend position, passed to [graphics::legend()].  Default
#'   `"topleft"`.
#' @param add Logical.  If `TRUE`, draw into the current plot instead of
#'   starting a new one.  Default `FALSE`.
#' @param ... Further arguments passed to [graphics::plot.default()].
#'
#' @return Invisibly, a `data.frame` with one row per nucleotide and columns
#'   `position`, `x`, `y` (layout coordinates), and -- when `orfs` was supplied
#'   -- `n_plus` and `n_minus` giving the number of ORFs covering each base on
#'   each strand.  Useful for adding further annotation to the plot.
#'
#' @details
#' ## Dependency
#' The 2D layout is computed by ViennaRNA's `RNAplot` program, which must be
#' installed and on the system `PATH`; see [fold_rna()] for installation notes.
#' Only the coordinates are taken from `RNAplot` -- all drawing is done with R
#' graphics, so the result composes with `par()`, `layout()` and any graphics
#' device.
#'
#' ## Strand conventions
#' `strand_colors` is keyed on the `strand` column of `orfs`, which for
#' [scan_orfs()] is `"+"` for the given sequence and `"-"` for its reverse
#' complement.  In vertebrate mitochondrial genomes as deposited (e.g. the rCRS)
#' the `"+"` strand corresponds to the heavy (H) strand and `"-"` to the light
#' (L) strand, so the default labels name both conventions.
#'
#' ORFs that wrap the origin of a circular genome (`wraps_around == TRUE`)
#' cannot be drawn as a single stretch of backbone and are dropped with a
#' warning.
#'
#' ## Naming ORFs
#' When `orfs` carries an `mdp_name` column -- as produced by
#' [annotate_mdp_orfs()] -- named ORFs are drawn in a per-MDP colour and
#' labelled on the figure, while ORFs that matched nothing stay in the muted
#' strand colours as context.  This is the distinction the plot usually needs to
#' make: in the human 16S every ORF above 60 nt is a known MDP, whereas in the
#' 12S only one of twenty-one is, and the reader has to be able to tell which.
#' Without an `mdp_name` column the function behaves exactly as before, colouring
#' purely by strand.
#'
#' @seealso [annotate_mdp_orfs()], [fold_rna()], [scan_orfs()],
#'   [collapse_nested_orfs()]
#'
#' @examples
#' \dontrun{
#' # Fold the human 16S rRNA (MT-RNR2) and highlight the ORFs inside it
#' rcrs <- fetch_sequences(accessions = "NC_012920.1")
#' r16s <- Biostrings::subseq(rcrs[[1]], 1671, 3229)
#'
#' fold <- fold_rna(as.character(r16s))
#' orfs <- scan_orfs(Biostrings::DNAStringSet(c(rnr2 = as.character(r16s))),
#'                   min_orf_length = 60,
#'                   genetic_code   = "SGC0",   # MDPs are read on 80S ribosomes
#'                   both_strands   = TRUE,
#'                   circular       = FALSE)
#'
#' plot_rna_structure(fold, orfs = collapse_nested_orfs(orfs),
#'                    main = "MT-RNR2 (rCRS)")
#'
#' # ORFs in genomic coordinates, marking a variant site
#' plot_rna_structure(fold, orfs = genomic_orfs, offset = 1671,
#'                    mark = 3206, mark_labels = "m.3206")
#' }
#'
#' @importFrom graphics plot lines segments points text title mtext legend par
#' @importFrom grDevices col2rgb rgb
#' @export
plot_rna_structure <- function(x,
                               orfs            = NULL,
                               sequence        = NULL,
                               seq_id          = NULL,
                               offset          = 1L,
                               strand_colors   = c("+" = "#0072B2",
                                                   "-" = "#D55E00"),
                               strand_labels   = c("+" = "+ strand (H) ORF",
                                                   "-" = "- strand (L) ORF"),
                               mdp_colors       = NULL,
                               mdp_label        = c("name", "name_pct", "none"),
                               mdp_label_cex    = 0.72,
                               mdp_label_offset = 0.075,
                               unnamed_alpha    = 0.45,
                               mark            = NULL,
                               mark_labels      = NULL,
                               highlight_lwd   = 4.2,
                               highlight_alpha = 0.55,
                               backbone_col    = "#9A9A9A",
                               pair_col        = "#D8D8D8",
                               layout          = c("naview", "simple",
                                                   "circular", "turtle",
                                                   "puzzler"),
                               main            = NULL,
                               sub             = NULL,
                               legend          = TRUE,
                               legend_pos      = "topleft",
                               add             = FALSE,
                               ...) {

  # ── resolve structure + sequence ────────────────────────────────────────────
  if (is.character(x) && length(x) == 1L) {
    structure <- x
    if (is.null(sequence))
      stop("When `x` is a dot-bracket string, `sequence` must be supplied.",
           call. = FALSE)
  } else if (is.list(x) && all(c("sequence", "structure") %in% names(x))) {
    structure <- x$structure
    if (is.null(sequence)) sequence <- x$sequence
  } else {
    stop("`x` must be an \"rna_fold\" object, a list with `sequence` and ",
         "`structure`, or a dot-bracket string.", call. = FALSE)
  }
  sequence  <- as.character(sequence)
  structure <- as.character(structure)
  n         <- nchar(structure)

  layout_type <- switch(match.arg(layout),
                        simple = 0L, naview = 1L, circular = 2L,
                        turtle = 3L, puzzler = 4L)
  mdp_label <- match.arg(mdp_label)

  # ── layout + base pairs ─────────────────────────────────────────────────────
  xy <- .rnaplot_layout(sequence, structure, layout_type)
  bp <- .parse_dot_bracket(structure)

  # ── ORFs -> local coordinates ───────────────────────────────────────────────
  cov_plus <- cov_minus <- integer(n)
  spans    <- NULL
  if (!is.null(orfs) && nrow(orfs) > 0L) {
    if (!all(c("start", "end") %in% names(orfs)))
      stop("`orfs` must have `start` and `end` columns.", call. = FALSE)

    if (!is.null(seq_id)) {
      if (!"seq_id" %in% names(orfs))
        stop("`seq_id` was given but `orfs` has no `seq_id` column.",
             call. = FALSE)
      orfs <- orfs[orfs$seq_id == seq_id, , drop = FALSE]
    } else if ("seq_id" %in% names(orfs) &&
               length(unique(orfs$seq_id)) > 1L) {
      stop("`orfs` spans ", length(unique(orfs$seq_id)),
           " sequences; supply `seq_id` to choose one.", call. = FALSE)
    }

    if ("wraps_around" %in% names(orfs)) {
      nw <- sum(orfs$wraps_around %in% TRUE)
      if (nw > 0L) {
        warning(nw, " origin-spanning ORF(s) dropped; they cannot be drawn ",
                "as a single backbone stretch.", call. = FALSE)
        orfs <- orfs[!(orfs$wraps_around %in% TRUE), , drop = FALSE]
      }
    }

    strand <- if ("strand" %in% names(orfs)) as.character(orfs$strand)
              else rep("+", nrow(orfs))
    mdp_nm <- if ("mdp_name" %in% names(orfs)) as.character(orfs$mdp_name)
              else rep(NA_character_, nrow(orfs))
    mdp_pc <- if ("mdp_pct_identity" %in% names(orfs)) orfs$mdp_pct_identity
              else rep(NA_real_, nrow(orfs))
    lo <- pmin(orfs$start, orfs$end) - offset + 1L
    hi <- pmax(orfs$start, orfs$end) - offset + 1L

    keep <- hi >= 1L & lo <= n
    if (any(!keep))
      message("plot_rna_structure(): ", sum(!keep),
              " ORF(s) fall outside the folded region and were dropped.")
    lo <- pmax(1L, lo[keep]); hi <- pmin(n, hi[keep]); strand <- strand[keep]
    mdp_nm <- mdp_nm[keep]; mdp_pc <- mdp_pc[keep]

    if (length(lo) > 0L) {
      spans <- data.frame(lo = lo, hi = hi, strand = strand,
                          mdp_name = mdp_nm, mdp_pct = mdp_pc,
                          stringsAsFactors = FALSE)
      for (k in seq_len(nrow(spans))) {
        idx <- spans$lo[k]:spans$hi[k]
        if (spans$strand[k] == "-") cov_minus[idx] <- cov_minus[idx] + 1L
        else                        cov_plus[idx]  <- cov_plus[idx]  + 1L
      }
    }
  }

  # ── draw ────────────────────────────────────────────────────────────────────
  if (!add) {
    pad <- 0.045 * max(diff(range(xy[, "x"])), diff(range(xy[, "y"])))
    plot(NA, asp = 1, axes = FALSE, xlab = "", ylab = "",
         xlim = range(xy[, "x"]) + c(-pad, pad),
         ylim = range(xy[, "y"]) + c(-pad, pad), ...)
  }

  if (nrow(bp) > 0L) {
    pcol <- if (is.function(pair_col)) pair_col(bp) else pair_col
    if (length(pcol) != 1L && length(pcol) != nrow(bp))
      stop("`pair_col` must be one colour, a function of the base-pair table, ",
           "or one colour per base pair (", nrow(bp), ").", call. = FALSE)
    segments(xy[bp$pos5, "x"], xy[bp$pos5, "y"],
             xy[bp$pos3, "x"], xy[bp$pos3, "y"], col = pcol, lwd = 0.7)
  }
  lines(xy[, "x"], xy[, "y"], col = backbone_col, lwd = 1.0)

  if (!is.null(spans)) {
    named    <- !is.na(spans$mdp_name)
    mdp_pal  <- .mdp_palette(spans$mdp_name[named], mdp_colors)

    for (k in seq_len(nrow(spans))) {
      idx <- spans$lo[k]:spans$hi[k]
      if (named[k]) {
        col <- mdp_pal[[spans$mdp_name[k]]]
        a   <- highlight_alpha
      } else {
        col <- strand_colors[[spans$strand[k]]]
        if (is.null(col)) col <- strand_colors[[1L]]
        a   <- highlight_alpha * unnamed_alpha
      }
      lines(xy[idx, "x"], xy[idx, "y"], col = .alpha_col(col, a),
            lwd = highlight_lwd, lend = 1)
    }

    if (mdp_label != "none" && any(named)) {
      lab <- spans$mdp_name[named]
      if (mdp_label == "name_pct" && any(!is.na(spans$mdp_pct[named])))
        lab <- ifelse(is.na(spans$mdp_pct[named]), lab,
                      sprintf("%s %.0f%%", lab, spans$mdp_pct[named]))
      mid <- floor((spans$lo[named] + spans$hi[named]) / 2)
      .draw_mdp_labels(xy, mid, lab,
                       col    = unname(mdp_pal[spans$mdp_name[named]]),
                       cex    = mdp_label_cex,
                       offset = mdp_label_offset)
    }
  }

  if (!is.null(mark) && length(mark) > 0L) {
    mi   <- as.integer(mark) - offset + 1L
    ok   <- mi >= 1L & mi <= n
    if (any(!ok))
      warning(sum(!ok), " `mark` position(s) fall outside the folded region.",
              call. = FALSE)
    mi   <- mi[ok]
    if (length(mi) > 0L) {
      points(xy[mi, "x"], xy[mi, "y"], pch = 21, bg = "white", col = "black",
             cex = 1.15, lwd = 1.4)
      if (!is.null(mark_labels))
        text(xy[mi, "x"], xy[mi, "y"], labels = mark_labels[ok],
             pos = 3, offset = 0.7, cex = 0.68, font = 2)
    }
  }

  if (!is.null(main)) title(main = main, font.main = 2)
  if (!is.null(sub))  mtext(sub, side = 3, line = 0.2, cex = 0.76,
                            col = "grey25")

  if (isTRUE(legend) && !is.null(spans)) {
    named   <- !is.na(spans$mdp_name)
    present <- unique(spans$strand[!named])
    present <- present[present %in% names(strand_colors)]

    leg_txt <- character(0); leg_col <- character(0)

    # When MDPs are labelled on the drawing itself, repeating every name in the
    # legend is redundant; show the strand key for the unnamed context ORFs and
    # name the MDPs here only if they were not drawn with text.
    if (any(named) && mdp_label == "none") {
      pal <- .mdp_palette(spans$mdp_name[named], mdp_colors)
      leg_txt <- names(pal)
      leg_col <- .alpha_col(unname(pal), highlight_alpha)
    }
    if (length(present) > 0L) {
      lbl <- unname(strand_labels[present])
      if (any(named)) lbl <- paste0(lbl, " (unassigned)")
      leg_txt <- c(leg_txt, lbl)
      leg_col <- c(leg_col,
                   .alpha_col(unname(strand_colors[present]),
                              highlight_alpha *
                                if (any(named)) unnamed_alpha else 1))
    }

    if (length(leg_txt) > 0L)
      legend(legend_pos, bty = "n", cex = 0.8, seg.len = 1.6, lty = 1,
             lwd = highlight_lwd, col = leg_col, legend = leg_txt)
  }

  out <- data.frame(position = seq_len(n), x = xy[, "x"], y = xy[, "y"],
                    stringsAsFactors = FALSE)
  if (!is.null(orfs)) {
    out$n_plus  <- cov_plus
    out$n_minus <- cov_minus
  }
  invisible(out)
}
