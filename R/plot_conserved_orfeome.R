# ---------------------------------------------------------------------------
# plot_conserved_orfeome() — one-call circos view of an ORFeome with conserved
# (shared) ORFs highlighted, driven by a find_conserved_windows() result.
# ---------------------------------------------------------------------------

#' Circos plot of an ORFeome with conserved ORFs highlighted
#'
#' Convenience wrapper that turns a [find_conserved_windows()] result straight
#' into a [plot_orfeome_circos()] figure for one genome.  It fetches the
#' genome, scans its full ORFeome, flags each ORF as **shared** (overlapping a
#' conserved window) or **unique**, and draws the six-frame circos plot with
#' the shared ORFs outlined and labelled.
#'
#' Reading frame is encoded by the concentric rings (as in
#' [plot_orfeome_circos()]); conservation is encoded by reusing that function's
#' MDP-highlight mechanism: the shared ORFs' own protein sequences are passed
#' as `custom_mdps`, so each shared ORF self-matches and is outlined in
#' `highlight_color`.  Unique ORFs appear as plain frame-coloured bars.
#'
#' By default the plot is drawn to the active graphics device (it "pops up").
#' Pass `file = "x.png"` (forwarded via `...`) to write a file instead.
#'
#' @param windows A `data.frame` returned by [find_conserved_windows()].
#' @param acc Accession of the genome to plot.  Conserved windows are matched
#'   wherever `acc` appears, whether as `acc_a` or `acc_b`, using the
#'   corresponding genome coordinates, so either member of a pair can be
#'   plotted directly.
#' @param genetic_code Genetic-code identifier passed to [scan_orfs()].
#'   Required (no default) — use `"SGC1"` for canonical mito genes or `"SGC0"`
#'   for a noncanonical/MDP scan.
#' @param min_orf_length Integer minimum ORF length in nt, passed to
#'   [scan_orfs()].  Default `30L`.
#' @param highlight_identity Minimum percent identity for the self-match that
#'   outlines a shared ORF, passed as `mdp_min_identity`.  Keep high
#'   (default `95`) so shared ORFs light up while the built-in human MDP
#'   reference set stays quiet.
#' @param highlight_color Outline/label colour for shared ORFs.  Default
#'   `"#7a0177"`.
#' @param title Plot title.  Default is derived from `acc`.
#' @param ... Further arguments forwarded to [plot_orfeome_circos()] (e.g.
#'   `file`, `collapse`, `frame_colors`, `track_height`).
#'
#' @return Invisibly, the scanned ORF `data.frame` with an added logical
#'   `shared` column.  Called mainly for the side effect (the plot).
#'
#' @seealso [find_conserved_windows()], [plot_orfeome_circos()],
#'   [scan_orfs()], [fetch_mito_genbank()].
#'
#' @examples
#' \dontrun{
#' genomes <- fetch_mito_genbank(c("MK848690.1", "KT289925.1"))
#' windows <- find_conserved_windows(genomes, genetic_codes = "SGC0")
#'
#' # Pops up on screen; shared ORFs outlined, frame = ring
#' plot_conserved_orfeome(windows, acc = "MK848690.1", genetic_code = "SGC0")
#'
#' # The other side of the pair, written to a file
#' plot_conserved_orfeome(windows, acc = "KT289925.1", genetic_code = "SGC0",
#'                        file = "kt289925_orfeome.png")
#' }
#'
#' @importFrom Biostrings DNAStringSet
#' @importFrom stats setNames
#' @export
plot_conserved_orfeome <- function(windows,
                                    acc,
                                    genetic_code      = .stop_no_code("genetic_code"),
                                    min_orf_length    = 30L,
                                    highlight_identity = 95,
                                    highlight_color   = "#7a0177",
                                    title             = NULL,
                                    ...) {
  req <- c("acc_a", "acc_b", "genome_start_a", "genome_end_a",
           "genome_start_b", "genome_end_b")
  if (!is.data.frame(windows) || !all(req %in% names(windows)))
    stop("'windows' must be a data.frame from find_conserved_windows() ",
         "(needs columns ", paste(req, collapse = ", "), ").", call. = FALSE)
  if (missing(acc) || length(acc) != 1L || is.na(acc))
    stop("'acc' must be a single accession to plot.", call. = FALSE)

  ## ---- genome record (sequence + features for the central gene ring) -------
  g <- fetch_mito_genbank(acc)[[acc]]
  if (is.null(g))
    stop("fetch_mito_genbank() returned nothing for '", acc, "'.", call. = FALSE)

  ## ---- full ORFeome of this genome -----------------------------------------
  orfs <- scan_orfs(Biostrings::DNAStringSet(g$sequence),
                    genetic_code   = genetic_code,
                    min_orf_length = min_orf_length,
                    both_strands   = TRUE,
                    circular       = TRUE)
  if (is.null(orfs) || nrow(orfs) == 0L) {
    message("No ORFs found in ", acc, "; nothing to plot.")
    return(invisible(orfs))
  }

  ## ---- flag shared ORFs (overlap a conserved window for this acc) ----------
  # find_conserved_windows() fixes an A/B orientation per row, so gather the
  # window spans for this accession from whichever side(s) it appears on.
  as_a <- windows[windows$acc_a == acc, , drop = FALSE]
  as_b <- windows[windows$acc_b == acc, , drop = FALSE]
  win <- rbind(
    data.frame(start = as_a$genome_start_a, end = as_a$genome_end_a),
    data.frame(start = as_b$genome_start_b, end = as_b$genome_end_b))
  if (nrow(win) == 0L)
    message("No conserved windows mention '", acc, "' (as acc_a or acc_b); ",
            "all ORFs will read as unique.")
  glen <- as.numeric(g$length)
  orfs$shared <- vapply(seq_len(nrow(orfs)), function(i) {
    if (nrow(win) == 0L) return(FALSE)
    s <- orfs$start[i]; e <- orfs$end[i]
    if (!is.na(e) && !is.na(s) && e < s) {
      # wrap-around ORF: covers [s, glen] U [1, e]
      any((win$end >= s) | (win$start <= e))
    } else {
      any(s <= win$end & e >= win$start)
    }
  }, logical(1))

  ## ---- shared ORFs become custom MDP refs -> self-match -> outlined --------
  shared_refs <- if (any(orfs$shared) && "protein_sequence" %in% names(orfs))
    stats::setNames(orfs$protein_sequence[orfs$shared],
                    paste0("shared_", which(orfs$shared)))
  else NULL

  if (is.null(title))
    title <- paste0(acc, " ORFeome — shared (outlined) vs unique")

  plot_orfeome_circos(
    orfs,
    genome_length    = glen,
    genes            = g$features,
    title            = title,
    mark_mdps        = !is.null(shared_refs),
    custom_mdps      = shared_refs,
    mdp_min_identity = highlight_identity,
    mdp_color        = highlight_color,
    ...)

  invisible(orfs)
}
