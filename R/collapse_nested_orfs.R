# ============================================================================
# collapse_nested_orfs — reduce nested start-codon ORFs to one locus per stop.
# ============================================================================

#' Collapse nested ORFs to the longest ORF per stop codon
#'
#' [scan_orfs()] / [find_orfs()] report *every* start-codon-to-stop ORF, so a
#' single stop-free region yields many nested ORFs (one per internal start
#' codon, all sharing the same stop).  This collapses each such family to the
#' single longest ORF terminating at a given stop, within each reading frame —
#' i.e. the distinct ORF *loci* rather than every nested start.  Use it to get
#' meaningful ORF counts and uncluttered plots (e.g. [plot_orfeome_circos()]).
#'
#' @param orfs A `data.frame` from [scan_orfs()] / [find_orfs()], requiring
#'   `frame`, `strand`, `start`, `end` columns.
#' @param length_col Name of the column to maximise when choosing the
#'   representative ORF per stop.  Default `"length_nt"`; if absent it is
#'   derived as `end - start + 1`.
#'
#' @details
#' The stop codon's genomic coordinate is the `end` on the `+` strand and the
#' `start` on the `-` strand (the table keeps `start <= end`, so on the minus
#' strand the stop sits in the `start` column).  ORFs are grouped by
#' `(frame, stop coordinate)` and the longest member of each group is kept.
#'
#' @return The input `data.frame` with nested ORFs removed — one row per
#'   `(frame, stop)` — preserving all original columns and row order by
#'   genomic `start`.  Returned unchanged (with a warning) if the required
#'   columns are absent.
#'
#' @seealso [scan_orfs()], [find_orfs()], [plot_orfeome_circos()]
#'
#' @examples
#' \dontrun{
#' orfs <- scan_orfs(seqs, genetic_code = "SGC1", min_orf_length = 150)
#' nrow(orfs)                       # many, incl. nested starts
#' loci <- collapse_nested_orfs(orfs)
#' nrow(loci)                       # one per (frame, stop)
#' }
#'
#' @export
collapse_nested_orfs <- function(orfs, length_col = "length_nt") {
  if (is.null(orfs) || nrow(orfs) == 0L) return(orfs)
  if (!all(c("frame", "strand", "start", "end") %in% names(orfs))) {
    warning("collapse_nested_orfs(): need columns frame/strand/start/end; ",
            "returning input unchanged.", call. = FALSE)
    return(orfs)
  }
  len <- if (length_col %in% names(orfs)) orfs[[length_col]]
         else abs(orfs$end - orfs$start) + 1L

  # Stop coordinate: 'end' on +, 'start' on - (start <= end convention).
  stop_coord <- ifelse(orfs$strand == "+", orfs$end, orfs$start)
  key <- paste(orfs$frame, stop_coord, sep = "|")

  ord <- order(key, -len)                       # longest first within each stop
  keep_ord <- ord[!duplicated(key[ord])]        # first (= longest) per stop
  out <- orfs[sort(keep_ord), , drop = FALSE]   # back to genomic-start order
  rownames(out) <- NULL
  out
}
