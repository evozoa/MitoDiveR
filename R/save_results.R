# ============================================================================
# save_results — write MitoDiveR pipeline outputs to disk
# ============================================================================

#' Save MitoDiveR Pipeline Results to Disk
#'
#' Writes one or more `data.frame` outputs from the MitoDiveR pipeline to
#' CSV files in a single call.  Any combination of the three main outputs
#' ([scan_orfs()], [find_conserved_orfs()], [score_mdp_similarity()]) can be
#' supplied; only non-`NULL` arguments are written.
#'
#' @param orfs `data.frame` returned by [scan_orfs()], or `NULL` to skip.
#' @param conserved `data.frame` returned by [find_conserved_orfs()], or
#'   `NULL` to skip.
#' @param mdp_hits `data.frame` returned by [score_mdp_similarity()], or
#'   `NULL` to skip.
#' @param outdir Character scalar.  Directory to write files into.  Created
#'   recursively if it does not exist.  Defaults to the current working
#'   directory (`"."`).
#' @param prefix Character scalar.  Prefix prepended to every output filename.
#'   Default `"MitoDiveR"`.  A timestamp is appended automatically so that
#'   repeated runs do not overwrite previous results.
#' @param timestamp Logical.  Append a `YYYYMMDD_HHMMSS` timestamp to
#'   filenames?  Default `TRUE`.  Set to `FALSE` for reproducible filenames
#'   (will overwrite previous results).
#'
#' @return A named character vector of the file paths written, invisibly.
#'
#' @examples
#' \dontrun{
#' seqs      <- fetch_sequences(c("NC_012920.1", "NC_001643.1"))
#' orfs      <- scan_orfs(seqs, min_orf_length = 60, max_orf_length = 900)
#' conserved <- find_conserved_orfs(orfs, min_identity = 50)
#' hits      <- score_mdp_similarity(conserved)
#'
#' # Save all three outputs to ~/results/
#' save_results(orfs, conserved, hits, outdir = "~/results")
#'
#' # Save with a fixed filename (no timestamp)
#' save_results(orfs, conserved, hits,
#'              outdir    = "~/results",
#'              prefix    = "NC012920_vs_NC001643",
#'              timestamp = FALSE)
#'
#' # Save only the MDP hits
#' save_results(mdp_hits = hits, outdir = "~/results")
#' }
#'
#' @export
save_results <- function(orfs      = NULL,
                         conserved = NULL,
                         mdp_hits  = NULL,
                         outdir    = ".",
                         prefix    = "MitoDiveR",
                         timestamp = TRUE) {

  if (is.null(orfs) && is.null(conserved) && is.null(mdp_hits))
    stop("Supply at least one of: 'orfs', 'conserved', 'mdp_hits'.",
         call. = FALSE)

  outdir <- path.expand(outdir)
  if (!dir.exists(outdir))
    dir.create(outdir, recursive = TRUE)

  ts <- if (timestamp) paste0("_", format(Sys.time(), "%Y%m%d_%H%M%S")) else ""

  outputs <- list(
    orfs      = orfs,
    conserved = conserved,
    mdp_hits  = mdp_hits
  )

  written <- character(0L)

  for (nm in names(outputs)) {
    df <- outputs[[nm]]
    if (is.null(df)) next

    if (!is.data.frame(df))
      stop("'", nm, "' must be a data.frame.", call. = FALSE)

    path <- file.path(outdir, paste0(prefix, "_", nm, ts, ".csv"))
    write.csv(df, path, row.names = FALSE)
    written[nm] <- path
    message("Wrote ", nrow(df), " rows -> ", path)
  }

  invisible(written)
}
