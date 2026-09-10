# ============================================================================
# cluster_conserved_orfs — DEPRECATED: thin shim over find_conserved_orfs()
# ============================================================================

#' Cluster ORFs Conserved Across Species (Deprecated)
#'
#' @description
#' **Deprecated.**  `cluster_conserved_orfs()` has been merged into
#' [find_conserved_orfs()], which now supports coverage gating, a configurable
#' percent-identity metric, and configurable alignment parameters.  This
#' function is retained as a thin wrapper so existing scripts keep working; it
#' simply forwards to [find_conserved_orfs()] with the equivalent settings and
#' will be removed in a future release.
#'
#' The mapping is:
#' \preformatted{
#' cluster_conserved_orfs(orfs, min_pct_identity = X, min_coverage = Y, ...)
#'   ->  find_conserved_orfs(orfs,
#'                           min_identity        = X,
#'                           min_coverage        = Y,
#'                           min_genomes         = 2,      # "shared by >= 2 genomes"
#'                           identity_metric     = "PID1", # cluster_*'s metric
#'                           substitution_matrix = ...,
#'                           gap_opening         = ...,
#'                           gap_extension       = ...)
#' }
#'
#' @param orfs A `data.frame` of ORFs with `seq_id` and `protein_sequence`
#'   columns (typically the combined output of [scan_orfs()] across genomes).
#' @param min_pct_identity Numeric.  Minimum percent identity (`PID1`) to
#'   connect two ORFs.  Forwarded as `min_identity`.  Default `40`.
#' @param min_coverage Numeric in \[0, 1\].  Minimum fraction of the *shorter*
#'   protein covered by the alignment.  Default `0.5`.
#' @param substitution_matrix Character.  Substitution matrix name.  Default
#'   `"BLOSUM62"`.
#' @param gap_opening Numeric.  Gap-opening penalty.  Default `10`.
#' @param gap_extension Numeric.  Gap-extension penalty.  Default `4`.
#' @param verbose Logical.  Print progress?  Default `TRUE`.
#'
#' @return The same conserved-ORF table returned by [find_conserved_orfs()]:
#'   the input rows that belong to a conserved cluster, with `cluster_id`,
#'   `n_sequences`, and `cluster_size` columns appended.  (Cluster *numbering*
#'   may differ from older releases, but the set of conserved ORFs is the
#'   same; `cluster_id` has always been an arbitrary within-call label.)
#'
#' @seealso [find_conserved_orfs()], which supersedes this function.
#'
#' @examples
#' \dontrun{
#' orfs_all <- scan_orfs(seqs, genetic_code = "SGC1", min_orf_length = 75)
#'
#' # Old (deprecated):
#' conserved <- cluster_conserved_orfs(orfs_all)
#'
#' # New equivalent:
#' conserved <- find_conserved_orfs(orfs_all, min_identity = 40,
#'                                  min_coverage = 0.5, min_genomes = 2,
#'                                  identity_metric = "PID1", gap_extension = 4)
#' }
#'
#' @export
cluster_conserved_orfs <- function(orfs,
                                    min_pct_identity    = 40,
                                    min_coverage        = 0.5,
                                    substitution_matrix = "BLOSUM62",
                                    gap_opening         = 10,
                                    gap_extension       = 4,
                                    verbose             = TRUE) {

  .Deprecated(
    new     = "find_conserved_orfs",
    package = "MitoDiveR",
    msg = paste0(
      "'cluster_conserved_orfs()' is deprecated and will be removed in a future ",
      "release.\n",
      "  Use 'find_conserved_orfs()', which now supports coverage gating and ",
      "configurable alignment.\n",
      "  Equivalent call:\n",
      "    find_conserved_orfs(orfs, min_identity = ", min_pct_identity,
      ", min_coverage = ", min_coverage, ",\n",
      "                        min_genomes = 2, identity_metric = \"PID1\", ",
      "gap_extension = ", gap_extension, ")"
    )
  )

  find_conserved_orfs(
    orfs                = orfs,
    min_identity        = min_pct_identity,
    alignment_type      = "local",
    min_coverage        = min_coverage,
    min_genomes         = 2L,
    identity_metric     = "PID1",
    substitution_matrix = substitution_matrix,
    gap_opening         = gap_opening,
    gap_extension       = gap_extension,
    verbose             = verbose
  )
}
