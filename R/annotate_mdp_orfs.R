# ============================================================================
# annotate_mdp_orfs — attach known-MDP identities to a table of ORFs
# ============================================================================

# No default reference, by design: which peptide set a screen was run against is
# part of the result, so it must be stated explicitly.  Mirrors .stop_no_code().
.stop_no_mdp_reference <- function(arg = "reference") {
  stop("'", arg, "' must be specified (no default, by design). Supply the ",
       "peptide sequences to match against: a named character vector, an ",
       "AAStringSet, or a path to a FASTA file — e.g.\n",
       "  annotate_mdp_orfs(orfs, reference = mdp_sequences())        # human MDPs\n",
       "  annotate_mdp_orfs(orfs, reference = \"my_peptides.fasta\")\n",
       "mdp_sequences() holds the human set (Swiss-Prot sequences, each verified ",
       "to be encoded in the rCRS); pass it explicitly if that is what you want.",
       call. = FALSE)
}

# Coerce the many accepted reference forms to a named character vector.
.resolve_mdp_reference <- function(reference) {
  if (methods::is(reference, "AAStringSet") || methods::is(reference, "XStringSet")) {
    out <- stats::setNames(as.character(reference), names(reference))
  } else if (is.character(reference) && length(reference) == 1L &&
             !grepl("^[A-Za-z*]+$", reference) && file.exists(reference)) {
    aa  <- Biostrings::readAAStringSet(reference)
    out <- stats::setNames(as.character(aa), names(aa))
  } else if (is.character(reference)) {
    out <- reference
  } else if (is.list(reference)) {
    out <- unlist(reference)
  } else {
    stop("`reference` must be a named character vector, an AAStringSet, or a ",
         "path to a FASTA file.", call. = FALSE)
  }

  if (is.null(names(out)) || any(!nzchar(names(out))))
    stop("`reference` must be named — the names are used as MDP identifiers.",
         call. = FALSE)

  # FASTA headers keep only the first whitespace-delimited token
  names(out) <- sub("^(\\S+).*$", "\\1", names(out))
  out <- toupper(gsub("[^A-Za-z*]", "", out))
  out <- sub("[*]+$", "", out)

  if (anyDuplicated(names(out)))
    stop("`reference` has duplicate names: ",
         paste(unique(names(out)[duplicated(names(out))]), collapse = ", "),
         call. = FALSE)
  out
}

#' Label ORFs With the Names of Known Mitochondrial-Derived Peptides
#'
#' Attaches the identity of a known mitochondrial-derived peptide (MDP) to each
#' row of an ORF table, so that downstream summaries and figures -- notably
#' [plot_rna_structure()] -- can name the ORFs they draw instead of showing them
#' as anonymous spans.
#'
#' Two matching methods are offered because the two situations they serve are
#' genuinely different.  Within the species the reference peptides came from,
#' an ORF either *is* a given MDP or it is not, and `method = "exact"` settles
#' that definitionally with no thresholds to tune.  Across species the peptides
#' have diverged and only a scored alignment can decide, which is what
#' `method = "homology"` provides.
#'
#' @param orfs `data.frame` of ORFs -- typically from [scan_orfs()],
#'   [find_orfs()] or [collapse_nested_orfs()].  Must contain a
#'   `protein_sequence` column.  A trailing stop-codon `*` is ignored.
#' @param reference Peptide sequences to match against.  A named `character`
#'   vector, an [Biostrings::AAStringSet], or a path to a FASTA file.  Names are
#'   used as MDP identifiers.  **Required -- there is no default**, because which
#'   peptide set a screen was run against belongs in the record.  [mdp_sequences()]
#'   supplies the human set (Swiss-Prot sequences, each verified to be encoded in
#'   the rCRS under the standard genetic code); pass it explicitly to use it.
#' @param method Either `"exact"` (default) or `"homology"`.  `"exact"` matches
#'   by string identity, including the containment cases where the ORF extends
#'   past the annotated peptide or vice versa.  `"homology"` runs
#'   [score_mdp_similarity()] and keeps hits passing `min_identity` and
#'   `min_coverage`.
#' @param min_identity,min_coverage Numeric (0--100).  Thresholds applied when
#'   `method = "homology"`; ignored for `"exact"`.  Defaults `60` and `60`.
#'   The permissive defaults of [score_mdp_similarity()] are deliberately *not*
#'   inherited: at 20% identity, unrelated peptides of MDP length routinely
#'   score as hits.
#' @param prefix Character scalar prepended to the names of the added columns.
#'   Default `"mdp"`.
#'
#' @return `orfs` with four columns added (named using `prefix`):
#' \describe{
#'   \item{`mdp_name`}{Name of the matched reference peptide, or `NA` when the
#'     ORF matched nothing.}
#'   \item{`mdp_match`}{How it matched: `"identical"`, `"orf_contains_mdp"`,
#'     `"mdp_contains_orf"`, `"homology"`, or `NA`.}
#'   \item{`mdp_pct_identity`}{Percent identity.  `100` for `"identical"` and
#'     for the containment cases (the shared region is exact); the alignment
#'     value for `"homology"`; `NA` for no match.}
#'   \item{`mdp_coverage`}{Percentage of the reference peptide accounted for by
#'     the ORF.  `NA` for no match.}
#' }
#' Row order and all existing columns are preserved.
#'
#' @details
#' When an ORF matches more than one reference peptide, the longest matching
#' reference wins under `"exact"` and the highest alignment score wins under
#' `"homology"`.  A single reference peptide may legitimately be assigned to
#' several ORFs -- overlapping and nested ORFs are common in rRNA genes -- so no
#' one-to-one constraint is imposed.
#'
#' Matching is on amino acid sequence alone and takes no account of coordinates,
#' so the caller remains responsible for having scanned the right molecule with
#' the right genetic code.  Human MDPs are translated on cytoplasmic 80S
#' ribosomes and therefore require the standard code (`"SGC0"`), not the
#' vertebrate mitochondrial code.
#'
#' @examples
#' \dontrun{
#' mdps <- c(Humanin = "MAPRGFSCLLLLTSEIDLPVKRRA",
#'           SHLP6   = "MLDQDIPMVQPLLKVRLFND")
#'
#' orfs <- scan_orfs(rnr2, min_orf_length = 60, genetic_code = "SGC0")
#' orfs <- annotate_mdp_orfs(orfs, reference = mdps)
#' orfs[!is.na(orfs$mdp_name), c("start", "end", "strand", "mdp_name")]
#'
#' # Named ORFs are drawn and labelled; the rest stay as context
#' plot_rna_structure(fold, orfs = orfs, offset = 1671)
#'
#' # Across species, where the peptides have diverged
#' annotate_mdp_orfs(killifish_orfs, reference = mdps, method = "homology",
#'                   min_identity = 50)
#' }
#'
#' @seealso [plot_rna_structure()], [score_mdp_similarity()], [scan_orfs()]
#' @importFrom Biostrings readAAStringSet
#' @importFrom stats setNames
#' @importFrom methods is
#' @export
annotate_mdp_orfs <- function(orfs,
                              reference    = .stop_no_mdp_reference(),
                              method       = c("exact", "homology"),
                              min_identity = 60,
                              min_coverage = 60,
                              prefix       = "mdp") {

  if (!is.data.frame(orfs))
    stop("`orfs` must be a data.frame.", call. = FALSE)
  if (!"protein_sequence" %in% names(orfs))
    stop("`orfs` must have a `protein_sequence` column.", call. = FALSE)

  method <- match.arg(method)
  ref    <- .resolve_mdp_reference(reference)

  col_name <- paste0(prefix, c("_name", "_match", "_pct_identity", "_coverage"))

  n <- nrow(orfs)
  hit_name <- rep(NA_character_, n)
  hit_how  <- rep(NA_character_, n)
  hit_pid  <- rep(NA_real_, n)
  hit_cov  <- rep(NA_real_, n)

  if (n == 0L || length(ref) == 0L) {
    orfs[[col_name[1L]]] <- hit_name; orfs[[col_name[2L]]] <- hit_how
    orfs[[col_name[3L]]] <- hit_pid;  orfs[[col_name[4L]]] <- hit_cov
    return(orfs)
  }

  prot <- toupper(sub("[*]+$", "", as.character(orfs$protein_sequence)))
  prot[is.na(prot)] <- ""

  if (method == "exact") {

    ref_len <- nchar(ref)
    ord     <- order(-ref_len)          # longest reference wins ties

    for (i in seq_len(n)) {
      p <- prot[i]
      if (!nzchar(p)) next
      for (k in ord) {
        s <- ref[[k]]
        if (!nzchar(s)) next
        if (identical(p, s)) {
          hit_name[i] <- names(ref)[k]; hit_how[i] <- "identical"
          hit_pid[i]  <- 100; hit_cov[i] <- 100
          break
        }
        if (grepl(s, p, fixed = TRUE)) {
          hit_name[i] <- names(ref)[k]; hit_how[i] <- "orf_contains_mdp"
          hit_pid[i]  <- 100; hit_cov[i] <- 100
          break
        }
        if (grepl(p, s, fixed = TRUE)) {
          hit_name[i] <- names(ref)[k]; hit_how[i] <- "mdp_contains_orf"
          hit_pid[i]  <- 100
          hit_cov[i]  <- round(100 * nchar(p) / nchar(s), 1)
          break
        }
      }
    }

  } else {

    min_identity <- as.numeric(min_identity)
    min_coverage <- as.numeric(min_coverage)
    if (is.na(min_identity) || min_identity < 0 || min_identity > 100)
      stop("`min_identity` must be between 0 and 100.", call. = FALSE)
    if (is.na(min_coverage) || min_coverage < 0 || min_coverage > 100)
      stop("`min_coverage` must be between 0 and 100.", call. = FALSE)

    keep <- which(nzchar(prot))
    if (length(keep) > 0L) {
      q <- data.frame(seq_id           = as.character(keep),
                      protein_sequence = prot[keep],
                      stringsAsFactors = FALSE)
      # score against `reference` only -- the built-in set must not leak in
      hits <- suppressMessages(
        score_mdp_similarity(q, min_identity = 0, custom_mdps = ref))
      hits <- hits[hits$mdp_name %in% names(ref), , drop = FALSE]
      hits <- hits[hits$pct_identity >= min_identity &
                   hits$mdp_coverage >= min_coverage, , drop = FALSE]

      if (nrow(hits) > 0L) {
        # score_mdp_similarity() builds query_id from seq_id when unclustered
        idx <- as.integer(sub("^.*?(\\d+)$", "\\1", hits$query_id))
        if (anyNA(idx)) idx <- match(hits$query_id, q$seq_id)
        for (r in order(hits$alignment_score)) {   # best score written last
          i <- keep[idx[r]]
          if (is.na(i)) next
          hit_name[i] <- hits$mdp_name[r]
          hit_how[i]  <- "homology"
          hit_pid[i]  <- hits$pct_identity[r]
          hit_cov[i]  <- hits$mdp_coverage[r]
        }
      }
    }
  }

  orfs[[col_name[1L]]] <- hit_name
  orfs[[col_name[2L]]] <- hit_how
  orfs[[col_name[3L]]] <- hit_pid
  orfs[[col_name[4L]]] <- hit_cov

  nm <- sum(!is.na(hit_name))
  message("annotate_mdp_orfs(): ", nm, "/", n, " ORF(s) matched a reference ",
          "peptide (", method, "); ",
          length(unique(stats::na.omit(hit_name))), " distinct MDP(s).")
  orfs
}
