# ============================================================================
# blast_orfs — remote BLASTP of ORF protein sequences via NCBI BLAST API
# ============================================================================

#' BLAST ORF Protein Sequences Against NCBI
#'
#' Submits protein sequences from a [find_conserved_orfs()] (or [scan_orfs()])
#' result to NCBI BLASTP and returns a tidy table of hits.  By default one
#' representative sequence per cluster is submitted (the longest), keeping
#' API usage minimal.
#'
#' @param orfs `data.frame` from [find_conserved_orfs()] or [scan_orfs()].
#'   Must contain a `protein_sequence` column.  If a `cluster_id` column is
#'   present, one representative per cluster is submitted.
#' @param database Character scalar.  NCBI BLAST database.  Default `"nr"`.
#'   Use `"swissprot"` for faster, curated hits.
#' @param max_hits Integer.  Maximum number of hits to return per query.
#'   Default `10L`.
#' @param evalue_threshold Numeric.  Maximum E-value to report.  Default `0.01`.
#' @param poll_interval Integer.  Seconds between status polls.  Default `15L`.
#' @param timeout Integer.  Maximum seconds to wait for results.  Default `300L`.
#'
#' @return A `data.frame` with columns:
#' \describe{
#'   \item{`query_id`}{Sequence identifier submitted (cluster_id + seq_id, or
#'     row index).}
#'   \item{`cluster_id`}{Cluster identifier, if present in input.}
#'   \item{`hit_accession`}{NCBI accession of the hit.}
#'   \item{`hit_description`}{Hit sequence description.}
#'   \item{`evalue`}{Alignment E-value.}
#'   \item{`pct_identity`}{Percent identity.}
#'   \item{`query_coverage`}{Fraction of query covered by the alignment.}
#'   \item{`hit_length`}{Length of the hit sequence.}
#' }
#'
#' @details
#' Uses the NCBI BLAST URL API.  An internet connection is required.  NCBI
#' requests that automated submissions include a contact email — set one with:
#' ```r
#' options(MitoDiveR.email = "you@example.com")
#' ```
#'
#' @seealso [find_conserved_orfs()], [scan_orfs()]
#'
#' @importFrom httr2 request req_url_query req_body_form req_perform
#'   resp_body_string
#' @export
blast_orfs <- function(orfs,
                       database        = "nr",
                       max_hits        = 10L,
                       evalue_threshold = 0.01,
                       poll_interval   = 15L,
                       timeout         = 300L) {

  if (!is.data.frame(orfs) || !"protein_sequence" %in% names(orfs))
    stop("'orfs' must be a data.frame with a 'protein_sequence' column.",
         call. = FALSE)

  has_clusters <- "cluster_id" %in% names(orfs)

  ## ---- select representatives ------------------------------------------------
  if (has_clusters) {
    # One representative per cluster: the row with the longest protein sequence
    orfs$..aa_len <- nchar(gsub("\\*", "", orfs$protein_sequence))
    reps <- do.call(rbind, lapply(split(orfs, orfs$cluster_id), function(cl) {
      cl[which.max(cl$..aa_len), , drop = FALSE]
    }))
    reps$..aa_len <- NULL
    orfs$..aa_len <- NULL
  } else {
    reps <- orfs
  }

  n <- nrow(reps)
  message("Submitting ", n, " sequence(s) to NCBI BLASTP (database: ", database, ")...")

  ## ---- build FASTA -----------------------------------------------------------
  seqs_clean <- gsub("\\*", "", reps$protein_sequence)

  if (has_clusters) {
    query_ids <- paste0("cluster_", reps$cluster_id, "_", reps$seq_id)
  } else {
    query_ids <- paste0("seq_", seq_len(n))
  }

  fasta_lines <- unlist(mapply(
    function(id, sq) c(paste0(">", id), sq),
    query_ids, seqs_clean,
    SIMPLIFY = FALSE, USE.NAMES = FALSE
  ))
  fasta_str <- paste(fasta_lines, collapse = "\n")

  email <- getOption("MitoDiveR.email", default = "mitodiver@example.com")

  ## ---- submit ----------------------------------------------------------------
  submit_resp <- httr2::req_perform(
    httr2::req_body_form(
      httr2::request("https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"),
      CMD          = "Put",
      PROGRAM      = "blastp",
      DATABASE     = database,
      QUERY        = fasta_str,
      HITLIST_SIZE = as.character(max_hits),
      EXPECT       = as.character(evalue_threshold),
      FORMAT_TYPE  = "XML",
      EMAIL        = email,
      TOOL         = "MitoDiveR"
    )
  )

  submit_body <- httr2::resp_body_string(submit_resp)

  rid <- regmatches(submit_body,
                    regexpr("(?<=RID = )[A-Z0-9]+", submit_body, perl = TRUE))
  if (length(rid) == 0L || nchar(rid) == 0L)
    stop("Failed to obtain a BLAST RID. NCBI may be unavailable.", call. = FALSE)

  message("BLAST RID: ", rid, " — polling every ", poll_interval, "s...")

  ## ---- poll ------------------------------------------------------------------
  elapsed <- 0L
  xml_result <- NULL

  repeat {
    Sys.sleep(poll_interval)
    elapsed <- elapsed + poll_interval

    status_resp <- httr2::req_perform(
      httr2::req_url_query(
        httr2::request("https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"),
        CMD         = "Get",
        RID         = rid,
        FORMAT_TYPE = "XML"
      )
    )
    status_body <- httr2::resp_body_string(status_resp)

    if (grepl("Status=WAITING", status_body)) {
      message("  ...waiting (", elapsed, "s elapsed)")
      if (elapsed >= timeout)
        stop("BLAST timed out after ", timeout, "s. Try again later or increase 'timeout'.",
             call. = FALSE)
      next
    }

    if (grepl("Status=FAILED", status_body))
      stop("BLAST job failed (RID: ", rid, ").", call. = FALSE)

    if (grepl("Status=UNKNOWN", status_body))
      stop("BLAST RID expired or unknown (RID: ", rid, ").", call. = FALSE)

    # Status=READY
    xml_result <- status_body
    break
  }

  message("Results received. Parsing...")

  ## ---- parse XML -------------------------------------------------------------
  .parse_blast_xml(xml_result, query_ids,
                   if (has_clusters) reps$cluster_id else NULL,
                   evalue_threshold)
}


# --------------------------------------------------------------------------- #
# Internal: parse BLAST XML output                                             #
# --------------------------------------------------------------------------- #
.parse_blast_xml <- function(xml_text, query_ids, cluster_ids, evalue_threshold) {

  iterations <- .xml_extract_all(xml_text, "Iteration")
  if (length(iterations) == 0L) {
    message("No BLAST hits found.")
    return(invisible(data.frame(
      query_id        = character(),
      cluster_id      = if (!is.null(cluster_ids)) integer() else NULL,
      hit_accession   = character(),
      hit_description = character(),
      evalue          = numeric(),
      pct_identity    = numeric(),
      query_coverage  = numeric(),
      hit_length      = integer(),
      stringsAsFactors = FALSE
    )))
  }

  rows <- list()

  for (i in seq_along(iterations)) {
    iter  <- iterations[[i]]
    qid   <- query_ids[min(i, length(query_ids))]
    clust <- if (!is.null(cluster_ids)) cluster_ids[min(i, length(cluster_ids))] else NA_integer_

    # Query length (for coverage calculation)
    qlen_str <- .xml_extract_first(iter, "Iteration_query-len")
    qlen <- if (nchar(qlen_str) > 0L) as.integer(qlen_str) else NA_integer_

    hits <- .xml_extract_all(iter, "Hit")
    if (length(hits) == 0L) next

    for (hit in hits) {
      accession   <- .xml_extract_first(hit, "Hit_accession")
      description <- .xml_extract_first(hit, "Hit_def")
      hit_len     <- as.integer(.xml_extract_first(hit, "Hit_len"))

      hsps <- .xml_extract_all(hit, "Hsp")
      if (length(hsps) == 0L) next
      hsp <- hsps[[1L]]   # top HSP only

      evalue   <- as.numeric(.xml_extract_first(hsp, "Hsp_evalue"))
      if (!is.na(evalue) && evalue > evalue_threshold) next

      identity <- as.integer(.xml_extract_first(hsp, "Hsp_identity"))
      aln_len  <- as.integer(.xml_extract_first(hsp, "Hsp_align-len"))
      q_from   <- as.integer(.xml_extract_first(hsp, "Hsp_query-from"))
      q_to     <- as.integer(.xml_extract_first(hsp, "Hsp_query-to"))

      pct_id  <- if (!is.na(identity) && !is.na(aln_len) && aln_len > 0L)
                   round(identity / aln_len * 100, 1) else NA_real_
      coverage <- if (!is.na(q_from) && !is.na(q_to) && !is.na(qlen) && qlen > 0L)
                    round((q_to - q_from + 1L) / qlen * 100, 1) else NA_real_

      rows[[length(rows) + 1L]] <- data.frame(
        query_id        = qid,
        cluster_id      = clust,
        hit_accession   = accession,
        hit_description = description,
        evalue          = evalue,
        pct_identity    = pct_id,
        query_coverage  = coverage,
        hit_length      = hit_len,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0L) {
    message("No hits below the E-value threshold.")
    return(invisible(data.frame(stringsAsFactors = FALSE)))
  }

  out <- do.call(rbind, rows)
  if (!is.null(cluster_ids)) out <- out[order(out$cluster_id, out$evalue), ]
  rownames(out) <- NULL
  out
}


# --------------------------------------------------------------------------- #
# Minimal XML helpers (no xml2 dependency)                                    #
# --------------------------------------------------------------------------- #
.xml_extract_all <- function(text, tag) {
  pattern <- sprintf("(?s)<%s[^>]*>(.*?)</%s>", tag, tag)
  m <- gregexpr(pattern, text, perl = TRUE)
  regmatches(text, m)[[1L]]
}

.xml_extract_first <- function(text, tag) {
  pattern <- sprintf("(?s)<%s[^>]*>(.*?)</%s>", tag, tag)
  m <- regexpr(pattern, text, perl = TRUE)
  if (m == -1L) return("")
  inner <- regmatches(text, m)
  # Strip outer tags
  gsub(sprintf("^<%s[^>]*>|</%s>$", tag, tag), "", inner)
}
