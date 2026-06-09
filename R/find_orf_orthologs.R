# ============================================================================
# find_orf_orthologs — locate the ortholog of a known noncanonical / MDP ORF
# across unannotated mitogenomes ("ORFologs")
# ============================================================================

#' Find ORF Orthologs ("ORFologs") of a Known Peptide Across Taxa
#'
#' Locates the ortholog of a known noncanonical open reading frame — typically a
#' mitochondria-derived peptide (MDP) — in other taxa whose mitogenomes are not
#' annotated for such ORFs.  Because these ORFs are never in the feature table,
#' the search is two-step: discover candidate coding stretches, then translate
#' them in all six frames ([scan_orfs()]) and match each ORF to the query
#' peptide ([score_mdp_similarity()]).
#'
#' Three discovery strategies (`method`) feed the same local confirmation core,
#' so the output table is identical regardless of how candidates were found:
#' \describe{
#'   \item{`"local"`}{No remote search.  Scan whole RefSeq mitogenomes for a
#'     clade (`taxa`) or a set you supply (`genomes`).  Each mitogenome is only
#'     ~16 kb, so scanning the entire genome is cheap; orthology is confirmed
#'     positionally with [annotate_genomic_regions()].  Best when the taxonomic
#'     scope is already known.}
#'   \item{`"tblastn"`}{Remote, peptide-level discovery.  The `query` peptide is
#'     submitted as a TBLASTN search (protein vs. translated nucleotide), the
#'     most sensitive option across deep divergence.  Hit subject ranges are
#'     harvested (± `flank`) and confirmed locally.}
#'   \item{`"blastn"`}{Remote, host-window discovery.  A conserved nucleotide
#'     `host_region` the ORF sits within (e.g. a stretch of 12S/16S rRNA) is
#'     submitted as a BLASTN search.  Anchoring on the conserved host locus is
#'     itself the positional/orthology filter.}
#' }
#'
#' For the remote methods the search space is constrained to the RefSeq
#' mitochondrion set (`refseq_only`), which collapses taxa with vast sequence
#' holdings (e.g. *Danio*) to a single curated ~16 kb mitogenome and dedupes to
#' one record per species — the same constraint that makes orthology meaningful.
#'
#' @param query The known ORF's amino-acid sequence as a character scalar.  May
#'   be named (the name labels the query in the output); otherwise labelled
#'   `"query_orf"`.
#' @param taxa Optional character vector of taxon names (e.g.
#'   `c("Sebastes", "Danio")`) or NCBI Taxonomy IDs.  For `method = "local"`
#'   their RefSeq mitogenomes are fetched and scanned; for the remote methods
#'   they scope the BLAST search via `ENTREZ_QUERY`.
#' @param genomes Optional sequences to scan directly under `method = "local"`,
#'   in any form accepted by the package input resolver (a `DNAStringSet`, NCBI
#'   accessions, a FASTA file, or raw text).  Bypasses NCBI discovery.
#' @param host_region Nucleotide sequence of the conserved host locus the ORF is
#'   embedded in.  Required for `method = "blastn"` (used as the BLAST query);
#'   ignored otherwise.
#' @param method Discovery strategy; one of `"local"`, `"tblastn"`, `"blastn"`.
#' @param genetic_code Character scalar genetic-code id passed to [scan_orfs()].
#'   Required (no default).  For MDPs use `"SGC0"` (the cytoplasmic 80S ribosome
#'   translates them), not the mitochondrial `"SGC1"`.
#' @param start_codons Character vector of recognised start codons.  Default is
#'   the permissive mitochondrial set `c("ATG","GTG","ATA","ATT","ATC")`.
#' @param min_orf_length Integer.  Minimum ORF length (nt) for [scan_orfs()].
#'   Default `30L` (~10 aa; MDPs are short).
#' @param min_identity Numeric (0–100).  Minimum percent identity of the
#'   ORF-to-query alignment to report a hit.  Default `40`.
#' @param min_coverage Numeric (0–100).  Minimum percent of the query peptide
#'   covered by the alignment.  Default `50`.  Suppresses short spurious motif
#'   matches.
#' @param refseq_only Logical.  Constrain remote searches (and `taxa` discovery)
#'   to RefSeq mitochondrion records.  Default `TRUE`.
#' @param max_genomes Integer.  Cap on genomes fetched (`taxa`) or hit subject
#'   ranges harvested (remote methods).  Default `50L`.
#' @param flank Integer.  Nucleotides added either side of a harvested BLAST hit
#'   range before scanning, so a start codon just outside the aligned segment is
#'   not missed.  Default `150L`.  Remote methods only.
#' @param annotate Logical.  For accession-backed genomes, annotate each hit's
#'   `genomic_region` via [fetch_mito_genbank()] + [annotate_genomic_regions()].
#'   Default `TRUE`.  Silently skipped when features are unavailable.
#' @param api_key Optional NCBI API key (raises the Entrez rate limit).
#' @param poll_interval,timeout Integers.  BLAST poll cadence and overall
#'   timeout (seconds) for the remote methods.  Defaults `15L` and `300L`.
#'
#' @return A `data.frame` (invisibly empty if nothing is found) with one row per
#'   confirmed ortholog ORF:
#' \describe{
#'   \item{`query_name`}{Label of the query peptide.}
#'   \item{`seq_id`}{Source genome / window identifier (accession when known).}
#'   \item{`start`,`end`,`strand`,`frame`}{ORF coordinates on `seq_id`.}
#'   \item{`length_nt`,`protein_sequence`}{ORF length and translation.}
#'   \item{`pct_identity`}{Percent identity of the ORF-to-query alignment.}
#'   \item{`query_coverage`}{Percent of the known query peptide covered by the
#'     alignment (the orthology-relevant metric `min_coverage` gates on).}
#'   \item{`orf_coverage`}{Percent of the discovered ORF covered by the
#'     alignment.}
#'   \item{`genomic_region`}{Host-feature label when annotation succeeded, else
#'     `NA` (remote windows carry their BLAST provenance instead).}
#'   \item{`method`,`blast_accession`,`blast_from`,`blast_to`,`blast_evalue`}{
#'     Discovery provenance; BLAST columns are `NA` under `method = "local"`.}
#' }
#'
#' @details
#' This separates *homology* (the peptide aligns) from *orthology* (it aligns at
#' the corresponding locus).  Under `"local"`/`"blastn"` the positional anchor is
#' explicit — the annotated host region, or the conserved host window searched.
#' Under `"tblastn"` confirm orthology by inspecting `blast_accession`/`blast_*`
#' against the expected host gene.  Remote methods require network access and
#' degrade with a clear error when NCBI is unreachable.
#'
#' @seealso [scan_orfs()], [score_mdp_similarity()], [annotate_genomic_regions()],
#'   [find_clade_conserved_orfs()], [blast_orfs()]
#'
#' @examples
#' \dontrun{
#' humanin <- c(Humanin = "MAPRGFSCLLLLTSEIDLPVKRRA")
#'
#' # Defined clade: scan whole RefSeq mitogenomes, confirm positionally
#' hits <- find_orf_orthologs(humanin, taxa = "Sebastes",
#'                            method = "local", genetic_code = "SGC0")
#'
#' # Open-ended discovery: peptide-level TBLASTN across RefSeq mitogenomes
#' hits <- find_orf_orthologs(humanin, method = "tblastn",
#'                            genetic_code = "SGC0")
#'
#' # Host-anchored discovery: BLASTN a conserved 16S window the ORF rides on
#' hits <- find_orf_orthologs(humanin, host_region = my_16s_window,
#'                            method = "blastn", genetic_code = "SGC0")
#' }
#'
#' @importFrom httr2 request req_url_query req_body_form req_perform
#'   resp_body_string
#' @importFrom rentrez entrez_search entrez_summary entrez_fetch set_entrez_key
#' @importFrom Biostrings DNAStringSet
#' @export
find_orf_orthologs <- function(query,
                               taxa           = NULL,
                               genomes        = NULL,
                               host_region    = NULL,
                               method         = c("local", "tblastn", "blastn"),
                               genetic_code   = .stop_no_code("genetic_code"),
                               start_codons   = c("ATG", "GTG", "ATA", "ATT", "ATC"),
                               min_orf_length = 30L,
                               min_identity   = 40,
                               min_coverage   = 50,
                               refseq_only    = TRUE,
                               max_genomes    = 50L,
                               flank          = 150L,
                               annotate       = TRUE,
                               api_key        = NULL,
                               poll_interval  = 15L,
                               timeout        = 300L) {

  method <- match.arg(method)

  ## ---- normalise the query peptide -------------------------------------------
  if (!is.character(query) || length(query) != 1L || is.na(query) ||
      nchar(query) == 0L)
    stop("'query' must be a single non-empty amino-acid string.", call. = FALSE)
  query_name <- if (!is.null(names(query)) && nzchar(names(query)))
                  names(query) else "query_orf"
  query_aa   <- toupper(unname(query))

  if (!is.null(api_key)) rentrez::set_entrez_key(api_key)
  max_genomes <- if (is.infinite(max_genomes)) 1e6L else as.integer(max_genomes)

  ## ---- gather candidate sequences per method ---------------------------------
  prov <- NULL          # per-seq BLAST provenance (NULL for local)

  if (method == "local") {
    if (is.null(genomes) && is.null(taxa))
      stop("method = \"local\" needs either 'genomes' or 'taxa'.", call. = FALSE)
    seqs <- if (!is.null(genomes)) {
      # .resolve_dna_input() yields a named list of DNAString; scan_orfs()
      # wants a DNAStringSet.
      lst <- .resolve_dna_input(genomes, api_key = api_key)
      Biostrings::DNAStringSet(vapply(lst, as.character, character(1L)))
    } else {
      .fetch_clade_mitogenomes(taxa, refseq_only, max_genomes)
    }
    circular <- TRUE

  } else {
    # Remote discovery -> harvest windows -> confirm locally (linear windows)
    blast_query <- if (method == "blastn") {
      if (is.null(host_region) || !nzchar(host_region))
        stop("method = \"blastn\" requires 'host_region' (nucleotide).",
             call. = FALSE)
      toupper(host_region)
    } else {
      query_aa
    }
    harvest <- .blast_harvest(
      program      = method,
      blast_query  = blast_query,
      taxa         = taxa,
      refseq_only  = refseq_only,
      max_hits     = max_genomes,
      flank        = flank,
      poll_interval = poll_interval,
      timeout      = timeout
    )
    if (is.null(harvest) || length(harvest$seqs) == 0L) {
      message("No BLAST hits to confirm.")
      return(invisible(.empty_ortholog_table()))
    }
    seqs     <- harvest$seqs
    prov     <- harvest$provenance
    circular <- FALSE
  }

  ## ---- confirm: scan ORFs, score against the query, annotate -----------------
  hits <- .confirm_orf_orthologs(
    query_name     = query_name,
    query_aa       = query_aa,
    seqs           = seqs,
    genetic_code   = genetic_code,
    start_codons   = start_codons,
    min_orf_length = min_orf_length,
    min_identity   = min_identity,
    min_coverage   = min_coverage,
    circular       = circular
  )
  if (nrow(hits) == 0L) {
    message("No ORF in the candidate set matched the query at the given thresholds.")
    return(invisible(.empty_ortholog_table()))
  }

  ## ---- positional orthology: annotate host region (accession-backed only) ----
  hits$genomic_region <- NA_character_
  if (annotate && method != "tblastn")
    hits <- .annotate_hits_by_accession(hits, api_key = api_key)

  ## ---- attach provenance -----------------------------------------------------
  hits$method <- method
  if (is.null(prov)) {
    hits$blast_accession <- NA_character_
    hits$blast_from      <- NA_integer_
    hits$blast_to        <- NA_integer_
    hits$blast_evalue    <- NA_real_
  } else {
    m <- match(hits$seq_id, prov$seq_id)
    hits$blast_accession <- prov$accession[m]
    hits$blast_from      <- prov$hit_from[m]
    hits$blast_to        <- prov$hit_to[m]
    hits$blast_evalue    <- prov$evalue[m]
  }

  col_order <- c("query_name", "seq_id", "start", "end", "strand", "frame",
                 "length_nt", "protein_sequence", "pct_identity",
                 "query_coverage", "orf_coverage", "genomic_region", "method",
                 "blast_accession", "blast_from", "blast_to", "blast_evalue")
  hits <- hits[, intersect(col_order, names(hits)), drop = FALSE]
  rownames(hits) <- NULL
  hits[order(-hits$pct_identity, hits$seq_id), , drop = FALSE]
}


# --------------------------------------------------------------------------- #
# Internal: local confirmation core (scan -> score -> tidy). Offline-testable. #
# --------------------------------------------------------------------------- #
.confirm_orf_orthologs <- function(query_name, query_aa, seqs, genetic_code,
                                   start_codons, min_orf_length, min_identity,
                                   min_coverage, circular) {

  orfs <- scan_orfs(seqs,
                    min_orf_length = min_orf_length,
                    genetic_code   = genetic_code,
                    start_codons   = start_codons,
                    both_strands   = TRUE,
                    circular       = circular)
  if (nrow(orfs) == 0L) return(.empty_ortholog_table())
  # Collapse nested ATG starts to the longest ORF per (frame, stop), so a locus
  # is reported once rather than as every internal-start fragment.
  orfs <- collapse_nested_orfs(orfs)

  scored <- score_mdp_similarity(
    orfs,
    min_identity = min_identity,
    custom_mdps  = stats::setNames(query_aa, query_name)
  )
  # score_mdp_similarity returns invisible empty data.frame when nothing matched
  if (!is.data.frame(scored) || nrow(scored) == 0L ||
      !"mdp_name" %in% names(scored))
    return(.empty_ortholog_table())

  # Gate on coverage of the KNOWN query peptide. In score_mdp_similarity's
  # frame the peptide is the reference, so that is `mdp_coverage`;
  # `query_coverage` there is coverage of the candidate ORF.
  keep <- scored$mdp_name == query_name & scored$mdp_coverage >= min_coverage
  scored <- scored[keep, , drop = FALSE]
  if (nrow(scored) == 0L) return(.empty_ortholog_table())

  # No cluster_id present, so score_mdp_similarity labels queries "seq_<row>",
  # positionally indexing the (unclustered) ORF table. Map back by that index.
  idx <- as.integer(sub("^seq_", "", scored$query_id))
  if (anyNA(idx) || max(idx) > nrow(orfs))
    stop("Internal: could not map scored hits back to ORF rows.", call. = FALSE)

  data.frame(
    query_name       = query_name,
    seq_id           = orfs$seq_id[idx],
    start            = orfs$start[idx],
    end              = orfs$end[idx],
    strand           = orfs$strand[idx],
    frame            = orfs$frame[idx],
    length_nt        = orfs$length_nt[idx],
    protein_sequence = orfs$protein_sequence[idx],
    pct_identity     = scored$pct_identity,
    query_coverage   = scored$mdp_coverage,    # % of the known query peptide
    orf_coverage     = scored$query_coverage,  # % of the discovered ORF
    stringsAsFactors = FALSE
  )
}


# --------------------------------------------------------------------------- #
# Internal: empty result skeleton                                             #
# --------------------------------------------------------------------------- #
.empty_ortholog_table <- function() {
  data.frame(
    query_name = character(), seq_id = character(),
    start = integer(), end = integer(), strand = character(),
    frame = integer(), length_nt = integer(), protein_sequence = character(),
    pct_identity = numeric(), query_coverage = numeric(),
    orf_coverage = numeric(), stringsAsFactors = FALSE
  )
}


# --------------------------------------------------------------------------- #
# Internal: fetch RefSeq mitogenomes for a set of taxa (one search per taxon)  #
# --------------------------------------------------------------------------- #
.fetch_clade_mitogenomes <- function(taxa, refseq_only, max_genomes) {
  accs <- character()
  for (tx in taxa) {
    term_parts <- c(sprintf("%s[ORGN]", tx), "mitochondrion[FILT]")
    if (refseq_only) term_parts <- c(term_parts, "srcdb_refseq[PROP]")
    term <- paste(term_parts, collapse = " AND ")
    res  <- rentrez::entrez_search(db = "nucleotide", term = term,
                                   retmax = max_genomes)
    if (length(res$ids) == 0L) {
      message("No RefSeq mitogenome hits for taxon '", tx, "'.")
      next
    }
    s <- rentrez::entrez_summary(db = "nucleotide", id = res$ids)
    if (length(res$ids) == 1L) s <- list(s)
    accs <- c(accs, vapply(s, function(x) x$accessionversion, character(1L)))
  }
  accs <- unique(accs)
  if (length(accs) == 0L)
    stop("No RefSeq mitogenomes found for the supplied taxa.", call. = FALSE)
  if (length(accs) > max_genomes) accs <- accs[seq_len(max_genomes)]
  fetch_sequences(accessions = accs)
}


# --------------------------------------------------------------------------- #
# Internal: remote BLAST discovery -> harvest windowed subsequences            #
# Reuses the QBlast submit/poll pattern and the XML helpers in blast_orfs.R.   #
# --------------------------------------------------------------------------- #
.blast_harvest <- function(program, blast_query, taxa, refseq_only,
                           max_hits, flank, poll_interval, timeout) {

  ## ---- ENTREZ_QUERY scoping --------------------------------------------------
  ent <- "mitochondrion[filter]"
  if (refseq_only) ent <- paste(ent, "AND refseq[filter]")
  if (!is.null(taxa)) {
    orgs <- paste(sprintf("%s[Organism]", taxa), collapse = " OR ")
    ent  <- paste0(ent, " AND (", orgs, ")")
  }

  email <- getOption("MitoDiveR.email", default = "mitodiver@example.com")

  ## ---- submit ----------------------------------------------------------------
  message("Submitting ", toupper(program), " query to NCBI (scope: ", ent, ")...")
  submit_resp <- httr2::req_perform(httr2::req_body_form(
    httr2::request("https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"),
    CMD          = "Put",
    PROGRAM      = program,
    DATABASE     = "nt",
    QUERY        = blast_query,
    ENTREZ_QUERY = ent,
    HITLIST_SIZE = as.character(max_hits),
    FORMAT_TYPE  = "XML",
    EMAIL        = email,
    TOOL         = "MitoDiveR"
  ))
  submit_body <- httr2::resp_body_string(submit_resp)
  rid <- regmatches(submit_body,
                    regexpr("(?<=RID = )[A-Z0-9]+", submit_body, perl = TRUE))
  if (length(rid) == 0L || nchar(rid) == 0L)
    stop("Failed to obtain a BLAST RID. NCBI may be unavailable.", call. = FALSE)
  message("BLAST RID: ", rid, " — polling every ", poll_interval, "s...")

  ## ---- poll ------------------------------------------------------------------
  elapsed <- 0L; xml_result <- NULL
  repeat {
    Sys.sleep(poll_interval)
    elapsed <- elapsed + poll_interval
    status_body <- httr2::resp_body_string(httr2::req_perform(httr2::req_url_query(
      httr2::request("https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"),
      CMD = "Get", RID = rid, FORMAT_TYPE = "XML")))
    if (grepl("Status=WAITING", status_body)) {
      if (elapsed >= timeout)
        stop("BLAST timed out after ", timeout, "s.", call. = FALSE)
      next
    }
    if (grepl("Status=FAILED", status_body))
      stop("BLAST job failed (RID: ", rid, ").", call. = FALSE)
    if (grepl("Status=UNKNOWN", status_body))
      stop("BLAST RID expired or unknown (RID: ", rid, ").", call. = FALSE)
    xml_result <- status_body; break
  }

  ## ---- parse best HSP per subject accession ----------------------------------
  hits <- .parse_blast_subjects(xml_result)
  if (nrow(hits) == 0L) return(NULL)
  if (nrow(hits) > max_hits) hits <- hits[seq_len(max_hits), , drop = FALSE]

  ## ---- harvest windowed subsequences via efetch ------------------------------
  seq_list <- list(); prov_rows <- list()
  for (k in seq_len(nrow(hits))) {
    acc <- hits$accession[k]
    lo  <- max(1L, min(hits$hit_from[k], hits$hit_to[k]) - flank)
    hi  <- max(hits$hit_from[k], hits$hit_to[k]) + flank
    fasta <- tryCatch(
      rentrez::entrez_fetch(db = "nucleotide", id = acc, rettype = "fasta",
                            retmode = "text", seq_start = lo, seq_stop = hi),
      error = function(e) NULL)
    if (is.null(fasta) || !nzchar(fasta)) next
    win <- tryCatch(fetch_sequences(fasta_text = fasta), error = function(e) NULL)
    if (is.null(win) || length(win) == 0L) next
    sid <- paste0(acc, ":", lo, "-", hi)
    seq_list[[sid]] <- as.character(win[[1]])
    prov_rows[[length(prov_rows) + 1L]] <- data.frame(
      seq_id = sid, accession = acc, hit_from = lo, hit_to = hi,
      evalue = hits$evalue[k], stringsAsFactors = FALSE)
  }
  if (length(seq_list) == 0L) return(NULL)

  list(
    seqs       = Biostrings::DNAStringSet(unlist(seq_list)),
    provenance = do.call(rbind, prov_rows)
  )
}


# --------------------------------------------------------------------------- #
# Internal: parse subject accession + best-HSP range from BLAST XML            #
# --------------------------------------------------------------------------- #
.parse_blast_subjects <- function(xml_text) {
  rows <- list()
  for (hit in .xml_extract_all(xml_text, "Hit")) {
    acc <- .xml_extract_first(hit, "Hit_accession")
    hsps <- .xml_extract_all(hit, "Hsp")
    if (length(hsps) == 0L) next
    hsp <- hsps[[1L]]   # top HSP
    rows[[length(rows) + 1L]] <- data.frame(
      accession = acc,
      hit_from  = as.integer(.xml_extract_first(hsp, "Hsp_hit-from")),
      hit_to    = as.integer(.xml_extract_first(hsp, "Hsp_hit-to")),
      evalue    = as.numeric(.xml_extract_first(hsp, "Hsp_evalue")),
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0L)
    return(data.frame(accession = character(), hit_from = integer(),
                      hit_to = integer(), evalue = numeric(),
                      stringsAsFactors = FALSE))
  out <- do.call(rbind, rows)
  # one window per accession: keep the lowest-evalue HSP
  out <- out[order(out$evalue), , drop = FALSE]
  out[!duplicated(out$accession), , drop = FALSE]
}


# --------------------------------------------------------------------------- #
# Internal: annotate each hit's genomic_region using its source accession      #
# --------------------------------------------------------------------------- #
.annotate_hits_by_accession <- function(hits, api_key) {
  # seq_id is an accession (local) or "acc:lo-hi" (blastn window); only the
  # bare-accession case can be annotated against the genome feature table.
  accs <- unique(hits$seq_id)
  for (acc in accs) {
    if (grepl(":", acc, fixed = TRUE)) next            # harvested window
    gb <- tryCatch(fetch_mito_genbank(acc, api_key = api_key),
                   error = function(e) NULL)
    if (is.null(gb) || is.null(gb[[acc]]$features)) next
    sel <- hits$seq_id == acc
    sub <- hits[sel, , drop = FALSE]
    ann <- tryCatch(
      annotate_genomic_regions(sub, gb[[acc]]$features),
      error = function(e) NULL)
    if (!is.null(ann) && "genomic_region" %in% names(ann))
      hits$genomic_region[sel] <- ann$genomic_region
  }
  hits
}
