# ============================================================================
# find_clade_conserved_orfs — taxon-scope conserved ORF discovery
# ============================================================================

#' Discover Conserved ORFs Across RefSeq Mitogenomes in a Taxon
#'
#' Pulls all RefSeq complete mitogenomes annotated to a user-supplied
#' taxon (genus / family / order), scans each for ORFs, and clusters
#' the result with [find_conserved_orfs()] to identify ORFs shared
#' across the clade.  Intended as a one-call entry point for clade-
#' level survey work.
#'
#' @param taxon Character scalar.  Taxon name (e.g. `"Sebastes"`,
#'   `"Nothobranchius"`, `"Primates"`) or an NCBI Taxonomy ID as a
#'   character (e.g. `"6790"`).  Passed verbatim to [rentrez] as
#'   `<taxon>[ORGN]`.
#' @param max_genomes Integer.  Cap on the number of mitogenomes to
#'   fetch and process.  Default `20L`.  Set to `Inf` for no cap.
#' @param refseq_only Logical.  Restrict the query to RefSeq deposits
#'   (`srcdb_refseq[PROP]`).  Default `TRUE`.  Set `FALSE` to include
#'   non-RefSeq GenBank mitogenomes — useful for clades with sparse
#'   RefSeq coverage but introduces redundancy from multiple
#'   submissions of the same species.
#' @param min_length Integer.  Minimum sequence length (nt) for a hit
#'   to be considered a complete mitogenome.  Default `13000L`
#'   (excludes partial vertebrate mtDNA deposits).
#' @param min_orf_length Integer.  Passed to [scan_orfs()].  Default
#'   `60L`.
#' @param genetic_code Character scalar.  Genetic code identifier
#'   passed to [scan_orfs()].  Required (no default) — e.g. `"SGC1"`
#'   (vertebrate mitochondrial).
#' @param start_codons Character vector of recognised start codons.
#'   Default `c("ATG","GTG","ATA","ATT","ATC")`.
#' @param min_identity Numeric in (0, 100].  Minimum percent amino-acid
#'   identity for [find_conserved_orfs()] clustering.  Default `50`.
#' @param min_sequences Numeric.  Minimum fraction of sequences (or
#'   absolute count if `>= 1`) in which a cluster must be present.
#'   Default `0.5`.
#' @param api_key Optional NCBI API key (character).  If supplied, set
#'   via [rentrez::set_entrez_key()] for the duration of the call to
#'   raise the request limit from 3 to 10 requests per second.
#'
#' @return A list with components:
#'   \describe{
#'     \item{`taxon`}{The taxon string supplied.}
#'     \item{`accessions`}{Character vector of accessions fetched.}
#'     \item{`metadata`}{`data.frame` with one row per accession
#'       (`acc`, `title`, `length`).}
#'     \item{`all_orfs`}{The combined [scan_orfs()] table.}
#'     \item{`conserved`}{The [find_conserved_orfs()] clustered table.}
#'   }
#'
#' @details
#' Network access is required.  NCBI rate-limits unauthenticated
#' Entrez queries to three requests per second; set an API key via
#' `api_key` for higher throughput on large clades.  When `max_genomes`
#' caps the search, the first `max_genomes` accessions returned by the
#' Entrez query are used (typically ordered by deposit date,
#' newest-first).  For deterministic taxon sampling, supply a curated
#' accession list to [scan_orfs()] + [find_conserved_orfs()] directly.
#'
#' @seealso [scan_orfs()], [find_conserved_orfs()], [fetch_sequences()]
#'
#' @examples
#' \dontrun{
#' # Conserved ORFs across Sebastes rockfishes
#' res <- find_clade_conserved_orfs("Sebastes", max_genomes = 20)
#' res$metadata
#' head(res$conserved)
#'
#' # Larger primate survey under permissive RefSeq-or-GenBank query
#' primates <- find_clade_conserved_orfs(
#'   "Primates",
#'   max_genomes = 30,
#'   refseq_only = FALSE
#' )
#' }
#'
#' @importFrom rentrez entrez_search entrez_summary set_entrez_key
#' @export
find_clade_conserved_orfs <- function(taxon,
                                      max_genomes    = 20L,
                                      refseq_only    = TRUE,
                                      min_length     = 13000L,
                                      min_orf_length = 60L,
                                      genetic_code   = .stop_no_code("genetic_code"),
                                      start_codons   = c("ATG","GTG","ATA","ATT","ATC"),
                                      min_identity   = 50,
                                      min_sequences  = 0.5,
                                      api_key        = NULL) {

  if (!is.character(taxon) || length(taxon) != 1L || nchar(taxon) == 0L)
    stop("'taxon' must be a single non-empty character string.", call. = FALSE)
  if (!is.null(api_key)) rentrez::set_entrez_key(api_key)
  max_genomes <- if (is.infinite(max_genomes)) 1e6L else as.integer(max_genomes)

  ## ---- 1. Search NCBI nucleotide for mitogenomes in the taxon --------------
  term_parts <- c(
    sprintf("%s[ORGN]", taxon),
    "mitochondrion[FILT]"
  )
  if (refseq_only) term_parts <- c(term_parts, "srcdb_refseq[PROP]")
  term <- paste(term_parts, collapse = " AND ")

  message(sprintf("Searching NCBI nucleotide: %s", term))
  search_res <- rentrez::entrez_search(
    db        = "nucleotide",
    term      = term,
    retmax    = max_genomes,
    use_history = FALSE
  )

  if (length(search_res$ids) == 0L)
    stop(sprintf("No NCBI nucleotide hits for taxon '%s' under term '%s'.",
                 taxon, term), call. = FALSE)
  message(sprintf("Found %d hit(s); fetching summaries for filtering.",
                  length(search_res$ids)))

  ## ---- 2. Fetch summaries, filter by length --------------------------------
  s <- rentrez::entrez_summary(db = "nucleotide", id = search_res$ids)
  if (length(search_res$ids) == 1L) s <- list(s)
  meta <- data.frame(
    acc    = vapply(s, function(x) x$accessionversion, character(1L)),
    title  = vapply(s, function(x) x$title,            character(1L)),
    length = vapply(s, function(x) as.integer(x$slen), integer(1L)),
    stringsAsFactors = FALSE
  )
  meta <- meta[meta$length >= min_length, , drop = FALSE]
  if (nrow(meta) == 0L)
    stop(sprintf("No mitogenomes >= %d nt for taxon '%s'.",
                 min_length, taxon), call. = FALSE)
  if (nrow(meta) > max_genomes) meta <- meta[seq_len(max_genomes), , drop = FALSE]

  message(sprintf("Retained %d mitogenome(s) after length filter (>= %d nt).",
                  nrow(meta), min_length))

  ## ---- 3. Fetch sequences --------------------------------------------------
  message("Fetching FASTA sequences...")
  seqs <- fetch_sequences(accessions = meta$acc)

  ## ---- 4. Scan ORFs --------------------------------------------------------
  message(sprintf("Scanning ORFs (code %s, min %d nt)...",
                  genetic_code, min_orf_length))
  all_orfs <- scan_orfs(seqs,
                        min_orf_length = min_orf_length,
                        genetic_code   = genetic_code,
                        start_codons   = start_codons,
                        both_strands   = TRUE,
                        circular       = TRUE)
  message(sprintf("  Found %d ORFs across %d genomes.",
                  nrow(all_orfs), nrow(meta)))

  ## ---- 5. Cluster conserved ORFs -------------------------------------------
  message(sprintf("Clustering conserved ORFs (>= %g%% identity, present in >= %s of sequences)...",
                  min_identity,
                  if (min_sequences < 1) sprintf("%g%%", 100 * min_sequences) else as.character(min_sequences)))
  conserved <- find_conserved_orfs(all_orfs,
                                   min_identity  = min_identity,
                                   min_sequences = min_sequences)
  n_clusters <- if (nrow(conserved) > 0L)
    length(unique(conserved$cluster_id)) else 0L
  message(sprintf("  Conserved clusters: %d", n_clusters))

  list(
    taxon       = taxon,
    accessions  = meta$acc,
    metadata    = meta,
    all_orfs    = all_orfs,
    conserved   = conserved
  )
}
