# ============================================================================
# fetch_gene — extract one mitochondrial gene across species from NCBI
# ============================================================================

#' Fetch a Mitochondrial Gene for One or More Species
#'
#' Retrieves a single mitochondrial gene (e.g. COI, CYTB, 16S) for a set of
#' species or accessions, cutting it out of each GenBank record using the
#' record's own feature annotation.  Species are looked up on NCBI by name, so
#' no accession numbers are needed.  The result is a [Biostrings::DNAStringSet]
#' ready for alignment or for [scan_orfs()].
#'
#' @param gene Character scalar.  Gene name; common synonyms are recognised
#'   case-insensitively (e.g. `"COI"`, `"CO1"`, `"COX1"` and `"MT-CO1"` all
#'   mean cytochrome c oxidase subunit I; `"16S"`, `"rrnL"` and `"MT-RNR2"`
#'   all mean the large-subunit rRNA).  Supported: COX1-3, CYTB, ND1-6,
#'   ND4L, ATP6, ATP8, 12S, 16S.  Other names are matched literally against
#'   the `/gene` qualifier.
#' @param species Character vector of species (or other taxon) names, e.g.
#'   `c("Sander vitreus", "Perca flavescens")`.  Each is searched on NCBI
#'   nuccore as `"<species>"[Organism]`.  `NULL` to skip.
#' @param accessions Character vector of NCBI accessions to extract the gene
#'   from directly, skipping the search.  `NULL` to skip.  At least one of
#'   `species` or `accessions` is required.
#' @param source Which records to search for when `species` is given:
#'   * `"mitogenome"` (default): complete mitochondrial genomes, preferring
#'     RefSeq (`NC_`) records and falling back to any complete mitogenome.
#'   * `"any"`: any mitochondrial record annotated with the gene, including
#'     single-gene barcode submissions.  Useful for building eDNA reference
#'     sets with several sequences per species (see `max_records`).
#' @param max_records Integer.  Maximum number of records per species.
#'   Default `1`.
#' @param term Optional extra Entrez query text ANDed onto every species
#'   search, e.g. `"2015:2025[PDAT]"`.
#' @param api_key Character scalar.  NCBI API key for higher rate limits.
#'   Defaults to the `ENTREZ_KEY` environment variable if set.
#' @param verbose Logical.  Print progress messages?  Default `TRUE`.
#'
#' @return A [Biostrings::DNAStringSet] named by accession, with the gene in
#'   its coding (5'->3') orientation, reverse-complemented when annotated on the
#'   minus strand.  Its `mcols()` hold `species`, `accession`, `gene`
#'   (as annotated in the record), `product`, `start`, `end`, `strand` and
#'   `definition`.  Species with no matching records, and records in which the
#'   gene is not annotated, are skipped with a warning.
#'
#' @details
#' The gene is located from the GenBank feature table (`CDS` features for
#' protein-coding genes, `rRNA` features for 12S/16S), matching the `/gene`
#' qualifier against known synonyms and falling back to the `/product`
#' qualifier.  Records without feature annotation (e.g. unannotated
#' assemblies) cannot be used.
#'
#' Species names are passed to NCBI as-is, so use current NCBI taxonomy (e.g.
#' *Barbodes semifasciolatus*, not *Puntius semifasciolatus*).  Records are
#' returned in NCBI relevance order; check `mcols(x)$species` for subspecies
#' or misidentified submissions.
#'
#' @seealso [fetch_mito_genbank()] for whole records with all features,
#'   [fetch_sequences()] for sequences by accession.
#'
#' @examples
#' \dontrun{
#' # COI from the RefSeq mitogenome of three percids
#' coi <- fetch_gene("COI", species = c("Sander vitreus", "Sander canadensis",
#'                                      "Perca flavescens"))
#' coi
#' S4Vectors::mcols(coi)
#'
#' # Up to 20 COI records per species, including barcode submissions
#' refs <- fetch_gene("COI", species = "Sander vitreus",
#'                    source = "any", max_records = 20)
#'
#' # From accessions you already have
#' cytb <- fetch_gene("CYTB", accessions = c("NC_028285", "NC_012920.1"))
#' }
#'
#' @importFrom rentrez entrez_fetch entrez_search
#' @export
fetch_gene <- function(gene,
                       species     = NULL,
                       accessions  = NULL,
                       source      = c("mitogenome", "any"),
                       max_records = 1L,
                       term        = NULL,
                       api_key     = NULL,
                       verbose     = TRUE) {

  if (!is.character(gene) || length(gene) != 1L || !nzchar(trimws(gene)))
    stop("'gene' must be a single gene name, e.g. \"COI\".", call. = FALSE)
  if (is.null(species) && is.null(accessions))
    stop("Supply 'species' and/or 'accessions'.", call. = FALSE)
  source <- match.arg(source)
  max_records <- as.integer(max_records)
  if (is.na(max_records) || max_records < 1L)
    stop("'max_records' must be a positive integer.", call. = FALSE)

  g <- .resolve_mito_gene(gene)

  ## ---- 1. Collect record IDs -------------------------------------------------
  ids <- character()
  if (!is.null(accessions)) {
    accs <- .parse_accession_arg(accessions)
    ids <- c(ids, accs)
  }

  if (!is.null(species)) {
    species <- unique(trimws(as.character(species)))
    species <- species[nzchar(species)]
    for (k in seq_along(species)) {
      sp <- species[[k]]
      if (k > 1L) .entrez_pause(api_key)
      if (verbose) message("Searching NCBI for ", g$label, " in ", sp, " ...")
      hits <- .search_gene_records(sp, g, source, max_records, term, api_key)
      if (!length(hits)) {
        warning(sp, ": no ", if (source == "mitogenome") "complete mitogenome"
                else paste(g$label, "records"),
                " found on NCBI; skipped.", call. = FALSE)
        next
      }
      ids <- c(ids, hits)
    }
  }
  ids <- unique(ids)
  if (!length(ids))
    stop("No NCBI records found for the requested species.", call. = FALSE)

  ## ---- 2. Fetch GenBank records in batches ------------------------------------
  if (verbose) message("Fetching ", length(ids), " GenBank record(s) ...")
  fetch_gb <- function(batch) {
    tryCatch(
      .entrez_retry(.efetch, db = "nuccore", id = batch,
                    rettype = "gbwithparts", retmode = "text",
                    api_key = api_key),
      error = function(e) {
        # A rejected batch falls back to one-by-one so good IDs survive.
        if (length(batch) > 1L)
          return(paste(vapply(batch, fetch_gb, character(1)), collapse = "\n"))
        warning(batch, ": could not be retrieved from NCBI (",
                .ncbi_error_reason(e), "); skipped.", call. = FALSE)
        ""
      })
  }
  records <- list()
  batches <- split(ids, ceiling(seq_along(ids) / 20L))
  for (b in seq_along(batches)) {
    if (b > 1L) .entrez_pause(api_key)
    records <- c(records, .split_genbank_records(fetch_gb(batches[[b]])))
  }

  ## ---- 3. Extract the gene from each record -----------------------------------
  seqs <- list()
  info <- list()
  for (rec_text in records) {
    ver <- .genbank_version(rec_text)
    rec <- tryCatch(.parse_genbank(rec_text, ver), error = function(e) NULL)
    if (is.null(rec)) next
    hit <- .find_gene_feature(rec$features, g)
    if (is.null(hit)) {
      warning(ver, " (", rec$organism, "): ", g$label,
              " not annotated in this record; skipped.", call. = FALSE)
      next
    }
    # The location parser collapses origin-spanning joins to min..max, which
    # would return most of the genome; no mitochondrial gene is > 5 kb.
    if (hit$end - hit$start + 1L > 5000L) {
      warning(ver, ": ", g$label, " spans the genome origin; skipped.",
              call. = FALSE)
      next
    }
    s <- Biostrings::subseq(rec$sequence, hit$start, hit$end)
    if (hit$strand == "-") s <- Biostrings::reverseComplement(s)
    seqs[[ver]] <- s
    info[[ver]] <- data.frame(species    = rec$organism,
                              accession  = ver,
                              gene       = hit$gene,
                              product    = hit$product,
                              start      = hit$start,
                              end        = hit$end,
                              strand     = hit$strand,
                              definition = rec$definition,
                              stringsAsFactors = FALSE)
  }

  # Accessions requested directly but not returned at all
  if (!is.null(accessions)) {
    ret  <- vapply(records, .genbank_version, character(1))
    miss <- accs[!vapply(accs, function(a)
      any(ret == a | sub("\\.\\d+$", "", ret) == a), logical(1))]
    if (length(miss))
      warning("Not returned by NCBI: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  if (!length(seqs))
    stop("Could not extract ", g$label, " from any record.", call. = FALSE)

  out <- Biostrings::DNAStringSet(seqs)
  S4Vectors::mcols(out) <- S4Vectors::DataFrame(do.call(rbind, unname(info)))
  if (verbose) message("Extracted ", g$label, " from ", length(out), " record(s).")
  out
}


# --------------------------------------------------------------------------- #
# Gene synonym table                                                          #
# --------------------------------------------------------------------------- #
.mito_gene_table <- function() {
  pcg <- function(label, syn, product)
    list(label = label, types = "CDS", syn = syn, product = product)
  rrna <- function(label, syn, product, title)
    list(label = label, types = "rRNA", syn = syn, product = product,
         title = title)
  list(
    COX1 = pcg("COX1", c("COX1", "COXI", "COI", "CO1", "MT-CO1", "COX-1", "CO-I"),
               "cytochrome c oxidase subunit (I|1)\\b"),
    COX2 = pcg("COX2", c("COX2", "COXII", "COII", "CO2", "MT-CO2", "COX-2", "CO-II"),
               "cytochrome c oxidase subunit (II|2)\\b"),
    COX3 = pcg("COX3", c("COX3", "COXIII", "COIII", "CO3", "MT-CO3", "COX-3", "CO-III"),
               "cytochrome c oxidase subunit (III|3)\\b"),
    CYTB = pcg("CYTB", c("CYTB", "CYT B", "CYT-B", "CYB", "COB", "MT-CYB"),
               "cytochrome b\\b"),
    ND1  = pcg("ND1", c("ND1", "NAD1", "NADH1", "MT-ND1", "NADH-1"),
               "NADH dehydrogenase subunit 1\\b"),
    ND2  = pcg("ND2", c("ND2", "NAD2", "NADH2", "MT-ND2", "NADH-2"),
               "NADH dehydrogenase subunit 2\\b"),
    ND3  = pcg("ND3", c("ND3", "NAD3", "NADH3", "MT-ND3", "NADH-3"),
               "NADH dehydrogenase subunit 3\\b"),
    ND4  = pcg("ND4", c("ND4", "NAD4", "NADH4", "MT-ND4", "NADH-4"),
               "NADH dehydrogenase subunit 4(?!L)\\b"),
    ND4L = pcg("ND4L", c("ND4L", "NAD4L", "NADH4L", "MT-ND4L", "NADH-4L"),
               "NADH dehydrogenase subunit 4L\\b"),
    ND5  = pcg("ND5", c("ND5", "NAD5", "NADH5", "MT-ND5", "NADH-5"),
               "NADH dehydrogenase subunit 5\\b"),
    ND6  = pcg("ND6", c("ND6", "NAD6", "NADH6", "MT-ND6", "NADH-6"),
               "NADH dehydrogenase subunit 6\\b"),
    ATP6 = pcg("ATP6", c("ATP6", "ATPASE6", "ATPASE 6", "MT-ATP6", "ATP-6"),
               "ATP synthase (F0 )?subunit (6|a)\\b"),
    ATP8 = pcg("ATP8", c("ATP8", "ATPASE8", "ATPASE 8", "MT-ATP8", "ATP-8"),
               "ATP synthase (F0 )?subunit 8\\b"),
    RNR1 = rrna("12S rRNA", c("12S", "12S RRNA", "RRNS", "RRN12", "MT-RNR1", "S-RRNA", "RNR1"),
                "12S ribosomal RNA|small subunit ribosomal RNA|\\bs-rRNA",
                c("12S ribosomal RNA", "12S rRNA")),
    RNR2 = rrna("16S rRNA", c("16S", "16S RRNA", "RRNL", "RRN16", "MT-RNR2", "L-RRNA", "RNR2"),
                "16S ribosomal RNA|large subunit ribosomal RNA|\\bl-rRNA",
                c("16S ribosomal RNA", "16S rRNA"))
  )
}

.resolve_mito_gene <- function(gene) {
  key <- toupper(trimws(gene))
  for (g in .mito_gene_table())
    if (key %in% g$syn) return(g)
  # Unknown gene: match it literally
  list(label = trimws(gene), types = c("CDS", "rRNA", "tRNA", "gene"),
       syn = key, product = NULL)
}


# --------------------------------------------------------------------------- #
# Search NCBI for records of one species; returns nuccore UIDs.               #
# --------------------------------------------------------------------------- #
.search_gene_records <- function(sp, g, source, max_records, term, api_key) {
  org   <- sprintf("\"%s\"[Organism]", gsub("\"", "", sp))
  extra <- if (!is.null(term) && nzchar(term)) paste0(" AND (", term, ")") else ""

  run <- function(q) {
    res <- .entrez_retry(.esearch, db = "nuccore", term = paste0(q, extra),
                         retmax = max_records, api_key = api_key)
    as.character(res$ids)
  }

  if (source == "mitogenome") {
    base <- paste(org, "AND mitochondrion[Filter] AND \"complete genome\"[Title]")
    ids  <- run(paste(base, "AND refseq[Filter]"))
    if (!length(ids)) { .entrez_pause(api_key); ids <- run(base) }
    return(ids)
  }

  gene_q <- paste0(sprintf("\"%s\"[Gene]", g$syn), collapse = " OR ")
  if (!is.null(g$title))
    gene_q <- paste(gene_q, "OR",
                    paste0(sprintf("\"%s\"[Title]", g$title), collapse = " OR "))
  run(sprintf("%s AND mitochondrion[Filter] AND (%s)", org, gene_q))
}


# --------------------------------------------------------------------------- #
# GenBank helpers                                                             #
# --------------------------------------------------------------------------- #
.split_genbank_records <- function(raw) {
  if (!nzchar(raw)) return(list())
  parts <- strsplit(raw, "\n//[ \t]*(\n|$)")[[1L]]
  parts <- parts[grepl("(^|\n)LOCUS", parts)]
  lapply(parts, function(p) paste0(sub("^[\\s\\S]*?(?=LOCUS)", "", p, perl = TRUE),
                                   "\n//\n"))
}

.genbank_version <- function(rec_text) {
  m <- regmatches(rec_text, regexpr("(?m)^VERSION\\s+\\S+", rec_text, perl = TRUE))
  if (length(m)) return(sub("^VERSION\\s+", "", m))
  m <- regmatches(rec_text, regexpr("(?m)^LOCUS\\s+\\S+", rec_text, perl = TRUE))
  if (length(m)) sub("^LOCUS\\s+", "", m) else NA_character_
}

.find_gene_feature <- function(features, g) {
  if (!nrow(features)) return(NULL)
  f <- features[features$type %in% g$types, , drop = FALSE]
  if (!nrow(f)) return(NULL)
  by_gene <- !is.na(f$gene) & toupper(trimws(f$gene)) %in% g$syn
  hit <- which(by_gene)
  if (!length(hit) && !is.null(g$product))
    hit <- which(!is.na(f$product) &
                 grepl(g$product, f$product, ignore.case = TRUE, perl = TRUE))
  if (!length(hit)) return(NULL)
  f[hit[[1L]], , drop = FALSE]
}
