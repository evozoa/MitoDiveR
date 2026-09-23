# Offline tests mock the NCBI layer (.efetch / .esearch); live tests at the
# bottom are skipped on CRAN and when offline.

gb_record <- function(acc, organism = "Testus testus") {
  paste0(
"LOCUS       ", sub("\\..*", "", acc), "                 60 bp    DNA     circular VRT 01-JAN-2020
DEFINITION  ", organism, " mitochondrion, complete genome.
ACCESSION   ", sub("\\..*", "", acc), "
VERSION     ", acc, "
SOURCE      mitochondrion ", organism, "
  ORGANISM  ", organism, "
            Eukaryota; Chordata.
FEATURES             Location/Qualifiers
     source          1..60
                     /organism=\"", organism, "\"
     gene            1..12
                     /gene=\"COI\"
     CDS             1..12
                     /gene=\"COI\"
                     /product=\"cytochrome c oxidase subunit I\"
     CDS             complement(21..29)
                     /product=\"NADH dehydrogenase subunit 6\"
ORIGIN
        1 atgaaacccg ggtaaaaaaa ttttccccgg ggaaaaccct tgggtaaatt tcccgggaat
//
")
}

fasta_record <- function(acc) paste0(">", acc, " Testus testus mitochondrion\nACGTACGTAC\n")

# Mock efetch: serves the IDs in `known`, silently drops the rest (as NCBI
# does in a mixed batch) and errors with HTTP 400 when none are known.
mock_efetch <- function(known) {
  function(db, id, rettype, retmode, ...) {
    base <- sub("\\.\\d+$", "", id)
    ok   <- base %in% sub("\\.\\d+$", "", known)
    if (!any(ok)) stop("HTTP failure: 400\nError+F+a+i+l+e+d")
    full <- known[match(base[ok], sub("\\.\\d+$", "", known))]
    fmt  <- if (rettype == "fasta") fasta_record else gb_record
    paste(vapply(full, fmt, character(1)), collapse = "")
  }
}


# ---- input cleaning ---------------------------------------------------------

test_that(".split_accessions trims, upper-cases and splits", {
  expect_equal(.split_accessions(c(" nc_028285, NC_012920.1 ", "AB123456")),
               c("NC_028285", "NC_012920.1", "AB123456"))
})

test_that("species names are rejected with a pointer to fetch_gene()", {
  expect_error(fetch_sequences("Homo sapiens"), "fetch_gene")
  expect_error(fetch_mito_genbank("Homo sapiens"), "fetch_gene")
})

test_that("a missing file path mentions the working directory", {
  expect_error(fetch_mito_genbank("data/walleye.fasta"), "working directory")
})

test_that("an accession passed as fasta_file gets a c() hint", {
  local_mocked_bindings(.efetch = mock_efetch("NC_028285.1"))
  expect_error(suppressMessages(fetch_sequences("NC_028285", "NC_012920.1")),
               'c\\("NC_028285", "NC_012920.1"\\)', fixed = FALSE)
})


# ---- retry ------------------------------------------------------------------

test_that(".entrez_retry retries transient errors, not permanent ones", {
  n <- 0L
  flaky <- function() { n <<- n + 1L; if (n < 3L) stop("HTTP failure: 429") else "ok" }
  expect_equal(.entrez_retry(flaky, .wait = 0), "ok")
  expect_equal(n, 3L)

  n <- 0L
  bad <- function() { n <<- n + 1L; stop("HTTP failure: 400") }
  expect_error(.entrez_retry(bad, .wait = 0), "400")
  expect_equal(n, 1L)
})

test_that(".ncbi_error_reason decodes HTTP status", {
  expect_match(.ncbi_error_reason(simpleError("HTTP failure: 400\n+F+a+i+l")),
               "not found or invalid")
  expect_match(.ncbi_error_reason(simpleError("HTTP failure: 429")), "rate-limited")
})


# ---- fetch_sequences ----------------------------------------------------------

test_that("fetch_sequences names results as requested and warns on missing", {
  local_mocked_bindings(.efetch = mock_efetch(c("NC_028285.1", "NC_012920.1")))
  expect_warning(
    s <- suppressMessages(fetch_sequences(c("nc_028285", "NC_999999999",
                                            "NC_012920.1"))),
    "NC_999999999")
  expect_equal(names(s), c("NC_028285", "NC_012920.1"))
  expect_s4_class(s[["NC_028285"]], "DNAString")
})

test_that("fetch_sequences splits a comma-separated string", {
  local_mocked_bindings(.efetch = mock_efetch(c("NC_028285.1", "NC_012920.1")))
  s <- suppressMessages(fetch_sequences("NC_028285, NC_012920.1"))
  expect_length(s, 2L)
})

test_that("fetch_sequences errors clearly when nothing comes back", {
  local_mocked_bindings(.efetch = mock_efetch(character()))
  expect_error(suppressMessages(fetch_sequences("NC_999999999")),
               "NCBI returned no sequences")
})


# ---- fetch_mito_genbank -----------------------------------------------------

test_that("fetch_mito_genbank skips a bad accession but keeps the rest", {
  local_mocked_bindings(.efetch = mock_efetch(c("NC_028285.1", "NC_012920.1")),
                        .entrez_pause = function(...) NULL)
  expect_warning(
    g <- fetch_mito_genbank(c("NC_028285.1", "NC_999999999", " NC_012920.1"),
                            verbose = FALSE),
    "NC_999999999: could not be retrieved.*not found or invalid")
  expect_equal(names(g), c("NC_028285.1", "NC_012920.1"))
  expect_equal(g[["NC_028285.1"]]$organism, "Testus testus")
})

test_that("fetch_mito_genbank errors when every accession fails", {
  local_mocked_bindings(.efetch = mock_efetch(character()))
  expect_error(suppressWarnings(fetch_mito_genbank("NC_999999999", verbose = FALSE)),
               "None of the requested accessions")
})

test_that("fetch_mito_genbank still accepts raw sequence", {
  g <- fetch_mito_genbank("ACGTACGTNNACGT", verbose = FALSE)
  expect_equal(g$seq_1$length, 14L)
})


# ---- fetch_gene ---------------------------------------------------------------

test_that("gene synonyms resolve", {
  expect_equal(.resolve_mito_gene("coi")$label, "COX1")
  expect_equal(.resolve_mito_gene("MT-CYB")$label, "CYTB")
  expect_equal(.resolve_mito_gene("rrnL")$label, "16S rRNA")
  expect_equal(.resolve_mito_gene("FOO1")$label, "FOO1")
})

test_that("fetch_gene extracts by /gene and by /product, on both strands", {
  local_mocked_bindings(.efetch = mock_efetch("NC_000001.1"),
                        .entrez_pause = function(...) NULL)
  coi <- fetch_gene("CO1", accessions = "NC_000001", verbose = FALSE)
  expect_equal(as.character(coi[[1]]), "ATGAAACCCGGG")
  expect_equal(S4Vectors::mcols(coi)$species, "Testus testus")

  # ND6 has no /gene qualifier (product match) and is on the minus strand
  nd6 <- fetch_gene("ND6", accessions = "NC_000001", verbose = FALSE)
  expect_equal(as.character(nd6[[1]]),
               as.character(Biostrings::reverseComplement(
                 Biostrings::DNAString("TTTTCCCCG"))))
})

test_that("fetch_gene searches by species and warns on no hits", {
  searched <- character()
  local_mocked_bindings(
    .esearch = function(db, term, ...) {
      searched <<- c(searched, term)
      list(ids = if (grepl("Testus testus", term)) "NC_000001.1" else character())
    },
    .efetch = mock_efetch("NC_000001.1"),
    .entrez_pause = function(...) NULL)
  expect_warning(
    x <- fetch_gene("COI", species = c("Testus testus", "Nullus nullus"),
                    verbose = FALSE),
    "Nullus nullus: no complete mitogenome")
  expect_equal(names(x), "NC_000001.1")
  expect_true(any(grepl("refseq\\[Filter\\]", searched)))
})

test_that("fetch_gene warns when the gene is not annotated", {
  local_mocked_bindings(.efetch = mock_efetch("NC_000001.1"))
  expect_error(expect_warning(
    fetch_gene("ATP8", accessions = "NC_000001", verbose = FALSE),
    "not annotated"), "Could not extract")
})


# ---- live NCBI ----------------------------------------------------------------

test_that("live: fetch_gene gets walleye COI from its mitogenome", {
  skip_on_cran()
  skip_if_offline()
  coi <- fetch_gene("COI", species = "Sander vitreus", verbose = FALSE)
  expect_length(coi, 1L)
  expect_gt(Biostrings::width(coi), 1500L)
})

test_that("live: bad accessions are skipped, not fatal", {
  skip_on_cran()
  skip_if_offline()
  expect_warning(
    s <- suppressMessages(fetch_sequences(c("NC_028285", "NC_999999999"))),
    "NC_999999999")
  expect_equal(names(s), "NC_028285")
})

test_that("a gene spanning a whole short (barcode) record is not flagged", {
  barcode <- sub("CDS             1..12", "CDS             1..60",
                 gb_record("MN000001.1"), fixed = TRUE)
  local_mocked_bindings(.efetch = function(...) barcode)
  x <- fetch_gene("COI", accessions = "MN000001", verbose = FALSE)
  expect_equal(Biostrings::width(x), 60L)
})
