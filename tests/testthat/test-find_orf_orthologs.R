library(testthat)
library(MitoDiveR)
library(Biostrings)

# ---- helpers -----------------------------------------------------------------
# A peptide and its EXACT standard-code (SGC0) reverse translation + TAA stop,
# so scan_orfs(genetic_code = "SGC0") recovers it at 100% identity.
QUERY_PEP <- c(MyMDP = "MFADRWLQGSTKVPEAILNCYDEHKR")
ORF_DNA   <- paste0(
  "ATG","TTT","GCT","GAT","CGT","TGG","CTT","CAA","GGT","TCT","ACT","AAA","GTT",
  "CCT","GAA","GCT","ATT","CTT","AAT","TGT","TAT","GAT","GAA","CAT","AAA","CGT","TAA")
FLANK <- strrep("C", 24)   # ATG-free filler

with_orf    <- paste0(FLANK, ORF_DNA, FLANK)
without_orf <- strrep("C", 90)

# ---- local discovery: recover a planted ortholog -----------------------------
test_that("recovers a planted ortholog ORF and excludes genomes without it", {
  genomes <- DNAStringSet(c(taxonA = with_orf, taxonB = without_orf))

  hits <- find_orf_orthologs(QUERY_PEP, genomes = genomes, method = "local",
                             genetic_code = "SGC0", annotate = FALSE,
                             min_orf_length = 30, min_identity = 80,
                             min_coverage = 50)

  expect_equal(nrow(hits), 1L)
  expect_equal(hits$seq_id, "taxonA")
  expect_equal(hits$query_name, "MyMDP")
  expect_equal(hits$strand, "+")
  expect_equal(hits$pct_identity, 100)
  expect_equal(hits$query_coverage, 100)
  expect_equal(hits$method, "local")
  expect_true(is.na(hits$blast_accession))   # local => no BLAST provenance
  expect_match(hits$protein_sequence, "^MFADRWLQGSTKVPEAILNCYDEHKR")
})

test_that("output carries the documented columns", {
  genomes <- DNAStringSet(c(taxonA = with_orf))
  hits <- find_orf_orthologs(QUERY_PEP, genomes = genomes, method = "local",
                             genetic_code = "SGC0", annotate = FALSE,
                             min_orf_length = 30, min_identity = 80)
  expect_true(all(c("query_name", "seq_id", "start", "end", "strand", "frame",
                    "length_nt", "protein_sequence", "pct_identity",
                    "query_coverage", "orf_coverage", "genomic_region",
                    "method", "blast_accession") %in% names(hits)))
})

# ---- coverage gate is on the KNOWN peptide -----------------------------------
test_that("min_coverage gates on coverage of the query peptide", {
  # ORF encoding only the first 10 aa of the 26-aa peptide (~38% coverage).
  partial <- paste0(FLANK,
    "ATG","TTT","GCT","GAT","CGT","TGG","CTT","CAA","GGT","TCT","TAA", FLANK)
  g <- DNAStringSet(c(taxonC = partial))

  strict <- find_orf_orthologs(QUERY_PEP, genomes = g, method = "local",
                               genetic_code = "SGC0", annotate = FALSE,
                               min_orf_length = 30, min_identity = 80,
                               min_coverage = 50)
  expect_equal(nrow(strict), 0L)

  loose <- find_orf_orthologs(QUERY_PEP, genomes = g, method = "local",
                              genetic_code = "SGC0", annotate = FALSE,
                              min_orf_length = 30, min_identity = 80,
                              min_coverage = 30)
  expect_equal(nrow(loose), 1L)
  expect_lt(loose$query_coverage, 50)   # partial peptide coverage retained
})

test_that("returns an empty (zero-row) table when nothing matches", {
  g <- DNAStringSet(c(taxonB = without_orf))
  hits <- find_orf_orthologs(QUERY_PEP, genomes = g, method = "local",
                             genetic_code = "SGC0", annotate = FALSE,
                             min_orf_length = 30, min_identity = 80)
  expect_s3_class(hits, "data.frame")
  expect_equal(nrow(hits), 0L)
})

# ---- input validation --------------------------------------------------------
test_that("validates inputs and required arguments", {
  genomes <- DNAStringSet(c(taxonA = with_orf))

  # genetic_code is mandatory (no default, by package convention)
  expect_error(find_orf_orthologs(QUERY_PEP, genomes = genomes, method = "local"),
               "genetic_code")
  # local needs a sequence source
  expect_error(
    find_orf_orthologs(QUERY_PEP, method = "local", genetic_code = "SGC0"),
    "'genomes' or 'taxa'")
  # blastn needs a host region
  expect_error(
    find_orf_orthologs(QUERY_PEP, method = "blastn", genetic_code = "SGC0"),
    "host_region")
  # query must be a non-empty amino-acid string
  expect_error(
    find_orf_orthologs("", genomes = genomes, method = "local",
                       genetic_code = "SGC0"),
    "non-empty amino-acid")
  # bad method
  expect_error(
    find_orf_orthologs(QUERY_PEP, genomes = genomes, method = "nope",
                       genetic_code = "SGC0"))
})
