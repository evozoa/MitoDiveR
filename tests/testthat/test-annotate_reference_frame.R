core <- "ATGCGTACGGATCCAAGCTTGGGCCCAAATTTGGGCCCTTTAAACCCGGG"  # 51 nt, non-repetitive

test_that("frame label agrees across genomes despite a local coordinate shift", {
  genomes <- Biostrings::DNAStringSet(c(
    anchor = core,
    query  = paste0("AA", core)   # 2-nt insertion before the shared block
  ))
  orfs <- data.frame(
    seq_id = c("anchor", "query"), strand = "+",
    start  = c(1L, 3L), end = c(9L, 11L),
    stringsAsFactors = FALSE
  )
  out <- annotate_reference_frame(orfs, genomes, anchor = "anchor", verbose = FALSE)

  # local starts differ (1 vs 3) but both map to the same anchor position
  expect_equal(out$ref_5p, c(1L, 1L))
  expect_equal(out$ref_start, c(1L, 1L))
  expect_equal(out$ref_end, c(9L, 9L))
  expect_equal(out$ref_frame, c("+_0", "+_0"))
})

test_that("minus-strand ORFs use the end coordinate as the 5' (start-codon) position", {
  genomes <- Biostrings::DNAStringSet(c(anchor = core))
  orfs <- data.frame(
    seq_id = "anchor", strand = "-", start = 5L, end = 13L,
    stringsAsFactors = FALSE
  )
  out <- annotate_reference_frame(orfs, genomes, anchor = "anchor", verbose = FALSE)
  expect_equal(out$ref_5p, 13L)
  expect_equal(out$ref_start, 5L)
  expect_equal(out$ref_end, 13L)
  expect_equal(out$ref_frame, paste0("-_", (13L - 1L) %% 3L))
})

test_that("positions aligning before the anchor's first base are NA", {
  genomes <- Biostrings::DNAStringSet(c(
    anchor = core,
    query  = paste0("GGGGG", core)  # unique 5-nt insertion, no counterpart in anchor
  ))
  orfs <- data.frame(
    seq_id = "query", strand = "+", start = 1L, end = 5L,  # entirely within the insertion
    stringsAsFactors = FALSE
  )
  out <- annotate_reference_frame(orfs, genomes, anchor = "anchor", verbose = FALSE)
  expect_true(is.na(out$ref_5p))
  expect_true(is.na(out$ref_frame))
})

test_that("defaults to the first genome as anchor and messages about it", {
  genomes <- Biostrings::DNAStringSet(c(anchor = core, query = core))
  orfs <- data.frame(
    seq_id = c("anchor", "query"), strand = "+",
    start = c(1L, 1L), end = c(9L, 9L),
    stringsAsFactors = FALSE
  )
  expect_message(
    out <- annotate_reference_frame(orfs, genomes),
    "anchor"
  )
  expect_equal(out$ref_5p, c(1L, 1L))
})

test_that("errors on missing required columns", {
  bad <- data.frame(start = 1L, end = 10L)
  expect_error(annotate_reference_frame(bad, Biostrings::DNAStringSet(c(a = core))),
               "seq_id")
})

test_that("errors when genomes is missing a sequence referenced in orfs", {
  genomes <- Biostrings::DNAStringSet(c(anchor = core))
  orfs <- data.frame(seq_id = "missing_genome", strand = "+", start = 1L, end = 9L,
                      stringsAsFactors = FALSE)
  expect_error(annotate_reference_frame(orfs, genomes, anchor = "anchor"),
               "missing_genome")
})

test_that("row order and existing columns are preserved", {
  genomes <- Biostrings::DNAStringSet(c(anchor = core, query = core))
  orfs <- data.frame(
    seq_id = c("query", "anchor"), strand = "+",
    start = c(20L, 1L), end = c(28L, 9L),
    length_nt = c(9L, 9L),
    stringsAsFactors = FALSE
  )
  out <- annotate_reference_frame(orfs, genomes, anchor = "anchor", verbose = FALSE)
  expect_equal(out$seq_id, c("query", "anchor"))
  expect_equal(out$length_nt, c(9L, 9L))
})
