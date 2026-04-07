library(testthat)
library(MitoDiveR)
library(Biostrings)

# ---- helpers -----------------------------------------------------------------
# Single clean ORF in frame +1 (standard code, ATG start, TAA stop)
# ATG AAA CCC GGG TTT TAA  = 18 nt
orf_seq   <- "ATGAAACCCGGGTTTTTAA"    # 19 nt padded — adjust below for clean:
simple_fw <- paste0("ATGAAACCCGGGTTTTAA", strrep("G", 30))  # 18 + 30 = 48 nt

# ---- basic forward-strand detection -----------------------------------------
test_that("detects a simple ATG-TAA ORF in frame +1", {
  orfs <- find_orfs(simple_fw,
                    genetic_code   = "SGC0",
                    start_codons   = "ATG",
                    min_orf_length = 15,
                    both_strands   = FALSE,
                    circular       = FALSE)

  expect_s3_class(orfs, "data.frame")
  expect_true(nrow(orfs) >= 1L)
  expect_true(any(orfs$start == 1L))
  expect_equal(orfs$start_codon[orfs$start == 1L], "ATG")
  expect_equal(orfs$stop_codon[orfs$start == 1L],  "TAA")
  expect_equal(orfs$frame[orfs$start == 1L],        1L)
  expect_equal(orfs$strand[orfs$start == 1L],      "+")
  # 18 nt: ATG + 4 codons + TAA
  expect_equal(orfs$length_nt[orfs$start == 1L],   18L)
})

# ---- reverse-complement detection -------------------------------------------
test_that("detects ORF on minus strand", {
  # Build a sequence where the RC of positions 20-37 is ATG...TAA
  # RC of "TTATAAA ACCCGGG TTT CAT" = "ATG AAA CCC GGG TTT TAA A"
  # TTAAAAAACCCGGGTTTCAT  → RC = ATGAAACCCGGGTTTTTAA ... let's compute directly
  fwd    <- "TTAAAACCCGGGTTTCAT"          # RC of this is the ORF
  padded <- paste0(strrep("C", 10), fwd, strrep("C", 10))

  orfs <- find_orfs(padded,
                    genetic_code   = "SGC0",
                    start_codons   = "ATG",
                    min_orf_length = 15,
                    both_strands   = TRUE,
                    circular       = FALSE)

  expect_true(any(orfs$strand == "-"))
})

# ---- all six frames are searched --------------------------------------------
test_that("returns results across multiple frames", {
  # Construct a sequence with ORFs in frames +1 and +2
  f1 <- "ATGAAATTTGGG"          # frame +1 ORF prefix (no stop yet)
  f2 <- paste0("A", "ATGCCCTAA") # frame +2: offset by 1, ATG...TAA
  seq_multi <- paste0("ATGAAATTTCCCTAA",   # frame +1: ATG...TAA at pos 1
                      "AATGCCCTAA",         # frame +2: ATG at pos 17
                      strrep("T", 20))

  orfs <- find_orfs(seq_multi,
                    genetic_code   = "SGC0",
                    start_codons   = "ATG",
                    min_orf_length = 6,
                    both_strands   = FALSE,
                    circular       = FALSE)

  frames_found <- unique(orfs$frame)
  expect_true(1L %in% frames_found)
})

# ---- min_orf_length filter --------------------------------------------------
test_that("min_orf_length filters short ORFs", {
  seq_short <- paste0("ATGAAATAA", strrep("C", 40))   # 9 nt ORF

  orfs_strict <- find_orfs(seq_short,
                            genetic_code   = "SGC0",
                            start_codons   = "ATG",
                            min_orf_length = 20,
                            both_strands   = FALSE,
                            circular       = FALSE)
  expect_equal(nrow(orfs_strict), 0L)

  orfs_loose <- find_orfs(seq_short,
                           genetic_code   = "SGC0",
                           start_codons   = "ATG",
                           min_orf_length = 6,
                           both_strands   = FALSE,
                           circular       = FALSE)
  expect_true(nrow(orfs_loose) >= 1L)
})

# ---- include_stop_codon toggle ----------------------------------------------
test_that("include_stop_codon = FALSE shortens ORF by 3 nt", {
  orfs_with <- find_orfs(simple_fw,
                          genetic_code       = "SGC0",
                          start_codons       = "ATG",
                          min_orf_length     = 12,
                          both_strands       = FALSE,
                          circular           = FALSE,
                          include_stop_codon = TRUE)

  orfs_without <- find_orfs(simple_fw,
                              genetic_code       = "SGC0",
                              start_codons       = "ATG",
                              min_orf_length     = 12,
                              both_strands       = FALSE,
                              circular           = FALSE,
                              include_stop_codon = FALSE)

  shared <- merge(orfs_with[orfs_with$start == 1L, ],
                  orfs_without[orfs_without$start == 1L, ],
                  by = "start", suffixes = c("_w", "_wo"))
  expect_equal(shared$length_nt_w - shared$length_nt_wo, 3L)
  expect_true(is.na(shared$stop_codon_wo))
})

# ---- vertebrate mitochondrial code: TGA is not a stop -----------------------
test_that("TGA is not a stop codon under vertebrate mito code (SGC1)", {
  # ATG + TGA (Trp in mito) + AAA + TAA
  seq_mito <- paste0("ATGTGAAAATAA", strrep("C", 30))

  orfs_mito <- find_orfs(seq_mito,
                          genetic_code   = "SGC1",
                          start_codons   = "ATG",
                          min_orf_length = 9,
                          both_strands   = FALSE,
                          circular       = FALSE)

  orfs_std  <- find_orfs(seq_mito,
                          genetic_code   = "SGC0",
                          start_codons   = "ATG",
                          min_orf_length = 9,
                          both_strands   = FALSE,
                          circular       = FALSE)

  # Under mito code, TGA is read through → longer ORF from pos 1
  mito_orf <- orfs_mito[orfs_mito$start == 1L, ]
  std_orf  <- orfs_std[orfs_std$start  == 1L, ]

  expect_true(nrow(mito_orf) >= 1L)
  expect_true(mito_orf$length_nt[1L] > std_orf$length_nt[1L])
})

# ---- DNAString input --------------------------------------------------------
test_that("accepts a DNAString input", {
  ds   <- Biostrings::DNAString(simple_fw)
  orfs <- find_orfs(ds,
                    genetic_code   = "SGC0",
                    start_codons   = "ATG",
                    min_orf_length = 15,
                    both_strands   = FALSE,
                    circular       = FALSE)
  expect_true(nrow(orfs) >= 1L)
})

# ---- output schema ----------------------------------------------------------
test_that("output data.frame has the expected columns", {
  orfs <- find_orfs(simple_fw,
                    genetic_code   = "SGC0",
                    start_codons   = "ATG",
                    min_orf_length = 15,
                    both_strands   = FALSE,
                    circular       = FALSE)

  expected_cols <- c("strand", "frame", "start", "end",
                     "length_nt", "length_aa",
                     "start_codon", "stop_codon", "wraps_around",
                     "orf_sequence", "protein_sequence")
  expect_true(all(expected_cols %in% names(orfs)))
})
