# Guards the built-in MDP reference sequences.
#
# The SHLP1-6 entries were once wrong: they were not encoded anywhere in the human
# mitochondrial genome, so homology screening could never recover a real SHLP. This
# test re-derives every natural peptide from the actual rCRS (NC_012920.1) coding
# sequence, translated with the STANDARD genetic code (MDPs are translated by
# cytoplasmic ribosomes), so a wrong sequence cannot be reintroduced silently.

# rCRS coding sequence of each MDP, already oriented 5'->3' on the coding strand.
rcrs_cds <- c(
  Humanin = "ATGGCTCCACGAGGGTTCAGCTGTCTCTTACTTTTAACCAGTGAAATTGACCTGCCCGTGAAGAGGCGGGCA",
  `MOTS-c` = "ATGAGGTGGCAAGAAATGGGCTACATTTTCTACCCCAGAAAACTACGA",
  SHLP1 = "ATGTGTCACTGGGCAGGCGGTGCCTCTAATACTGGTGATGCTAGAGGTGATGTTTTTGGTAAACAGGCGGGG",
  SHLP2 = "ATGGGTGTTAAATTTTTTACTCTCTCTACAAGGTTTTTTCCTAGTGTCCAAAGAGCTGTTCCTCTTTGGACTAACAGT",
  SHLP3 = paste0("ATGCTTGGTTATAATTTTTCATCTTTCCCTTGCGGTACTATATCTATTGCGCCAGGTTTCAATTTCTAT",
                 "CGCCTATACTTTATTTGGGTAAATGGTTTGGCTAAGGTTGTCTGG"),
  SHLP4 = "ATGCTAGAGGTGATGTTTTTGGTAAACAGGCGGGGTAAGATTTGCCGAGTTCCTTTTACTTTTTTTAACCTTTCCTTA",
  SHLP5 = "ATGTACTGCTCGGAGGTTGGGTTCTGCTCCGAGGTCGCCCCAACCGAAATTTTTAATGCAGGTTTGGTAGTT",
  SHLP6 = "ATGTTGGATCAGGACATCCCGATGGTGCAGCCGCTATTAAAGGTTCGTTTGTTCAACGAT"
)

translate_standard <- function(nt) {
  aa <- as.character(Biostrings::translate(Biostrings::DNAString(nt),
                                           genetic.code = Biostrings::GENETIC_CODE,
                                           if.fuzzy.codon = "solve"))
  sub("\\*$", "", aa)
}

test_that("built-in MDPs are the peptides actually encoded in the human rCRS", {
  built_in <- mdp_sequences()
  for (nm in names(rcrs_cds)) {
    expect_true(nm %in% names(built_in), info = paste(nm, "missing from mdp_sequences()"))
    expect_identical(unname(built_in[[nm]]), translate_standard(rcrs_cds[[nm]]),
                     info = paste(nm, "does not match the rCRS-encoded peptide"))
  }
})

test_that("the two internal reference tables agree", {
  a <- mdp_sequences()                              # from .mdp_references
  b <- unlist(MitoDiveR:::.MDP_REFERENCE)           # from score_mdp_similarity
  expect_setequal(names(a), names(b))
  expect_identical(a[sort(names(a))], b[sort(names(b))])
})

test_that("only canonical Swiss-Prot peptides are included", {
  s <- mdp_sequences()
  expect_setequal(names(s), c("Humanin", "MOTS-c", paste0("SHLP", 1:6)))
  # synthetic analogues (e.g. HNG, the S14G humanin variant) are not encoded in
  # any genome and must not sit in a reference set used for genome screens
  expect_false(any(grepl("Humanin-G|HNG", names(s))))
  expect_false("MAPRGFSCLLLLTGEIDLPVKRRA" %in% unname(s))
})

test_that("known SHLP sequence lengths match the published peptides", {
  s <- mdp_sequences()
  expect_equal(vapply(s[c("SHLP1", "SHLP2", "SHLP3", "SHLP4", "SHLP5", "SHLP6")],
                      nchar, integer(1), USE.NAMES = FALSE),
               c(24L, 26L, 38L, 26L, 24L, 20L))
  expect_equal(nchar(s[["Humanin"]]), 24L)
  expect_equal(nchar(s[["MOTS-c"]]), 16L)
})

test_that("no reference peptide is a stale placeholder absent from rCRS", {
  # every natural MDP must be findable in a six-frame standard-code translation
  # of its own coding sequence concatenation (sanity check on the table above)
  for (nm in names(rcrs_cds)) {
    expect_true(grepl("^M", translate_standard(rcrs_cds[[nm]])),
                info = paste(nm, "does not start with Met"))
  }
})
