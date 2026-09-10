ref <- c(Humanin = "MAPRGFSCLLLLTSEIDLPVKRRA",
         SHLP6   = "MLDQDIPMVQPLLKVRLFND")

mk <- function(prot, ...) {
  data.frame(start = seq_along(prot), end = seq_along(prot) + 10L,
             protein_sequence = prot, stringsAsFactors = FALSE, ...)
}

test_that("reference is required, with guidance naming the built-in trap", {
  o <- mk("MAPRGFSCLLLLTSEIDLPVKRRA")
  expect_error(annotate_mdp_orfs(o), "must be specified")
  expect_error(annotate_mdp_orfs(o), "mdp_sequences\\(\\) holds the human set")
})

test_that("input validation rejects malformed arguments", {
  expect_error(annotate_mdp_orfs(42, reference = ref), "data\\.frame")
  expect_error(annotate_mdp_orfs(data.frame(a = 1), reference = ref),
               "protein_sequence")
  expect_error(annotate_mdp_orfs(mk("MAP"), reference = c("MAPRG")),
               "must be named")
  expect_error(annotate_mdp_orfs(mk("MAP"), reference = c(A = "MA", A = "MB")),
               "duplicate names")
  expect_error(
    annotate_mdp_orfs(mk("MAP"), reference = ref, method = "homology",
                      min_identity = 140),
    "between 0 and 100"
  )
})

test_that("exact matching identifies peptides and leaves others NA", {
  o <- mk(c("MAPRGFSCLLLLTSEIDLPVKRRA",     # Humanin, identical
            "MLDQDIPMVQPLLKVRLFND",         # SHLP6, identical
            "WWWWWWWWWWWWWWWW"))            # nothing
  a <- suppressMessages(annotate_mdp_orfs(o, reference = ref))

  expect_equal(a$mdp_name, c("Humanin", "SHLP6", NA))
  expect_equal(a$mdp_match, c("identical", "identical", NA))
  expect_equal(a$mdp_pct_identity, c(100, 100, NA))
  expect_equal(a$mdp_coverage, c(100, 100, NA))
})

test_that("a trailing stop codon does not defeat matching", {
  a <- suppressMessages(
    annotate_mdp_orfs(mk("MAPRGFSCLLLLTSEIDLPVKRRA*"), reference = ref))
  expect_equal(a$mdp_name, "Humanin")
  expect_equal(a$mdp_match, "identical")
})

test_that("containment is reported in the direction it occurred", {
  # ORF runs past the annotated peptide
  a <- suppressMessages(
    annotate_mdp_orfs(mk("QQMAPRGFSCLLLLTSEIDLPVKRRAQQ"), reference = ref))
  expect_equal(a$mdp_name, "Humanin")
  expect_equal(a$mdp_match, "orf_contains_mdp")
  expect_equal(a$mdp_coverage, 100)

  # ORF is a fragment of the peptide
  b <- suppressMessages(
    annotate_mdp_orfs(mk("MAPRGFSCLLLL"), reference = ref))
  expect_equal(b$mdp_name, "Humanin")
  expect_equal(b$mdp_match, "mdp_contains_orf")
  expect_lt(b$mdp_coverage, 100)
})

test_that("the longest reference wins when several match", {
  r2 <- c(Short = "MAPRG", Long = "MAPRGFSCLLLL")
  a  <- suppressMessages(annotate_mdp_orfs(mk("MAPRGFSCLLLL"), reference = r2))
  expect_equal(a$mdp_name, "Long")
})

test_that("existing columns and row order are preserved", {
  o <- mk(c("MAPRGFSCLLLLTSEIDLPVKRRA", "WWWW"))
  o$strand <- c("+", "-")
  o$note   <- c("a", "b")
  a <- suppressMessages(annotate_mdp_orfs(o, reference = ref))

  expect_equal(a[, names(o)], o)
  expect_true(all(c("mdp_name", "mdp_match", "mdp_pct_identity",
                    "mdp_coverage") %in% names(a)))
})

test_that("prefix renames the added columns", {
  a <- suppressMessages(
    annotate_mdp_orfs(mk("MAPRGFSCLLLLTSEIDLPVKRRA"), reference = ref,
                      prefix = "pep"))
  expect_true("pep_name" %in% names(a))
  expect_false("mdp_name" %in% names(a))
  expect_equal(a$pep_name, "Humanin")
})

test_that("an empty ORF table still gains the columns", {
  o <- mk(character(0))
  a <- annotate_mdp_orfs(o, reference = ref)
  expect_equal(nrow(a), 0L)
  expect_true("mdp_name" %in% names(a))
})

test_that("a FASTA file is accepted as the reference", {
  f <- tempfile(fileext = ".fasta")
  on.exit(unlink(f), add = TRUE)
  writeLines(c(">Humanin desc here", "MAPRGFSCLLLLTSEIDLPVKRRA"), f)

  a <- suppressMessages(
    annotate_mdp_orfs(mk("MAPRGFSCLLLLTSEIDLPVKRRA"), reference = f))
  expect_equal(a$mdp_name, "Humanin")   # header truncated at first space
})

test_that("homology matching applies its thresholds", {
  # one substitution: still clearly the same peptide
  near <- "MAPRGFSCLLLLTSEIDLPVKRRW"
  hit  <- suppressMessages(
    annotate_mdp_orfs(mk(near), reference = ref, method = "homology",
                      min_identity = 60, min_coverage = 60))
  expect_equal(hit$mdp_name, "Humanin")
  expect_equal(hit$mdp_match, "homology")
  expect_gt(hit$mdp_pct_identity, 60)

  # the same ORF fails an identity threshold it cannot meet
  miss <- suppressMessages(
    annotate_mdp_orfs(mk(near), reference = ref, method = "homology",
                      min_identity = 99.5, min_coverage = 60))
  expect_true(is.na(miss$mdp_name))
})

test_that("homology does not fall back to the built-in MDP set", {
  # SHLP2's real sequence is absent from `ref`, so nothing may be assigned
  a <- suppressMessages(
    annotate_mdp_orfs(mk("MGVKFFTLSTRFFPSVQRAVPLWTNS"), reference = ref,
                      method = "homology", min_identity = 60,
                      min_coverage = 60))
  expect_true(is.na(a$mdp_name))
})
