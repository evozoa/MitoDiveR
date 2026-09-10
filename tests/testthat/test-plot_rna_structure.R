skip_if_no_rnaplot <- function() {
  testthat::skip_if(nchar(Sys.which("RNAplot")) == 0L,
                    "RNAplot (ViennaRNA) not available")
}

fold_stub <- list(
  #        1234567890123456789012345678901234567890
  sequence  = "AUGGCUCCACGAGGGUUCAGCUGUCUCUUACUUUUAACCAG",
  structure = "((((((....))))))...((((((......))))))...."
)

test_that("input validation rejects malformed arguments", {
  expect_error(plot_rna_structure(42), "rna_fold")
  expect_error(plot_rna_structure("((..))"), "sequence")
  expect_error(
    plot_rna_structure(list(sequence = "AUGC", structure = "((..))")),
    "same length"
  )
})

test_that("a structure with no ORFs draws and returns layout coordinates", {
  skip_if_no_rnaplot()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  res <- plot_rna_structure(fold_stub)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), nchar(fold_stub$structure))
  expect_named(res, c("position", "x", "y"))
  expect_true(all(is.finite(res$x)) && all(is.finite(res$y)))
})

test_that("ORF coverage is reported per strand and respects offset", {
  skip_if_no_rnaplot()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  orfs <- data.frame(
    seq_id = "a", strand = c("+", "-"),
    start = c(1L, 20L), end = c(10L, 30L),
    stringsAsFactors = FALSE
  )
  res <- plot_rna_structure(fold_stub, orfs = orfs)
  expect_equal(sum(res$n_plus  > 0L), 10L)
  expect_equal(sum(res$n_minus > 0L), 11L)

  # same ORFs expressed in genomic coordinates
  shifted <- transform(orfs, start = start + 999L, end = end + 999L)
  res2 <- plot_rna_structure(fold_stub, orfs = shifted, offset = 1000L)
  expect_equal(res2$n_plus,  res$n_plus)
  expect_equal(res2$n_minus, res$n_minus)
})

test_that("ORFs outside the folded region are clipped or dropped", {
  skip_if_no_rnaplot()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  outside <- data.frame(seq_id = "a", strand = "+", start = 500L, end = 600L,
                        stringsAsFactors = FALSE)
  expect_message(res <- plot_rna_structure(fold_stub, orfs = outside),
                 "outside the folded region")
  expect_equal(sum(res$n_plus), 0L)

  overhang <- data.frame(seq_id = "a", strand = "+", start = 35L, end = 900L,
                         stringsAsFactors = FALSE)
  res2 <- plot_rna_structure(fold_stub, orfs = overhang)
  expect_equal(sum(res2$n_plus > 0L), nchar(fold_stub$structure) - 34L)
})

test_that("multi-sequence ORF tables need a seq_id, and wrapping ORFs warn", {
  skip_if_no_rnaplot()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  two <- data.frame(seq_id = c("a", "b"), strand = "+",
                    start = c(1L, 5L), end = c(10L, 15L),
                    stringsAsFactors = FALSE)
  expect_error(plot_rna_structure(fold_stub, orfs = two), "seq_id")
  res <- plot_rna_structure(fold_stub, orfs = two, seq_id = "a")
  expect_equal(sum(res$n_plus > 0L), 10L)

  wrapped <- data.frame(seq_id = "a", strand = "+", start = 1L, end = 10L,
                        wraps_around = TRUE, stringsAsFactors = FALSE)
  expect_warning(plot_rna_structure(fold_stub, orfs = wrapped),
                 "origin-spanning")
})

test_that("pair_col accepts a scalar, a vector, and a function", {
  skip_if_no_rnaplot()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  npair <- nrow(MitoDiveR:::.parse_dot_bracket(fold_stub$structure))
  expect_silent(plot_rna_structure(fold_stub, pair_col = "grey80"))
  expect_silent(plot_rna_structure(fold_stub,
                                   pair_col = rep("red", npair)))
  expect_silent(plot_rna_structure(
    fold_stub,
    pair_col = function(bp) ifelse(bp$pos3 - bp$pos5 > 10, "red", "grey80")))
  expect_error(plot_rna_structure(fold_stub, pair_col = c("red", "blue")),
               "per base pair")
})

test_that("an mdp_name column is accepted and does not disturb coverage", {
  skip_if_no_rnaplot()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  orfs <- data.frame(
    seq_id = "a", strand = c("+", "-"),
    start = c(1L, 20L), end = c(10L, 30L),
    mdp_name = c("Humanin", NA),
    mdp_pct_identity = c(100, NA),
    stringsAsFactors = FALSE
  )
  res <- plot_rna_structure(fold_stub, orfs = orfs)
  expect_equal(sum(res$n_plus  > 0L), 10L)
  expect_equal(sum(res$n_minus > 0L), 11L)

  # every label style draws
  for (style in c("name", "name_pct", "none"))
    expect_silent(plot_rna_structure(fold_stub, orfs = orfs,
                                     mdp_label = style))
  expect_error(plot_rna_structure(fold_stub, orfs = orfs, mdp_label = "nope"),
               "'arg'")
})

test_that("mdp_colors accepts a named vector or a function, and reports gaps", {
  skip_if_no_rnaplot()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  orfs <- data.frame(seq_id = "a", strand = "+", start = c(1L, 20L),
                     end = c(10L, 30L), mdp_name = c("Humanin", "SHLP6"),
                     stringsAsFactors = FALSE)

  expect_silent(plot_rna_structure(
    fold_stub, orfs = orfs,
    mdp_colors = c(Humanin = "red", SHLP6 = "blue")))
  expect_silent(plot_rna_structure(
    fold_stub, orfs = orfs,
    mdp_colors = function(n) rep("red", n)))

  expect_error(plot_rna_structure(fold_stub, orfs = orfs,
                                  mdp_colors = c(Humanin = "red")),
               "no entry for: SHLP6")
  expect_error(plot_rna_structure(fold_stub, orfs = orfs,
                                  mdp_colors = c("red", "blue")),
               "must be named")
})

test_that("ORFs with no MDP still draw when others are named", {
  skip_if_no_rnaplot()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  orfs <- data.frame(seq_id = "a", strand = c("+", "-"),
                     start = c(1L, 20L), end = c(10L, 30L),
                     mdp_name = c("Humanin", NA), stringsAsFactors = FALSE)
  res <- plot_rna_structure(fold_stub, orfs = orfs, unnamed_alpha = 1)
  expect_equal(sum(res$n_minus > 0L), 11L)
})
