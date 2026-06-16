library(testthat)
library(MitoDiveR)

# The reference-anchored depth math lives in the un-exported core
# .orf_depth_profile(); these tests exercise it directly with synthetic
# alignments + ORF tables (no network, no MAFFT).
depth_profile <- MitoDiveR:::.orf_depth_profile

# Two ungapped 30-nt genomes, reference = g1, so genome bp == reference bp.
flat_aln <- c(g1 = strrep("A", 30L), g2 = strrep("A", 30L))

test_that("species unit counts each species once; orf unit stacks every ORF", {
  # g1 carries two identical +1 ORFs at 5-14 (same species), g2 one at 5-14.
  orfs <- data.frame(
    seq_id = c("g1", "g1", "g2"),
    frame  = c(1L, 1L, 1L),
    start  = c(5L, 5L, 5L),
    end    = c(14L, 14L, 14L),
    stringsAsFactors = FALSE)
  grp <- c(g1 = "sp1", g2 = "sp2")

  d_sp <- depth_profile(flat_aln, orfs, reference = "g1",
                        depth_unit = "species", group_of = grp)
  d_or <- depth_profile(flat_aln, orfs, reference = "g1",
                        depth_unit = "orf")

  expect_equal(dim(d_sp), c(30L, 7L))
  expect_equal(colnames(d_sp),
               c("+1", "+2", "+3", "-1", "-2", "-3", "sum"))
  # species: two species present -> depth 2 across the locus
  expect_equal(max(d_sp[, "+1"]), 2)
  expect_true(all(d_sp[5:14, "+1"] == 2))
  expect_equal(unname(d_sp[1, "+1"]), 0)  # outside the locus
  # orf: three ORFs stack -> depth 3
  expect_equal(max(d_or[, "+1"]), 3)
  expect_true(all(d_or[5:14, "+1"] == 3))
})

test_that("ORFs land only in their own frame track and sum is the row total", {
  orfs <- data.frame(
    seq_id = c("g1", "g2"),
    frame  = c(1L, -2L),
    start  = c(5L, 20L),
    end    = c(14L, 25L),
    stringsAsFactors = FALSE)
  grp <- c(g1 = "sp1", g2 = "sp2")

  d <- depth_profile(flat_aln, orfs, reference = "g1",
                     depth_unit = "species", group_of = grp)

  expect_true(all(d[5:14, "+1"] == 1))
  expect_true(all(d[5:14, c("+2","+3","-1","-2","-3")] == 0))
  expect_true(all(d[20:25, "-2"] == 1))
  expect_true(all(d[20:25, c("+1","+2","+3","-1","-3")] == 0))
  # sum column equals the sum of the six frame columns at every position
  expect_equal(unname(d[, "sum"]),
               unname(rowSums(d[, c("+1","+2","+3","-1","-2","-3")])))
})

test_that("wrap-around ORFs cover both ends of the genome", {
  orfs <- data.frame(
    seq_id       = "g1",
    frame        = 3L,
    start        = 28L,
    end          = 3L,
    wraps_around = TRUE,
    stringsAsFactors = FALSE)

  d <- depth_profile(flat_aln, orfs, reference = "g1", depth_unit = "orf")

  expect_true(all(d[c(28:30, 1:3), "+3"] == 1))   # both sides of the origin
  expect_equal(unname(d[15, "+3"]), 0)             # untouched middle
})

test_that("columns inserted relative to the reference fold onto the preceding bp", {
  # g1 (reference) has a 5-col gap; g2 fills it. ref ungapped length = 10.
  aln <- c(g1 = "AAAAA-----AAAAA",      # ref: bp 1..10
           g2 = "AAAAACCCCCAAAAA")      # bp  : 1..15
  # g2 ORF over its bp 6-10 sits entirely in g1's insertion -> ref bp 5.
  orfs <- data.frame(seq_id = "g2", frame = 1L, start = 6L, end = 10L,
                     stringsAsFactors = FALSE)

  d <- depth_profile(aln, orfs, reference = "g1", depth_unit = "orf")

  expect_equal(nrow(d), 10L)             # reference length, not 15
  expect_equal(unname(d[5, "+1"]), 1)    # insertion folds onto preceding ref bp
  expect_true(all(d[6:10, "+1"] == 0))
})
