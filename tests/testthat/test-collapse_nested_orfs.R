test_that("nested starts sharing a stop collapse to the longest", {
  orfs <- data.frame(
    seq_id = "a", strand = "+", frame = 1L,
    start  = c(1L, 10L, 19L), end = c(36L, 36L, 36L),
    length_nt = c(36L, 27L, 18L),
    stringsAsFactors = FALSE
  )
  out <- collapse_nested_orfs(orfs)
  expect_equal(nrow(out), 1L)
  expect_equal(out$start, 1L)
  expect_equal(out$length_nt, 36L)
})

test_that("identical coordinates in different sequences are kept separately", {
  # regression: grouping once omitted seq_id, so an ORF present at the same
  # coordinates in two genomes was collapsed to a single row
  orfs <- data.frame(
    seq_id = c("a", "b"), strand = "+", frame = 3L,
    start  = c(963L, 963L), end = c(1037L, 1037L),
    length_nt = c(75L, 75L),
    stringsAsFactors = FALSE
  )
  out <- collapse_nested_orfs(orfs)
  expect_equal(nrow(out), 2L)
  expect_setequal(out$seq_id, c("a", "b"))
})

test_that("non-overlapping ORFs in the same frame are all retained", {
  orfs <- data.frame(
    seq_id = "a", strand = "+", frame = 3L,
    start  = c(963L, 1320L), end = c(1037L, 1382L),
    length_nt = c(75L, 63L),
    stringsAsFactors = FALSE
  )
  expect_equal(nrow(collapse_nested_orfs(orfs)), 2L)
})

test_that("opposite strands at the same coordinates do not collapse", {
  orfs <- data.frame(
    seq_id = "a", strand = c("+", "-"), frame = c(1L, -1L),
    start  = c(100L, 100L), end = c(160L, 160L),
    length_nt = c(61L, 61L),
    stringsAsFactors = FALSE
  )
  expect_equal(nrow(collapse_nested_orfs(orfs)), 2L)
})

test_that("minus-strand ORFs group on their stop (the start column)", {
  # start <= end convention: on '-' the stop sits in `start`
  orfs <- data.frame(
    seq_id = "a", strand = "-", frame = -2L,
    start  = c(50L, 50L), end = c(140L, 110L),
    length_nt = c(91L, 61L),
    stringsAsFactors = FALSE
  )
  out <- collapse_nested_orfs(orfs)
  expect_equal(nrow(out), 1L)
  expect_equal(out$end, 140L)
})

test_that("empty and malformed input are handled", {
  expect_equal(nrow(collapse_nested_orfs(data.frame())), 0L)
  bad <- data.frame(start = 1L, end = 10L)
  expect_warning(res <- collapse_nested_orfs(bad), "frame")
  expect_identical(res, bad)
})
