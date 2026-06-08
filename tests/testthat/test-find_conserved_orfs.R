library(testthat)
library(MitoDiveR)

# ---- helpers -----------------------------------------------------------------
# Three genomes sharing one clearly-homologous protein (~96% identity), each
# also carrying a genome-private ORF that should never cluster.
shared_a <- "MFADRWLQGSTKVPEAILNCYDEHKR"
shared_b <- "MFADRWLQGSTKVPEAILNCYDEHKH"   # 1 aa diff vs a
shared_c <- "MFADRWLQGSTKVPEAILNCYDERKR"   # 1 aa diff vs a

three_genome_orfs <- function() {
  data.frame(
    seq_id           = c("g1", "g1", "g2", "g2", "g3"),
    start            = c(10, 500, 12, 600, 15),
    end              = c(88, 540, 90, 640, 93),
    protein_sequence = c(paste0(shared_a, "*"), "MQQQWWWEEE",
                         paste0(shared_b, "*"), "MKKKLLLPPP",
                         paste0(shared_c, "*")),
    stringsAsFactors = FALSE
  )
}

# A pair that shares only a 6-aa local motif (MFADRW). Under PID1 the local
# alignment is 100% identical, but it covers only 6/20 of the shorter protein.
partial_motif_orfs <- function() {
  data.frame(
    seq_id           = c("g1", "g2"),
    start            = c(10, 20),
    end              = c(70, 80),
    protein_sequence = c("MFADRWLLLLLLLLLLLLLL", "QQQQQQQQQQQQQQMFADRW"),
    stringsAsFactors = FALSE
  )
}

# ---- core clustering ---------------------------------------------------------
test_that("clusters a shared protein across all genomes", {
  res <- find_conserved_orfs(three_genome_orfs(), min_identity = 60,
                             verbose = FALSE)

  expect_s3_class(res, "data.frame")
  expect_true(all(c("cluster_id", "n_sequences", "cluster_size") %in% names(res)))
  # Only the shared protein survives; the three private ORFs do not.
  expect_equal(nrow(res), 3L)
  expect_equal(unique(res$cluster_id), 1L)
  expect_equal(unique(res$n_sequences), 3L)
  expect_equal(unique(res$cluster_size), 3L)
  expect_setequal(res$seq_id, c("g1", "g2", "g3"))
})

test_that("returns an empty annotated frame when nothing clusters", {
  res <- suppressMessages(
    find_conserved_orfs(three_genome_orfs(), min_identity = 100,
                        min_coverage = 1, verbose = FALSE)
  )
  expect_equal(nrow(res), 0L)
  expect_true(all(c("cluster_id", "n_sequences", "cluster_size") %in% names(res)))
})

# ---- min_genomes override ----------------------------------------------------
test_that("min_genomes overrides the min_sequences fraction", {
  orfs <- three_genome_orfs()
  # Drop g3's shared copy: the cluster now spans only 2 genomes.
  orfs <- orfs[orfs$seq_id != "g3", ]

  kept2 <- find_conserved_orfs(orfs, min_identity = 60, min_genomes = 2,
                               verbose = FALSE)
  expect_equal(nrow(kept2), 2L)

  kept3 <- find_conserved_orfs(orfs, min_identity = 60, min_genomes = 3,
                               verbose = FALSE)
  expect_equal(nrow(kept3), 0L)
})

# ---- coverage gate -----------------------------------------------------------
test_that("default PID3 metric rejects a short high-identity motif", {
  # PID3 = matches / alignment length; 6 matches over a 20-residue protein = 30%.
  res <- find_conserved_orfs(partial_motif_orfs(), min_identity = 80,
                             min_genomes = 2, verbose = FALSE)
  expect_equal(nrow(res), 0L)
})

test_that("min_coverage gates a short motif that PID1 would otherwise link", {
  no_gate <- find_conserved_orfs(partial_motif_orfs(), min_identity = 80,
                                 identity_metric = "PID1", min_genomes = 2,
                                 verbose = FALSE)
  expect_equal(nrow(no_gate), 2L)   # PID1 = 100 links them without a gate

  with_gate <- find_conserved_orfs(partial_motif_orfs(), min_identity = 80,
                                   identity_metric = "PID1", min_genomes = 2,
                                   min_coverage = 0.5, verbose = FALSE)
  expect_equal(nrow(with_gate), 0L) # coverage 0.3 < 0.5 drops the link
})

# ---- input validation --------------------------------------------------------
test_that("validates inputs", {
  good <- three_genome_orfs()

  expect_error(find_conserved_orfs(good[, "seq_id", drop = FALSE]),
               "protein_sequence")
  expect_error(find_conserved_orfs(good[, "protein_sequence", drop = FALSE]),
               "seq_id")
  expect_error(find_conserved_orfs(good[good$seq_id == "g1", ]),
               "at least two sequences")
  expect_error(find_conserved_orfs(good, min_identity = 150),
               "between 0 and 100")
  expect_error(find_conserved_orfs(good, min_coverage = 2),
               "min_coverage")
  expect_error(find_conserved_orfs(good, min_genomes = 0),
               "min_genomes")
  expect_error(find_conserved_orfs(good, identity_metric = "PIDX"))
})

# ---- deprecation shim --------------------------------------------------------
test_that("cluster_conserved_orfs() is deprecated", {
  expect_warning(
    suppressMessages(cluster_conserved_orfs(three_genome_orfs(), verbose = FALSE)),
    "deprecated"
  )
})

test_that("cluster_conserved_orfs() matches its find_conserved_orfs() equivalent", {
  orfs <- three_genome_orfs()

  shim <- suppressWarnings(suppressMessages(
    cluster_conserved_orfs(orfs, min_pct_identity = 40, min_coverage = 0.5,
                           gap_extension = 4, verbose = FALSE)
  ))
  equiv <- find_conserved_orfs(orfs, min_identity = 40, min_coverage = 0.5,
                               min_genomes = 2, identity_metric = "PID1",
                               gap_extension = 4, verbose = FALSE)

  expect_equal(shim, equiv)
  expect_gt(nrow(shim), 0L)   # guard against trivially-equal empty frames
})
