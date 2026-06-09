#!/usr/bin/env Rscript
# ============================================================================
# validate_find_orf_orthologs.R  —  ONLINE-ONLY smoke test for the network
# paths of find_orf_orthologs() (method = "tblastn" / "blastn").
#
# The local confirmation core is covered by tests/testthat/test-find_orf_orthologs.R.
# This script exercises the parts that need NCBI and therefore cannot run in
# unit tests: ENTREZ_QUERY scoping, windowed efetch harvest, and BLAST-XML
# subject-coordinate parsing.  Run it by hand when you are online:
#
#     Rscript dev/validate_find_orf_orthologs.R
#
# NCBI BLAST jobs take minutes; the whole script may run 5–15 min.  Set a
# contact email (NCBI asks for one) and, optionally, an API key:
#
#     export ENTREZ_KEY=xxxx           # optional, raises the rate limit
#     export MITODIVER_EMAIL=you@x.edu # recommended
# ============================================================================

suppressMessages({
  library(MitoDiveR)
  library(Biostrings)
})

options(MitoDiveR.email = Sys.getenv("MITODIVER_EMAIL",
                                     unset = "mitodiver@example.com"))
api_key <- Sys.getenv("ENTREZ_KEY", unset = "")
if (!nzchar(api_key)) api_key <- NULL

# Humanin: a real MDP encoded within MT-RNR2 (16S rRNA) on the human heavy
# strand; translated by the cytoplasmic ribosome -> standard code (SGC0).
HUMANIN <- c(Humanin = "MAPRGFSCLLLLTSEIDLPVKRRA")
HUMAN_ACC <- "NC_012920.1"

# ---- tiny reporting harness --------------------------------------------------
.results <- list()
check <- function(label, expr) {
  t0 <- Sys.time()
  val <- tryCatch(expr, error = function(e) structure(list(msg = conditionMessage(e)),
                                                      class = "validation_error"))
  secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  ok <- !inherits(val, "validation_error")
  cat(sprintf("[%s] %-46s (%ss)\n", if (ok) "PASS" else "FAIL", label, secs))
  if (!ok) cat("       ", val$msg, "\n")
  .results[[label]] <<- ok
  invisible(val)
}

cat("== find_orf_orthologs() network validation ==\n\n")

# ---------------------------------------------------------------------------
# Baseline (local, real fetch): recover humanin in the human mitogenome.
# Confirms fetch -> scan -> score -> annotate on a real genome before BLAST.
# ---------------------------------------------------------------------------
human <- check("fetch human mitogenome", {
  fetch_sequences(accessions = HUMAN_ACC)
})

local_hit <- check("local: recover humanin in human mtDNA", {
  h <- find_orf_orthologs(HUMANIN, genomes = human, method = "local",
                          genetic_code = "SGC0", min_orf_length = 30,
                          min_identity = 80, min_coverage = 60,
                          annotate = TRUE, api_key = api_key)
  stopifnot(nrow(h) >= 1L, max(h$pct_identity) >= 80)
  cat("       best identity:", max(h$pct_identity),
      "| region:", h$genomic_region[which.max(h$pct_identity)], "\n")
  cat("       (expect a 12S/16S rRNA region label — humanin sits in MT-RNR2)\n")
  h
})

# Build a conserved host window from the recovered humanin locus (± 60 nt of
# 16S) to use as the BLASTN query in the host-anchored test.
host_region <- check("derive 16S host window for BLASTN", {
  if (inherits(local_hit, "validation_error") || nrow(local_hit) == 0L)
    stop("no local hit to derive a host window from")
  best <- local_hit[which.max(local_hit$pct_identity), ]
  lo <- max(1L, best$start - 60L); hi <- best$end + 60L
  win <- as.character(Biostrings::subseq(human[[1]], start = lo, end = hi))
  cat("       host window:", nchar(win), "nt around", best$start, "-", best$end, "\n")
  win
})

# ---------------------------------------------------------------------------
# Remote path A — TBLASTN with the peptide, scoped to RefSeq mitochondrion.
# Keep the scope small (a few hominids) so the job returns quickly.
# ---------------------------------------------------------------------------
tblastn_hits <- check("tblastn: peptide search across Hominidae", {
  h <- find_orf_orthologs(HUMANIN, method = "tblastn",
                          taxa = "Hominidae", genetic_code = "SGC0",
                          min_orf_length = 30, min_identity = 40,
                          min_coverage = 50, max_genomes = 5,
                          poll_interval = 20L, timeout = 600L,
                          api_key = api_key)
  cat("       hits:", nrow(h), "\n")
  if (nrow(h) > 0L) {
    stopifnot(all(c("blast_accession", "blast_from", "blast_to") %in% names(h)))
    stopifnot(!all(is.na(h$blast_accession)))   # provenance must populate
    cat("       accessions:", paste(unique(h$blast_accession), collapse = ", "), "\n")
    cat("       window seq_id example:", h$seq_id[1], "(expect 'acc:lo-hi')\n")
  } else {
    cat("       WARNING: zero hits — check ENTREZ_QUERY scope / thresholds\n")
  }
  h
})

# ---------------------------------------------------------------------------
# Remote path B — BLASTN with the conserved 16S host window.
# ---------------------------------------------------------------------------
blastn_hits <- check("blastn: host-window search across Hominidae", {
  if (inherits(host_region, "validation_error"))
    stop("no host_region available")
  h <- find_orf_orthologs(HUMANIN, method = "blastn",
                          host_region = host_region, taxa = "Hominidae",
                          genetic_code = "SGC0", min_orf_length = 30,
                          min_identity = 40, min_coverage = 50,
                          max_genomes = 5, poll_interval = 20L,
                          timeout = 600L, api_key = api_key)
  cat("       hits:", nrow(h), "\n")
  if (nrow(h) > 0L) {
    stopifnot(!all(is.na(h$blast_accession)))
    cat("       accessions:", paste(unique(h$blast_accession), collapse = ", "), "\n")
  } else {
    cat("       WARNING: zero hits — check host window / scope\n")
  }
  h
})

# ---- summary -----------------------------------------------------------------
cat("\n== summary ==\n")
for (nm in names(.results))
  cat(sprintf("  %-46s %s\n", nm, if (.results[[nm]]) "PASS" else "FAIL"))
cat("\nRisk points this exercises: ENTREZ_QUERY refseq+mitochondrion scoping,\n",
    "windowed efetch (seq_start/seq_stop), Hsp_hit-from/to XML parsing,\n",
    "DATABASE='nt' reachability, and 'acc:lo-hi' window provenance.\n", sep = "")

if (any(!unlist(.results))) quit(status = 1L)
