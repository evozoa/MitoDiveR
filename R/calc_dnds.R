# ============================================================================
# calc_dnds — pairwise dN/dS (Ka/Ks) for conserved ORF clusters
# ============================================================================

#' Calculate pairwise dN/dS for Conserved ORF Clusters
#'
#' Estimates the ratio of nonsynonymous to synonymous substitution rates
#' (dN/dS, also written Ka/Ks or ω) for each orthologous ORF cluster
#' identified by [find_conserved_orfs()].  Protein sequences are globally
#' aligned to guide a codon-level alignment, and substitution counts are
#' computed using the Nei-Gojobori (1986) method with Jukes-Cantor distance
#' correction.
#'
#' @param orfs `data.frame` from [find_conserved_orfs()].  Must contain
#'   `cluster_id`, `seq_id`, `orf_sequence`, and `protein_sequence` columns.
#' @param genetic_code Character scalar.  Biostrings genetic-code identifier.
#'   Default `"SGC1"` (vertebrate mitochondrial code).  Must match the code
#'   used when scanning ORFs.
#'
#' @return A `data.frame` with one row per sequence pair per cluster:
#' \describe{
#'   \item{`cluster_id`}{Cluster identifier from [find_conserved_orfs()].}
#'   \item{`seq_id_1`, `seq_id_2`}{Identifiers of the two sequences
#'     compared.}
#'   \item{`n_codons`}{Number of codon positions in the alignment (gap-free
#'     columns only).}
#'   \item{`S`}{Total synonymous sites (averaged across both sequences).}
#'   \item{`N`}{Total nonsynonymous sites.}
#'   \item{`Sd`}{Observed synonymous differences.}
#'   \item{`Nd`}{Observed nonsynonymous differences.}
#'   \item{`pS`}{Proportion of synonymous differences (`Sd / S`).}
#'   \item{`pN`}{Proportion of nonsynonymous differences (`Nd / N`).}
#'   \item{`dS`}{Synonymous distance after Jukes-Cantor correction.  `NA`
#'     when `pS >= 0.75` (saturation).}
#'   \item{`dN`}{Nonsynonymous distance after Jukes-Cantor correction.  `NA`
#'     when `pN >= 0.75`.}
#'   \item{`omega`}{dN/dS ratio.  `NA` when either distance is `NA`.
#'     `Inf` when `dS = 0` and `dN > 0` (no synonymous change detected).
#'     Values < 1 indicate purifying selection; ~1 neutral evolution;
#'     > 1 positive selection.}
#' }
#' Returns an empty `data.frame` invisibly (with a message) if no valid codon
#' alignments are found.
#'
#' @details
#' ## Method
#' The Nei-Gojobori (1986) method counts synonymous and nonsynonymous sites
#' per codon by considering all single-nucleotide changes at each codon
#' position.  For codon pairs differing at two or three positions, observed
#' differences are averaged over all possible mutational pathways, excluding
#' paths that pass through stop codons.  The raw proportions pS and pN are
#' corrected with the Jukes-Cantor formula:
#' \deqn{d = -\tfrac{3}{4}\ln\!\left(1 - \tfrac{4}{3}\,p\right)}
#' Codons containing ambiguous bases (non-ACGT) and aligned stop codons are
#' skipped.
#'
#' ## Interpretation caveats
#' dN/dS is estimated per cluster from as few as two sequences, so individual
#' values can be noisy.  Very short ORFs (< 100 nt) yield wide confidence
#' intervals.  A saturated synonymous distance (`pS >= 0.75`) returns `NA`
#' for dS and omega; this is common for highly diverged sequences or
#' ORFs in fast-evolving regions.
#'
#' @references
#' Nei M, Gojobori T (1986) Simple methods for estimating the numbers of
#' synonymous and nonsynonymous nucleotide substitutions.
#' *Molecular Biology and Evolution* **3**(5):418–426.
#'
#' @seealso [find_conserved_orfs()], [calc_codon_usage()]
#'
#' @examples
#' \dontrun{
#' seqs      <- fetch_sequences(c("NC_012920.1", "NC_011137.1"))
#' orfs      <- scan_orfs(seqs, min_orf_length = 60)
#' conserved <- find_conserved_orfs(orfs, min_identity = 50)
#' dnds      <- calc_dnds(conserved)
#'
#' # Clusters under putative positive selection
#' dnds[!is.na(dnds$omega) & dnds$omega > 1, ]
#'
#' # Clusters with strong purifying selection
#' dnds[!is.na(dnds$omega) & dnds$omega < 0.3, ]
#' }
#'
#' @importFrom Biostrings getGeneticCode
#' @importFrom pwalign pairwiseAlignment pattern subject
#' @importFrom utils combn
#' @export
calc_dnds <- function(orfs, genetic_code = "SGC1") {

  ## ---- validate ---------------------------------------------------------------
  if (!is.data.frame(orfs))
    stop("'orfs' must be a data.frame (output of find_conserved_orfs()).",
         call. = FALSE)
  required <- c("cluster_id", "seq_id", "orf_sequence", "protein_sequence")
  missing  <- setdiff(required, names(orfs))
  if (length(missing) > 0L)
    stop("'orfs' is missing columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  if (nrow(orfs) == 0L)
    stop("'orfs' is empty.", call. = FALSE)

  code     <- Biostrings::getGeneticCode(genetic_code)
  stop_cds <- names(code)[code == "*"]

  ## ---- process each cluster ---------------------------------------------------
  cluster_ids <- unique(orfs$cluster_id)
  message("Computing dN/dS for ", length(cluster_ids), " cluster(s)...")

  results <- lapply(cluster_ids, function(cid) {
    cl      <- orfs[orfs$cluster_id == cid, , drop = FALSE]
    seq_ids <- unique(cl$seq_id)
    if (length(seq_ids) < 2L) return(NULL)

    pairs <- utils::combn(seq_ids, 2L, simplify = FALSE)

    lapply(pairs, function(pair) {
      r1 <- cl[cl$seq_id == pair[1L], , drop = FALSE][1L, ]
      r2 <- cl[cl$seq_id == pair[2L], , drop = FALSE][1L, ]

      nt1   <- .strip_stop_nt(r1$orf_sequence,   r1$protein_sequence)
      nt2   <- .strip_stop_nt(r2$orf_sequence,   r2$protein_sequence)
      prot1 <- gsub("\\*", "", r1$protein_sequence)
      prot2 <- gsub("\\*", "", r2$protein_sequence)

      if (nchar(nt1) < 3L || nchar(nt2) < 3L) return(NULL)

      ## global protein alignment to guide codon alignment
      aln   <- pwalign::pairwiseAlignment(
        pattern            = prot2,
        subject            = prot1,
        type               = "global",
        substitutionMatrix = "BLOSUM62",
        gapOpening         = 10,
        gapExtension       = 0.5
      )
      aln1 <- as.character(pwalign::subject(aln))
      aln2 <- as.character(pwalign::pattern(aln))

      ## thread nucleotide sequences through protein alignment
      cdn1 <- .thread_codons(nt1, aln1)
      cdn2 <- .thread_codons(nt2, aln2)

      ## keep only fully aligned codon pairs
      keep <- cdn1 != "---" & cdn2 != "---"
      cdn1 <- cdn1[keep];  cdn2 <- cdn2[keep]
      n_aligned <- length(cdn1)
      if (n_aligned == 0L) return(NULL)

      ## accumulate S, N, Sd, Nd
      S_tot <- 0; N_tot <- 0; Sd_tot <- 0; Nd_tot <- 0

      for (i in seq_len(n_aligned)) {
        c1 <- cdn1[i]; c2 <- cdn2[i]
        if (grepl("[^ACGT]", c1, perl = TRUE) ||
            grepl("[^ACGT]", c2, perl = TRUE)) next
        if (c1 %in% stop_cds || c2 %in% stop_cds) next

        s1 <- .ng_syn_sites(c1, code, stop_cds)
        s2 <- .ng_syn_sites(c2, code, stop_cds)
        avg_s  <- (s1 + s2) / 2
        S_tot  <- S_tot + avg_s
        N_tot  <- N_tot + (3 - avg_s)

        diff <- .ng_codon_diff(c1, c2, code, stop_cds)
        if (!anyNA(diff)) {
          Sd_tot <- Sd_tot + diff[1L]
          Nd_tot <- Nd_tot + diff[2L]
        }
      }

      pS <- if (S_tot > 0) Sd_tot / S_tot else NA_real_
      pN <- if (N_tot > 0) Nd_tot / N_tot else NA_real_
      dS <- .jc_correct(pS)
      dN <- .jc_correct(pN)

      omega <- if (anyNA(c(dS, dN))) NA_real_
               else if (dS == 0 && dN == 0)  NA_real_
               else if (dS == 0)              Inf
               else                           dN / dS

      data.frame(
        cluster_id = cid,
        seq_id_1   = pair[1L],
        seq_id_2   = pair[2L],
        n_codons   = n_aligned,
        S          = round(S_tot,  3L),
        N          = round(N_tot,  3L),
        Sd         = round(Sd_tot, 3L),
        Nd         = round(Nd_tot, 3L),
        pS         = round(pS,     4L),
        pN         = round(pN,     4L),
        dS         = round(dS,     4L),
        dN         = round(dN,     4L),
        omega      = round(omega,  4L),
        stringsAsFactors = FALSE
      )
    })
  })

  ## ---- flatten and return -----------------------------------------------------
  flat <- Filter(Negate(is.null), unlist(results, recursive = FALSE))
  if (length(flat) == 0L) {
    message("No valid codon alignments found.")
    return(invisible(data.frame()))
  }
  out <- do.call(rbind, flat)
  rownames(out) <- NULL
  out
}


# --------------------------------------------------------------------------- #
# Internal helpers                                                             #
# --------------------------------------------------------------------------- #

# Strip the terminal stop codon (3 nt) from a nucleotide ORF sequence when
# the corresponding protein ends with "*".
.strip_stop_nt <- function(nt_seq, prot_seq) {
  if (grepl("\\*$", prot_seq))
    nt_seq <- substr(nt_seq, 1L, nchar(nt_seq) - 3L)
  nt_seq
}

# Thread a nucleotide sequence through a gapped amino-acid alignment string,
# returning a character vector of codons with "---" at gap positions.
.thread_codons <- function(nt_seq, gapped_aa) {
  aa_chars <- strsplit(gapped_aa, "")[[1L]]
  codons   <- character(length(aa_chars))
  nt_pos   <- 1L
  for (i in seq_along(aa_chars)) {
    if (aa_chars[i] == "-") {
      codons[i] <- "---"
    } else {
      codons[i] <- substr(nt_seq, nt_pos, nt_pos + 2L)
      nt_pos    <- nt_pos + 3L
    }
  }
  codons
}

# Count synonymous sites in a single codon (Nei-Gojobori 1986).
# Returns S for the codon; N = 3 - S.
.ng_syn_sites <- function(codon, code, stop_cds) {
  bases      <- c("A", "C", "G", "T")
  ref_aa     <- code[[codon]]
  codon_chars <- strsplit(codon, "")[[1L]]
  s <- 0
  for (pos in 1:3) {
    n_syn <- 0L
    for (b in bases) {
      if (b == codon_chars[pos]) next
      mut_chars      <- codon_chars
      mut_chars[pos] <- b
      mut            <- paste(mut_chars, collapse = "")
      if (!mut %in% names(code) || mut %in% stop_cds) next
      if (code[[mut]] == ref_aa) n_syn <- n_syn + 1L
    }
    s <- s + n_syn / 3
  }
  s
}

# Count synonymous (Sd) and nonsynonymous (Nd) differences for a codon pair
# (Nei-Gojobori 1986).  For codons differing at >1 position, averages over
# all valid mutational pathways (excluding those through stop codons).
# Returns c(Sd, Nd); returns c(NA, NA) if all pathways pass through stops.
.ng_codon_diff <- function(c1, c2, code, stop_cds) {
  if (c1 == c2) return(c(0, 0))

  chars1    <- strsplit(c1, "")[[1L]]
  chars2    <- strsplit(c2, "")[[1L]]
  diff_pos  <- which(chars1 != chars2)

  if (length(diff_pos) == 1L) {
    syn <- as.integer(code[[c1]] == code[[c2]])
    return(c(syn, 1L - syn))
  }

  perms    <- .permutations(diff_pos)
  total_sd <- 0;  total_nd <- 0;  valid <- 0L

  for (perm in perms) {
    path_sd <- 0;  path_nd <- 0;  ok <- TRUE
    cur     <- chars1

    for (pos in perm) {
      nxt      <- cur
      nxt[pos] <- chars2[pos]
      cur_cdn  <- paste(cur, collapse = "")
      nxt_cdn  <- paste(nxt, collapse = "")

      if (nxt_cdn %in% stop_cds) { ok <- FALSE; break }

      if (code[[cur_cdn]] == code[[nxt_cdn]]) {
        path_sd <- path_sd + 1L
      } else {
        path_nd <- path_nd + 1L
      }
      cur <- nxt
    }

    if (ok) {
      total_sd <- total_sd + path_sd
      total_nd <- total_nd + path_nd
      valid    <- valid + 1L
    }
  }

  if (valid == 0L) return(c(NA_real_, NA_real_))
  c(total_sd / valid, total_nd / valid)
}

# All permutations of a vector (used for multi-difference codon pathways).
.permutations <- function(x) {
  if (length(x) <= 1L) return(list(x))
  result <- list()
  for (i in seq_along(x)) {
    for (p in .permutations(x[-i]))
      result[[length(result) + 1L]] <- c(x[i], p)
  }
  result
}

# Jukes-Cantor distance correction: d = -3/4 * ln(1 - 4/3 * p).
# Returns NA when p >= 0.75 (saturation) or p is NA.
.jc_correct <- function(p) {
  if (is.na(p) || p >= 0.75) return(NA_real_)
  if (p == 0)                 return(0)
  -0.75 * log(1 - (4 / 3) * p)
}
