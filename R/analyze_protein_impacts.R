# ============================================================================
# analyze_protein_impacts — protein-level impact table from compare_mitogenomes
# ============================================================================

# ── Lookup tables ─────────────────────────────────────────────────────────────

# Chou-Fasman helix (Pa) and sheet (Pb) propensities (Chou & Fasman 1978)
.CF_HELIX <- c(
  A=1.42, R=0.98, N=0.67, D=1.01, C=0.70, Q=1.11, E=1.51, G=0.57,
  H=1.00, I=1.08, L=1.21, K=1.16, M=1.45, F=1.13, P=0.57, S=0.77,
  T=0.83, W=1.08, Y=0.69, V=1.06
)
.CF_SHEET <- c(
  A=0.83, R=0.93, N=0.89, D=0.54, C=1.19, Q=1.10, E=0.37, G=0.75,
  H=0.87, I=1.60, L=1.30, K=0.74, M=1.05, F=1.38, P=0.55, S=0.75,
  T=1.19, W=1.37, Y=1.47, V=1.70
)

# Kyte-Doolittle hydrophobicity scale
.KD_HYDRO <- c(
  A= 1.8, R=-4.5, N=-3.5, D=-3.5, C= 2.5, Q=-3.5, E=-3.5, G=-0.4,
  H=-3.2, I= 4.5, L= 3.8, K=-3.9, M= 1.9, F= 2.8, P=-1.6, S=-0.8,
  T=-0.7, W=-0.9, Y=-1.3, V= 4.2
)

# Physicochemical class per amino acid
.PHYSCHEM <- c(
  A="nonpolar", C="polar",   D="acidic",   E="acidic",   F="nonpolar",
  G="nonpolar", H="basic",   I="nonpolar", K="basic",    L="nonpolar",
  M="nonpolar", N="polar",   P="nonpolar", Q="polar",    R="basic",
  S="polar",    T="polar",   V="nonpolar", W="nonpolar", Y="polar",
  `*`="stop"
)

# ── Internal helpers ──────────────────────────────────────────────────────────

# Strip trailing stop from a protein string
.strip_stop <- function(s) sub("[*]$", "", s)

# Mean Chou-Fasman propensity for a protein string (ignoring unknowns/stops)
.cf_mean <- function(prot, table) {
  aas <- strsplit(.strip_stop(prot), "")[[1L]]
  vals <- table[aas[aas %in% names(table)]]
  if (length(vals) == 0L) return(NA_real_)
  round(mean(vals), 4L)
}

# Sliding-window Kyte-Doolittle; returns 1-based positions of TM-predicted
# residues (regions where mean KD > threshold over a window).
.tm_positions <- function(prot, window = 19L, threshold = 1.6) {
  aas  <- strsplit(.strip_stop(prot), "")[[1L]]
  n    <- length(aas)
  if (n < window) return(integer(0L))
  vals <- .KD_HYDRO[aas]
  vals[is.na(vals)] <- 0
  tm_pos <- integer(0L)
  for (i in seq_len(n - window + 1L)) {
    if (mean(vals[i:(i + window - 1L)]) >= threshold)
      tm_pos <- union(tm_pos, i:(i + window - 1L))
  }
  tm_pos
}

# BLOSUM62 score for a single substitution (returns NA for gaps/unknowns)
.blosum62 <- local({
  mat <- NULL
  function(aa1, aa2) {
    if (is.null(mat)) {
      e <- new.env(parent = emptyenv())
      tryCatch(
        utils::data("BLOSUM62", package = "Biostrings", envir = e),
        error = function(err) NULL
      )
      mat <<- if (exists("BLOSUM62", envir = e)) e$BLOSUM62 else NULL
    }
    if (is.null(mat) || aa1 == "-" || aa2 == "-") return(NA_real_)
    aa1 <- toupper(aa1); aa2 <- toupper(aa2)
    if (!(aa1 %in% rownames(mat)) || !(aa2 %in% colnames(mat))) return(NA_real_)
    mat[aa1, aa2]
  }
})

# Align two protein strings globally; return data.frame of changed positions
# (1-based in ref protein, skipping gaps in ref).
.protein_changes <- function(prot_ref, prot_alt) {
  p1 <- .strip_stop(prot_ref)
  p2 <- .strip_stop(prot_alt)
  if (p1 == p2) return(data.frame(pos=integer(), ref=character(), alt=character(),
                                   stringsAsFactors=FALSE))
  aln <- tryCatch(
    pwalign::pairwiseAlignment(
      Biostrings::AAString(p1), Biostrings::AAString(p2),
      substitutionMatrix = "BLOSUM62",
      gapOpening = 10, gapExtension = 0.5,
      type = "global"
    ),
    error = function(e) NULL
  )
  if (is.null(aln)) return(NULL)

  p  <- strsplit(as.character(pwalign::pattern(aln)),  "")[[1L]]
  s  <- strsplit(as.character(pwalign::subject(aln)), "")[[1L]]
  ref_pos <- 0L
  rows <- list()
  for (i in seq_along(p)) {
    if (p[i] != "-") ref_pos <- ref_pos + 1L
    if (p[i] == s[i]) next
    rows[[length(rows) + 1L]] <- data.frame(
      pos = ref_pos,
      ref = p[i],
      alt = s[i],
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0L)
    return(data.frame(pos=integer(), ref=character(), alt=character(),
                      stringsAsFactors=FALSE))
  do.call(rbind, rows)
}

# Build a one-row summary for one affected ORF
.orf_impact_row <- function(orf_type, gene, strand, frame, orf_start, orf_end,
                             prot_ref, prot_alt, snp_positions) {
  changes <- .protein_changes(prot_ref, prot_alt)

  if (is.null(changes)) {
    # alignment failed
    return(data.frame(
      orf_type=orf_type, gene=gene, strand=strand, frame=frame,
      orf_start=orf_start, orf_end=orf_end,
      orf_length_aa=nchar(.strip_stop(prot_ref)),
      n_snps=length(snp_positions),
      snp_positions=paste(snp_positions, collapse=","),
      protein_ref=.strip_stop(prot_ref),
      protein_alt=.strip_stop(prot_alt),
      n_aa_changes=NA_integer_,
      aa_changes=NA_character_,
      blosum_sum=NA_real_, blosum_mean=NA_real_,
      physchem_changes=NA_character_, n_charge_changes=NA_integer_,
      cf_helix_ref=NA_real_, cf_helix_alt=NA_real_, cf_helix_delta=NA_real_,
      cf_sheet_ref=NA_real_, cf_sheet_alt=NA_real_, cf_sheet_delta=NA_real_,
      tm_ref=NA_integer_, tm_alt=NA_integer_, tm_change_overlap=NA,
      stringsAsFactors=FALSE
    ))
  }

  p_ref_clean <- .strip_stop(prot_ref)
  p_alt_clean <- .strip_stop(prot_alt)

  ## AA changes string
  n_ch <- nrow(changes)
  aa_change_str <- if (n_ch == 0L) NA_character_ else
    paste(sprintf("%s%d%s", changes$ref, changes$pos, changes$alt), collapse="; ")

  ## BLOSUM62
  if (n_ch > 0L) {
    bl_scores <- vapply(seq_len(n_ch), function(j) {
      r <- changes$ref[j]; a <- changes$alt[j]
      if (r == "-" || a == "-") NA_real_ else .blosum62(r, a)
    }, numeric(1L))
    blosum_sum  <- round(sum(bl_scores, na.rm=TRUE), 2L)
    blosum_mean <- round(mean(bl_scores, na.rm=TRUE), 2L)
  } else {
    blosum_sum <- blosum_mean <- 0
  }

  ## Physicochemical changes
  physchem_str <- NA_character_
  n_charge_ch  <- 0L
  if (n_ch > 0L) {
    pc_events <- character(0L)
    for (j in seq_len(n_ch)) {
      r <- changes$ref[j]; a <- changes$alt[j]
      if (r == "-" || a == "-") next
      cl_r <- .PHYSCHEM[r]; cl_a <- .PHYSCHEM[a]
      if (!is.na(cl_r) && !is.na(cl_a) && cl_r != cl_a) {
        pc_events <- c(pc_events, sprintf("%s→%s@%d", cl_r, cl_a, changes$pos[j]))
        if ((cl_r == "acidic" || cl_r == "basic" || cl_a == "acidic" || cl_a == "basic") &&
            cl_r != cl_a)
          n_charge_ch <- n_charge_ch + 1L
      }
    }
    if (length(pc_events) > 0L) physchem_str <- paste(pc_events, collapse="; ")
  }

  ## Chou-Fasman
  cf_h_ref <- .cf_mean(p_ref_clean, .CF_HELIX)
  cf_h_alt <- .cf_mean(p_alt_clean, .CF_HELIX)
  cf_s_ref <- .cf_mean(p_ref_clean, .CF_SHEET)
  cf_s_alt <- .cf_mean(p_alt_clean, .CF_SHEET)

  ## TM prediction (only meaningful for proteins > 30 aa)
  len_ref <- nchar(p_ref_clean)
  tm_ref_pos <- if (len_ref >= 30L) .tm_positions(p_ref_clean) else integer(0L)
  tm_alt_pos <- if (nchar(p_alt_clean) >= 30L) .tm_positions(p_alt_clean) else integer(0L)
  tm_overlap  <- if (n_ch > 0L && length(tm_ref_pos) > 0L)
                   any(changes$pos[changes$ref != "-"] %in% tm_ref_pos)
                 else FALSE

  data.frame(
    orf_type          = orf_type,
    gene              = gene,
    strand            = strand,
    frame             = frame,
    orf_start         = as.integer(orf_start),
    orf_end           = as.integer(orf_end),
    orf_length_aa     = nchar(p_ref_clean),
    n_snps            = length(snp_positions),
    snp_positions     = paste(sort(unique(snp_positions)), collapse=","),
    protein_ref       = p_ref_clean,
    protein_alt       = p_alt_clean,
    n_aa_changes      = as.integer(n_ch),
    aa_changes        = aa_change_str,
    blosum_sum        = blosum_sum,
    blosum_mean       = blosum_mean,
    physchem_changes  = physchem_str,
    n_charge_changes  = as.integer(n_charge_ch),
    cf_helix_ref      = cf_h_ref,
    cf_helix_alt      = cf_h_alt,
    cf_helix_delta    = if (!is.na(cf_h_ref) && !is.na(cf_h_alt))
                          round(cf_h_alt - cf_h_ref, 4L) else NA_real_,
    cf_sheet_ref      = cf_s_ref,
    cf_sheet_alt      = cf_s_alt,
    cf_sheet_delta    = if (!is.na(cf_s_ref) && !is.na(cf_s_alt))
                          round(cf_s_alt - cf_s_ref, 4L) else NA_real_,
    tm_ref            = as.integer(length(tm_ref_pos) > 0L),
    tm_alt            = as.integer(length(tm_alt_pos) > 0L),
    tm_change_overlap = tm_overlap,
    stringsAsFactors  = FALSE
  )
}

# ── Main function ─────────────────────────────────────────────────────────────

#' Protein-Level Impact Table for a Mitogenome Comparison
#'
#' Takes the output of [compare_mitogenomes()] and produces a companion table
#' with one row per affected protein (canonical CDS or MDP ORF).  For each
#' protein, the function aligns the reference and alternate amino acid sequences
#' and computes amino-acid-level structural metrics.
#'
#' @param comparison A `mito_comparison` object returned by
#'   [compare_mitogenomes()].
#' @param canonical Logical.  Include canonical CDS proteins (default `TRUE`).
#' @param mdp Logical.  Include MDP ORFs (default `TRUE`).  Requires that
#'   `scan_mdps = TRUE` was used in [compare_mitogenomes()].
#'
#' @return A `data.frame` with one row per affected ORF and columns:
#' \describe{
#'   \item{`orf_type`}{`"CDS"` or `"MDP"`.}
#'   \item{`gene`}{Gene name for CDS; `NA` for MDP.}
#'   \item{`strand`, `frame`}{Orientation of the ORF.}
#'   \item{`orf_start`, `orf_end`}{Genomic coordinates in the reference.}
#'   \item{`orf_length_aa`}{Reference protein length (aa, excluding stop).}
#'   \item{`n_snps`, `snp_positions`}{Number and positions of contributing
#'     mutations from the comparison.}
#'   \item{`protein_ref`, `protein_alt`}{Amino acid sequences (stop stripped).}
#'   \item{`n_aa_changes`}{Number of differing residue positions.}
#'   \item{`aa_changes`}{Semicolon-separated list of substitutions, e.g.
#'     `"T23A; Q45K"`.}
#'   \item{`blosum_sum`}{Sum of BLOSUM62 scores across all substitutions
#'     (more negative = more structurally disruptive).}
#'   \item{`blosum_mean`}{Mean BLOSUM62 score per substituted position.}
#'   \item{`physchem_changes`}{Semicolon-separated physicochemical class
#'     changes, e.g. `"nonpolar→polar@23"`.}
#'   \item{`n_charge_changes`}{Number of positions where charge class changed
#'     (basic/acidic/neutral).}
#'   \item{`cf_helix_ref`, `cf_helix_alt`, `cf_helix_delta`}{Mean Chou-Fasman
#'     helix propensity for reference and alternate proteins, and their
#'     difference.}
#'   \item{`cf_sheet_ref`, `cf_sheet_alt`, `cf_sheet_delta`}{As above for
#'     beta-sheet propensity.}
#'   \item{`tm_ref`, `tm_alt`}{Whether a transmembrane helix is predicted
#'     (Kyte-Doolittle window 19, threshold 1.6) in reference/alternate.}
#'   \item{`tm_change_overlap`}{Logical: does any substitution fall within a
#'     predicted TM region in the reference protein?}
#' }
#'
#' @details
#' ## Metrics
#' **BLOSUM62** scores reflect evolutionary acceptability of each substitution;
#' large negative sums indicate radical changes unlikely to be structurally
#' tolerated.
#'
#' **Chou-Fasman propensities** (Chou & Fasman 1978) provide a residue-level
#' estimate of secondary-structure tendency.  `cf_helix_delta > 0` means the
#' alternate sequence has higher mean helix propensity.
#'
#' **Transmembrane prediction** uses a Kyte-Doolittle sliding window of 19
#' residues (threshold 1.6) and is only informative for proteins longer than
#' 30 aa.
#'
#' ## MDP matching
#' Reference MDP ORFs that overlap a mutation position are matched to the
#' closest alternate ORF by strand, frame, and maximum protein-sequence
#' identity.  ORFs with no match above 40% identity in the alternate genome
#' are reported as abolished (protein_alt = NA).
#'
#' @seealso [compare_mitogenomes()]
#'
#' @examples
#' \dontrun{
#' comp    <- compare_mitogenomes("NC_012920.1", "AP008907.1")
#' impacts <- analyze_protein_impacts(comp)
#' impacts[impacts$n_charge_changes > 0, ]
#' }
#'
#' @importFrom Biostrings AAString translate DNAString getGeneticCode
#' @importFrom pwalign pairwiseAlignment pattern subject pid
#' @importFrom utils data
#' @export
analyze_protein_impacts <- function(comparison,
                                     canonical = TRUE,
                                     mdp       = TRUE) {

  if (!inherits(comparison, "mito_comparison"))
    stop("'comparison' must be a mito_comparison object from compare_mitogenomes().",
         call. = FALSE)

  seq1_str  <- attr(comparison, "seq1_str")
  seq2_str  <- attr(comparison, "seq2_str")
  features  <- attr(comparison, "features")
  features2 <- attr(comparison, "features2")
  mdp_ref   <- attr(comparison, "mdp_ref")
  mdp_alt   <- attr(comparison, "mdp_alt")

  if (is.null(seq1_str) || is.null(seq2_str))
    stop("Comparison object is missing stored sequences. ",
         "Re-run compare_mitogenomes() with the current version of MitoDiveR.",
         call. = FALSE)

  rows <- list()

  ## =========================================================================
  ## 1. Canonical CDS impacts
  ## Build the alt protein by applying the AA changes already computed in the
  ## comparison table (avoids coordinate-shift issues from upstream indels).
  ## =========================================================================
  if (isTRUE(canonical) && !is.null(features)) {
    cds_muts <- comparison[grepl("^CDS:", comparison$region) &
                             !is.na(comparison$pos_ref), , drop = FALSE]
    affected_genes <- unique(cds_muts$canonical_gene)
    affected_genes <- affected_genes[!is.na(affected_genes)]

    for (gene in affected_genes) {
      feat_row <- features[features$feature == "CDS" &
                             features$gene == gene, , drop = FALSE]
      if (nrow(feat_row) == 0L) next

      ## Translate reference CDS
      prot_ref <- tryCatch({
        cds_seq <- .extract_cds_seq(seq1_str, feat_row)
        as.character(Biostrings::translate(
          Biostrings::DNAString(cds_seq),
          genetic.code   = Biostrings::getGeneticCode("SGC1"),
          if.fuzzy.codon = "solve"
        ))
      }, error = function(e) NULL)
      if (is.null(prot_ref)) next

      gene_muts <- cds_muts[!is.na(cds_muts$canonical_gene) &
                               cds_muts$canonical_gene == gene, , drop = FALSE]

      ## Primary: translate alt protein using the alt-genome annotation (features2).
      ## This correctly handles cases where the two genomes differ in gene order,
      ## strand orientation, or absolute gene position (e.g. family-level divergence).
      prot_alt <- tryCatch({
        feat2_row <- features2[features2$feature == "CDS" &
                                 features2$gene == gene, , drop = FALSE]
        if (nrow(feat2_row) == 0L) stop("no alt feature")
        cds_alt_seq <- .extract_cds_seq(seq2_str, feat2_row)
        as.character(Biostrings::translate(
          Biostrings::DNAString(cds_alt_seq),
          genetic.code   = Biostrings::getGeneticCode("SGC1"),
          if.fuzzy.codon = "solve"
        ))
      }, error = function(e) NULL)

      ## Fallback: apply missense/nonsense changes from the comparison table.
      ## Used when features2 is unavailable (e.g. sequences supplied directly).
      if (is.null(prot_alt)) {
        p_chars <- strsplit(.strip_stop(prot_ref), "")[[1L]]
        for (mi in seq_len(nrow(gene_muts))) {
          effect <- gene_muts$aa_effect[mi]
          aa_a   <- gene_muts$aa_alt[mi]
          pos_r  <- gene_muts$pos_ref[mi]
          if (is.na(effect) || is.na(pos_r)) next
          s      <- as.integer(feat_row$start[1L])
          e      <- as.integer(feat_row$end[1L])
          aa_pos <- if (feat_row$strand[1L] == "+") (pos_r - s) %/% 3L + 1L
                    else                             (e - pos_r) %/% 3L + 1L
          if (effect == "missense" && !is.na(aa_a) &&
              aa_pos >= 1L && aa_pos <= length(p_chars)) {
            p_chars[aa_pos] <- aa_a
          } else if (effect == "nonsense" &&
                     aa_pos >= 1L && aa_pos <= length(p_chars)) {
            p_chars <- p_chars[seq_len(aa_pos - 1L)]
          }
        }
        prot_alt <- paste(p_chars, collapse = "")
      }

      if (.strip_stop(prot_ref) == prot_alt) next   # only synonymous changes

      snp_pos <- gene_muts$pos_ref[!is.na(gene_muts$pos_ref)]

      rows[[length(rows) + 1L]] <- .orf_impact_row(
        orf_type      = "CDS",
        gene          = gene,
        strand        = feat_row$strand[1L],
        frame         = "+1",
        orf_start     = feat_row$start[1L],
        orf_end       = feat_row$end[1L],
        prot_ref      = prot_ref,
        prot_alt      = prot_alt,
        snp_positions = snp_pos
      )
    }
  }

  ## =========================================================================
  ## 2. MDP ORF impacts
  ## =========================================================================
  if (isTRUE(mdp)) {
    if (is.null(mdp_ref) || is.null(mdp_alt)) {
      warning("MDP ORF data not found in comparison object. ",
              "Re-run compare_mitogenomes() with scan_mdps = TRUE.",
              call. = FALSE)
    } else {

      ## Collect all ref ORFs that overlap any mutation position
      mut_pos <- comparison$pos_ref[!is.na(comparison$pos_ref)]

      overlap_idx <- which(vapply(seq_len(nrow(mdp_ref)), function(i) {
        any(mdp_ref$start[i] <= mut_pos & mdp_ref$end[i] >= mut_pos)
      }, logical(1L)))

      if (length(overlap_idx) > 0L) {
        ref_affected <- mdp_ref[overlap_idx, , drop = FALSE]

        ## Deduplicate by protein sequence + strand + frame
        orf_key <- paste0(ref_affected$strand, "_", ref_affected$frame, "_",
                          ref_affected$protein_sequence)
        ref_affected <- ref_affected[!duplicated(orf_key), , drop = FALSE]

        for (i in seq_len(nrow(ref_affected))) {
          r <- ref_affected[i, , drop = FALSE]

          ## Which mutations overlap this ORF?
          orf_muts <- comparison[!is.na(comparison$pos_ref) &
                                   comparison$pos_ref >= r$start &
                                   comparison$pos_ref <= r$end, , drop = FALSE]
          snp_pos  <- orf_muts$pos_ref

          ## Find best matching alt ORF by strand + frame + protein identity
          cand <- mdp_alt[mdp_alt$strand == r$strand &
                            mdp_alt$frame  == r$frame,  , drop = FALSE]
          prot_alt_seq <- NA_character_
          if (nrow(cand) > 0L) {
            ## Score each candidate by sequence identity to ref protein
            ref_p <- .strip_stop(r$protein_sequence)
            pid_vals <- vapply(seq_len(nrow(cand)), function(j) {
              alt_p <- .strip_stop(cand$protein_sequence[j])
              if (nchar(ref_p) == 0L || nchar(alt_p) == 0L) return(0)
              aln <- tryCatch(
                pwalign::pairwiseAlignment(
                  Biostrings::AAString(ref_p), Biostrings::AAString(alt_p),
                  substitutionMatrix="BLOSUM62",
                  gapOpening=10, gapExtension=0.5, type="global"),
                error = function(e) NULL)
              if (is.null(aln)) return(0)
              pwalign::pid(aln, type="PID2")
            }, numeric(1L))

            best_j <- which.max(pid_vals)
            if (pid_vals[best_j] >= 40) {
              prot_alt_seq <- .strip_stop(cand$protein_sequence[best_j])
            }
          }

          ref_prot_clean <- .strip_stop(r$protein_sequence)

          ## Skip if proteins are identical
          if (!is.na(prot_alt_seq) && ref_prot_clean == prot_alt_seq) next

          rows[[length(rows) + 1L]] <- .orf_impact_row(
            orf_type      = "MDP",
            gene          = NA_character_,
            strand        = r$strand,
            frame         = as.character(r$frame),
            orf_start     = r$start,
            orf_end       = r$end,
            prot_ref      = ref_prot_clean,
            prot_alt      = if (is.na(prot_alt_seq)) "" else prot_alt_seq,
            snp_positions = snp_pos
          )
        }
      }
    }
  }

  ## =========================================================================
  ## 3. Assemble and return
  ## =========================================================================
  if (length(rows) == 0L) {
    message("No protein-level impacts found.")
    return(invisible(data.frame(stringsAsFactors = FALSE)))
  }

  out <- do.call(rbind, rows)
  out <- out[order(out$orf_type, out$orf_start), ]
  rownames(out) <- NULL
  out
}
