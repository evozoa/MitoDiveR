# MitoDiveR

**Comparative mitochondrial genomics in R**

MitoDiveR is an R toolkit for comparative analysis, annotation parsing, and visualization of mitochondrial genomes across species. The package is designed to facilitate evolutionary and functional investigations of mitogenome structure, gene content, and sequence variation in a reproducible workflow.

---

## Overview

Mitochondrial genomes play central roles in metabolism, aging, and evolutionary diversification. Comparative analysis of mitogenomes across taxa can reveal patterns of gene retention, rearrangement, sequence divergence, and adaptation.

MitoDiveR provides tools to:

* Import and parse mitochondrial genome annotations
* Extract gene sequences and features
* Compare gene content and order across taxa
* Analyze sequence divergence and variation
* Prepare data for phylogenetic and evolutionary analyses
* Generate publication-ready visualizations

The package is designed for researchers in evolutionary biology, genomics, fisheries science, and comparative physiology.

---

## Installation

MitoDiveR is currently available from GitHub.

```r
install.packages("remotes")
remotes::install_github("YOUR-USERNAME/MitoDiveR")
```

Replace `YOUR-USERNAME` with the repository owner.

---

## Quick Start

### Load the package

```r
library(MitoDiveR)
```

### Import mitochondrial genomes

```r
mt_list <- mito_import_genbank("data/mitogenomes/")
```

### Extract gene annotations

```r
gene_table <- mito_extract_features(mt_list, feature = "CDS")
```

### Compare gene content across species

```r
content_summary <- mito_compare_gene_content(mt_list)
```

### Visualize gene order differences

```r
mito_plot_gene_order(mt_list)
```

---

## Key Functions (planned core API)

### Data Import

* `mito_import_genbank()` — Import annotated mitochondrial genomes
* `mito_import_fasta()` — Import sequences without annotations

### Annotation & Parsing

* `mito_extract_features()` — Extract genes, rRNAs, tRNAs, and control regions
* `mito_get_gene_sequences()` — Retrieve sequences by gene name

### Comparative Analysis

* `mito_compare_gene_content()` — Compare presence/absence across taxa
* `mito_compare_gene_order()` — Detect rearrangements
* `mito_sequence_divergence()` — Compute divergence metrics

### Visualization

* `mito_plot_gene_order()` — Visualize genome architecture
* `mito_plot_content_heatmap()` — Compare gene content
* `mito_plot_divergence()` — Display sequence differences

---

## Use Cases

MitoDiveR is particularly suited for:

* Comparative mitogenomics across clades
* Evolution of mitochondrial gene content
* Structural rearrangement analysis
* Phylogenetic dataset preparation
* Integrative studies of mitonuclear evolution

---

## Relationship to SNPmineR

MitoDiveR complements the SNPmineR package, which focuses on cross-species mapping of nuclear variants. Together, these tools enable integrated analyses of nuclear and mitochondrial genomic evolution.

---

## Development Status

MitoDiveR is under active development. Features and function names may change prior to the first stable release.

Planned milestones:

* v0.1 — Core import and parsing
* v0.5 — Comparative analyses and visualization
* v1.0 — Stable release

---

## Contributing

Contributions, feature requests, and bug reports are welcome via GitHub Issues and Pull Requests.

---

## Citation

If you use MitoDiveR in your research, please cite the package (citation details forthcoming).

---

## License

MIT License

---

## Contact

Maintained by Michael W. Sandel
Mississippi State University
Laboratory of Aquatic Evolution
