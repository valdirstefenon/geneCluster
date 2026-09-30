# geneCluster

A tool for detecting resistance gene clusters and biosynthesis gene clusters from GFF3 annotations.

geneCluster is a Python application that parses a GFF3 annotation file, cross-references gene identifiers against resistance and pharmaceutical annotation tables, and identifies spatial clusters of functionally related genes along scaffolds. It detects resistance gene clusters, core-validated biosynthesis clusters, and tandem duplicate pairs, then produces CSV tables and publication-ready figures.

The tool combines multiple evidence sources:

- GFF3 gene coordinates and transcript-to-gene relationships
- Resistance annotation tables with class and PFAM information
- Pharmaceutical annotation tables with PFAM domain information
- Curated PFAM dictionaries for biosynthesis core and accessory domains
- Sliding-window clustering with configurable distance and intervening-gene limits

This multi-layered strategy ensures maximum recovery of biologically meaningful gene clusters, including those in non-model plant species.

## Table of Contents

- [Overview](#overview)
- [What geneCluster Does](#what-genecluster-does)
- [Pipeline Workflow](#pipeline-workflow)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Input Files](#input-files)
- [Output Files](#output-files)
- [Using the Graphical Interface](#using-the-graphical-interface)
- [Understanding the Results](#understanding-the-results)
- [Resistance Cluster Analysis](#resistance-cluster-analysis)
- [Biosynthesis Cluster Analysis](#biosynthesis-cluster-analysis)
- [Tandem Duplicate Detection](#tandem-duplicate-detection)
- [GFF3 Validation](#gff3-validation)
- [Troubleshooting](#troubleshooting)
- [FAQ](#faq)
- [License](#license)

## Overview

geneCluster identifies and characterizes gene clusters from annotated plant genomes using a multi-evidence strategy:

- GFF3 gene coordinate parsing
- Transcript-to-gene mapping from mRNA features
- Resistance annotation table integration
- Pharmaceutical annotation table integration
- PFAM-based biosynthesis core and accessory classification
- Sliding-window cluster detection with intervening-gene limits

This approach maximizes recovery of biologically meaningful clusters and works with any species for which a GFF3 annotation and functional annotation tables are available.

## What geneCluster Does

### Core Functions

**Resistance Gene Cluster Detection**

Detects spatial clusters of resistance genes along scaffolds using a sliding window. A cluster is reported when at least two resistance-flagged genes fall within the configured distance and the number of intervening genes stays below the configured limit. Resistance classes recognized include NLR, TNL, TNJ, TNJ_like, CC_NLR, NLR_Jacalin, TIR_NLR, RPW8, Jacalin, LRR_repeat, NB_ARC, TIR, and CC.

**Biosynthesis Cluster Detection**

Detects biosynthesis gene clusters with a two-tier validation logic. A cluster is reported when it contains at least N core genes (default 2), or at least one core gene plus enough accessory genes to reach the core-plus-accessory threshold (default 4). Core and accessory domains are drawn from curated PFAM dictionaries covering terpene synthases, chalcone synthases, methyltransferases, cytochrome P450s, glycosyltransferases, and related families.

**Tandem Duplicate Detection**

Identifies tandem duplicate pairs among resistance genes, defined as consecutive resistance-flagged genes whose intergenic distance falls below the configured tandem window (default 50 kb).

**GFF3 Validation**

A dedicated validation mode checks the structural integrity of a GFF3 file, reporting field-count errors, invalid coordinates, missing gene IDs, duplicate gene IDs, and features lacking a Parent attribute.

**Comprehensive Visualization**

Cluster map across scaffolds, cluster size distributions, and per-cluster summaries. Figures are exported in PDF, PNG, and SVG.

**Tabular Outputs**

Gene table, resistance clusters, biosynthesis clusters, tandem duplicates, and a plain-text summary report.

## Pipeline Workflow

```text
┌─────────────────────────────────────────────────────────────────────────────┐
│                              INPUT FILES                                    │
├──────────────────────┬──────────────────────┬───────────────────────────────┤
│      GFF3 file       │  Resistance table    │   Pharmaceutical table        │
│   gene coordinates   │  annoMining output   │   annoMining output           │
└──────────┬───────────┴──────────┬───────────┴───────────────┬───────────────┘
           │                      │                           │
           ▼                      ▼                           ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                         PARSE AND MERGE                                     │
│  • Extract gene features and coordinates                                    │
│  • Build transcript-to-gene map from mRNA features                          │
│  • Load resistance and pharma annotation tables                             │
│  • Normalize gene identifiers across sources                                │
└─────────────────────────────────────────────────────────────────────────────┘
           │
           ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                         GENE CLASSIFICATION                                 │
│  • Flag resistance genes by class and PFAM                                  │
│  • Flag biosynthesis core genes by curated PFAM                             │
│  • Flag biosynthesis accessory genes by curated PFAM                        │
└─────────────────────────────────────────────────────────────────────────────┘
           │
           ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                    RESISTANCE CLUSTER DETECTION                             │
│  • Sliding window along each scaffold                                       │
│  • Distance threshold and intervening-gene limit                            │
│  • Cluster summaries with classes and gene lists                            │
└─────────────────────────────────────────────────────────────────────────────┘
           │
           ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                    BIOSYNTHESIS CLUSTER DETECTION                           │
│  • Sliding window over biosynthesis-flagged genes                           │
│  • Core-validated filtering                                                 │
│  • Core and accessory domain summaries                                      │
└─────────────────────────────────────────────────────────────────────────────┘
           │
           ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                    TANDEM DUPLICATE DETECTION                               │
│  • Consecutive resistance genes within tandem window                        │
│  • Distance reporting per pair                                              │
└─────────────────────────────────────────────────────────────────────────────┘
           │
           ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                    VISUALIZATION AND REPORTING                              │
│  • Cluster map across scaffolds                                             │
│  • Cluster size distributions                                               │
│  • CSV tables and plain-text summary                                        │
└─────────────────────────────────────────────────────────────────────────────┘
```

## System Requirements

### Minimum Hardware

| Component | Recommendation |
|-----------|----------------|
| CPU | 2+ cores recommended |
| RAM | 4 GB minimum; 8 GB+ recommended for large genomes |
| Storage | 2 GB free space for input files and outputs |

### Software Dependencies

- **Operating System:** Linux, macOS, or Windows
- **Python:** 3.8 or higher
- **Conda:** for environment management (recommended)

### Required Python Libraries

| Library | Purpose |
|---------|---------|
| pandas | Tabular data handling |
| numpy | Numerical operations |
| matplotlib | Figure generation |
| tkinter | GUI (python3-tk) |

## Installation Guide

### Step 1 — Create a Conda Environment (recommended)

```bash
conda create -n genecluster python=3.10 -y
conda activate genecluster
```

### Step 2 — Install core dependencies

```bash
conda install -c conda-forge pandas numpy matplotlib openpyxl -y
conda install -c conda-forge tk -y
```

### Step 3 — Download the script

Save the script as `genecluster.py` in your working directory.

### Step 4 — Verify installation

```bash
python3 -c "import pandas, numpy, matplotlib, tkinter; print('All libraries installed')"
```

## Input Files

### Mandatory Inputs

| File | Format | Description | Source |
|------|--------|-------------|--------|
| GFF3 annotation | `.gff3` | Gene predictions with coordinates and IDs | BRAKER3 or similar |
| Resistance table | `.tsv`, `.csv`, `.xlsx` | Resistance gene annotations with classes and PFAMs | annoMining or similar |

### Optional but Recommended Inputs

| File | Format | Description | Source |
|------|--------|-------------|--------|
| Pharmaceutical table | `.tsv`, `.csv`, `.xlsx` | PFAM domain annotations for pharmaceutical potential | annoMining or similar |

### File Preparation

**For the GFF3**

The GFF3 should contain gene features with `ID=` attributes and mRNA or transcript features with `ID=` and `Parent=` attributes. Transcript-to-gene mapping is built automatically from these features.

**For the resistance table**

The table should contain a `Gene_ID` column. Recognized optional columns include `Resistance_Classes`, `Resistance_PFAMs`, and `Resistance_Score`. Multiple values within a cell are separated by semicolons.

**For the pharmaceutical table**

The table should contain a `Gene_ID` column and a `PFAM_Domains` column. Multiple PFAM identifiers within a cell are separated by semicolons.

## Output Files

### Core Reports

| File | Description |
|------|-------------|
| `gene_table.csv` | Complete gene table with coordinates, classes, PFAMs, and flags |
| `resistance_clusters.csv` | Resistance clusters with coordinates, span, gene count, and classes |
| `biosynthesis_clusters.csv` | Biosynthesis clusters with core and accessory domain summaries |
| `tandem_duplicates.csv` | Tandem duplicate pairs among resistance genes |
| `cluster_summary.txt` | Plain-text summary report with global and per-cluster statistics |

### Visualizations

| File | Description |
|------|-------------|
| `cluster_cluster_map.pdf` | Cluster map across scaffolds |
| `cluster_cluster_map.png` | Cluster map across scaffolds |
| `cluster_cluster_map.svg` | Cluster map across scaffolds |
| `cluster_cluster_sizes.pdf` | Histograms of resistance and biosynthesis cluster sizes |
| `cluster_cluster_sizes.png` | Histograms of resistance and biosynthesis cluster sizes |
| `cluster_cluster_sizes.svg` | Histograms of resistance and biosynthesis cluster sizes |

## Using the Graphical Interface

### Launching the GUI

```bash
conda activate genecluster
python3 genecluster.py
```

### Main Tabs

#### 1. Cluster Analysis Tab

Runs resistance cluster detection, biosynthesis cluster detection, and tandem duplicate detection.

Steps:

1. Select input files.
   - GFF3 file: click Browse and select your `.gff3` file.
   - Resistance table: click Browse and select your resistance annotation file.
   - Pharma table (optional): click Browse and select your pharmaceutical annotation file.
2. Configure cluster parameters.
   - Resistance window (kb): default 500.
   - Max intervening genes: default 10.
   - BGC window (kb): default 200.
   - BGC min core genes: default 2.
   - BGC core+accessory threshold: default 4.
   - Tandem window (kb): default 50.
3. Set output directory. Default: `cluster_analysis`.
4. Click **Run Cluster Analysis**. Progress appears in the log window.

#### 2. GFF Validation Tab

Validates the structural integrity of a GFF3 file.

Steps:

1. Select the GFF3 file.
2. Click **Validate**.
3. Review errors and warnings in the output panel.

#### 3. About Tab

Describes the steps performed, the biosynthesis validation logic, the core and accessory PFAM dictionaries, the default parameters, and the required libraries.

## Understanding the Results

### Resistance Cluster Parameters

A resistance cluster is reported when:

- At least two resistance-flagged genes fall within the resistance window (default 500 kb), measured from the end of the previous gene to the start of the next.
- The number of intervening genes between consecutive flagged genes does not exceed the configured limit (default 10).

### Biosynthesis Cluster Validation

A biosynthesis cluster is reported when:

- It contains at least the configured number of core genes (default 2), or
- It contains at least one core gene and enough accessory genes to reach the core-plus-accessory threshold (default 4).

Core PFAMs include terpene synthases, squalene synthase, prenyltransferase, chalcone synthase, flavonoid 3-hydroxylase, SAM methyltransferase, O-methyltransferase, and HMGL-like domains. Accessory PFAMs include cytochrome P450, UDP glucosyltransferase, aminotransferase, pyridoxal enzyme, AMP binding, peroxidase, glycosyl hydrolase, FAD binding, NAD binding, and enoyl reductase.

### Tandem Duplicate Detection

A tandem duplicate pair is reported when two consecutive resistance-flagged genes on the same scaffold have an intergenic distance of at most the tandem window (default 50 kb).

### Interpreting the Summary File

```text
GFF3 CLUSTER ANALYSIS SUMMARY
Total genes parsed: 29,698
Scaffolds: 412
Resistance genes: 1,247
Biosynthesis core genes: 186
Biosynthesis accessory genes: 421
Biosynthesis total genes: 512

RESISTANCE CLUSTERS
Total clusters: 48
Total genes in clusters: 312
Mean cluster size: 6.50
Largest cluster: 24

BIOSYNTHESIS CLUSTERS (CORE-VALIDATED)
Total clusters: 17
Total genes in clusters: 94
Mean cluster size: 5.53
Total core genes: 38
Total accessory genes: 56

TANDEM DUPLICATES (RESISTANCE)
Total pairs: 63
```

## Resistance Cluster Analysis

### Resistance Classes Detected

| Class | Domains | Description |
|-------|---------|-------------|
| NLR | NB-ARC + LRR | Classic NLR resistance genes |
| TIR_NLR | TIR + NB-ARC + LRR | TIR-domain containing NLRs |
| TNJ | TIR + NB-ARC + Jacalin | Myrtaceae-specific resistance genes |
| TNJ_like | TIR + NB-ARC + Jacalin + LRR | Extended TNJ architecture |
| CC_NLR | CC + NB-ARC + LRR | Coiled-coil NLRs |
| TNL | TIR + NB-ARC + LRR | TIR-NBS-LRR |
| RPW8 | RPW8 domains | RPW8-type resistance |
| NB_ARC | NB-ARC | Nucleotide-binding domain |
| TIR | TIR | Toll / Interleukin-1 receptor |
| LRR_repeat | LRR | Leucine-rich repeats |
| Jacalin | Jacalin | Jacalin lectin domain |
| NLR_Jacalin | NB-ARC + Jacalin | NLR with Jacalin, no TIR |

### Resistance PFAMs Recognized

| PFAM | Domain |
|------|--------|
| PF00931 | NB_ARC |
| PF01582 | TIR |
| PF13676 | TIR |
| PF00560 | LRR |
| PF07723 | LRR |
| PF13855 | LRR |
| PF01419 | Jacalin |
| PF05659 | RPW8 |
| PF05660 | RPW8 |

## Biosynthesis Cluster Analysis

### Core PFAMs

| PFAM | Domain |
|------|--------|
| PF01397 | Terpene_synthase |
| PF03936 | Terpene_synthase_N |
| PF13243 | Terpene_synthase_C |
| PF13249 | Terpene_synthase_C2 |
| PF08491 | Squalene_synthase |
| PF00494 | Prenyltransferase |
| PF00195 | Chalcone_synthase |
| PF02797 | Flavonoid_3_hydroxylase |
| PF01596 | SAM_methyltransferase |
| PF00891 | O_methyltransferase |
| PF00682 | HMGL_like |

### Accessory PFAMs

| PFAM | Domain |
|------|--------|
| PF00067 | Cytochrome_P450 |
| PF00201 | UDP_glucosyltransferase |
| PF00155 | Aminotransferase |
| PF00282 | Pyridoxal_enzyme |
| PF00501 | AMP_binding |
| PF00141 | Peroxidase |
| PF00933 | Glycosyl_hydrolase |
| PF01565 | FAD_binding_4 |
| PF13460 | NAD_binding |
| PF13561 | Enoyl_reductase |

## Tandem Duplicate Detection

Tandem duplicates among resistance genes are reported with the scaffold identifier, both gene identifiers, and the intergenic distance in kilobases. This output is useful for identifying recent local expansions of resistance gene families, a common pattern in plant NLR clusters.

## GFF3 Validation

The validation module checks:

- Each non-comment line has exactly nine tab-separated fields.
- Start and end coordinates are integers.
- Gene features carry an `ID=` attribute.
- Gene IDs are unique within the file.
- mRNA, CDS, exon, and transcript features carry a `Parent=` attribute.

Errors and warnings are reported with line numbers. Up to fifty of each are displayed in the GUI.

## Troubleshooting

| Issue | Solution |
|-------|----------|
| GUI does not start | Install tkinter: `conda install -c conda-forge tk` |
| No resistance clusters detected | Check that resistance classes or PFAMs are present in your annotation table and that gene IDs match the GFF3 |
| No biosynthesis clusters detected | Confirm that PFAM identifiers in the pharma table match the curated core and accessory dictionaries |
| Gene IDs do not match | The tool normalizes common suffixes such as `.t1`, `.mrna1`, and `.1`, and strips common prefixes. Check the log for the normalized intersection count |
| ValueError on a numeric parameter | Ensure all parameter fields contain numeric values |
| MemoryError | Reduce the number of genes analyzed or increase available RAM |
| Figures not generated | Ensure matplotlib is installed and the output directory is writable |

### Debug Mode

Run from the terminal to see verbose logging:

```bash
python3 genecluster.py 2>&1 | tee output.log
```

## FAQ

**Q: Why do I need both a GFF3 and annotation tables?**

The GFF3 provides gene coordinates and the physical layout along scaffolds, which is what makes cluster detection possible. The annotation tables provide the functional information, resistance classes and PFAM domains, that determines which genes are flagged.

**Q: What if my annotation table uses transcript IDs instead of gene IDs?**

The tool builds a transcript-to-gene map from the GFF3 mRNA features and normalizes common suffixes. As long as the transcript IDs in your annotation table match the mRNA IDs in the GFF3, the mapping is resolved automatically.

**Q: Can I run the tool without the pharmaceutical table?**

Yes. Resistance cluster detection and tandem duplicate detection work with the GFF3 and resistance table alone. Biosynthesis cluster detection requires the pharmaceutical table because it depends on PFAM domain annotations.

**Q: How do I tune cluster detection sensitivity?**

Increase the resistance window or the max intervening genes to detect larger, more dispersed clusters. Decrease them to detect tighter clusters. For biosynthesis clusters, adjust the BGC window, the minimum core genes, and the core-plus-accessory threshold.

**Q: What does the "span_kb" column mean?**

It is the distance in kilobases between the start of the first gene and the end of the last gene in the cluster, giving the physical footprint of the cluster on the scaffold.

**Q: Can I add my own PFAM domains to the biosynthesis dictionaries?**

Yes. Edit the `BIOSYNTHESIS_CORE_PFAMS` and `BIOSYNTHESIS_ACCESSORY_PFAMS` dictionaries at the top of the script. Adding a PFAM to the core dictionary raises the likelihood that a cluster containing it will pass validation.

**Q: Does the tool handle gzipped GFF3 files?**

Yes. Files ending in `.gz` are read transparently with gzip.

**Q: Can I run the analysis on multiple genomes?**

Run the tool once per genome with separate output directories, then compare the resulting CSV tables externally. There is no built-in multi-genome mode in this version.

## License

geneCluster is distributed under the MIT License.
