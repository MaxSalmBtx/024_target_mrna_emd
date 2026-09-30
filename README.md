# 024_target_mrna_emd: Tumour mRNA Earth Mover's Distance (EMD) Pipeline

> **Quantitative characterisation and comparative distributional analysis of target mRNA expression across tumour profiles using Earth Mover's Distance and tensor decomposition.**

[![Pipeline Status](https://img.shields.io/badge/targets-reproducible-blue.svg)](#)
[![R Version](https://img.shields.io/badge/R-%3E%3D4.1.0-276DC3.svg)](#)
[![License](https://img.shields.io/badge/License-Proprietary-red.svg)](#)
[![Status](https://img.shields.io/badge/Status-Active-success.svg)](#)

---

## Value Proposition

The **024_target_mrna_emd** repository provides an end-to-end analytical workflow designed to profile, compare, and rank therapeutic targets by evaluating mRNA expression distribution across tumour models. Leveraging the Earth Mover's Distance (EMD / Wasserstein metric) alongside PARAFAC tensor decomposition, this pipeline quantifies multi-dimensional expression divergence across cell lineages and primary tumour profiles, equipping translational biology and discovery teams with robust statistical metrics for target prioritization.

---

## Key Features

* **Reproducible `{targets}` Orchestration:** Built upon the R `{targets}` ecosystem for automated dependency tracking, selective recomputation, and transparent caching of intermediate analytical stages.
* **Distributional Dissimilarity via EMD:** Quantifies subtle distributional shifts in mRNA expression beyond standard mean/median summaries by calculating the Earth Mover's Distance across tumour cohorts.
* **Multi-way PARAFAC Tensor Decomposition:** Decomposes higher-order target-by-tissue-by-condition tensors to identify latent expression signatures and underlying biological variations.
* **Configurable Target Evaluation:** Dynamically filters, stratifies, and scores targets defined via configuration manifests (`target_list.json`).
* **Publication-Quality Vector Artifacts:** Generates high-resolution visualization deliverables (`output/tumour_mrna_emd.svg`) seamlessly integrated with institutional visual assets.

---

## 🛠️ Tech Stack & Architecture

### Technology Overview

| Layer | Component / Package | Purpose |
| :--- | :--- | :--- |
| **Language** | R (>= 4.1.0) | Core data science and statistical computing runtime |
| **Workflow Management** | `{targets}`, `{tarchetypes}` | Directed acyclic graph (DAG) pipeline orchestration and cache management |
| **Data Manipulation** | `{tidyverse}`, `{data.table}`, `{jsonlite}` | High-throughput data wrangling and target configuration ingestion |
| **Visualization & Reporting** | `{ggplot2}`, `{svglite}` | Generation of vector graphics and distributional diagnostics |
| **Execution Scripting** | Bash / Windows Shell | Automated batch rendering and headless execution |

### Directory Breakdown

```text
024_target_mrna_emd/
├── code/
│   └── R/
│       ├── _targets.R                       # Central pipeline workflow definition and DAG configuration
│       ├── _targets/                        # Cache metadata, workspace snapshots, and pipeline objects
│       ├── 01a_install_packages.R           # Dependency bootstrapping and package installation script
│       ├── batch_render_windows.sh          # Batch script for headless execution on Windows environments
│       ├── target_list.json                 # Target list definitions and configuration parameters
│       ├── function_lib/                    # Modular analysis routines and helper functions
│       │   ├── target_blocks.R              # Pipeline target blocks, modular transforms, and DAG stages
│       │   └── tumour_mrna_deg_module_lib.R # distribution modeling & EMD distance routines
│       └── IMG/                             # Branding assets and visualization overlays
├── input/                                   # Source datasets, transcriptomic profiles, and metadata
└── output/                                  # Exported analysis artifacts, summaries, and SVG figures
    └── tumour_mrna_emd.svg                  # Rendered target mRNA EMD comparative figure
```

---

## Getting Started & Installation

### Prerequisites

* **R (>= 4.1.0)** installed on your machine.
* **Rtools** (for Windows users) configured on your system path to compile C/C++ dependencies.
* **Git** installed and authenticated.

### 1. Clone the Repository

```bash
git clone https://github.com/BicycleTx/024_target_mrna_emd.git
cd 024_target_mrna_emd
```

### 2. Environment Configuration & Data Preparation

```bash
# Verify directory structure
mkdir -p input output
```

Inspect and modify `code/R/target_list.json` to define the target genes of interest.

### 3. Dependency Installation

Launch an R session or execute the installation script from your terminal:

```bash
Rscript code/R/01a_install_packages.R
```

Alternatively, from within R:

```R
source("code/R/01a_install_packages.R")
```

---

## Usage Examples

### Running the Complete Workflow via R Console

To execute the computational pipeline with full dependency resolution:

```R
# Set working directory to code/R or repository root as appropriate
setwd("code/R")

library(targets)

# Inspect pending targets and dependencies
tar_manifest()

# Visualise the pipeline execution graph
tar_visnetwork()

# Run the complete pipeline
tar_make()
```

### Headless Batch Execution (Windows / Shell)

To run the pipeline via the automated batch script:

```bash
cd code/R
bash batch_render_windows.sh
```

### Inspecting Pipeline Outputs in R

Once the pipeline has completed, you can interactively inspect generated objects without recomputing:

```R
library(targets)
setwd("code/R")

# Read computed Earth Mover's Distance matrix
gene_emd_matrix <- tar_read(gene_emd_matrix)
head(gene_emd_matrix)

# Inspect PARAFAC model results
parafac_mdl <- tar_read(parafac_mdl)
summary(parafac_mdl)
```

Generated publication-grade figures will be exported to `output/`:

---

## Pipeline Validation & Diagnostics

To verify pipeline integrity and audit intermediate stages:

```R
library(targets)
setwd("code/R")

# Validate target definition syntax and dependencies
tar_validate()

# Check pipeline progress and outdated targets
tar_outdated()

# Audit compute performance and execution times
tar_meta(fields = c("name", "seconds", "bytes"))
```

---

## Contributing

Contributions, extensions, and target additions are welcome. Please adhere to the following workflow:

1. **Create a Feature Branch:**
   ```bash
   git checkout -b feature/target-expansion
   ```
2. **Follow Coding Standards:**
   - Adhere to the [Tidyverse style guide](https://style.tidyverse.org/).
   - Ensure all analysis steps are wrapped into modular functions within `code/R/function_lib/`.
   - Update `_targets.R` with deterministic targets and valid dependency references.
3. **Validate Pipeline:**
   - Run `targets::tar_make()` to confirm the pipeline builds cleanly without broken dependencies.
4. **Submit a Pull Request:** Provide a concise summary of new features, parameter additions, or dataset updates.

---

