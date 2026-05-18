# SynergyAnalysis

R framework for identifying **synergistic genes** from transcriptomics data — genes where the combination treatment **A+B** produces an effect that exceeds the sum of A alone and B alone.

Ships with a single-file Shiny web app for interactive exploration, ORA enrichment (GO BP / KEGG), and one-click Excel export including per-sample FPKM.

---

## Synergy criteria

For each gene, the framework merges 5 pairwise differential-expression comparisons (C vs NT, A vs NT, B vs NT, C vs A, C vs B) and applies one of two filters.

### Strict mode (4 criteria)

A gene is called **synergistic UP** when all four hold:

1. `C vs NT` significant (p or q < cutoff) and log2FC > 0
2. `C vs A`  significant and log2FC > 0
3. `C vs B`  significant and log2FC > 0
4. `Increase(C vs NT) > Increase(A vs NT) + Increase(B vs NT)`
   where `Increase = FC − 1` and `FC = 2^log2FC`

For **synergistic DOWN**, replace `> 0` with `< 0` and `Increase` with `Decrease = 1 − FC`.

### Relaxed mode

Only requires criterion 1 (C vs NT significant + direction) plus the magnitude-additivity criterion (4). Useful when A or B alone is already strongly significant, making C vs A / C vs B hard to detect at FDR thresholds even when supra-additivity is real.

---

## Installation

> **Tip for users in Mainland China** — set CRAN and Bioconductor to a domestic mirror (e.g. Westlake University) before installing to avoid slow / failing downloads. Run this **once per R session**, before the `install.packages()` / `BiocManager::install()` calls below:
>
> ```r
> options(repos       = c(CRAN = "https://mirrors.westlake.edu.cn/CRAN/"))
> options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")
> ```
>
> Other widely used mirrors: Tsinghua (`https://mirrors.tuna.tsinghua.edu.cn/CRAN/` + `https://mirrors.tuna.tsinghua.edu.cn/bioconductor`), USTC, etc. To make the setting permanent, add the two lines to `~/.Rprofile`.

### Option A — install with `devtools` (recommended)

```r
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Bioconductor dependencies for ORA (install these first)
BiocManager::install(c(
  "clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE"
))

# Install SynergyAnalysis from GitHub
devtools::install_github("puweilin/SynergyAnalysis")
```

Then load and launch:

```r
library(SynergyAnalysis)
run_synergy_app()        # opens the Shiny app in your browser
```

To install from a local clone instead:

```r
# After: git clone https://github.com/puweilin/SynergyAnalysis.git
devtools::install_local("SynergyAnalysis")
```

### Option B — source the scripts directly (no install)

```bash
git clone https://github.com/puweilin/SynergyAnalysis.git
cd SynergyAnalysis
```

```r
install.packages(c(
  "dplyr", "tidyr", "ggplot2", "ggrepel", "openxlsx",
  "pheatmap", "rmarkdown", "knitr", "stringr",
  "shiny", "bslib", "DT", "plotly"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c(
  "clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE"
))

# Then source the files directly
source("R/synergy_io.R")
source("R/synergy_core.R")
source("R/synergy_plot.R")
source("R/synergy_report.R")
```

---

## Quick start (command line)

After installing via devtools:

```r
library(SynergyAnalysis)

files <- list(
  c_vs_nt = "path/to/C_vs_Control.xls",
  a_vs_nt = "path/to/A_vs_Control.xls",
  b_vs_nt = "path/to/B_vs_Control.xls",
  c_vs_a  = "path/to/C_vs_A.xls",
  c_vs_b  = "path/to/C_vs_B.xls"
)

res <- calculate_synergy(
  files,
  p_cutoff   = 0.05,
  use_qvalue = FALSE,           # use raw P (TRUE = use FDR)
  mode       = "strict",        # or "relaxed"
  labels     = c(nt = "NT", a = "A", b = "B", c = "A+B")
)

print(res)
plot_synergy_volcano(res)
export_synergy_excel(res, "synergy_results.xlsx")
render_synergy_report(res, "synergy_report.html")
```

Without installing (source-only workflow), replace `library(SynergyAnalysis)` with the four `source()` calls shown above.

---

## Shiny app

```r
# If installed as a package:
SynergyAnalysis::run_synergy_app()

# Or from a clone, without installing:
shiny::runApp("shiny/", launch.browser = TRUE)
```

Tabs:

- **Data Input** — upload 5 files or specify local paths; configure column mapping
- **Results** — sortable table of synergy genes, switch UP/DOWN, P/Q-value columns follow the sidebar selection
- **Visualization** — interactive volcano (WebGL), top-gene effect contribution bars, per-gene inspector with Increase/Decrease vs A+B reference line and metric panel
- **Enrichment** — ORA via `clusterProfiler` on synergy UP/DOWN sets, GO BP or KEGG, p-value or p.adjust filter, downloadable .xlsx
- **Help** — criteria definitions, input format reference

---

## Input format

Tab-delimited file per pairwise comparison (compatible with edgeR `_all.xls` and standard DESeq2 output). Required columns:

| Column | Description |
|---|---|
| `gene_id` | Stable identifier (Ensembl, etc.) |
| `gene_name` | Display name; falls back to `gene_id` if missing |
| `log2FC` | log2 fold change |
| `Pvalue` | Raw p-value |
| `Qvalue` | FDR-adjusted p-value |
| `updown` | `"UP"` / `"DOWN"` / `"-"` (informational only; direction is taken from log2FC sign) |

If the file also contains `<sample>_FPKM` columns (as in edgeR `_all.xls`), per-sample FPKM is automatically extracted across all 5 files, deduplicated, and exported on the `FPKM_Synergy_UP` / `FPKM_Synergy_DOWN` sheets.

---

## Output

`export_synergy_excel()` produces a workbook with:

- `Summary` — counts and parameters
- `Synergy_UP` — UP synergy genes with log2FC, Q-values, FC, Increase, Increase_sum_ab
- `Synergy_DOWN` — same for DOWN with Decrease columns
- `FPKM_Synergy_UP` — per-sample FPKM matrix for UP genes
- `FPKM_Synergy_DOWN` — per-sample FPKM matrix for DOWN genes

`render_synergy_report()` produces a self-contained HTML report with the volcano, contribution bar, heatmap, and tables.

---

## License

MIT
