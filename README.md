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

### 1. Clone

```bash
git clone https://github.com/puweilin/SynergyAnalysis.git
cd SynergyAnalysis
```

### 2. R dependencies

Core analysis:

```r
install.packages(c(
  "dplyr", "tidyr", "ggplot2", "ggrepel", "openxlsx",
  "pheatmap", "rmarkdown", "knitr", "stringr"
))
```

Shiny app and enrichment:

```r
install.packages(c(
  "shiny", "bslib", "DT", "plotly"
))

# Bioconductor packages for ORA
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c(
  "clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE"
))
```

### 3. Sanity check

```r
source("R/synergy_io.R")
source("R/synergy_core.R")
# Should print no errors
```

---

## Quick start (command line)

```r
source("R/synergy_io.R")
source("R/synergy_core.R")
source("R/synergy_plot.R")
source("R/synergy_report.R")

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

---

## Shiny app

```r
# Launch in external browser
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
