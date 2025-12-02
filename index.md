# bccAnalysis

[![R-CMD-check](https://github.com/YOURUSERNAME/bccAnalysis/workflows/R-CMD-check/badge.svg)](https://github.com/Jacob1106/bccAnalysis/actions)
[![CRAN
status](https://www.r-pkg.org/badges/version/bccAnalysis)](https://CRAN.R-project.org/package=bccAnalysis)

**Interactive ChIP-Seq peak analysis with Shiny: overlap detection, peak
annotation, width comparison, and visualization.**

## ✨ Features

- **Peak overlap analysis** with custom Venn diagrams
- **Peak annotation** using ChIPseeker + TxDb.Hsapiens.UCSC.hg38
- **Width comparison** density plots
- **Export tables** (Excel) and plots (PDF)
- **Multi-sample filtering** by p-value and overlap count

## 🚀 Installation

Install: devtools::install_github(“YOURUSERNAME/bccAnalysis”)

Load package: library(bccAnalysis)

## 🎛️ Shiny App Interface

| Tab                   | Features                                                     |
|-----------------------|--------------------------------------------------------------|
| **Load Files**        | Upload .rds GRangesList, select samples, set overlap filters |
| **Venn Diagram**      | Visual peak overlap between conditions                       |
| **Table of Overlaps** | Downloadable proportions + peak counts                       |
| **Width Comparison**  | Log-scale density plot                                       |
| **Peak Analysis**     | Annotation barplots (main/unique/overlap peaks)              |

## 🔗 Dependencies

- Bioconductor: `GenomicRanges`, `rtracklayer`,
  `TxDb.Hsapiens.UCSC.hg38.knownGene`, `ChIPseeker`
- CRAN: `shiny`, `bslib`, `ggplot2`, `writexl`, `tibble`, `gridExtra`

## 📄 License

MIT © \[Jacob Martin\] 2025 \[memory:21\]

## 🙏 Acknowledgments

Built with using R package development best practices from [Hadley
Wickham](https://r-pkgs.org/).
