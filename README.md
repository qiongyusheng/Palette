# Palette
Robust integration of single-cell datasets with imbalanced modality composition

## Description

Palette is a computational framework for mosaic integration of single-cell multimodal data. Palette supports a wide range of integration scenarios, including mosaic integration with (1) supervised (when cell type annotations are available) and (2) unsupervised mode, and (3) reference-based integration between the integrated reference and unintegrated query dataset.

![Overview](https://github.com/qiongyusheng/Palette/blob/main/tutorials/Palette_overview.png)

## Installation
To run `Palette`, install from GitHub through ``devtools`` directly:

```R
install.packages('devtools')
library(devtools)
devtools::install_github("qiongyusheng/Palette")
```

## Dependencies
Palette has been successfully installed and used on Windows, Linux and Mac OS (R version >= 4.2.3). The dependencies including:
- Rcpp (v1.0.13)
- RcppArmadillo (v14.0.2.1)
- methods (v4.2.3)
- Matrix (v1.6.5)
- RSpectra (v0.16.1)
- Seurat (v4.4.0)
- bigstatsr (v1.5.12)
- igraph (v2.0.3)
- irlba (v2.3.5.1)
- preprocessCore (v1.60.2)
- rlang (v1.1.6)
- MASS (v7.3.60)
- Signac (v1.11.0)
- bachelor (v1.14.1)
- class (v7.3.22)
- harmony (v1.2.0)

## Tutorials

Example 1: [Unsupervised mosaic integration](https://github.com/qiongyusheng/Palette/blob/main/tutorials/Unsupervised%20mosaic%20integration%20using%20Palette.ipynb)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Data is available from [here.](https://drive.google.com/drive/folders/1odz_MkWqfNY-MpzPYju6JGN1K-1AxtWj?usp=drive_link)

Example 2: [Supervised mosaic integration](https://github.com/qiongyusheng/Palette/blob/main/tutorials/Supervised%20mosaic%20integration%20using%20Palette.ipynb)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Data is available from [here.](https://drive.google.com/drive/folders/1odz_MkWqfNY-MpzPYju6JGN1K-1AxtWj?usp=drive_link)

Example 3: [Reference-based integration](https://github.com/qiongyusheng/Palette/blob/main/tutorials/Reference-based%20integration%20using%20Palette.ipynb)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Data is available from [here.](https://cellxgene.cziscience.com/collections/9b02383a-9358-4f0f-9795-a891ec523bcc)

Example 4: [Diagonal integration](https://github.com/qiongyusheng/Palette/blob/main/tutorials/Diagonal%20integration%20using%20Palette.ipynb)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Data is available from [here.](https://drive.google.com/drive/folders/1p8zRSE6hcQIM4giQRY6k3Ijw6wUUrxeD?usp=drive_link)

Example 5: [Low-resolution spatial multimodal data integration](https://github.com/qiongyusheng/Palette/blob/main/tutorials/10x_Visium_Human_Tonsil.ipynb)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Data is available from [here.](https://drive.google.com/drive/folders/1AD7rS5muCiq1cINcDs_lLP48PAjO3_Iv?usp=drive_link)
