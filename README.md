# Palette
Robust integration of single-cell datasets with imbalanced modality composition

## Description

**Palette** is a computational framework for mosaic integration of single-cell multimodal data. It supports a wide range of integration scenarios, including (1) supervised mosaic integration when cell type annotations are available, (2) unsupervised mosaic integration, and (3) reference-based integration between an integrated reference and an unintegrated query dataset.

![Overview](https://github.com/qiongyusheng/Palette/blob/main/tutorials/Palette_overview.png)

## Installation
The Palette method is implemented in the R package `PaletteSC`.

### Option 1 — Install from GitHub

To use Palette, install the package directly from GitHub using `devtools`:

```R
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 600)

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("qiongyusheng/Palette")
```

### Option 2 — Install from GitHub Releases

```R
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_url(
  "https://github.com/qiongyusheng/Palette/releases/download/v0.1.0/PaletteSC_0.1.0.tar.gz"
)
```

### Option 3 — Install from a locally downloaded release tarball
1. Download the release file ([PaletteSC_0.1.0.tar.gz](https://github.com/qiongyusheng/Palette/releases/download/v0.1.0/PaletteSC_0.1.0.tar.gz)) to your local machine.

2. Install from the local file path:

```R
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_local(
  "/full/path/to/PaletteSC_0.1.0.tar.gz"
)
```
### Quick sanity check
After installation, verify that the package loads correctly:
```R
library(PaletteSC)
packageVersion("PaletteSC")
```
   

> Optional dependency: The reference-based integration module rely on additional packages.
> If needed, please install the corresponding optional dependencies:
> 
> ```R
> if (!requireNamespace("BiocManager", quietly = TRUE)) {
> install.packages("BiocManager")
> }
> options(repos = BiocManager::repositories())
> BiocManager::install("batchelor")
> install.packages("harmony")
> install.packages("Signac")
>  ```

## Quick Start

Here, we provide a simulated dataset to demonstrate the basic workflow of Palette. For more detailed explanations of different integration scenarios and advanced usage, please refer to the **`Tutorials`** section below. Here, users can download the [dataset](https://github.com/qiongyusheng/Palette/blob/main/tutorials/Demo.rds) and corresponding [metadata](https://github.com/qiongyusheng/Palette/blob/main/tutorials/Demo_meta.rds) and run the following code:

```R
library(PaletteSC)

# Data loading and preprocessing
obj <- readRDS('./Demo.rds') # replace the path to your real data path
obj <- Normalize.Data(obj,modal = "rna",normal.method = "LogNormalize")
obj <- Normalize.Data(obj,modal = "adt",normal.method = "CLR",margin = 2)
obj <- Add.HVFs(obj, modal = "rna")
obj <- Add.HVFs(obj, modal = "adt")

# Clustering and representative cells sampling
obj <- Find.Cluster(object = obj, 
                    modal = c("rna",'adt'),
                    method = c("PCA",'PCA'))

# Intra-modal joint dimensionality reduction using Bi-sPCA
obj <- Find.Subspace(obj,modal = c("rna",'adt'),
                     joint = TRUE,
                     sub.dims = list(1:20,1:15))

# MBG-guided inferring of missing modality matrices and Cross-batch alignment
obj <- Run.Palette(obj)

# Visualization
meta <- readRDS('./Demo_meta.rds') # replace the path to your real data path

library(uwot)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(scattermore)
library(ggsci)

Palette_emb <- as.matrix(obj@Int.result[["bind"]][rownames(meta),])
Palette_umap = as.data.frame(umap(Palette_emb))
colnames(Palette_umap) = c("UMAP1", "UMAP2")
Palette_umap = cbind.data.frame(meta, Palette_umap)
fig.size = function(height, width) {
  options(repr.plot.height = height, repr.plot.width = width, repr.plot.res = 300)
}
p1 <- ggplot(Palette_umap, aes(UMAP1, UMAP2, color = Batch)) +
  geom_scattermore(pointsize = 1) +
  scale_color_npg() +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme_classic()+
  theme(legend.position = "right")

p2 <- ggplot(Palette_umap, aes(UMAP1, UMAP2, color = Modality)) +
  geom_scattermore(pointsize = 1) +
  scale_color_npg() +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme_classic()+
  theme(legend.position = "right")


Colors <- colorRampPalette(brewer.pal(12, "Accent"))(16)
p3 <- ggplot(Palette_umap, aes(x = UMAP1, y = UMAP2, color = CellType)) +
  geom_scattermore(pointsize = 1) +   
  scale_color_manual(values = Colors) +  
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme_classic() +
  theme(legend.position = "right")


fig.size(4,15)
plot_grid(p1,p2,p3,align = 'h', axis = "b",nrow = 1,rel_widths = c(5,5,5))
```
![result](https://github.com/qiongyusheng/Palette/blob/main/tutorials/Demo_umap.png)

## Tutorials

The following tutorials provide step-by-step examples demonstrating how to apply Palette to a variety of single-cell integration scenarios. 

 - Example 1: [Unsupervised mosaic integration](https://github.com/qiongyusheng/Palette/blob/main/tutorials/Unsupervised%20mosaic%20integration%20using%20Palette.ipynb) | [Data](https://drive.google.com/drive/folders/1odz_MkWqfNY-MpzPYju6JGN1K-1AxtWj?usp=sharing)

- Example 2: [Supervised mosaic integration](https://github.com/qiongyusheng/Palette/blob/main/tutorials/Supervised%20mosaic%20integration%20using%20Palette.ipynb) | [Data](https://drive.google.com/drive/folders/1odz_MkWqfNY-MpzPYju6JGN1K-1AxtWj?usp=sharing)

- Example 3: [Reference-based integration](https://github.com/qiongyusheng/Palette/blob/main/tutorials/Reference-based%20integration%20using%20Palette.ipynb) | [Data](https://drive.google.com/drive/folders/1p8zRSE6hcQIM4giQRY6k3Ijw6wUUrxeD?usp=drive_link)

- Example 4: [Diagonal integration](https://github.com/qiongyusheng/Palette/blob/main/tutorials/Diagonal%20integration%20using%20Palette.ipynb) | [Data](https://drive.google.com/drive/folders/1Oy-PQRXn7BFfr1MQ6n7fL_XsQmn0_aSK?usp=sharing)

- Example 5: [Low-resolution spatial multimodal data integration](https://github.com/qiongyusheng/Palette/blob/main/tutorials/10x_Visium_Human_Tonsil.ipynb) | [Data](https://drive.google.com/drive/folders/1AD7rS5muCiq1cINcDs_lLP48PAjO3_Iv?usp=sharing)

## Dependencies
Palette has been successfully installed and tested on Windows, Linux, and macOS (R version ≥ 4.2.0). The main dependencies include:
- Rcpp (v1.0.13)
- RcppArmadillo (v14.0.2.1)
- methods (v4.2.3)
- Matrix (v1.6.5)
- RSpectra (v0.16.1)
- Seurat (v4.4.0)
- bigstatsr (v1.5.12)
- igraph (v2.0.3)
- irlba (v2.3.5.1)
- rlang (v1.1.6)
- MASS (v7.3.60)
- class (v7.3.22)
  
