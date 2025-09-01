# Palette
Robust integration of single-cell datasets with imbalanced modality composition

Palette is a computational framework for mosaic integration of single-cell multimodal data. Palette supports a wide range of integration scenarios, including mosaic integration with (1) supervised (when cell type annotations are available) and (2) unsupervised mode, and (3) reference-based integration between the integrated reference and unintegrated query dataset.

![Overview](https://github.com/qiongyusheng/Palette/blob/main/tutorials/Palette_overview.png)

# Installation
To run `Palette`, install from GitHub through ``devtools`` directly:

```R
install.packages('devtools')
library(devtools)
devtools::install_github("qiongyusheng/Palette")
```

# Tutorials

1.[Unsupervised mosaic integration](https://github.com/qiongyusheng/Palette/blob/main/tutorials/Unsupervised%20mosaic%20integration%20using%20Palette.ipynb)

2.[Supervised mosaic integration](https://github.com/qiongyusheng/Palette/blob/main/tutorials/Supervised%20mosaic%20integration%20using%20Palette.ipynb)

3.[Reference-based integration](https://github.com/qiongyusheng/Palette/blob/main/tutorials/Reference-based%20integration%20using%20Palette.ipynb)

4.[Diagonal integration](https://github.com/qiongyusheng/Palette/blob/main/tutorials/Diagonal%20integration%20using%20Palette.ipynb)
