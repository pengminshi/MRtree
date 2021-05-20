# mrtree 

R package for MRtree.

MRtree for Multi-resolution Reconciled tree, is a post clustering procedure to build cluster hierarchy based on flat clustering obtained for multi-resolutions. The tool is developed for single-cell RNAseq clustering analysis to identify sub-cell types. It can couple with many popular clustering tools such as Seurat. 


## Citation

```
@article{peng2020cfit,
  title={Cell Type Hierarchy Reconstruction via Reconciliation of Multi-resolution Cluster Tree},
  author={Minshi Peng, Brie Wamsley, Andrew G Elkins, Daniel M Geschwind, Yuting Wei and Kathryn Roeder},
  journal={bioRxiv},
  year={2021},
  doi={10.1101/430067}
}
```


## Installation
First, make sure to install all the dependent packages indicated in the `DESCRIPTION` file.
`mrtree` package can then be installed through `devtools` in R with the following commands,
```{r}
library("devtools")
devtools::install_github("pengminshi/mrtree")
```

`mrtree` provides wrappers for multiple clustering functions (see [vignette 1](https://htmlpreview.github.io/?https://github.com/pengminshi/MRtree/blob/master/vignettes/MRtree-tutorial.html)). The required packages for respective clustering methods are only required if the method is needed. This avoids the overloading dependencies.

## Examples
Please follow the [vignette 1](https://htmlpreview.github.io/?https://github.com/pengminshi/MRtree/blob/master/vignettes/MRtree-tutorial.html) for an example of using this R package with Seurat and SC3 on a simulated dataset, and [vignette 2](https://htmlpreview.github.io/?https://github.com/pengminshi/MRtree/blob/master/vignettes/zeisel_example.html) for the analysis on a mouse brain data.
