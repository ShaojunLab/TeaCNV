# TeaCNV: Tumor Epigenomic Absolute Copy Number Variation


## Description
TeaCNV is a computational framework designed to estimate tumor absolute copy number variations from epigenomics data by aggregating the sparse signal of multi-individual cells to cell populations.
TeaCNV allows the inference of absolute copy number profiles from low coverage sparse data, including single-cell ATAC-seq and Multiome-seq data. 

## System requirements and dependency
This framework runs on R (version > 4.2.0)

## Installing TeaCNV
#### Option A: Install TeaCNV within R using devtools
If installing from directly within R, you can instead use the following command from within R.
```
library("devtools")
devtools::install_github("ShaojunLab/TeaCNV")
```
#### Option B: Install TeaCNV by pulling the code using git followed by source installation:
Alternatively, you can pull the code from github and install it:
```
git clone https://github.com/ShaojunLab/TeaCNV
cd TeaCNV
R
install.packages("./", repos=NULL, type="source")
```

## Usage

If you have installed TeaCNV, you can set the working path, load the sample data, and specify the reference normal cell annotation by setting the `ref_group_names` parameter, then run TeaCNV as follows:

```
library(infercnv)
setwd("./TeaCNV")
load("./example/atac_count.RData")
cell_meta <- read.csv("./example/cell_meta.csv",row.name=1)
ref_group_names <- c("Myeloid","Tcell","Bcell")
```

We use the `count_lim` parameter to limit abnormally high peak counts, and it is recommended to use the 99% quantile count value.
The `seu_resolution` parameter is the `resolution` to cluster cells in **Seurat**, use a value above 1.2 if you want to obtain a larger number of clusters (default 1.0).
We set the minimum size of each clone through the `min_cells_in_group` parameter (default value is 20). `delt_lim` represents the relative CN ratio spacing of one absolute copy number change (default value is 0.3). The larger the delt_lim, the lower the clonal-level ploidy estimate.
```
cnv_obj <- CreateTeaCNVObject(input = mtx,
                             annotationFile = cell_meta,
                             ref_group_names = ref_group_names,
                             ChrRemove = c('chrX', 'chrY', 'chrM'),
                             genome = "hg38",
                             count_lim = 4
                            )
res <- runTeaCNV(
	        input_obj = cnv_obj,
	        outdir = "./example",
	        delt_lim = 0.3,
		min_cells_in_group = 10,
	        seu_resolution = 1.2
	    )
```

