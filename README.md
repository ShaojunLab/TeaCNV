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
devtools::install_github("ShaojunLab/TeaCNV",INSTALL_opts = c("--no-staged-install"))
```
#### Option B: Install TeaCNV by pulling the code using git followed by source installation:
Alternatively, you can pull the code from github and install it:
```
git clone https://github.com/ShaojunLab/TeaCNV
cd TeaCNV
R
install.packages("./", repos=NULL, type="source",INSTALL_opts = c("--no-staged-install"))
```

## Usage

If TeaCNV is already installed, follow the steps below to set up the working directory. This step is essential if the user has not provided the necessary reference files, as TeaCNV's built-in reference data —defaulting to the HG38 genome version— is stored in the working directory. After setting up, you can load the sample data and run TeaCNV:

```
library(TeaCNV)
#Download the 'example' folder to the current  work path
load("./example/atac_count.RData")
mtx <- as.matrix(mtx)
cell_meta <- read.csv("./example/cell_meta.csv",row.name=1)
```
TeaCNV takes a peak-cell count matrix as input. A data table is provided via the `annotationFile`, where the first column contains annotations for the cells represented in the columns of the input matrix. These annotations typically distinguish between reference and observed cell groups. The `ref_group_names` parameter identifies the group name associated with the reference cells.

To filter out abnormally high peak counts, the `count_lim` parameter can be used, and it is recommended to set it to the 99th percentile of peak counts. 
The `seu_resolution` parameter controls the clustering resolution in **Seurat**; a value above 1.2 is suggested if a higher number of clusters is desired (default: 1.0).
The `min_cells_in_group` parameter sets the minimum required size for each clone, with a default value of 20 cells. `delt_lim` defines the relative copy number (CN) ratio interval for a single absolute CN change (default: 0.3). Increasing `delt_lim` results in a lower clonal-level ploidy estimate.

We set `Correct_by_length = TRUE` to normalize the counts of peaks with different lengths to counts per kilobase. If the input matrix features (bins) are of equal length, we set it to `FALSE`.
```

cnv_obj <- CreateTeaCNVObject(input = mtx,
                             annotationFile = cell_meta,
                             ref_group_names = "reference",
                             ChrRemove = c('chrX', 'chrY', 'chrM'),
                             genome = "hg38",
                             count_lim = 4,
			     Correct_by_length = TRUE
)
res <- runTeaCNV(input_obj = cnv_obj,
	        outdir = "./example",
	        delt_lim = 0.3,
		min_cells_in_group = 10,
	        seu_resolution = 1.2)
```

The output of TeaCNV is saved in the 'final.CNVres.rds' file. If you need to modify parameters and rerun the analysis, you can load the 'TeaCNV.obj' file as the `cnv_obj` object and continue from there.

For TeaCNV-related questions, refer to the README, FAQ, and publication. If something is unclear or missing, submit a Documentation Request or Feature Request on GitHub. For further assistance, contact us via email: Ying Wang [yingwang0727@outlook.com].


