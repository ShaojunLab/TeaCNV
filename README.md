# TeaCNV: Tumor Epigenomic Absolute Copy Number Variation


## Description
TeaCNV is a computational framework designed to estimate tumor absolute copy number variations from single-cell chromatin accessibility data by aggregating the sparse signal of multi-individual cells to cell populations.
TeaCNV allows the inference of absolute copy number profiles from low coverage sparse data, including scATAC-seq and scATAC&RNA-seq co-assay data. 

![Overview](example/Figure1.png)

## System requirements and dependency
This framework runs on R (version > 4.2.0)

## Installing
#### Option A: Install TeaCNV within R using devtools
If installing from directly within R, you can instead use the following command from within R.
```
library("devtools")
devtools::install_github("ShaojunLab/TeaCNV",build_opts = c("--no-staged-install"))
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

### Input
TeaCNV takes a peak-by-cell count matrix and a data table of cell annotation as input. 

Input files (demo: `/example/`):
- `atac_count.RData`: raw scATAC count matrix

						  cell1  cell2  cell3   ......
		chr1-9906-10568    0      1      2    ......
		chr1-16010-16366   2      0      1    ......


- `cell_meta.csv`: cell type annotation (e.g., observed, reference)

				     Cluster
		cell1      observed
		cell2      observed
		cell3      reference
		......     ...... 
  

#### Run TeaCNV on the demo dataset

```
library(TeaCNV)
#Download the 'example' folder to the current  work path
load("./example/atac_count.RData")
mtx <- as.matrix(mtx)
cell_meta <- read.csv("./example/cell_meta.csv",row.name=1)
```
Once the input data are prepared, the next step is to initialize a TeaCNV object, which handles data preprocessing and quality control. After initialization, run the main TeaCNV pipeline to perform CNV and subclonal inference.

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

- `annotationFile`: annotation data.
- `ref_group_names`: one or more group labels in annotation data, and cells assigned to these label(s) are treated as the normal/reference population for baseline normalization in downstream CNV inference.
- `count_lim`: filters out peaks with abnormally high counts; it is recommended to set this threshold to the 99th percentile of peak counts.
- `seu_resolution`: controls clustering granularity in **Seurat**; higher values (e.g., > 1.2) produce more clusters, while the default is 1.0.
- `min_cells_in_group`: sets the minimum number of cells required per clone; the default is 20. 
- `delt_lim`: defines the relative copy number (CN) ratio interval corresponding to a one-copy change in absolute CN (default: 0.3). Increasing `delt_lim` results in a lower estimated clonal-level ploidy.
- `Correct_by_length`: when set to `TRUE`, normalizes peak counts to counts per kilobase to account for varying peak lengths; set it to `FALSE` if the input matrix features (bins) are of equal length.


### Output and Visualization

The final output of TeaCNV is saved in the 'final.CNVres.rds', which contains the complete analysis results.

If you wish to adjust parameters and rerun the analysis, you can simply reload the intermediate file 'TeaCNV.obj' as the `cnv_obj` object and continue the workflow from that point without repeating earlier steps.

All visualization outputs are stored in the `/Figure/` directory:
- 'heatmap_CNratio.pdf' — heatmap showing copy-ratio profiles with subclone annotations.

- 'heatmap_cloneCNV.pdf' — heatmap displaying the inferred integer copy number states for each subclone.

- 'clonalCN_final_noDots.pdf' — composite figure illustrating genome-wide segmentation and corresponding copy-ratio and CN profiles for each subclone.

For TeaCNV-related questions, refer to the README, FAQ, and publication. If something is unclear or missing, submit a Documentation Request or Feature Request on GitHub. For further assistance, contact us via email: Ying Wang [yingwang0727@outlook.com].


## SessionInfo
```r
R version 4.4.3 (2025-02-28)
Platform: x86_64-apple-darwin20
Running under: macOS Sequoia 15.4.1

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] TeaCNV_0.1.0    ggplot2_3.5.2        data.table_1.17.6    GenomicRanges_1.58.0 GenomeInfoDb_1.42.3 
[6] IRanges_2.40.1       S4Vectors_0.44.0     BiocGenerics_0.52.0 

loaded via a namespace (and not attached):
 [1] vctrs_0.6.5             httr_1.4.7              cli_3.6.5               rlang_1.1.6             UCSC.utils_1.2.0       
 [6] generics_0.1.4          jsonlite_2.0.0          glue_1.8.0              scales_1.4.0            grid_4.4.3             
[11] tibble_3.2.1            lifecycle_1.0.4         compiler_4.4.3          dplyr_1.1.4             RColorBrewer_1.1-3     
[16] pkgconfig_2.0.3         XVector_0.46.0          rstudioapi_0.17.1       farver_2.1.2            R6_2.6.1               
[21] tidyselect_1.2.1        dichromat_2.0-0.1       pillar_1.10.2           GenomeInfoDbData_1.2.13 magrittr_2.0.3         
[26] withr_3.0.2             tools_4.4.3             gtable_0.3.6            zlibbioc_1.52.0   


