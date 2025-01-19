# CellSpectra

Thank you for checking out our R package CellSpectra

Please find additional information in our three Tutorials with example datasets 
and R markdown files (/CellSpectra/tutorials/).

## 1. System requirements:

    dplyr,
    magrittr,
    stats,
    utils,
    Matrix,
    KEGGREST, 
    methods, 
    tidyr

    CellSpectra "on-the-fly" references depend on Seurat and Seuratobjects. 
    While we tested it with Seurat 5, this version needs a Seurat 4 object 
    structure (see tutorials for details).
    
    No non-standard hardware is needed. 
    
Versions the software has been tested on:
(for complete sessionInfo() check the tutorials!)

KEGGREST_1.38.0 CellSpectra_0.1.0 Seurat_5.0.3 SeuratObject_5.0.1
sp_2.1-3 dplyr_1.1.4 Matrix_1.6-5

## 2. Installation
(If not installed already, we recommend to first install the latest version of Seurat 
and SeuratObject)

download the .tar.gz file and change to the respective path
install.packages("/.../CellSpectra_0.1.0.tar.gz", repos = NULL, type = "source")

See the Installation Tutorial for more details!

## 3. Demo
See our two tutorials for demo data, output, and code (R markdown)

Basics: 
/CellSpectra/tutorials/Query_Reference_Seurat_onthefly

Advanced: 
/CellSpectra/tutorials/CellSpectra_advanced

## 4. Instructions 
Code used to generate the analysis from the manuscript: 
/CellSpectra/tutorials/Analysis_NG

or

https://github.com/susztaklab/SISKA

    
  

