# System Requirements  
## OS Requirements  
The package has been tested on the following systems:  
Linux：Linux version 3.10.0-1160.el7.x86_64    
R：4.2.3  
## Packages Version  
DoubletFinder: v2.0.3  
Seurat: v4.4.0  
SeuratWrappers: v0.3.1  
dplyr: v1.1.4  
ggplot2: v3.4.4  
RColorBrewer: v1.1.3  
ComplexHeatmap: v2.14.0  
scRNAtoolVis: v0.0.7  

# Installation Guide:
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')(v2.0.4)
remotes::install_github('https://github.com/ekernf01/DoubletFinder', force = T)(v2.0.3)  

**notes**:use this command line to install DoubletFinder(v2.0.3) which adapts to Seurat(v4).If you use Seurat(v5),
please install DoubletFinder(v2.0.4).

# Demo and Instrctions for use
to run demo data from Quality Control,you can start at the first step in code file:**Load endometrium dataset and QC**  
to plot figures in paper,you can start from the third step in code file:**Intergate data form three samples and cluster**

