This repository contains code supporting the single-cell RNA-sequencing analyses presented in:

# Integrated transcriptional and spatial profiling identifies CTHRC1+ fibroblasts as fibrogenic signaling hubs in human lung fibrosis.
Natalia Vanegas-Avendano<sup>1,†</sup>, Hu Chen<sup>1,2,3,†</sup>, Jack H. Wellmerling<sup>4</sup>, Mahnoor Hamid<sup>1</sup>, Ahmed Ghobashi<sup>1,2,3</sup>, Victor Peters<sup>1</sup>, Chandani Sen<sup>5</sup>, Brenda F. Reader<sup>6</sup>, Konstantin Shilo<sup>7</sup>, Brigitte N. Gomperts<sup>5</sup>, Qin Ma<sup>2,3</sup>, Ana L. Mora<sup>1</sup>, Daniel J. Tschumperlin<sup>4</sup>, Mauricio Rojas<sup>1*</sup>

<sup>1</sup>Division of Pulmonary, Critical Care, & Sleep Medicine, Department of Internal Medicine, The Ohio State University, Wexner Medical Center, Columbus, OH, USA. 
<sup>2</sup>Department of Biomedical Informatics, College of Medicine, The Ohio State University, Columbus, OH, USA. 
<sup>3</sup>Pelotonia Institute for Immuno-Oncology, The James Comprehensive Cancer Center, The Ohio State University, Columbus, OH, USA. 
<sup>4</sup>Department of Physiology & Biomedical Engineering. Mayo Clinic College of Medicine and Science, Rochester, MN, USA. 
<sup>5</sup>Children’s Discovery and Innovation Institute, Broad StemCell Research Center, Jonsson Comprehensive Cancer Center, Departments of Pediatrics and Pulmonary Medicine, University of California, Los Angeles, CA, USA 
<sup>6</sup>Human Tissue Biorepository. OSU Wexner Medical Center Comprehensive Transplant Center, The Ohio State University, Columbus, OH, USA 
<sup>7</sup>Department of Pathology, The Ohio State University Wexner Medical Center, Columbus, OH, USA. 
<sup>†</sup>These authors contributed equally   
*Corresponding author: Mauricio.Rojas@osumc.edu 

## Repository overview
Code for the fibroblast annotation validation and the adventitial vs. inflammatory fibroblast comparison added during peer review (Supplementary Fig. 1f–h).
The main analyses (QC, integration, annotation, DGE, CellChat, trajectories, WOT) are in https://github.com/Tigerrr07/IPF-mesenchymal-heterogeneity.

## Data
scRNA-seq data are available from GEO:
- Lung tissue (age-matched healthy controls, IPF upper and lower lobes): GSE283885, GSE331286, GSE331144

## Requirements
- Python 3.10.8, Scanpy 1.11.3 
- R 4.5.2, Seurat 5.5.1

## Files and figure panels
| File | Purpose | Figure panels |

| `Adventitial_inflammatory_fibroblasts_annotation_validation.ipynb` | Annotation validation, adventitial vs. inflammatory differential expression comparisons (Wilcoxon, BH correction, ≥20-point detection difference), marker refinement, visualize candidate-marker expression, module scores | Supplementary Fig. 1f, g|
| `adventitial_vs_inflammatory_markers.R` | Specificity of the top adventitial and inflammatory genes in whole lung data | Supplementary Fig. 1h |

## How to run
1. Download the processed objects from GEO.
2. Run `Adventitial_inflammatory_fibroblasts_annotation_validation.ipynb`.
3. Run `adventitial_vs_inflammatory_markers.R`.

