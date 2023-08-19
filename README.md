# STABox: Large-Scale Spatial Transcriptome Computational Analysis and Visualization Platform

----------------------------------------------------

## Introduction

----------------------------------------------------

### **STABox**, an integrated large-scale analysis and visualization framework for spatial transcriptomics data. It is a Python-based integrated development framework that provides interfaces for multiple high-performance algorithms, allowing for convenient integration and utilization in comprehensive downstream analysis of spatial transcriptomics data


![STABox Frame](/STABox.png)


## Installation

----------------------------------------------------

### STABox requires Python version >= 3.7 to run. Recommend version: 3.7.<br>The STABox package is developed based on the Python libraries Scanpy, PyTorch and PyG (PyTorch Geometric) framework, and can be run on GPU (recommend) or CPU.


### First clone the repository.

````
git clone https://github.com/Tidebear/STABOX.git
cd STABOX
````

### It's recommended to create a separate conda environment for running STABox:

````
#create an environment called STABox_env
conda create -n env_STAligner python=3.7

#activate your environment
conda activate STABox_env
````

### Install all the required packgers
````
pip install -r requiement.txt
````

### The use of the mclust algorithm requires the rpy2 package (Python) and the mclust package (R). See https://pypi.org/project/rpy2/ and https://cran.r-project.org/web/packages/mclust/index.html for detail.<p>The torch-geometric library is also required, please see the installation steps in https://github.com/pyg-team/pytorch_geometric#installation

### Running STABox
````
python VIEW/Interface.py
````



## Data collection

| Spatial transcriptome sequencing technology | Tissue type                                    |
|---------------------------------------------|------------------------------------------------|
| 10X Visium                                  | Brain, Breast                                  |
| Slide-seq V1                                | Brain, Hippocampus, Embroys                    |
| Slide-seq V2                                | Brain, Cerebellum, Hippocampus, Olfactory bulb |
| Stereo-seq                                  | Embroys, Olfactory bulb, Drosophila            |
| SeqFISH                                     | Brain                                          |
| MERFISH                                     | Embroys                                        |
| ST                                          | Heart                                          |
| STARmap                                     | Visual cortex                                  |
| HDST                                        | Brain                                          |


## Getting started

--------------------------------------------------------

### Default parameter Settings of Six algorithm

| Algorithm | parameter Settings                                                         |
|-----------|----------------------------------------------------------------------------|
| STABox    | (rad_cutoff:150-300,alpha:0-1,n_cluster: random int, Denoising: gene name) |
| STAligner | (rad_cutoff:150-300, alpha:0-3, n_cluster: random int, Margin: 0-2)        |
| STAMarker | (rad_cutoff:150-300, alpha:0-3, n_cluster: random int, SVGs: gene name)    |
| STAGE     | (alpha:0-1, Gene: gene name, Details: Random)                              |
| SpaGCN    | (weight:49, alpha: 1, n_cluster: random int, percentage: 0-1)              |
| SEDR      | (k_value:10, alpha:0-1, n_cluster: random int, p_drop:0-1)                 |
| SCANPY    | (random)                                                                   |


## Running!!!!

-----------------------------------
