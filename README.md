# 🧬 Single-Cell RNA-Seq Analysis of Mouse Brain Using Seurat & Monocle3

This repository contains a comprehensive downstream analysis pipeline for single-cell RNA-sequencing (scRNA-seq) data generated using 10X Genomics. The project focuses on clustering, visualization, trajectory inference, and cell type annotation of mouse brain cells.

---

## 📌 Objectives

- Perform quality control and filtering of single-cell data
- Normalize, scale, and reduce dimensionality of gene expression matrix
- Identify clusters and marker genes
- Reconstruct cellular trajectories using Monocle3
- Annotate cell types using reference-based methods
- Evaluate clustering quality using silhouette scores

---

## 📂 Dataset

- **Source**: [10X Genomics – Mouse Brain](https://www.10xgenomics.com/resources/datasets)
- **Format**: 3 files from Cell Ranger output:
  - `matrix.mtx.gz`
  - `features.tsv.gz`
  - `barcodes.tsv.gz`

---

## ⚙️ Tools & Packages

- [`Seurat`](https://satijalab.org/seurat/)
- `monocle3`
- `SeuratWrappers`
- `SingleR` + `celldex`
- `ggplot2`, `dplyr`, `Matrix`, `cluster`, `randomForest`

---

## 🚀 Workflow

1. **Quality Control**  
   - Mitochondrial % calculation  
   - Violin plots  
   - Filtering low-quality cells

2. **Normalization & Feature Selection**  
   - Log normalization  
   - Identification of highly variable genes  
   - Scaling of gene expression

3. **Dimensionality Reduction & Clustering**  
   - PCA + Elbow plot  
   - UMAP visualization  
   - Louvain clustering

4. **Marker Gene Identification**  
   - Differential expression using `FindAllMarkers`

5. **Trajectory Analysis**  
   - Seurat → Monocle3 conversion  
   - UMAP + trajectory graph  
   - Pseudotime inference

6. **Cell Type Annotation**  
   - Reference-based prediction using `SingleR` and `MouseRNAseqData`

7. **Clustering Validation**  
   - Silhouette analysis to evaluate clustering structure

---

## 🧠 Author

**Gunjan Sarode**  
MSc Biotechnology | Bioinformatics & Genomics | scRNA-Seq | Machine Learning  
🔗 🌐 [LinkedIn](https://www.linkedin.com/in/gunjan-sarode/) | 📫 gunjansarode.bioinfo@gmail.com

---

## 📜 License

This project is for educational & academic use. Please cite original data sources (10X Genomics) if reused.
