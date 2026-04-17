[![DOI](https://img.shields.io/badge/DOI-10.1038%2Fs43856--026--01557--y-blue)](https://doi.org/10.1038/s43856-026-01557-y)
[![zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.19628260.svg)](https://doi.org/10.5281/zenodo.19628260)

# Geometric morphometrics based diagnostic model for Skeletal Class III patients 
This project aims to develop a data-driven approach to standardize the diagnosis of Class III skeletal malocclusion (SCIII). We use landmark coordinates annotated from lateral cephalometric radiographs to cluster patients based on craniofacial morphology.

We preprocess the data using geometric morphometric techniques, specifically Generalized Procrustes Analysis (GPA), and apply k-means clustering (with k=6) using a reduced set of 12 key anatomical landmarks.. Additionally, we explore other clustering algorithms and different values of k.

We provide a comprehensive evaluation of clustering robustness through a cross-validation framework, and implement a K-Nearest Neighbors (KNN) classifier to assign new patients to clusters; And interpretation of the obtained subphenotypes carried out in close collaboration with clinical experts in the field.

Further details on the methodology and findings are described in the associated publication: https://doi.org/10.1038/s43856-026-01557-y

To access the datasets—comprising a cohort of 655 Class III patients of White origin (used for training) and an external dataset of 186 patients of Korean origin and 85 White origin - please request access from the corresponding authors via the link below.

## Project overview

We start by  exploring our dataset of 655 SCIII (``step1_dataset_exploration.ipynb``) and further analyse the reliability of the annotated anatomical landmarks (``step2_annotation_reliability.ipynb``). Following a geometric morphometrics approach, we used Generalized Procrustes Analysis (``step3_generalized_procrustes_analysis.ipynb``). Finally, we present the results of the clustering analysis and introduce six subphenotypes of skeletal Class III malocclusion (``step4_clustering_analysis.ipynb``). Other notebooks also include the analyses conducted for algorithm selection (``step4.1_clustering_algorithms.ipynb``), as well as, robustness assessment and subphenotype assignment for new samples (``step4.2_cross_validation.ipynb``, ``step4.3_clustering_prediction_korean_population.ipynb``, ``step4.4_clustering_prediction_white_external_population.ipynb``).


## Directory structure
```
├── README.md          <- The top-level README for developers using this project.
│
│
├── notebooks          <- Jupyter notebooks.
|                       step1_dataset_exploration.ipynb
|                       step2_annotation_reliability.ipynb
|                       step3_generalized_procrustes_analysis.ipynb
|                       step4_clustering_analysis.ipynb
|                       step4.1_clustering_algorithms.ipynb
|                       step4.2_cross_validation.ipynb
|                       step4.3_clustering_prediction_korean_population.ipynb
|                       step4.4_clustering_prediction_white_external_population.ipynb    
│
├── outputs            <- Results from models/analysis go here, like figures, metrics, or predictions
│
├── config.yaml        <- config file with parameters for running the pipelines
│
└─── scripts	       <- scripts developed that are related to this project 


```
---

## Step-by-step Instructions

### 🔧 1. Install & Configure Environment

* Make sure you have **R** installed. Recommended version:
  `R version 4.4.3 (2025-02-28)`
* Install **Conda** (via [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/))
* Install **Jupyter** and the **R kernel**:

```r
Rscript -e "install.packages('IRkernel'); IRkernel::installspec()"
```

---

### 📂 2. Clone the Repository

Use the terminal to clone the GitHub repository:

```verbatim
git clone https://github.com/istars-fmul/proj-sciii-diagnostic.git
cd proj-sciii-diagnostic
```

---

### 📦 3. Set Up the Conda Environment

Move to the root of the project directory (where `environment.yml` is located), then run:

```verbatim
conda env create -f environment.yml
```

This will create a new Conda environment named:
**`r_env_sciii_diagnosis`**

---

### ▶️ 4. Activate the Environment

```verbatim
conda activate r_env_sciii_diagnosis
```

---

### 📁 5. Include the Dataset

To reproduce the results in this repository, a `data/` folder must be present in the project root. This folder should contain:

* **Cephalometric measurements dataset**
* **Landmark coordinate dataset**
* **Patient info dataset**

If need change the path in the `config.yml` file. 


## SCIII Diagnostic Tool Web Interface

To compute and obtain the subphenotype, you can use our [online application](https://tools.istars.pt/sciii/). Please ensure that your input Excel file containing the annotated landmarks is in the correct format.

## Citation 

```text
Faria-Teixeira, M.C., Carvalho, I.M.N., Dehesa-Santos, A. et al. Geometric morphometrics based diagnostic model for Skeletal Class III patients. Commun Med (2026). https://doi.org/10.1038/s43856-026-01557-y
```

## Contact

For questions, feedback, or collaboration, feel free to reach out:

**iStars Team:**

* Inês M. N. Carvalho: [ines.c@edu.ulisboa.pt](mailto:ines.c@edu.ulisboa.pt)
* João C. Guimarães: [joao.guimaraes@medicina.ulisboa.pt](mailto:joao.guimaraes@medicina.ulisboa.pt)

**Website:** [https://istars.pt/](https://istars.pt/)


**Clinical team:**

* Maria Cristina Faria-Teixeira: [cristina.vft@gmail.com](mailto:cristina.vft@gmail.com)
* Alejandro Iglesias-Linares: [aleigl01@ucm.es](mailto:aleigl01@ucm.es)


## License

This software is available under the [European Union Public Licence (EUPL) v1.2 or later](https://eupl.eu/1.2/en/). For proprietary use, commercial support, or alternative licensing terms, please contact the authors.

