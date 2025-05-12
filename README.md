# proj-sciii-diagnostic
This project aims to develop a data-driven approach to standardize the diagnosis of Class III skeletal malocclusion (SCIII). We use landmark coordinates annotated from lateral cephalometric radiographs to cluster patients based on craniofacial morphology.

We preprocess the data using geometric morphometric techniques, specifically Generalized Procrustes Analysis (GPA), and apply k-means clustering (with k=6) using a reduced set of 12 key anatomical landmarks.. Additionally, we explore other clustering algorithms and different values of k.

We provide a comprehensive evaluation of clustering robustness through a cross-validation framework, and implement a K-Nearest Neighbors (KNN) classifier to assign new patients to clusters; And interpretation of the obtained subphenotypes carried out in close collaboration with clinical experts in the field.

Further details on the methodology and findings are described in the associated publication: **LINK**

To access the datasets—comprising a cohort of 655 Class III patients of Caucasian origin (used for training) and an external dataset of 96 patients of Asian origin—please visit the following link: **LINK**

## Project overview

We start by  exploring our dataset of 655 SCIII (``step1_dataset_exploration.ipynb``) and further analyse the reliability of the annotated anatomical landmarks (``step2_annotation_reliability.ipynb``). Following a geometric morphometrics approach, we used Generalized Procrustes Analysis (``step3_generalized_procrustes_analysis.ipynb``)


## Directory structure
```
├── README.md          <- The top-level README for developers using this project.
│
├── docker             <- Dockerfile and VM configs to run the project
│
├── logs               <- Log files go here
│
├── notebooks          <- Jupyter notebooks.
|                       step1_dataset_exploration.ipynb
|                       step2_annotation_reliability.ipynb
|                       step3_generalized_procrustes_analysis.ipynb
│
├── outputs            <- Results from models/analysis go here, like figures, metrics, or
│			  predictions
│
├── config.yaml        <- config file with parameters for running the pipelines
│
├── pipelines          <- pipelines for the project
│
├── resources          <- other resources used by the project, such as SQL queries
│
├── scripts	       <- scripts developed that are related to this project 
│
└──src                <- Source code for use in this project.

```

## Prerequisites

## Install dependencies

## Usage

## Credentials

## Contact
