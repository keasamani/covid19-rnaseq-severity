# COVID-19 RNA-seq Severity Prediction

A reproducible bioinformatics pipeline predicting COVID-19 disease severity 
from bulk RNA-seq gene expression data using differential expression analysis 
and machine learning classification.

---

## Background

COVID-19 presents with a wide spectrum of clinical severity — from asymptomatic 
infection to critical illness requiring ICU admission. Identifying transcriptomic 
signatures that distinguish mild from severe disease could inform early clinical 
decision-making and reveal targetable host response pathways.

This project leverages publicly available bulk RNA-seq data from COVID-19 patients 
to build a reproducible pipeline that identifies differentially expressed genes and 
trains a machine learning classifier to predict disease severity.

---

## Objectives

- Identify differentially expressed genes between mild and severe COVID-19 patients
- Characterize host immune response pathways driving severity
- Build and validate a machine learning model predicting severity from gene expression
- Deliver a fully reproducible, documented pipeline

---

## Dataset

**GEO Accession:** GSE152641  
**Source:** Thair et al. (2021), *Cell Reports Medicine*  
**Samples:** ~100 whole blood RNA-seq samples (COVID-19 patients + healthy controls)  
**Severity labels:** Mild, Moderate, Severe, Critical  

---

## Methods Overview
