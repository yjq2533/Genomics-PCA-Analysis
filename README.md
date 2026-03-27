# Genomic PCA Analysis in R

## Overview
This project explores high-dimensional genomic data using Principal Component Analysis (PCA) to uncover patterns in population structure and genetic variation.

The dataset contains individuals from multiple populations, with nucleobase information across thousands of genomic positions. To make the data suitable for PCA, categorical nucleobase values are transformed into a binary mutation matrix.

## Approach

The analysis begins by preprocessing the genomic data and separating metadata such as individual ID, sex, and population group. The nucleobase data is then converted into a binary matrix, where each position indicates whether an individual's nucleobase differs from the most common value (mode) at that position.

PCA is applied to the resulting matrix without scaling, allowing the principal components to capture meaningful variation across individuals.



## Analysis and Results

The first two principal components are visualized to reveal clustering patterns across populations, showing clear structure in the data. 

Additional visualization using the first and third principal components highlights further variation, which can be interpreted using metadata such as sex.

To better understand the contribution of genomic positions, the absolute values of the third principal component loadings are plotted against nucleobase indices. This helps identify positions that have a stronger influence on the observed variation.



## Mathematical Validation

To reinforce understanding of PCA, the relationship between PCA, Singular Value Decomposition (SVD), and the eigen-decomposition of the covariance matrix is verified using a smaller synthetic dataset.



## Tools
- R
- data.table
- ggplot2



## Files
- `genomic_pca_analysis.R`: main analysis script

## Data
The dataset used in this project was provided as part of course materials and is not included in this repository.


## Timeline
Completed: Feburary 2026  
