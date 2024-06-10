# Microplastics-PBTK-Mice

June.10.2024

This repository contains all raw data and R code files for the project entitled "A Physiologically Based Toxicokinetic Model for Microplastics and Nanoplastics in Mice After Oral Exposure and Its Implications to Human Dietary Exposure Assessment". Datasets and Codes are provided. In 'Code' folder, there are four subfolders: Model Calibration, Model Evaluation, Model Optimization, and PCA.

Specifically,
1. Model calibration: 
	i. R file for 'MrgCode_PBTK model' introduces the PBTK model structure, parameters with initial values, and ODEs.
	ii. The 'Size-specific ModFit' folder includes four R files, each for one specific particle size. 
	iii. R file for 'NSC analysis' performs the normalized sensitivity analysis for AUCs to calibrated parameters.
2. Model evaluation: Four R files, each for one specific particle size, offer simulation results by directly applying calibrated models based on exposure scenarios of evaluation datasets.
3. Model optimization: There are four folders, each for one specific particle size. Each folder includes two R files: 'MrgCode-xx' for introducing PBTK models with size-specific calibrated parameters, and 'Model adaptation_xx' for optimizing PBTK models tailored to each specific evaluation dataset. 
4. PCA: R file for 'Correlation, PCA' performs the principal component analysis.

The manuscript describing this model is under review.
