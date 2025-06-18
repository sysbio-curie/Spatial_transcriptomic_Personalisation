# Spatial Transcriptomic Personalisation

Personalising Multiscale Models with Spatial Transcriptomics to Investigate T Cell Exclusion Mediated by Cancer-Associated Fibroblasts in Non-Small Cell Lung Cancer

This repository contains code and resources developed during a research project at Institut Curie that aimed to personalise multiscale models (MSMs) of tumour-immune interactions using spatial transcriptomic data. The focus was on understanding mechanisms of T cell exclusion mediated by cancer-associated fibroblasts (CAFs) in non-small cell lung cancer (NSCLC), particularly lung squamous cell carcinoma (LUSC).

## Project Overview

Spatial organisation within the tumour microenvironment (TME) significantly influences immune cell infiltration. CAFs are known to form barriers around tumour nests, excluding T cells and contributing to resistance to immunotherapy. This project explored CAF-mediated T cell exclusion using in silico multiscale simulations combining agent-based and Boolean intracellular models.

The MSMs were personalised using spatial transcriptomic (ST) data from LUSC samples, deconvoluted to infer cell types and gene expression at single-cell resolution.

## Key Components

    Multiscale Modelling (MSM): Built using the PhysiBoSS framework, integrating:

        Agent-based modelling with PhysiCell

        Boolean modelling with MaBoSS

    Spatial Transcriptomic Data:

        Visium v2 data from NSCLC patient samples (not avaialable on GitHub)

        Deconvolution using Cell2location and SpatialScope

    Model Personalisation Pipeline:

        Cell localisation and annotation via BioInformatics WalkThrough (BIWT)

## Repository Structure


## Getting Started

To reproduce or run the personalised simulations:

    Set up the environment:
    Use the provided Apptainer/Singularity definition files to ensure reproducibility across systems.

    Preprocess ST and scRNAseq data:

        Run deconvolution

    Personalise the model:

        Initialise agent-based models with inferred cell types and spatial locations

    Run simulations:

        Launch PhysiBoSS simulations from the model configurations

## Results Summary

    Cell2location and SpatialScope deconvolutions matched known histological structures, validated against immunohistochemistry (IHC).

    Personalised MSMs reproduced distinct tumour architectures and immune exclusion phenotypes.

    Simulations explored how CAF-derived ECM can mechanically hinder T cell infiltration.

## References

Key references include:

    Grout et al., Cancer Discovery 2022 – CAF states in NSCLC

    Sizek et al., Cell Cycle Model 2019 – Boolean modelling

    Sobkowicz A., Internship Report 2025

## Contact and Contributions

This work was developed by Agathe Sobkowicz during her Master’s internship at Institut Curie, under the supervision of Dr. Vincent Noël, Dr. Laurence Calzone, and Pr. Emmanuel Barillot.

Please open an issue or submit a pull request if you wish to contribute or have questions.
