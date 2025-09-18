# Strain-Analysis

This repository contains the methodology, calculations, and
visualizations for performing **geodetic strain analysis** on the
Iranian Plateau. The project uses GPS velocity observations and finite
difference methods to compute strain tensors, strain invariants, and
principal strains, ultimately producing maps of crustal deformation and
strain ellipses.

## Overview

The Iranian Plateau lies at the collision zone of the Arabian and
Eurasian tectonic plates. Continuous crustal deformation makes this
region an important case study for: - Understanding tectonic processes
- Assessing earthquake hazards
- Mapping active deformation patterns

This project quantifies crustal deformation by calculating strain
tensors on a regular grid, using GPS velocity data from 399 stations.

## Features

-   Computation of **strain tensors** from geodetic data\
-   Weighting scheme combining **distance** and **station clustering**
    (Voronoi-based)
-   Calculation of **strain invariants**, **principal strains**, and
    **rotation rates**
-   Visualization of:
    -   Normal strain maps
    -   Shear strain maps
    -   Strain invariants maps
    -   Principal strain maps
    -   Strain ellipses

## Methodology

1.  **Grid Definition**: 24°--43°N, 40°--64°E with 0.25° spacing
2.  **Finite Difference Method**: Solves for strain tensor components
    (εxx, εyy, εxy) and rotation (ω)
3.  **Weighting**:
    -   **Distance weighting** using a Gaussian function
    -   **Clustering weighting** using Voronoi diagrams
4.  **Strain Invariants & Ellipses**: Principal strains and ellipse
    orientations derived using eigenvalue analysis

## Results

-   Normal and shear strain distribution across the Iranian Plateau
-   Principal strain orientations and magnitudes
-   Strain ellipses visualizing crustal deformation
