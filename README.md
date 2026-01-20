# Image-Based Process Monitoring: Robust Tensorial Data Approaches

This repository contains the code, data processing pipelines, and simulation results associated with the MSc thesis **“Image-Based Process Monitoring: Robust Tensorial Data Approaches.”** The project focuses on the development and evaluation of robust monitoring schemes for image-based industrial processes, leveraging tensor representations to efficiently handle high-dimensional and structured data.

The proposed methodology is designed for Phase I statistical process monitoring and aims to improve robustness against noise, outliers, and atypical observations commonly encountered in real-world image data. The repository includes both simulation-based studies and applications to real image datasets in order to assess the statistical and computational performance of the proposed approaches. Real imgages can be download at the MVTec research datasets webpage.

## Repository Overview

The repository is organized to facilitate reproducibility and extensibility. It includes implementations of tensor-based decompositions for process monitoring, scripts for data preprocessing and tensor construction, simulation frameworks for performance evaluation, and real-data case studies. The code structure clearly separates model estimation, monitoring statistics, and evaluation metrics, allowing for transparent analysis and easy methodological extensions.

## Reproducibility

All simulation studies and real-data analyses presented in the thesis can be fully reproduced using the scripts provided in this repository. The code is documented to facilitate adaptation to new datasets, alternative contamination schemes, or related monitoring frameworks.

## Requirements

The code is written in **R** (version 4.5.2). Some simulation studies included in this repository are computationally intensive, particularly those involving large image tensors, robust estimation procedures, and extensive Monte Carlo replications. As a result, full replication of certain experiments may require considerable execution time, depending on the available hardware and computational resources.

## Contact

For questions, suggestions, or collaboration opportunities, feel free to contact the author through GitHub or via academic email.


---

