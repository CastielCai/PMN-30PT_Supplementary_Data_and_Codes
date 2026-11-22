# Supplementary Data and Codes

**Three-dimensional gradients of polar order in a relaxor ferroelectric**

Chen Cai<sup>1,2†</sup>, Pan Chen<sup>1†</sup>, Yuxuan Liao<sup>3†</sup>, Hongbo Han<sup>1,2†</sup>, Yuan-Jinsheng Liu<sup>4</sup>, Jianhua Ren<sup>4</sup>, Minh Pham<sup>3</sup>, Li Zhong<sup>5</sup>, Boyu Li<sup>5</sup>, Xiangying Liu<sup>5</sup>, Zitao Chen<sup>1</sup>, Aomiao Zhi<sup>1,2</sup>, Lifen Wang<sup>1</sup>, Zhi Xu<sup>1</sup>, Weijin Chen<sup>5</sup>, Peng Gao<sup>6</sup>, Shi Liu<sup>4&#42;</sup>, Xuezeng Tian<sup>1,5&#42;</sup>, Jianwei Miao<sup>3&#42;</sup>, and Xuedong Bai<sup>1,2&#42;</sup>

<sup>1</sup> Beijing National Laboratory for Condensed Matter Physics, Institute of Physics, Chinese Academy of Sciences, Beijing 100190, China.<br>
<sup>2</sup> School of Physics, University of Chinese Academy of Sciences, Beijing 100049, China.<br>
<sup>3</sup> Department of Physics & Astronomy and California NanoSystems Institute, University of California, Los Angeles, CA 90095, USA.<br>
<sup>4</sup> Westlake University, Hangzhou 310000, China.<br>
<sup>5</sup> Guangdong Provincial Key Laboratory of Magnetoelectric Physics and Devices, School of Physics, Sun Yat-sen University, Guangzhou 510275, China.<br>
<sup>6</sup> Electron Microscopy Laboratory and International Center for Quantum Materials, School of Physics, Peking University, Beijing 100871, China.

<sup>†</sup> These authors contributed equally to this work.<br>
<sup>&#42;</sup> Corresponding authors are Shi Liu (<liushi@westlake.edu.cn>), Xuezeng Tian (<tianxuezeng@mail.sysu.edu.cn>), Jianwei Miao (<j.miao@ucla.edu>), and Xuedong Bai (<xdbai@iphy.ac.cn>).

## Contents

- [Overview](#overview)
- [System requirements](#system-requirements)
- [Repository contents](#repository-contents)

## Overview

We used atomic electron tomography to determine the three-dimensional positions of individual cations in a lead magnesium niobate–lead titanate (PMN–0.3PT) single crystal. For each of 36,979 unit cells, we estimated local polarization from the displacement of the B-site cation relative to its eight Pb neighbors. Across the 60 nanodomains, polar order decreases on average from core to periphery. In 220 of 234 neighboring nanodomain pairs, the polarization order parameter (POP) decreases and then recovers along the path connecting the two cores. These observations form the basis of the polar order gradient model described in the manuscript.

This repository contains the experimental tilt series, reconstructed volume, atomic coordinates, multislice validation data, and MATLAB code for analyzing the polar displacement field. The original 2025 data and code release is archived on [Zenodo](https://doi.org/10.5281/zenodo.17661027).

## System requirements

- **MATLAB.** The analyses in folders 08–11 were tested with MATLAB R2023b on Windows. They use Statistics and Machine Learning Toolbox, Optimization Toolbox, Curve Fitting Toolbox, and Signal Processing Toolbox. Image Processing Toolbox is used for image morphology in folder 05 and for optional seed detection in folder 08.
- **Hardware.** The RESIRE reconstruction and projection-angle refinement routines in folders 02, 03, and 07 use CUDA MEX files and require a compatible NVIDIA GPU. The analyses in folders 08–11 run on a CPU using the deposited atomic model.
- **Data download.** Clone the repository with Git Large File Storage (Git LFS) enabled to retrieve the MATLAB data files. In an existing clone, run `git lfs pull` from the repository root.

## Repository contents

To reproduce the polar-order analyses from the coordinates in folder 06, follow steps 08–11 in order. In MATLAB, set the current folder to the directory listed for each step.

### [01. Experimental tilt series](./01_Measured_data/)

`proj_PMN30PT.mat` contains the 71 processed annular dark-field STEM projections. `angle_PMN30PT.mat` contains the corresponding projection angles.

### [02. RESIRE reconstruction](./02_RESIRE_package/)

`Main_RESIREcuda.m` reads the projections and angles from folder 01 and applies the supplied defocus parameters to reconstruct the experimental volume.

### [03. Projection-angle refinement](./03_Refinement_package/)

`RESIRE_multiGPU_AngleRefine_XYZTilt.m` refines the projection angles against a reconstruction. It calls `RESIRE_multiGPU_AngleRefine_GetOpAng.m` to evaluate candidate angles by comparison with the measured projections.

### [04. Reconstructed experimental volume](./04_final_reconstructoin_volume/)

`PMN30PT_Recon_Volume.mat` contains the three-dimensional volume used for atom tracing.

### [05. Atom tracing and classification](./05_Tracing_and_classification/)

The three `Modeling_step` scripts form an interactive sequence. Run steps 00, 01, and 02 in order to locate intensity maxima, trace the Pb sublattice, and locate the B-site cations.

### [06. Experimental atomic coordinates](./06_final_coordinates/)

`Model_output_02_Final_model.mat` contains the final coordinates of 41,969 Pb atoms and 36,979 B-site cations, their lattice indices, and the eight Pb neighbors of each B site. The B-site row order is retained throughout the polar-field analysis.

### [07. Multislice validation](./07_multislice_AET_simulation/)

Three overlapping regions are supplied with simulated tilt series, reconstructions, and traced coordinates. In `3_results`, `validation_of_tracing_result_whole_model.m` compares the recovered coordinates and polar displacements with the experimental model. The full-model comparison uses the supplied regional data. Reproducing the accepted-trace subset for the original panels also requires `Rind0208.mat`, a selection mask absent from the original deposit.

### 08. Polar displacement, polar order, and nanodomains

**Folder** [08_polar&domain&POP](./08_polar%26domain%26POP/)

Run [main.m](./08_polar%26domain%26POP/main.m) to calculate the B-site displacements relative to the centroids of their eight Pb neighbors, apply Gaussian averaging, and calculate the polarization order parameter (POP). The script saves the raw and averaged displacements in reconstruction pixels (0.3079 Å per pixel), POP, and the centers and mean polarization directions of the 60 cores in `output/polar_field.mat`. The reviewed domain and core labels are in `input/reviewed_domains.mat`, with zero marking unassigned sites. `grow_domains.m` implements the seed-based growth step that precedes review.

### [09. Stability analysis](./09_Stability_analysis/)

#### 09.1 Coordinate uncertainty

**Folder** [Coordinate_uncertainty](./09_Stability_analysis/Coordinate_uncertainty/)

Run [main.m](./09_Stability_analysis/Coordinate_uncertainty/main.m) to propagate 1,000 independent perturbations of the atomic coordinates through the displacement and averaging calculations (fig. S10). The supplied random streams reproduce the individual trials. The script saves the per-trial RMS deviations and the component-wise comparison for the first trial in `output/coordinate_uncertainty.mat`.

#### 09.2 Gaussian kernel sensitivity

**Folder** [Gaussian_kernel_sensitivity](./09_Stability_analysis/Gaussian_kernel_sensitivity/)

Run [main.m](./09_Stability_analysis/Gaussian_kernel_sensitivity/main.m) to compare Gaussian averaging kernels 1.20 and 2.66 nm in diameter (Fig. 2D and fig. S12). The radial correlation functions, fitted correlation lengths, and POP for the wider kernel are saved in `output/gaussian_comparison.mat`.

### [10. Polar order and spatial correlations](./10_Polar_order_and_spatial_correlations/)

#### 10.1 Polar-vector correlations

**Folder** [Polar_vector_correlations](./10_Polar_order_and_spatial_correlations/Polar_vector_correlations/)

Run [main.m](./10_Polar_order_and_spatial_correlations/Polar_vector_correlations/main.m) to calculate the distribution of angles between polar displacement vectors as a function of their separation (Fig. 2C). The script reads the reference correlation length from step 09.2 and saves the angular distributions in `output/angle_distributions.mat`.

#### 10.2 Intradomain polar order

**Folder** [Intradomain_polar_order](./10_Polar_order_and_spatial_correlations/Intradomain_polar_order/)

Run [main.m](./10_Polar_order_and_spatial_correlations/Intradomain_polar_order/main.m) to calculate domain sizes and radial POP profiles. The script saves both in `output/domain_order.mat`. For the equivalent spherical diameters in Fig. 2E, the code fits the experimental B-site positions to the lattice. The resulting unit-cell volume is 67.8477 Å³, and the calculation uses the unrounded value.

#### 10.3 Domain adjacency

**Folder** [Domain_adjacency](./10_Polar_order_and_spatial_correlations/Domain_adjacency/)

Run [main.m](./10_Polar_order_and_spatial_correlations/Domain_adjacency/main.m) to identify neighboring nanodomains whose unit cells share at least one face. The 234 pairs in `output/adjacent_pairs.mat` set the row order for the interdomain analyses in folder 11.

### [11. Interdomain transitions](./11_Interdomain_transition_analysis/)

#### 11.1 Core-to-core profiles

**Folder** [Core_to_core_profiles](./11_Interdomain_transition_analysis/Core_to_core_profiles/)

Run [main.m](./11_Interdomain_transition_analysis/Core_to_core_profiles/main.m) to sample each core-to-core path in the pair list from step 10.3 using a flat-ended cylinder of radius 8 Å and 2 Å axial bins. Each matrix in `output/profiles.mat` has five columns giving distance from the starting core (Å), mean POP, mean polarization magnitude (µC/cm²), and mean angles to the starting and ending core directions (°). Empty bins retain their distance coordinate and have `NaN` measurements.

#### 11.2 Transition widths

**Folder** [Transition_widths](./11_Interdomain_transition_analysis/Transition_widths/)

Run [main.m](./11_Interdomain_transition_analysis/Transition_widths/main.m) to measure the transition widths from the profiles in step 11.1 (Fig. 4C–E and fig. S15). For POP and polarization magnitude, the script measures the width at half depth relative to a line joining the first and last occupied-bin means. The angular width spans 12.5% to 87.5% of the cumulative absolute change in the two fitted angular profiles. In `output/transition_widths.mat`, the scalar-method flags identify observed crossings (1), an accepted count-weighted single-valley fit (2), or an unresolved width (0). The optional scripts [plot_population.m](./11_Interdomain_transition_analysis/Transition_widths/plot_population.m) and [plot_example.m](./11_Interdomain_transition_analysis/Transition_widths/plot_example.m) display the Fig. 4C distribution and the fig. S15 example, respectively.

#### 11.3 Angular monotonicity

**Folder** [Angular_monotonicity](./11_Interdomain_transition_analysis/Angular_monotonicity/)

Run [main.m](./11_Interdomain_transition_analysis/Angular_monotonicity/main.m) to calculate angular monotonicity from the net changes and total variation of the two angular fits in step 11.2. In `output/angular_monotonicity.mat`, `in_figure` selects the 233 pairs used in Fig. 4F, excluding D11–D14. The optional script [plot_results.m](./11_Interdomain_transition_analysis/Angular_monotonicity/plot_results.m) displays Fig. 4F.
