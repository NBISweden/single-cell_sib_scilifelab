Single-cell trajectories
================

## Theory session

[Slides](https://docs.google.com/presentation/d/1t_0yD7DxsMTK3fJPngNm9RN2CTgNJeeoBnxkd7c1mMc)

## Practical session

Download the data:
<https://drive.google.com/file/d/1eYDUry72RmQTz1t56DpXVUgkd7oEDaoq>

Unzip and put into the ‘data’ folder

① Quantification (Including RNA velocity):
[1_quantify.Rmd](session-trajectories/1_quantify.Rmd) ([html](session-trajectories/1_quantify.md))

② Preprocessing including RNA velocity):
[2_preprocessing.Rmd](session-trajectories/2_preprocessing.Rmd) ([html](session-trajectories/2_preprocessing.md))

③ Trajectory analysis using slingshot: [3_slingshot.Rmd](session-trajectories/3_slingshot.Rmd) ([html](session-trajectories/3_slingshot.md))

④ Velocity analysis using velocyto.R: [4_velocity.Rmd](session-trajectories/4_velocity.Rmd) ([html](session-trajectories/4_velocity.md))

We make use of several packages which should have already been installed
using the `install_packages.R` script. Several of these packages are
located on github only (but most will be available in the next
Bioconductor release). These packages are:

  - BUStools/**BUSpaRse**
  - velocyto-team/**velocyto.R**
  - satijalab/**seurat-wrappers**
  - dynverse/**tradeSeq**
