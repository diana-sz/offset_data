The script `get_T_sector.R` extracts protein mass fractions for the T sector of E. coli from published proteomics datasets and combines them with growth-rate data. The T sector consists of EF-Tu genes and aminoacyl-tRNA synthetase genes.
The script uses proteomics data from:
* Mori et al. (2021) [https://doi.org/10.15252/msb.20209536] - carbon-limited chemostat and minimal-media conditions
* Wu et al. (2023) [https://doi.org/10.1038/s41564-022-01310-w] - steady-state growth conditions
The resulting T-sector estimates are saved as phi_T_data.csv.

The script `get_Rb_data.R` combines ribosomal mass fraction data from multiple datasets, including
* Chure & Cremer 2023 [https://doi.org/10.7554/eLife.84878] - directly from publication + missing dataset for LB condition from Bremer 2008 added
* Mori et al. (2021) [https://doi.org/10.15252/msb.20209536] - calculated from proteomics data using a list of ribosomal protein genes
* Wu et al. (2023)  [https://doi.org/10.1038/s41564-022-01310-w] - calculated from proteomics data using a list of ribosomal protein genes
* Zhu (2025) [https://doi.org/10.1073/pnas.2427091122] - R/P ratio of E. coli taken from publication and multiplied by a factor of 0.5 to get phi_Rb
The resulting R-sector estimates are saved as phi_Rb_data.csv.