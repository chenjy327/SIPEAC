# Extended Analyses

| File | Analysis |
|---|---|
| 01_Spatial_pairing_validation.R | Repeated spatial events, random pairing and model validation |
| 02_RCTD_pairing.R | Single-cell reference, RCTD, singlet pairing and comparison |
| 03_Bin1_Bin50_pairing.R | Bin1 vs Bin50 pairing with a fixed ground truth |
| 04_Clonotype_normalization.R | Depth matching and repeated rarefaction |
| 05_Repair_comparison.R | SIPEAC vs repair comparison rebuilt from existing predictions (repair is not rerun) |
| 06_Sequencing_QC.R | Raw input preprocessing, coverage and gene saturation |
| 07_Spatial_BCR_GEX.R | Spatial concordance between BCR and GEX |
| 08_New_antibody_selection.R | Selection of new candidate antibodies and extraction of sequences for synthesis |
| 09_Bin50_pairing_examples.R | Bin50 pairing examples with and without AI |
| 10_Isotype_expansion.R | Isotypes and clonal expansion across tissues |
| 11_TCGA_antigen.R | TCGA TLS, tumor vs adjacent tissue |
| 12_Antibody_TLS.R | Proportion of antibody bins inside and outside TLS |
| 13_IFI30_distance.R | IFI30 expression around antibodies |
| 14_Antigen_neighborhood.R | Antigen and antibody neighborhoods |
| 15_Antibody_phylogeny.R | Antibody families and spatial lineage trees |
| 16_ssDNA_prepare.R | h5ad conversion and coordinate export |
| 17_ssDNA_pairing.R | SIPEAC with ssDNA vs scRNA reference |
| 18_ssDNA_reference.R | Direct ssDNA pairing vs scRNA reference |
| 19_ssDNA_comparison.R | Summary of ssDNA results across resolutions |
| 01_spatial_validation_functions.R … 09_bin50_example_functions.R | Functions used only by the script with the same number |
| 07_bcr_gex_concordance.R | BCR–GEX concordance metrics and block permutation test, used by 07 |
| utils_io.R | Project paths, reading and writing RDS files and tables |
| utils_plot.R | Plot themes, figure saving and spatial plot wrappers |
| utils_pairing.R | Ground truth construction and SIPEAC heavy/light chain pairing |
| utils_spatial_expand.R | Expanding spatial regions by n layers (`st_expand_n`) |
| utils_spatial_plot.py | Python spatial plotting called by utils_plot.R |
