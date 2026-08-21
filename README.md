# Analysis code for the transwell study

This repository contains the R scripts used for the targeted and untargeted metabolomics, cell-viability, TEER, and transwell analyses in this study. The analyses were performed with **R 4.2.2**.

## Download the input data

The data are stored as flat files in this [unpublished Zenodo preview](https://zenodo.org/records/22044327?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjBjOWNkZmU0LTQ0YWYtNGM1NC1iMmVjLWQ5ZjFkNWUyNTBmNSIsImRhdGEiOnt9LCJyYW5kb20iOiI1YTk4MjI3NjQ5NjcxOTM3ZTc2NGVjODE5ZDA1NjFmOCJ9.jEVYwx-yi_9xN8WLCdwnBe7zvyHErlPygNNLPLk9NNYubLpb0I4sLPb6Dn6ZIplN8oY1ZPQy9COVF0m8BMXDzg). The preview token is included because the record has not yet been published.

All raw mass-spectrometry files and associated study metadata are available from [MetaboLights study MTBLS12528](https://www.ebi.ac.uk/metabolights/MTBLS12528). The Zenodo record contains the processed input tables used by the analysis scripts in this repository.

Install the one package needed by the downloader, then run it from the repository root:

```r
install.packages("jsonlite")
```

```bash
Rscript download_zenodo_inputs.R
```

This command uses `download_zenodo_inputs.R` and `zenodo_utils.R` to:

1. resolve the unpublished Zenodo draft using its preview token;
2. create `input_folder/` if it does not exist; and
3. download all Zenodo files into that folder without recreating subdirectories.

Existing files are skipped. To replace them, use:

```bash
Rscript download_zenodo_inputs.R input_folder --overwrite
```

The analysis scripts expect this layout:

```text
.
├── README.md
├── download_zenodo_inputs.R
├── zenodo_utils.R
├── *.R
├── Fig2B_E_untar_analysis.Rmd
└── input_folder/
    ├── AB012I_TEER_24well.xlsx
    ├── AB012I_cellviab_24well.txt
    ├── AB012_drugs_mets_table.csv
    ├── Caco2_TEER.csv
    ├── Caco2_viability_data.csv
    ├── Data_AB012I_drugassay3.csv
    ├── Data_AB012I_markers_transwell1.csv
    ├── Data_AB012I_transwell_drugs.csv
    ├── FT5225_ion.csv
    ├── FT7682_ion.csv
    ├── Fig2A_drugs_all_raw_with_metadata.rds
    ├── MZMINE_*__iimn_gnps_quant.csv
    └── mzmine_dmso_api_baso__quant.csv
```

> **Pre-publication note:** the preview URL contains an access token. Anyone with access to this GitHub repository can use it to view and download the unpublished data. After the Zenodo record is published, replace the preview URL in this README and in `download_zenodo_inputs.R` with the permanent record URL or DOI.

## R packages

The scripts require the following CRAN and Bioconductor packages:

```r
install.packages(c(
  "tidyverse", "data.table", "RColorBrewer", "ggforce", "readxl",
  "xlsx", "hrbrthemes", "gplots", "fuzzyjoin", "gridExtra",
  "reshape2", "showtext", "plyr", "cowplot", "patchwork", "here",
  "rmarkdown", "jsonlite"
))

install.packages("BiocManager")
BiocManager::install(version = "3.16") # Bioconductor release for R 4.2
BiocManager::install(c(
  "ComplexHeatmap", "xcms", "MSnbase", "mzR", "BiocParallel"
))
```

Create the output directories used by the scripts:

```r
dir.create("Plot", showWarnings = FALSE)
dir.create("Plots", showWarnings = FALSE)
dir.create("Plots/drugassays", recursive = TRUE, showWarnings = FALSE)
dir.create("Plots/MZmine_plots", recursive = TRUE, showWarnings = FALSE)
dir.create("Env", showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)
dir.create("EIC_output", showWarnings = FALSE)
```

## Script-to-data map

All paths below are relative to the repository root.

| Script | Analysis or figure | File(s) read from `input_folder/` |
|---|---|---|
| `Fig1_TEER_viability_script.R` | Caco-2 viability and TEER | `Caco2_viability_data.csv`; `Caco2_TEER.csv` |
| `Fig1_marker_analysis.R` | Marker-compound analysis | `markers_AB012G1_AB012G5_AB012G6.csv`; `markers_AB012G2_AB012G3.csv` |
| `Fig1E_Cell_viability&TEER.R` | AB012I viability and TEER | `AB012I_cellviab_24well.txt`; `AB012I_TEER_24well.xlsx` |
| `Fig1F_markers_transwell1.R` | Markers in transwell experiment 1 | `Data_AB012I_markers_transwell1.csv` |
| `Fig2A_tardrug_analysis_pipeline.R` | Targeted screen and Figure 2A | `Fig2A_drugs_all_raw_with_metadata.rds`; `NEW_pooling_scheme_info_with_IS.csv`; `Mapping_BugDrug_all_binary.csv` |
| `Fig2B_E_untar_analysis.Rmd` | Untargeted MZmine analysis and Figures 2B–E | All 11 `MZMINE_*__iimn_gnps_quant.csv` files; `NEW_pooling_scheme_info_with_IS.csv`; `TP8_ttest_results.RDS`; `AB005_results_p02_10percent_241126.csv` |
| `Fig2EF_Drugassay.R` | Follow-up drug assays and Figures 2E–F | `Data_AB012I_drugassay3.csv` |
| `Fig3_plots.R` | Drug/metabolite time courses | `AB012_drugs_mets_table.csv` |
| `Fig3_bisacodyl.R` | Bisacodyl time course | `AB012_drugs_mets_table.csv` |
| `Fig3D_transwell_drugs.R` | Drug transport and metabolism, including dexamethasone | `Data_AB012I_transwell_drugs.csv` |
| `Fig4_analysis.R` | Untargeted apical/basolateral feature analysis | `mzmine_dmso_api_baso__quant.csv`; `mimedb_metabolites_v20240319.csv`; raw `.mzML` files downloaded from [MTBLS12528](https://www.ebi.ac.uk/metabolights/MTBLS12528) |
| `Fig4CD_plots.R` | Annotated-hit plots for Figures 4C–D | `FT7682_ion.csv`; `FT5225_ion.csv`; `AB007_baso_ttest_candidates.csv`; `AB007_api_ttest_candidates.csv`; raw `.mzML` files downloaded from [MTBLS12528](https://www.ebi.ac.uk/metabolights/MTBLS12528) |
| `compare_tonnobug_mzmine.R` | Helper for feature-wise bacterial-versus-control testing | No direct input; sourced by `Fig4_analysis.R` |
| `EIC_all_mzmine_function.R` | Helper for EIC extraction | Raw `.mzML` files supplied by `Fig4_analysis.R` |
| `EIC_all_annot_function.R` | Helper for annotated EIC extraction | Raw `.mzML` files supplied by `Fig4CD_plots.R` |

The 11 MZmine tables correspond to *A. naeslundii*, *A. omnicolens*, *B. thetaiotaomicron*, *B. uniformis*, *C. comes*, *C. ramosum*, *C. scindens*, *D. formicigenerans*, *E. faecalis*, *G. haemolysans*, and the no-bacteria control.

## Run the analyses

Run commands from the repository root. Examples:

```bash
Rscript Fig1_marker_analysis.R
Rscript 'Fig1E_Cell_viability&TEER.R'
Rscript -e 'rmarkdown::render("Fig2B_E_untar_analysis.Rmd")'
Rscript Fig2EF_Drugassay.R
```

The main dependent analyses should be run in this order:

1. `Fig2A_tardrug_analysis_pipeline.R`
2. `Fig2B_E_untar_analysis.Rmd`
3. `Fig4_analysis.R`
4. `Fig4CD_plots.R`

## Files still missing from the current Zenodo preview

The raw AB007 `.mzML` files are intentionally hosted in [MetaboLights study MTBLS12528](https://www.ebi.ac.uk/metabolights/MTBLS12528), rather than in the Zenodo record. Download the required raw files from MetaboLights and place them in `input_folder/` before running the EIC-extraction sections of `Fig4_analysis.R` or `Fig4CD_plots.R`.

## Citation

If you use this code or its associated data, please cite the accompanying publication and the final Zenodo dataset DOI.
