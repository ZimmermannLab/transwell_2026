# Suggested text for the Zenodo entry

## Suggested title

Processed data supporting "Establishing an in vitro pipeline for the high-throughput quantification of epithelial permeability of gut bacterial metabolites"

## Suggested description

This dataset contains the processed data files used for the analyses presented in the manuscript **"Establishing an in vitro pipeline for the high-throughput quantification of epithelial permeability of gut bacterial metabolites"** by Brauer-Nikonow, Voogdt, Ribeiro, and Zimmermann.

The deposited files support the targeted and untargeted metabolomics analyses, Caco-2 cell-viability and transepithelial electrical resistance (TEER) measurements, bacterial drug-metabolism experiments, and transwell transport experiments described in the manuscript. They include processed targeted-metabolomics tables, MZmine feature-quantification tables, preprocessed R data, and tables used to generate the reported figures.

The corresponding analysis scripts are available from the associated GitHub repository: **[ADD GITHUB REPOSITORY URL]**. The repository includes `download_zenodo_inputs.R` and `zenodo_utils.R`, which download all files from this Zenodo record into the `input_folder/` directory expected by the analysis scripts. The analyses were performed with **R version 4.2.2**; required R packages and the mapping between scripts and input files are documented in the repository README.

All raw mass-spectrometry files and associated study metadata are deposited in the EMBL-EBI MetaboLights repository under accession **[MTBLS12528](https://www.ebi.ac.uk/metabolights/MTBLS12528)**. The Zenodo deposit provides the processed inputs needed to reproduce the analyses, while MetaboLights provides the underlying raw data.

The associated manuscript is available as a preprint:

Brauer-Nikonow A, Voogdt CGP, Ribeiro NV, Zimmermann M. *Establishing an in vitro pipeline for the high-throughput quantification of epithelial permeability of gut bacterial metabolites.* bioRxiv. [https://doi.org/10.1101/2025.06.10.658935](https://doi.org/10.1101/2025.06.10.658935)

## Suggested keywords

- gut microbiome
- bacterial metabolites
- epithelial permeability
- transwell assay
- Caco-2
- metabolomics
- LC-MS
- drug metabolism
- TEER
- MZmine

## Suggested related identifiers

- Manuscript/preprint: `https://doi.org/10.1101/2025.06.10.658935` — relation: **is supplement to**
- MetaboLights: `https://www.ebi.ac.uk/metabolights/MTBLS12528` — relation: **is derived from** or **is supplement to**
- GitHub repository: `[ADD GITHUB REPOSITORY URL]` — relation: **is supplemented by**
