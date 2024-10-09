# TCR-pMHC Interface Analysis

The following project contains an analysis of the unbound (*apo*) and bound (*holo*) conformations of TCR CDR loops and pMHCs. The work is described in detail in our pre-print: [Quantifying conformational changes in the TCR:pMHC-I binding interface](https://www.biorxiv.org/content/10.1101/2024.08.13.607715v1).

## Key Links

* Source code: [https://github.com/benjiemc/tcr-pmhc-interface-analysis](https://github.com/benjiemc/tcr-pmhc-interface-analysis)
* Documentation: [https://benjiemc.github.io/tcr-pmhc-interface-analysis/](https://benjiemc.github.io/tcr-pmhc-interface-analysis/)
* DOI: [https://doi.org/10.1101/2024.08.13.607715](https://doi.org/10.1101/2024.08.13.607715)

## Overview

The project has the following structure for the code, data, and analysis notebooks.

```
.
├── data/                                   ----->  Structure of the data in the project
│   ├── external/                             ----->  Other data needed for some annotations and analysis in the project
│   ├── interim/                              ----->  Place for intermediate data processing steps
│   ├── processed/                            ----->  Finalised cleaned data and results used for the analysis
│   └── raw/                                  ----->  Baseline data used for the analysis
├── docs/                                   ----->  Configuration for documentation
├── notebooks/                              ----->  Results and conclusions of the analysis (more below)
├── scripts/                                ----->  Utility scripts
├── src/                                    ----->  Python package of the code used to process and analyse the data project
│   └── tcr_pmhc_interface_analysis/
└── tests/                                  ----->  Tests for the Python package
    ├── apps/                                 ----->  Tests for command line applications used throughout
    └── unit/                                 ----->  Tests for Python modules
```

The analysis was conducted in the following set of notebooks:

| Notebook Name | Description | Associated Figures in Manuscript |
| ------------- | ----------- | -------------------------------- |
| [Ascertaining_the_generalisability_of_the_structure_data.ipynb](notebooks/Ascertaining_the_generalisability_of_the_structure_data.ipynb) | Comparison of *apo*-*holo* structure data to other general TCR data sources. | Figure 1D-F, Figure S1 |
| [Centre_of_mass_analysis.ipynb](notebooks/Centre_of_mass_analysis.ipynb) | Analysis of where the centre of mass of each chain lies and the changes in relative angles of these domains between *apo* and *holo* states. | |
| [Compare_d_scores](notebooks/Compare_d_scores.ipynb)| Analysis of the changes in backbone dihedral angles between *apo* and *holo* conformations using D-scores. | Figure S5, Figure S6D |
| [Comparing_apo_and_holo_CDR_loop_clustering.ipynb](notebooks/Comparing_apo_and_holo_CDR_loop_clustering.ipynb) | Analysis of how CDR loops change clusters between *apo* and *holo* states. | Figure 3D, Table 1 |
| [Comparison_of_apo_and_holo_CDR_loops.ipynb](notebooks/Comparison_of_apo_and_holo_CDR_loops.ipynb) | Analysis of loop movement between *apo* and *holo* states | Figure 2B-D, Figure 3B, Figure S2 |
| [Correlating_conformational_changes_to_affinity.ipynb](notebooks/Correlating_conformational_changes_to_affinity.ipynb) | Analysis of how the movement of CDR loops correlates to the affinity of TCR-pMHC interactions where data is available. | Figure S6 , Figure S7 |
| [Identify_contact_residues_on_MHC_Class_I_molecules.ipynb](notebooks/Identify_contact_residues_on_MHC_Class_I_molecules.ipynb) | Mapping of the TCR contacts onto pMHC molecules. | Figure 4, Figure S4 |
| [Length_dependency_of_conformational_change.ipynb](notebooks/Length_dependency_of_conformational_change.ipynb) | Analysis of the correlation between conformational changes and length of CDR loops or peptides. | Figure S8 |
| [Modalities_of_TCR_interactions_with_pMHC.ipynb](notebooks/Modalities_of_TCR_interactions_with_pMHC.ipynb) | Investigation of whether TCR/peptides can be flexible and rigid depending on the interaction. | |
| [Per_residue_changes_of_TCR_CDR_between_apo_and_holo_structures.ipynb](notebooks/Per_residue_changes_of_TCR_CDR_between_apo_and_holo_structures.ipynb) | Analysis of how much each residue moves between *apo* and *holo* states for CDR loops. | Figure 3C, Figure S3 |
| [pMHC_movement_based_on_peptide_anchoring.ipynb](notebooks/pMHC_movement_based_on_peptide_anchoring.ipynb) | Analysis of how the anchoring of peptides in the MHC binding groove affects the conformational change of peptides between *apo* and *holo* states. | Figure 5B |
| [pMHC_movement_between_apo_and_holo_conformations.ipynb](notebooks/pMHC_movement_between_apo_and_holo_conformations.ipynb) | Comparison of how each part of the pMHC molecule moves between apo and holo states. | Figure 5A, Figure S5 |
| [Summary_of_apo_holo_data.ipynb](notebooks/Summary_of_apo_holo_data.ipynb) | Summary of the dataset used throughout the main analysis incuding TCR gene usage, MHC allele, and peptide similarities. | Figure 1A-C |
| [Visualising_CDR_loop_clustering.ipynb](notebooks/Visualising_CDR_loop_clustering.ipynb) | Visualisations of CDR loop clustering and the forms of canonical clusters. |  |

The analysis itself can be run as a pipeline using the provided [Makefile](Makefile).
Details on how to setup the environment and what commands to use can be found in the following sections.

## Installing and Setup

The code for this project can be installed from github using either of the following methods:
```
git clone git@github.com:benjiemc/tcr-pmhc-interface-analysis.git
cd tcr-pmhc-interface-analysis/

# Or using HTTPS
git clone https://github.com/benjiemc/tcr-pmhc-interface-analysis.git
cd tcr-pmhc-interface-analysis/
```

[Conda](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html) is used to manage the dependencies of the project.
Please ensure you have it downloaded and available as a command before running the setup commands.
To quickly build the environment, use the following command:
```
make environment
```

Which is a wrapper around:

```
conda env create -f environment.yml
conda run -n tcr-pmhc-interface-analysis python -m pip install .
```

This will install all of the required dependencies into a new environment and ensure the distributed code is installed and available.
The environment can then be activated using the following command:

```
conda activate tcr-pmhc-interface-analysis
```

The testing pipeline can be run to ensure everything has been installed correctly using the following command (you will need the additional testing tools that can be installed with `pip install '.[develop]'`):

```
make test
```

## Running the Analysis

The analysis can be run using the provided Makefile.
Once the environment is set, the whole analysis can be run using the following command:

```
make all
```

This provides a wrapper around the following steps `make data`, `make analysis`, and `make notebooks`.
However, running all stages in one command will be highly resource intensive and therefore it may be more desirable to run each stage individually depending on the target system.

the `make <COMMAND> --recon` may help ascertain what commands are run in each stage of the analysis and these can be run individually.

> **_IMPORTANT NOTE_**: The processed data used for the results reported in the manuscript has been provided for reproducibility.
> If you want to run the analysis with updated data, the provided data must be renamed or deleted (or each command can be run individually) as the make workflow will not run commands with existing outputs.

## Citing this Work

The results of this analysis are described in the article [here](https://www.biorxiv.org/content/10.1101/2024.08.13.607715v1). If you use the code, please cite:

```
@article{mcmasterQuantifyingConformationalChanges2024,
  title = {Quantifying Conformational Changes in the {{TCR}}:{{pMHC-I}} Binding Interface},
  author = {McMaster, Benjamin and Thorpe, Christopher and Rossjohn, Jamie and Deane, Charlotte and Koohy, Hashem},
  journal = {bioRxiv},
  doi = {10.1101/2024.08.13.607715},
  url = {https://www.biorxiv.org/content/10.1101/2024.08.13.607715v1},
  year = {2024},
}
```
