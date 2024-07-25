# TCR-pMHC Interface Analysis

The following project contains an analysis between the unbound (*apo*) and bound (*holo*) confomations of TCR CDR loops and pMHCs.

## Overview

The project has the following structure for the code, data, and analysis notebooks.

```
.
├── data/                                   ----->  Structure of the data in the project
│   ├── external/                             ----->  Other data needed for some annotations and analysis in the project
│   ├── interim/                              ----->  Place for intermediate data processing steps
│   ├── processed/                            ----->  Finalised cleaned data and results used for the analysis
│   └── raw/                                  ----->  Baseline data used for the analysis
├── notebooks/                              ----->  Results and conclusions of the analysis (more below)
├── scripts/                                ----->  Utility scripts
├── src/                                    ----->  Python package of the code used to process and analyse the data project
│   └── tcr_pmhc_interface_analysis/
└── tests/                                  ----->  Tests for the python package
    ├── apps/                                 ----->  Tests for command line applications used throughout
    └── unit/                                 ----->  Tests for python modules
```

The analysis was conducted in the following set of notebooks:

| Notebook Name | Description | Associated Figures in Manuscript |
| ------------- | ----------- | -------------------------------- |
| [Ascertaining_the_generalisability_of_the_structure_data.ipynb](notebooks/Ascertaining_the_generalisability_of_the_structure_data.ipynb) | Comparison of *apo*-*holo* structure data to other general TCR data sources. | |
| [Centre_of_mass_analysis.ipynb](notebooks/Centre_of_mass_analysis.ipynb) | Analysis of where the centre of mass of each chain lies and the changes in relative angles of these domains between *apo* and *holo* states. | |
| [Comparing_of_apo_and_holo_CDR_loop_clustering.ipynb](notebooks/Comparing_of_apo_and_holo_CDR_loop_clustering.ipynb) | Analysis of how CDR loops change clusters between *apo* and *holo* states. | Figure 2D, Table 1 |
| [Comparison_of_apo_and_holo_CDR_loops.ipynb](notebooks/Comparison_of_apo_and_holo_CDR_loops.ipynb) | Analysis of loop movement between *apo* and *holo* states | Figure 1B-D, Figure 2B, Figure S1 |
| [Correlating_conformational_changes_to_affinity.ipynb](notebooks/Correlating_conformational_changes_to_affinity.ipynb) | Analysis of how the movement of CDR loops correlates to the affinity of TCR-pMHC interactions where data is available. | Figure S5 , Figure S6 |
| [Identify_contact_residues_on_MHC_Class_I_molecules.ipynb](notebooks/Identify_contact_residues_on_MHC_Class_I_molecules.ipynb) | Mapping of the TCR contacts onto pMHC molecules. | Figure 3, Figure S3 |
| [Length_dependency_of_conformational_change.ipynb](notebooks/Length_dependency_of_conformational_change.ipynb) | Analysis of the correlation between conformational changes and length of CDR loops. | Figure S7 |
| [Per_residue_changes_of_TCR_CDR_between_apo_and_holo_structures.ipynb](notebooks/Per_residue_changes_of_TCR_CDR_between_apo_and_holo_structures.ipynb) | Analysis of how much each residue moves between *apo* and *holo* states for CDR loops. | Figure 2C, Figure S2 |
| [pMHC_movement_based_on_peptide_anchoring.ipynb](notebooks/pMHC_movement_based_on_peptide_anchoring.ipynb) | Analysis of how the anchoring of peptides in the MHC binding groove affects the conformational change of peptides between *apo* and *holo* states. | Figure 4B |
| [pMHC_movement_between_apo_and_holo_conformations.ipynb](notebooks/pMHC_movement_between_apo_and_holo_conformations.ipynb) | Comparison of how each part of the pMHC molecule moves between apo and holo states. | Figure 4A, Figure S4 |
| [Summary_of_apo_holo_data.ipynb](notebooks/Summary_of_apo_holo_data.ipynb) | Summary of the dataset used throughout the main analysis incuding TCR gene usage, MHC allele, and peptide similarities. | Figure 5 |
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
pip install .
```

This will install all of the required dependencies into a new environment and ensure the distributed code is installed and available as well.
The environment can then be activated using the following command:

```
conda activate tcr-pmhc-interface-analysis
```

Run the test pipeline to be sure everything has installed correctly:

```
make test
```

## Running the Analysis

The analysis can be run using the provide Makefile.
Once the environment is setup, the whole analysis can be run using the following command:

```
make all
```

This provides a wrapper around the following steps `make data`, `make analysis`, and `make notebooks`.
However, running all stages in one command will be highly resource intensive and therefore it may be more desirable to run each stage individually depending on the target system.

Alternatively, the `make <COMMAND> --recon` may be helpful in ascertaining what commands are run in each stage of the analysis and these can be run individually.

> **_IMPORTANT NOTE_**: The processed data used for the results reported in the manuscript has been provided for reproducibility.
> If you would like to run the analysis with updated data, the provided data will need to be renamed or deleted (or each command can be run individually) as the make workflow will not run commands with existing outputs.
