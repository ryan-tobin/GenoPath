# GenoPath

**Developers:** Ryan Tobin, Sayaka Miura, Sudhir Kumar  
**Institute:** Institute for Genomics and Evolutionary Medicine, Temple University  

## Introduction

GenoPath is a genomic analysis pipeline that combines efficiency and user-friendliness to ease tumor evolution analysis. It integrates cutting-edge bioinformatics tools for predicting clones, inferring cell migration routes, and more.

## Installation
### System Requirements

- Operating System: Windows
- Python: 3.11
- R: 4.3.0
- Java
- MegaCC (Download available from http://www.megasoftware.net)
- Anaconda or Miniconda (recommended)

### Steps

1. Clone the repository:
    ```sh
    git clone https://github.com/ryan-tobin/GenoPath
    cd GenoPath
    ```

2. Set up the environment using Anaconda or Miniconda:
    ```sh
    conda create -n genopath python=3.11
    conda activate genopath
    ```

3. Install required libraries and dependencies:
    ```sh
    pip install -r requirements.txt
    ```

## Getting Started

### Running GenoPath

To run GenoPath with the provided sample dataset:
```sh
python genopath.py --run_process All --target_dir [path/to/dir] snv Sample_Datasets/input.tsv --max_graphs_per_tree 50 --control_file Sample_Datasets/control.txt --abundance_weighted True --sv_file Sample_Datasets/sv.txt
```

## Pipeline Components

### CloneFinder

- **Purpose:** Infers tumor cell population genotypes from bulk sequencing data.
- **Usage:** 
    ```sh
    python genopath.py --run_process CloneFinder --target_dir [path/to/dir] snv [path/to/input.tsv]
    ```

### PathFinder

- **Purpose:** Reconstructs cancer cell migration routes.
- **Usage:**
    ```sh
    python genopath.py --run_process CloneFinder PathFinder --target_dir [path/to/dir] snv [path/to/input.tsv] --primary [tumor] --max_graphs_per_tree [int]
    ```

### PhyloSignare

- **Purpose:** Detects branch-specific mutation signatures.
- **Usage:**
    ```sh
    # With Known Driver Mutations
    python genopath.py --run_process CloneFinder PhyloSignare --target_dir [path/to/dir] --driver_mutation_file [path/to/file] snv [path/to/input.tsv] [path/to/control_file]

    # With Driver Mutation Calculation
    python genopath.py --run_process CloneFinder PhyloSignare --target_dir [path/to/dir] --ref_alt_file [path/to/file] --tool CGI --email [email] --token [token] --cancer_type_input [cancer_type] snv [path/to/input.tsv] [path/to/control_file]

    # No Driver Mutation Analysis
    python genopath.py --run_process CloneFinder PhyloSignare --target_dir [path/to/dir] snv [path/to/input.tsv] [path/to/control_file]
    ```

### Picante

- **Purpose:** Analyzes phylogenetic diversity and community structure.
- **Usage:**
    ```sh
    # Picante (All)
    python genopath.py --run_process CloneFinder Picante_all --target_dir [path/to/dir] snv [path/to/input.tsv] --abundance_weighted [True/False]
    ```

### Meltos

- **Purpose:** Constructs tumor phylogeny trees based on structural variants.
- **Usage:**
    ```sh
    python genopath.py --run_process CloneFinder Meltos --target_dir [path/to/dir] snv [path/to/input.tsv] --sv_file [path/to/sv_file]
    ```

### Driver Mutation Analysis

- **Purpose:** Identifies driver mutations in tumor sites.
- **Usage:**
    ```sh
    # Calculating Driver Mutations
    python genopath.py --run_process All --target_dir [path/to/dir] --ref_alt_file [path/to/file] --tool CGI --email [email] --token [token] --cancer_type_input [cancer_type] snv [path/to/input.tsv] --sv_file [path/to/sv_file]
    
    # Known Driver Mutations
    python genopath.py --run_process All --target_dir [path/to/dir] --driver_mutation_file [path/to/file] snv [path/to/input.tsv] --sv_file [path/to/sv_file]

    # No Driver Mutation Analysis
    python genopath.py --run_process All --target_dir [path/to/dir] snv [path/to/input.tsv] --sv_file [path/to/sv_file]
    ```
- **With OpenCRAVAT**
  ```sh
  python driver_mutation_with_cravat.py -i [/path/to/input] -o [/path/to/output] -c [cancer_type ]
  ```

## Output Files

Output files will be generated in the specified target directory.

## License

This project is licensed under the BSD-3 License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

GenoPath development is supported by the Institute for Genomics and Evolutionary Medicine at Temple University.
