﻿This project complements the [fused_tmap](https://github.com/afloresep/fused_target_tmap) project, where we visualize nearly 40,000 compounds obtained from ChEMBL. In that project, we merged molecular fingerprints of compounds sharing the same target and explored their relationships.

In contrast, this project focuses on visualizing the same dataset, but instead of using molecular fingerprints (i.e., the vector resulting from encoding the SMILES of each compound), we use the target information to create our topological map (TMAP). For this, we obtain target information (e.g., Target Taxonomy, Target Protein Class, etc.), tokenize the words, and generate a fingerprint associated with that information.

The repository is organized as follows:

`src/` Contains the main source code.
    - fingerprint_generator.py: Main class for generating fingerprints.
    - data_loader.py: Functions for loading and preprocessing data.

`notebooks/`: Jupyter notebooks for exploration and visualization.

`data/`: Directory for storing raw and processed data.

`tests/`: Unit tests to ensure code quality.

`requirements.txt`: List of project dependencies.

`setup.py`: Configuration for installing the project as a package.
