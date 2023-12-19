
# nuclease_design: ML-Guided Directed Evolution for Engineering a Better Nuclease Enzyme

This project accompanies the paper "Engineering highly active and diverse nuclease enzymes by combining machine learning and ultra-high-throughput screening". Here, we provide library code and notebooks for reproducing all of the paper's data analysis as well as a 'landscape.csv' file that provides enzyme activity for 56K enzyme variants.


## Usage

There are a few suggested ways to interact with this project:

* Use the landscape file to develop machine learning models or to generate new insights about NucB enzymes. Quick start: [TODO: open_in_colab link for plot_landscape_analysis.ipynb].


* Use the data analysis notebooks to learn more about the results in the paper or about our approaches for processing count data, assigning variants discrete enzyme activity labels, and analyzing hit rates. Quick start: [TODO: open_in_colab link for plot_hit_rates.ipynb].


## Code

Notebooks for data processing
and reproducing all of the figures and analysis in the paper's main text and
supplement are found in `notebooks/`.

Utilities for extracting enrichment factor data from raw NGS counts, analyzing
results, and making plots are available in `nuclease_design/`.

## Data
Data is available in a Google cloud storage (GCS) bucket:
`https://storage.googleapis.com/nuclease-design-data`.

It is organized in the following sub-directories:

*   `raw_count_data`: raw NGS count data for pre-sort and post-sort populations.

*   `processed_fiducial_data`: enrichment factors for synonyms of various
    'fiducial' sequences. Each row represents a distinct DNA sequence that
    translates to the same amino acid sequence.

*   `processed_data`: enrichment factors computed from the raw count data and
    the processed fiducial data. Each row represents a unique amino acid
    sequence. For each row and each fiducial, the row is assigned a p-value for
    observing its enrichment factor under the null distribution of enrichment
    factors from the fiducial.

*   `processed_data/landscape.csv`: A single file that merges data from all 4
    rounds of experiments and provides a multi-class catalytic activity
    labels for 56K distinct amino acid sequences.

*   `plate_data`: Data from the low-throughput purified protein experiments used
    to confirm hits.

*   `library_designs`: A mapping from amino acid sequences to the list of the
    names of the sub-libraries (corresponding to different sequence design
    methods) that proposed it. Note that some sequences were proposed by
    multiple methods.

*   `analysis`: Data used for creating certain tables and results in the paper
    that require expensive computations, such as clustering hits in order to
    quantify diversity.

*   `alignments`: A multiple sequence alignment used to fit our VAE model.


## Reproducing the paper's analysis

All data analysis in the paper is produced by the above notebooks. To fully
reproduce this analysis from the raw count data, run
`get_enrichment_factor_data.ipynb` with a local value of `DATA_DIR` and then
change `DATA_DIR` in the analysis notebooks. The analysis notebooks can also be
re-run without changing `DATA_DIR``, in which case they will use the enrichment
factor data that is already on GCS.

## Installation

The notebooks directly install this package from GitHub, so no installation is
necessary. However, you may want to locally install this package in order to run
tests. It can be installed as follows:

```
venv=/tmp/nuclease_design_venv
python3 -m venv $venv
source $venv/bin/activate
pip install -e .
```

After this, run the tests:

```
python -m pytest nuclease_design/*test.py
```


## Citing this work

TODO: add citation details here, usually a pastable BibTeX snippet, when the
paper appears on bioarxiv.

## License and disclaimer

Copyright 2023 DeepMind Technologies Limited

All software is licensed under the Apache License, Version 2.0 (Apache 2.0);
you may not use this file except in compliance with the Apache 2.0 license.
You may obtain a copy of the Apache 2.0 license at:
https://www.apache.org/licenses/LICENSE-2.0

All other materials are licensed under the Creative Commons Attribution 4.0
International License (CC-BY). You may obtain a copy of the CC-BY license at:
https://creativecommons.org/licenses/by/4.0/legalcode

Unless required by applicable law or agreed to in writing, all software and
materials distributed here under the Apache 2.0 or CC-BY licenses are
distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
either express or implied. See the licenses for the specific language governing
permissions and limitations under those licenses.

This is not an official Google product.



