# Mutation signatures visualization

![TravisCIBuildBadge](https://travis-ci.com/lrgr/mutation-signatures-viz.svg?token=xpopk4qvQzXty9qXHH3S&branch=master)

A repository with scripts and functions for visualizations related to mutation signatures.

## Setup

### Dependencies
The methods and experiments are written in Python 3. We recommend using Conda to manage dependencies, which you can do directly using the provided [`environment.yml`](environment.yml) file:

    conda env create -f environment.yml
    source activate mutation-signatures-viz-env

We use [`snakemake`](https://snakemake.readthedocs.io/en/stable/index.html) to manage workflows.

## Usage
The main expected usage of this repo is importing the various plotting functions into your own Python programs.

The list of available plotting functions is below, with further document in the Wiki (TODO).
* `sbs_signature_plot()`: plot one or more distributions over SBS categories separated by substitution.

### Example

We do provide one example in [`example/`](example/) to generate a single base substitution plot from the COSMIC mutation signatures for glioblastoma: COSMIC SBS Signatures 1, 5, and 11. To download the data and run the example:

    cd example
    snakemake all
