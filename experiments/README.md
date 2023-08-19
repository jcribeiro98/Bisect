# Experiments

This directory contais all of the necessary functions and files to reproduce the experiments from [2]. We also include the files to reproduce any of the appendix experiments.

## Files

We now include a list with all the files and its corresponding experiments:
    `new_execution_time_exp-synth.R`: Performs the time experiments with the synthetic datasets. The synthetic datasets are generated with a function during the experiment.
    `new_execution_time_exp-real.R`: Performs the real dataset experiments. The real datasets should be included in a file called dataset from the main directory (a tutorial is included in this very same file to download all of the datasets).
    `new_balancing_exp.R`: Functions to perform the Supervised Outlier Detection experiments. To run experiments with cWGAN, one should pull the repo and take the function from `src_extra_detectors/` and include it inside this file.
    `one_class_classification.R`: Functions to perform the One-class Classification experiments.
    `launch_balancing.R`: Runs the SOD experiments.
    `launch_occ.R`: Runs the OCC experiments.
    `lof_based_generator.R`: Code including the Hyperbox generation method. 

## Installing the datasets

To install all of the datasets create a new directory in this repo's root called `databases/`. Then, download all datasets from https://www.dbs.ifi.lmu.de/research/outlier-evaluation/DAMI/ and put them directly in `databases/` (i.e., do not include them inside any extra folder in `dataset/`).
