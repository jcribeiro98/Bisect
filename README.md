# R-Linear-search

Implementation of the BISECT algorithm for Hidden Outlier generation.

## Source Code (src)

Source Code can be found in the `src/` directory.
   Directory `HOGen/` contains all of the necessary code to run the Hidden Outlier generators. We implemented both BISECT and HIDDEN [1].
   Directory `ODM/` contains all the necessary code about the adversaries. This includes fitting, inference and handling
   Directory `registry_extras/` contais all the auxiliary files used for the implementation. It includes some new methods for `set` class objects and the dataset handler for the method.
  
## Instalation Guide

  To install, it suffices to add all of the packages from the 'R' files included via Rstudio package handler. If one wants to run python scripts (necessary to use any adversary beyond the Mahalanobis distance) one needs to install the following dependencies inside a virtual environment named `hidden_out`:
    Python 3.10
    NumPY 1.23
    Pandas 1.5.2
    Pillow 9.3
    PyOD 1.0.7
    TensorFLOW 2.11
    Scikit-Learn  1.2.0
    Scipy 1.9.3
    Seaborn 0.12.2
    Numba 0.56.4
    Jupyter 1.0.0
  We included the file `requirements.txt` including all dependencies from our virtual environment. After cloning the repository, simply execute:
  ```> pip install -r requirements.txt```

## User Guide

  An example of use to generate hidden outliers using this package is included in the file `example.R`. The directory `experiments/` contains all of the necesary files to reproduce the experiments from [2].

## References

[1] Steinbuss, G., Böhm, K. _Hiding outliers in high-dimensional data spaces_. Int J Data Sci Anal 4, 173–189 (2017).
[2] Anonimized, _Efficient Generation of Hidden Outliers for Improved Outlier Detection._
