# Source Code (`src/`)

This directory contains all the necessary files to generate hidden outliers utilizing different generators. The files contains the following algorithms:
    BISECT: Main algorithm from [2],
    HIDDEN: Implementation from [1]
    simpleBISECT: a BISECT implementation without the "cut trick" (See [2] for more information).

## HOGen/

The following files are contained in this directory:
- `bisect_strat.R`: Contains all the necessary code for performing the simpleBISECT algorithm,
- `hidden_strat.R`: Contains all the necessary code for performing the HIDDEN algorithm as in [1],
- `multi_bisect_strat.R`: Contains all the necessary code for performing the BISECT algorithm as in [2].

## ODM/

The following files are contained in this directory:
- `fit_method.R`: Functions for fitting the adversary as a full space model M and ensemble E_M,
- `inference_methods.R`: Functions for performing inference from M or E_M fitted with `fit_method.R`.

## registry_extras/

The following files are contained in this directory:
- `classes.R`: Contains the formulation as a `hog_method` class and all of its methods.
- `get_origin.R`: Function for obtaining the origin in Step 1 of BISECT [2]. Includes multiple other methods besides the weighted method (see the function description).
- `utils.R`: Extra functions for other object classes and various purposes.
