# BC2025 recipes: Updates of the 1973 Becker-Coppens extinction recipes

Code and documents related to an update of the 1973 Becker-Coppens extinction recipes from the paper:

"Extinction within the limit of validity of the Darwin transfer equations. I. General formalism for primary and secondary extinction and their applications to spherical crystals"
P. J. Becker and P. Coppens
Acta Cryst. (1974). A30, 129-147
https://doi.org/10.1107/S0567739474000337

## How to reference this work:

The paper describing the work is currently in preparation, and reading it is most likely a prerequisite before anything in this repository makes sense. For the time being you can reference the work as:

* "Revisiting Becker-Coppens (1974): Updated Recipes for Estimating Extinction Factors in Spherical Crystallites", T. Kittelmann et al., 2025 (in preparation)

## Download the BC2025 recipes for C, C++ or Python

In case you came here after having read the paper and now wants to download the recipes in C, C++ or Python format, they are located in the [recipes](recipes/) subdirectory:

* Implementation in C: [bc2025.c](recipes/bc2025.c)
* Implementation in C++: [bc2025.cpp](recipes/bc2025.cpp)
* Implementation in Python: [bc2025.py](recipes/bc2025.py)

The files above also contains a small test function, which contains a number of reference y(x,theta) values. This test function could be included in the unit tests of your framework if desired.

But in any case, please consult the following conditions of usage.

## Download the BC2025 reference data

The various bcdata_*.json files in In the [data/](data/) subdirectory contains dictionaries with high precision reference values of y(x,theta) for the various models.

## Conditions of usage:

You can use these files freely under the APACHE-2.0 license, but we naturally assume good scientific conducts on behalf of any users beyond mere software license issues.

That means that you should remember cite the paper (currently in preparation) if you use them directly or indirectly in your scientific work. And if you are integrating them into a piece of software on behalf of end-users, please make sure to properly reference the origin of them in a way that your users will notice (using the paper above as the reference) and reference themselves..


## Internal notes for how to use analysis code in this repository

These instructions are mostly for internal usage, but all the code is mainly kept here in the open as a future reference and in the spirit of openness. It is, however, not intended to support it as a regular software product going forward.

To run the code, first create a conda environment based on the `conda.yml` file found in the root of the repository.

Then launch the various scripts from the `bin/` subdirectory. It is in principle believed that the code will work on any platform supported by conda-forge, but note that this has only been tested on Ubuntu 24.

Compilation of the LaTeX sources in the `paper/` subdirectory is to be done with `pdflatex` and requires a suitable LaTeX installation (probably `sudo apt install texlive` is enough on Ubuntu 24 but that is not verified). For convenience the compilation can be carried out by stepping into the `paper/` subdirectory and running pdflatex. A convenience BASH script (`compile.x`) can be used to build it (not tested or supported outside of ubuntu of course).
