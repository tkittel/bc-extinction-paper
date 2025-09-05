# bc-extinction-paper

Code and documents related to an update of the 1973 Becker-Coppens extinction recipes from the paper:

"Extinction within the limit of validity of the Darwin transfer equations. I. General formalism for primary and secondary extinction and their applications to spherical crystals"
P. J. Becker and P. Coppens
Acta Cryst. (1974). A30, 129-147
https://doi.org/10.1107/S0567739474000337

To run the code, first create a conda environment based on the `conda.yml` file found in the root of the repository.

Then launch the various scripts from the `bin/` subdirectory. It is in principle believed that the code will work on any platform supported by conda-forge, but note that this has only been tested and intended for on Ubuntu 24.

Compilation of the LaTeX sources in the `paper/` subdirectory is to be done with `pdflatex` and requires a suitable LaTeX installation (probably `sudo apt install texlive` is enough on Ubuntu 24 but that is not verified). For convenience the compilation can be carried out by stepping into the `paper/` subdirectory and running pdflatex. A convenience BASH script (`compile.x`) can be used to build it.
