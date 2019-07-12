# Supporting scripts for Bacterial contribution to the genesis of the novel germ line determinant oskar
This repository contains all the scripts used in the analyses performed for the preprint Bacterial contribution to the genesis of the novel germ line determinant oskar. 
The preprint can be accessed on BiorXiv: https://doi.org/10.1101/453514

## Composition of the repository
This repository contains 3 folders:
* Iterative-HMMER: Git Submodule link to the Iterative HMMER script created to created Figure 1 of the paper. 
* utils: Scripts that were generated during the course of this study to help in the analysis of sequence data. Those scripts are not all necessary to reproduce the paper but are reported here for whomever might find them useful.
* notebooks: The ipython notebook used to perform the different analysis for the paper. All notebooks are referred to in the main text of the paper, see Method section.

## Requirements
The following libraries are used throughout the code:
* pandas
* numpy
* matplotlib
* scipy
* seaborn
* ETE3
* progressbar2
* scikit-learn
* BioPython
* prettyplotlib
* bs4
* dendropy
* networkx
