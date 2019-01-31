# Supporting scripts for Bacterial contribution to genesis of the novel germ line determinant oskar
This repository conatins all the scripts used in the analyses performed for the paper Bacterial contribution to genesis of the novel germ line determinant oskar. 
The paper can be access on BioArxiv: https://doi.org/10.1101/453514

## Composition of the repository
This repository contains 3 folders:
* Iterative-HMMER: Git Submodule link to the Itterative HMMER script created to created Figure 1 of the paper. 
* utils: Scripts that were generated during the course of this study to help in the analysis of sequence data. Thoose scripts are not all necessary to reproduce the paper but are repported here for the use by people who might find their need. 
* notebooks: The ipython notebook used to perform the different analysis for the paper. All notebooks are refered to in the main text of the paper, see Method section.

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
