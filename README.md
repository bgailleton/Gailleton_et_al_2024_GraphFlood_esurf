# Gailleton_et_al_2024_GraphFlood_esurf

Repository associated with "GraphFlood 1.0: an efficient algorithm to approximate 2D hydrodynamics for Landscape Evolution Models" submitted to esurf. 

Contains archived codes used for the analysis and may not represent the up-to-date current version of GraphFlood.

## Structure

`dagger` contains all the `c++` code required to run the analysis and its python binding thanks to `pybind11`. `scabbard` is a python-only wrapper on the top of it, helping with more general routines (e.g. reading/saving DEMs, visualisation).

## Installation

The easy way: within a conda environment run `conda install -c conda-forge daggerpy matplotlib rasterio ipympl jupyterlab` and `pip install pyscabbard`

## How to run the code

Comprehensive examples can be found in `Notebooks`


## Contact

Boris Gailleton - University of Rennes

boris.gailleton@univ-rennes.fr