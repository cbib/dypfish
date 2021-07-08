DypFISH
===

# Introduction

DypFISH is a python library designed for the analysis of mRNA and protein distributions in single cell confocal images.
Using this library it is possible to analyse patterns of clustering of mRNA and proteins, dependencies of mRNA-protein localization on MTOC orientation as well as interdependent dynamics globally and at specific subcellular locations. 

DypFISH library functions take as input the preprocessed images stored in HDF5 files. Examples of analysis scripts that use our library for different spatial and temporal statistics on mRNA FISH and protein IF data are also provided.

# Provided datasets

Datasets corresponding to the results presented in the manuscript "DypFISH: Dynamic Patterned FISH to Interrogate RNA and Protein Spatial and Temporal Subcellular Distribution" by A. Savulescu et al. are available on the [accompanying website](http://dypfish.org).
They should be appropriately placed under the `data` subdirectory, e.g. the `chx` dataset should uncompressed in the `data/savulescu/chx` subfolder.  
Each subdirectory  in `data` contains a README, configuration files (see explanation below) and will eventually stores all analysis results.


-----------------------------------------------------------------------------------------------------


# Installing DypFISH 

## System Requirements
DypFISH installation requires an Unix environment with [python 3.7](http://www.python.org/). 
It was tested under Linux and MacOS environments.

In order to run DypFISH your installation should include [pip](https://pypi.org/project/pip/).
         
## Installation 

The full installation process should take less than 15 minutes on a standard computer.

Clone the repository from [github](https://github.com/cbib/DypFISH)


`git clone https://github.com/cbib/dypfish.git`

Go to the DypFISH root folder

`cd dypfish/`

# Install dependencies via CONDA 
If you have Conda install, you can create the dypfish conda environment using dypfish.yml
First you need to set the prefix of your conda installation (path to envs directory) in dypfish.yml.
Then to create the environment: 

`conda env create -f tools/envs/dypfish.yml`

Then activate the environment:

`conda activate dypfish`

# Install dependencies via PIP (No conda installation) 

If you don't have Conda installed in your system, you can install python dependencies via pip program:

`pip install -r requirements.txt`

Add the src directory to the Python path:

`export PYTHONPATH=${PYTHONPATH}:$(pwd)/src`

Matplotlib targets many different use cases and output formats.
DypFISH uses non-interactive backend AGG (PNG output - high quality images using the Anti-Grain Geometry engine)

Set the MPLBACKEND environment variable:

`export MPLBACKEND="Agg"`


Helpers scripts to (i) setup a virtualenv with all requirements and (ii) run any script with enviroments variables configured are provided in the `tools` subdirectory.

## Code organization

* `src` directory contains the python library
* `src/analysis` directory contains the implemented high-level analysis scripts that produced the figures in the DypFISH paper
* `src/tests` directory contains unit tests
* `data` directory contains placeholder subdirectories for (i) an example HDF5 dataset  `data/example_hdf5` and (ii) all the datasets analysed for the paper
* `tools` directory contains venv setup scripts

## Running unit tests 

* Place yoursefl in the root directory, then execute: ```sh download_data_test.sh ```
* When the download is complete execute `export PYTHONPATH=${PYTHONPATH}:$(pwd)/src`
* To run the test, execute `python -m unittest`
* Expect ~130 unit tests to be run in ~15 minutes.   


# Using DypFISH 

DypFISH runs in a command line environment. The runtime is dependent on the hardware, certain analysis can be time consuming.


## Running available analysis on the provided data

The analysis scripts are in the `src/analysis` folder. 

To download the data from the paper: 
```sh
sh download_data_paper.sh
```

To run the available analysis on the provided data first execute at DypFISH root folder (see above):

```sh
export PYTHONPATH=${PYTHONPATH}:$(pwd)/src;export MPLBACKEND="Agg"
```

Then to execute an analysis such as cytoplasmic spread, execute: 
```sh 
python src/analysis/cytoplasmic_spread/cytoplasmic_spread.py`
```
| WARNING: **The clustering and colocalization analyses cannot be run using the same h5 secondary file as for MTOC analysis. For example, if you have already performed MTOC analysis, you need to temporarily move or rename all of your h5 secondary files (from data/savulescu/original/, data/savulescu/cytod/, etc.) prior to run the clustering and colocalization analyses!!**!|
| --- |

## Running available analysis on your data

To run the analysis on your own data, first the images have to be preprocessed and put in the HDF5 file format and second, you need to create the corresponding configuration files in dataset and analysis folders (see below) and make the analysis script refer to your analysis configuration file.

## Available analysis

Analysis scripts and corresponding `analysis configuration` files that were used to generate the figures in the DypFISH paper are provided as examples of how the DypFISH library can be used. They are contained in the `src/analysis` directory and are supposed to be run on the data in the provided data archive.

| WARNING: **DypFISH HDF5 files not allow write parallel access!! Please do not run more than one analysis at a time**!|
| --- |

## Your own analysis

You can write your own analysis scripts by placing them in the `analysis` directory and using classes and functions provided in the library (`src`).
    
## Input data format

DypFISH takes *(i) preprocessed images in HDF5 format* and *(ii) configuration files* as input. [Click here](https://www.hdfgroup.org/solutions/hdf5/) for further informations on the HDF5 format. 

### HDF5 files

Each HDF5 file contains serie(s) of images to be analysed. Each set of images is characterized according to the molecule type (mRNA or protein), name of the molecule and possibly additional information such as e.g. timepoints. This information is organized within the HDF5 hierarchy as shown below; the lowest level in the HDF5 file corresponds to individual images.

```
        -molecule_type (mrna, protein, etc.. )\
	--molecule name (arghdia, beat-actin, etc.)\
	---molecule acquisiton time (2h, 3h, 4h, etc..)\
	----molecule acquisition image number (image\_1, image\_2, 1, 2, etc...)\
```

Each image (leafs in the HDF5 file) contains the required basic image descriptors obtained by image preprocessing. These descriptors are the following:
* `cell_mask`, `nucleus_mask` and `nucleus_centroid` **required**
* `spots` coordinated for FISH images (mRNA) **required**
* `IF` for summed intensities for immunofluorescence images (protein) **required**
* `height_map` and `zero_level` for 3D images
* `mtoc_position` for the MTOC's coordinates

An example of an HDF5 file `basic.h5` is available for download and should be placed in `data/example_hdf5/`. This file is provided *both* as an example of data formatting and for unit testing of the code. It is part of the downloadable data avalaible on the [dypfish.org](http://dypfish.org).

In order to run DypFISH, you need to download the HDF5 representation of images from the website [(data.zip file)](http://dypfish.org/file/zip/data_and_conf.zip), place the zip file in the root directory `dypfish/` and unzip it there. This populates the directory `dypfish/data/` and all the subfolders with HDF5 files and enables running and testing the pipeline. 

If you wish to run the pipeline on your own data, place the HDF5 file in the `dypfish/data/name_of_your_data` directory and modify the analysis scripts to refer to your dataset files.
        
### Configuration files

There are two types of configuration files: `dataset configuration` files and `analysis configuration` files, both in json format. A given `dataset configuration` file refers to the corresponding HDF file by specifying its name, for example `basic.h5`:
```
"PRIMARY_FILE_NAME": "basic.h5"
```
It should be located in the same directory as the corresponding HDF5 file, and describes its content (molecules, their names, certain fixed image parameters such as height and width etc). It also provides the file name to store derived (secondary) image descriptors computed by various functions of the library in order to avoid recomputation when running the scripts; for example:
```
"SECONDARY_FILE_NAME": "secondary.h5"
```
An full example `example_datasets_config.json` is provided on [dypfish.org](http://dypfish.org) and `configuration files` for the datasets used in the publications are located in the same archive as the HDF5 files.

The `analysis configuration` file  should be located in the same folder as the python analysis script. This file is used to provide the parameters for the analysis scripts. An example of an `analysis confguration` file is provided in the `src/tests` directory, `test_config.json`, and covers the needs of unit testing; `srs/analysis` subfolders contain `analysis configuration` files for the corresponding analysis.

Thus, an `analysis configuration` file indicates where is located the `dataset configuration` file, which in its turn points to the HDF5 file contianing the preprocessed images and to the HDF5 file that will be created durung the analysis and conain the intermediary results.

# Getting help

For any information or help running DypFISH, you can get in touch with: 

* [Macha Nikolski](mailto:macha.nikolski[AT]u-bordeaux.fr)
* [Benjamin Dartigues](mailto:benjamin.dartigues[AT]u-bordeaux.fr)
* [Emmanuel Bouilhol](mailto:emmanuel.bouilhol[AT]u-bordeaux.fr)

# LICENSE MIT

    Copyright (c) 2020 
    Benjamin Dartigues (1)  (benjamin.dartigues@u-bordeaux.fr)
    Emmanuel Bouilhol (1,2) (emmanuel.bouilhol@u-bordeaux.fr 
    Hayssam Soueidan (1)    (massyah@gmail.com)
    Macha Nikolski (1,2)    (macha.nikolski@u-bordeaux.fr) 
    
        (1) CBiB - University of Bordeaux,
        146, rue Leo Saignat, 33076 Bordeaux, France

        (2) CNRS, IBGC - University of Bordeaux,
        1, rue Camille Saint-Saens, 33077 Bordeaux, France

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
