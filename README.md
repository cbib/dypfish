DypFISH : Dynamic Patterned FISH
===

# Introduction

DypFISH is a python library designed for the analysis of mRNA and protein distributions in single cell confocal images.
Using this library it is possible to analyse patterns of clustering of mRNA and proteins, dependencies of mRNA-protein localization on MTOC orientation as well as interdependent dynamics globally and at specific subcellular locations. 

DypFISH library functions take as input the preprocessed images stored in HDF5 files. Examples of analysis scripts that use our library for different spatial and temporal statistics on mRNA FISH and protein IF data are also provided.

# Provided datasets

Datasets corresponding to the results presented in the manuscript "DypFISH: Dynamic Patterned FISH to Interrogate RNA and Protein Spatial and Temporal Subcellular Distribution" by A. Savulescu et al. are available on the [accompanying website](http://dypfish.org). `data` directory both contains the corresponding `datset configuration` files (see explanation below) as well as plays the role of a place holder for the available data after download. 

# Installing DypFISH 

## System Requirements
DypFISH installation requires an Unix environment with [python 3.7](http://www.python.org/). 
DypFISH was implemented in Python 3 and tested under Linux environment.

In order to run DypFISH your installation should include [pip](https://pypi.org/project/pip/).
         
## Installation 

The full installation process should take less than 15 minutes on a standard computer.

Clone the repository from [github](https://github.com/cbib/DypFISH)

`git clone https://github.com/cbib/dypfish.git`

Go to the DypFISH root folder

`cd dypfish/`

Then install python dependencies :

`pip install -r requirements.txt`

Add the current directory to the Python path:

`export PYTHONPATH=${PYTHONPATH}:$(pwd)`

Helpers scripts to (i) setup a virtualenv with all requirements and (ii) run any script with enviroments variables configured are provided in the `tools` subdirectory.


# Using DypFISH 

DypFISH runs in a command line environment. The runtime is dependent on the hardware, certain analysis can be time consuming.

## Code organization

* `src` directory contains the python library
* `src/analysis` directory contains the implemented high-level analysis scripts that produced the figures in the DypFISH paper
* `src/tests` directory contains unit tests
* `data` directory contains `dataset configuration` files for the (i) example of a dataset in `data/example_hdf5` and all the analysis performed for the paper
    
## Input data

DypFISH takes (i) preprocessed images in HDF5 format and (ii) configuration files as input. [Click here](https://www.hdfgroup.org/solutions/hdf5/) for further informations on the HDF5 format. 

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

The `analysis configuration` file  should be located located in the same folder as the python analysis script. This file is used to provide the parameters for the analysis scripts. An example of an `analysis confguration` file is provided in the `src/tests` directory, `test_config.json`, and covers the needs of unit testing; `srs/analysis` subfolders contain `analysis configuration` files for the corresponding analysis.

Thus, an `analysis configuration` file indicates where is located the `dataset configuration` file, which in its turn points to the HDF5 file contianing the preprocessed images and to the HDF5 file that will be created durung the analysis and conain the intermediary results.

## Prepping the input data

If you wish to test the DypFISH analysis, you can download the data from [our website](http://dypfish.org/file/zip/data.zip) 
and the try running the analysis scripts in the `src/analysis` folder. 

To run the analysis on your own data, first the images have to be preprocessed and put in the HDF5 file format and second, you need to create the corresponding configuration files in dataset and analysis folders. 

## Available analysis

Analysis scripts and corresponding `analysis configuration` files that were used to generate the figures in the DypFISH paper are provided as examples of how the DypFISH library can be used. They are contained in the `src/analysis` directory and are supposed to be run on the data contained in the provided data archive.

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
        351 cours de la Libération, 33400 Talence, France

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
