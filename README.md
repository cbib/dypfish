DypFISH : Dynamic Patterned FISH
===

# Introduction

A large number of mRNAs exhibit tightly controlled localization patterns specific to particular cell types, which can act as determinants of protein localization patterns through processes such as local translation. Here we describe DypFISH, an approach to investigate the relationship between mRNA and corresponding protein distributions over time by combining micropatterning of cells with Nipkow disk confocal microscopy at high resolution. 
We introduce a range of analytical techniques for quantitatively interrogating single molecule RNA FISH data in combination with protein immunolabeling. Strikingly, our results show that micropatterning of cells reduces variation in subcellular mRNA and protein distributions, allowing the characterization of their localization and dynamics with high reproducibility.
The method reveals patterns of clustering, strong dependency of mRNA-protein localization on MTOC orientation and interdependent dynamics globally and at specific subcellular locations. We establish a robust approach that is broadly applicable to a wide range of systems.

DypFISH is a python library designed for the analysis of confocal microscopy images.
Using HDF5 files, users can run spatial and temporal statistics on mRNA FISH and protein IF data.

# Datasets

Datasets corresponding to the data analyzed in the manuscript "DypFISH: Dynamic Patterned FISH to Interrogate RNA and Protein Spatial and Temporal Subcellular Distribution" by A. Savulescu et al. are available on the [accompanying website](http://dypfish.org)


# Installation

## System Requirements
DypFISH installation requires an Unix environment with [python 2.7](python 2.7 (http://www.python.org/))
DypFISH was implemented in Python and tested under Linux environment.

In order to run DypFISH your installation should include [pip](https://pypi.org/project/pip/)
         
## Installing DypFISH 

The full installation process should take less than 15 minutes on a standard computer.

Clone the repository from [github](https://github.com/cbib/DypFISH)

`git clone https://github.com/cbib/dypfish.git`

Go to the DypFISH root folder

`cd dypfish/`

Install pip and system dependicies (debian-like system, else replace apt-get with appropriate package manager):

`chmod +x apt-install.sh ; ./apt-install.sh`


Then install python dependencies :

`sudo pip install -r requirements.txt`

Add the current directory to the Python path:

`export PYTHONPATH=${PYTHONPATH}:$(pwd)`
    
## Data input:
DypFISH takes HDF5 files as input. [Click here](https://www.hdfgroup.org/solutions/hdf5/) for further informations on the HDF5 format.

```
	-molecule_type (mrna, protein, etc.. )\
	--molecule name (arghdia, beat-actin, etc.)\
	---molecule acquisiton time (2h, 3h, 4h, etc..)\
	----molecule acquisition image number (image\_1, image\_2, 1, 2, etc...)\
```

![Arhgdia H2](http://image.noelshack.com/fichiers/2018/50/3/1544607676-capture-du-2018-12-12-10-31-06.png)

Examples of HDF5 encoding of image acquisitions can be found in `dypfish/analysis/example_hdf5/`. These files are provided *only* as an example of data formatting. To download data and run the code please download HDF5 files avalaible on the [dypfish.org](http://dypfish.org).

In order to run DypFISH, you need to download the HDF5 representation of images from the website [(data.zip file)](http://dypfish.org/file/zip/data_and_conf.zip), place the zip file in the root directory `dypfish/` and unzip it there. This creates the directory `dypfish/data/` and all the  subfolders containing data in HDF5 files and enables running the pipeline. 

If you wish to run the pipeline on your own data, place the HDF5 file in the `dypfish/data/name_of_your_data` directory and modify the source code (see below).

## Data settings

If you wish to run DypFISH analysis, you can download the data from [our website](http://dypfish.org/file/zip/data.zip) 
and the code will be executable as is. 

Here is the content of the data repository :

* data
    * original
        * basic.h5
        * secondary.h5
        * hstar.h5
        * mtoc.h5
        * config.json
    * nocodazole
        * basic.h5
        * secondary.h5
        * mtoc.h5
        * config.json
    * cytod
        * basic.h5
        * secondary.h5
        * config.json    
    * chx
        * basic.h5
        * secondary.h5
        * hstar.h5
        * config.json
    * muscle
        * basic.h5
        * secondary.h5
        * config.json
    * cultured
        * basic.h5
        * secondary.h5
        * config.json
        
## Configuration file        
        
Here is a example of json config file
```
{
    "PNG_IMAGES_MIME_TYPE" : "image/png",
    "BOOTSTRAP_MPI" : 100,
    "SIZE_COEFFICIENT" : 9.75,
    "VOLUME_OFFSET" : 1,
    "PERIPHERAL_FRACTION_THRESHOLD" : 30,
    "NUM_CONTOURS" : 100,
    "MAX_CELL_RADIUS" : 300,
    "DET_TOLERANCE" : 0.00000001,
    "PIXEL_PER_VOXEL" : 16,
    "VOXELS_PER_IMAGE" : 512,
    "RIPLEY_K_SIMULATION_NUMBER" : 20,
    "IMG_DIMENSION" : [512,512],
    "STRIPE_NUM" : 3,
    "IMAGE_WIDTH" : 512,
    "IMAGE_HEIGHT" : 512,
    "SLICE_THICKNESS" : 0.3,
    "Z_LINE_SPACING" : 20,
    "PIXELS_IN_SLICE" : "0.3 * 9.75",
    "VOLUME_COEFFICIENT" : "((1 / SIZE_COEFFICIENT)**2) * 0.3",
    
    "MOLECULE_TYPES": ["mrna", "protein"],
    "GENES" : ["beta_actin", "arhgdia", "gapdh", "pard3", "pkp4", "rab13"],
    "PROTEINS" : ["beta_actin", "arhgdia", "gapdh", "pard3"],
    "TIMEPOINTS_MRNA" : ["2h", "3h", "4h", "5h"],
    "TIMEPOINTS_PROTEIN" : ["2h", "3h", "5h", "7h"],
    "TIMEPOINTS_NUM_MRNA":[2, 3, 4, 5],
    "TIMEPOINTS_NUM_PROTEIN":[2, 3, 5, 7],
    "BASIC_FILE_NAME":"basic.h5",
    "SECONDARY_FILE_NAME":"secondary.h5",
    "HSTAR_FILE_NAME":"h_star.h5",
    "MTOC_FILE_NAME":"mtoc.h5"
}

```

To run it on your own data, first the data has to be compiled in the HDF5 file format and second, you need to create the corresponding configuration files called config.json. 

Then you have to create a folder in the data repository called "Your_data_label" and put all your h5 files as the corresponding config files



# Using the Pipeline

DypFISH runs in a command line environment. 

The runtime is dependent on the hardware.


Version 2.0 will include automatic parsing of gene names.

## Available analysis

The package contains one main script by general analysis called main.py that coordinates the execution of the whole analysis.
Here is the list of availbale analysis

* Cytoplasmic total count
* Peripheral fraction analysis 
* Volume corrected noise measure (analysis_density)
* Cytoplasmic spread analysis
* Stability analysis
* MTOC Polarity Index
* mRNA / Protein distribution profile (Correlation profile analysis)
* Temporal interaction analysis
* Degree of clustering analysis 
 
 
For the following analyses, you have to run the script corresponding to the sub analysis you needed (peripheral_fraction.py, cytoplasmic_spread.py, etc..). 
 * nocodazole 
 * cytod 
 * chx
 * muscle
 
For more details see section specific analysis:

### General analysis

DypFish implements several analysis:


This section describes how the user must call it from the DypFISH root folder `dypfish/`.


###### Cytoplasmic total count analysis
The cytoplasmic total count descriptor was calculated as the number of transcripts or protein within the cytoplasm.
```
        python analysis/analysis_cytoplasmic_total_count/main.py -i original
```

###### Peripheral fraction analysis 
Based on the cell masks, we calculated the peripheral fraction of mRNA and proteins at a given percent p of the radial distance.
```
        python analysis/analysis_peripheral_fraction_profile/main.py -i original
```

###### Volume corrected noise measure (Spots density analysis) 
In order to measure gene expression noise while accounting for cell volume, we computed the volume corrected noise measure Nm for micropatterned and standardly cultured cells.
```
        python analysis/analysis_spots_density/main.py -i original
```

###### Cytoplasmic spread analysis
Measures how evenly a molecule is spread across the cell
```
        python analysis/analysis_cytoplasmic_spread/main.py -i original
```

###### Stability analysis 
Compares the reproducibility of distributions in standard cultured and micropatterned cells
```
        python analysis/analysis_stability/main.py -i original
        python analysis/analysis_stability/compute_TIS_by_quad_df.py -i original
        python analysis/analysis_stability/search_enriched_quad.py -i original
        python analysis/analysis_stability/stability_analysis_by_fraction_profile.py -i original
        python analysis/analysis_stability/stability_analysis_cyt_spread.py -i original
```

###### MTOC Polarity Index 
Defines a polarity index that measures the enrichment of mRNA or protein signal for a given image acquisition series
```
        python analysis/analysis_MTOC/search_enriched_quad.py -i original

        python analysis/analysis_MTOC/plot_figures.py -i original
```

###### mRNA / Protein distribution profile (Correlation profile analysis) 
Defines a spatial distribution profile of mRNAs and proteins for images acquired at a given time point
```
        python analysis/analysis_correlation_profile/main.py -i original
```

###### Temporal interaction analysis
Measures the interdependence between the mRNA and protein dynamics. **This analysis needs the results of the MTOC Polarity Index analysis to be performed.**
For the cytoplasmique analysis: 
```
        python analysis/analysis_degree_of_clustering/compute_TIS_by_quad_df.py -i original

        python analysis/analysis_temporal_interactions/compute_TIS_analysis.py -i original
```

For the peripheral analysis:
```
        python analysis/analysis_degree_of_clustering/compute_TIS_by_periph_quad_df.py -i original
        
        python analysis/analysis_degree_of_clustering/compute_TIS_periph_analysis.py -i original
```

###### Degree of clustering analysis 
The degree of clustering is a unitless measure that can be used to compare clustering between different molecules and conditions.
```
        python analysis/analysis_degree_of_clustering/compute_TIS_analysis.py -i original
```

### specific analysis


###### Nocodazole analysis
Adaptations of the previous methods to the case of nocodazole treated cells.
```
        python analysis/analysis_nocodazole/cytoplasmic_spread.py -i nocodazole
        python analysis/analysis_nocodazole/cytoplasmic_total_count.py -i nocodazole
        python analysis/analysis_nocodazole/peripheral_fraction.py -i nocodazole
        python analysis/analysis_nocodazole/search_enriched_quad.py -i nocodazole
        python analysis/analysis_nocodazole/compute_TIS_by_quad_df.py -i nocodazole
        python analysis/analysis_nocodazole/compute_TIS_by_periph_quad_df.py -i nocodazole
        python analysis/analysis_nocodazole/compute_TIS_analysis.py -i nocodazole
        python analysis/analysis_nocodazole/compute_TIS_periph_analysis.py -i nocodazole
        python analysis/analysis_nocodazole/plot_figures_MTOC.py -i nocodazole
```
###### CytoD analysis
Adaptations of the previous methods to the case of nocodazole treated cells.
```
        python analysis/analysis_cytoD/preprocess_secondary_descriptors.py -i cytod
        python analysis/analysis_cytoD/cytoplasmic_spread.py -i cytod
        python analysis/analysis_cytoD/cytoplasmic_total.py -i cytod
        python analysis/analysis_cytoD/peripheral_fraction.py -i cytod
```

###### CHX analysis
Adaptations of the previous methods to the case of nocodazole treated cells.
```
        python analysis/analysis_CHX/preprocess_secondary_descriptors.py -i chx
        python analysis/analysis_CHX/degree_of_clustering.py -i chx
        python analysis/analysis_CHX/peripheral_fraction.py -i chx
```

###### Muscle analysis
Adaptations of the previous methods to the case of muscle cells.
```
        python analysis/analysis_muscle_data/compute_z_line_distance.py -i muscle
        python analysis/analysis_muscle_data/quadrat_analysis_graph.py-i muscle
        python analysis/analysis_muscle_data/quadrat_analysis_heatmap.py-i muscle
```

## Outputs:

The output files produced by DypFish will be stored in the corresponding analysis figures and dataframe folders. 
As an example here are the outputs produced by the cytoplasmic spread analysis, based on the full data from the website, as described in the Data Input section :
```
        analysis/analysis_cytoplasmic_spread/figures/cyt_spread_arhgdia.png
        analysis/analysis_cytoplasmic_spread/figures/cyt_spread_beta_actin_.png
        analysis/analysis_cytoplasmic_spread/figures/cyt_spread_gapdh.png
        analysis/analysis_cytoplasmic_spread/figures/cyt_spread_pard3.png
        analysis/analysis_cytoplasmic_spread/figures/cyt_spread_pkp4.png
        analysis/analysis_cytoplasmic_spread/figures/cyt_spread_rab13.png
        analysis/analysis_cytoplasmic_spread/figures/micropatterned/mrna_cytoplamsic_spread.png
        analysis/analysis_cytoplasmic_spread/figures/micropatterned/protein_cytoplamsic_spread.png
```

# Getting help

For any information or help running DypFISH, you can get in touch with: 

* [Macha Nikolski](mailto:macha.nikolski[AT]u-bordeaux.fr)
* [Benjamin Dartigues](mailto:benjamin.dartigues[AT]u-bordeaux.fr)
* [Emmanuel Bouilhol](mailto:emmanuel.bouilhol[AT]u-bordeaux.fr)

# LICENSE MIT

    Copyright (c) 2018 
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
