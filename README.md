DypFISH : Dynamic Patterned FISH
===

# Introduction

By : Macha Nikolski, Benjamin Dartigues and Emmanuel Bouilhol

A large number of mRNAs exhibit tightly controlled localization patterns specific to particular cell types, which can act as determinants of protein localization patterns through processes such as local translation. Here we describe DypFISH, an approach to investigate the relationship between mRNA and corresponding protein distributions over time by combining micropatterning of cells with Nipkow disk confocal microscopy at high resolution. 
We introduce a range of analytical techniques for quantitatively interrogating single molecule RNA FISH data in combination with protein immunolabeling. Strikingly, our results show that micropatterning of cells reduces variation in subcellular mRNA and protein distributions, allowing the characterization of their localization and dynamics with high reproducibility.
The method reveals patterns of clustering, strong dependency of mRNA-protein localization on MTOC orientation and interdependent dynamics globally and at specific subcellular locations. We establish a robust approach that is broadly applicable to a wide range of systems.

DypFISH is a python library designed to analysis confocal microscopy image.
Using hdf5 files, users can run spatial and temporal statistics on mRNA and protein signals


# Getting help

For any information or help running DypFISH, you can get in touch with: 

* [Macha Nikolski](mailto:macha.nikolski[AT]u-bordeaux.fr)
* [Benjamin Dartigues](mailto:benjamin.dartigues[AT]u-bordeaux.fr)
* [Emmanuel Bouilhol](mailto:emmanuel.bouilhol[AT]u-bordeaux.fr)


# Datasets

Datasets corresponding to the data analyzed in the manuscript are available on the [accompanying website](http://dypfish.org)


# Installation

## System Requirements
DypFISH installation requires an Unix environment with [python 2.7](python 2.7 (http://www.python.org/))
DypFISH was implemented in Python and tested under Linux environment.

In order to run DypFISH your installation should include [pip](https://pypi.org/project/pip/)
         
## Installing DypFISH 

Clone the repository from [github](https://github.com/cbib/DypFISH)

`git clone https://github.com/cbib/dypfish.git`

Go to the DypFISH root folder

`cd dypfish/`
    
Then install python dependencies :

`sudo pip install -r requirements.txt`

Add the current directory to the Python path:

`export PYTHONPATH=$(pwd)`
    
## Data input:
DypFISH takes HDF5 files as input. [Click here](https://www.hdfgroup.org/solutions/hdf5/) for further informations on the HDF5 format.

```
	-molecule_type (mrna, protein, etc.. )\
	--molecule name (arghdia, beat-actin, etc.)\
	---molecule acquisiton time (2h, 3h, 4h, etc..)\
	----molecule acquisition image number (image\_1, image\_2, 1, 2, etc...)\
```

![Arhgdia H2](http://image.noelshack.com/fichiers/2018/50/3/1544607676-capture-du-2018-12-12-10-31-06.png)

Examples of HDF5 encoding of image acquisitions can be found in DypFISH/analysis/data/. These files are provided *only* as an example of data formatting. To download data and run the code please download HDF5 files avalaible on the [dypfish.org](http://dypfish.org)

In order to run DypFISH, HDF5 representation of images should be stored in the directory DypFISH/analysis/data/


# Using the Pipeline

DypFish runs in a command line environment. 

If you wish to run DypFISH analysis, you can download the data from [our website](http://dypfish.org/file/zip/all_dypfish_data.zip) and the code will be executable as is. To run it on your own data, first the data as to be compiled in the HDF5 file format and second, the code has to be modified to include your own gene names. Version 2.0 will include automatic parsing of gene names.

The package contains one main script by analysis called main.py that coordinates the execution of the whole analysis. 
This section describes how the user must call it from the DypFISH root folder `dypfish/`.

DypFish implements several analysis:

###### Cytoplasmic total count analysis
The cytoplasmic total count descriptor was calculated as the number of transcripts or protein within the cytoplasm.
```
        python analysis/analysis_cytoplasmic_total_count/main.py
```

###### Peripheral fraction analysis 
Based on the cell masks, we calculated the peripheral fraction of mRNA and proteins at a given percent p of the radial distance.
```
        python analysis/analysis_peripheral_fraction_profile/main.py
```

###### Volume corrected noise measure (Spots density analysis) 
In order to measure gene expression noise while accounting for cell volume, we computed the volume corrected noise measure Nm for micropatterned and standardly cultured cells.
```
        python analysis/analysis_spots_density/main.py
```

###### Cytoplasmic spread analysis
Measures how evenly a molecule is spread across the cell
```
        python analysis/analysis_cytoplasmic_spread/main.py
```

###### Stability analysis 
Compares the reproducibility of distributions in standard cultured and micropatterned cells
```
        python analysis/analysis_stability/main.py
        python analysis/analysis_stability/compute_TIS_by_quad_df.py
        python analysis/analysis_stability/search_enriched_quad.py
        python analysis/analysis_stability/stability_analysis_by_fraction_profile.py
        python analysis/analysis_stability/stability_analysis_cyt_spread.py
```

###### MTOC Polarity Index 
Defines a polarity index that measures the enrichment of mRNA or protein signal for a given image acquisition series
```
        python analysis/analysis_MTOC/search_enriched_quad.py

        python analysis/analysis_MTOC/plot_figures.py
```

###### mRNA / Protein distribution profile (Correlation profile analysis) 
Defines a spatial distribution profile of mRNAs and proteins for images acquired at a given time point
```
        python analysis/analysis_correlation_profile/main.py
```

###### Temporal interaction analysis
Measures the interdependence between the mRNA and protein dynamics. **This analysis needs the results of the MTOC Polarity Index analysis to be performed.**
For the cytoplasmique analysis: 
```
        python analysis/analysis_degree_of_clustering/compute_TIS_by_quad_df.py

        python analysis/analysis_temporal_interactions/compute_TIS_analysis.py
```

For the peripheral analysis:
```
        python analysis/analysis_degree_of_clustering/compute_TIS_by_periph_quad_df.py
        
        python analysis/analysis_degree_of_clustering/compute_TIS_periph_analysis.py
```

###### Degree of clustering analysis 
The degree of clustering is a unitless measure that can be used to compare clustering between different molecules and conditions.
```
        python analysis/analysis_degree_of_clustering/compute_TIS_analysis.py
```

###### Muscle analysis
Adaptations of the previous methods to the case of muscle cells.
```
        python analysis/analysis_muscle_data/main.py
```

###### Nocodazole analysis
Adaptations of the previous methods to the case of nocodazole treated cells.
```
        python analysis/analysis_nocodazole/cytoplasmic_spread.py
        python analysis/analysis_nocodazole/cytoplasmic_total_count.py
        python analysis/analysis_nocodazole/peripheral_fraction.py
        python analysis/analysis_nocodazole/search_enriched_quad.py
        python analysis/analysis_nocodazole/compute_TIS_by_quad_df.py
        python analysis/analysis_nocodazole/compute_TIS_by_periph_quad_df.py
        python analysis/analysis_nocodazole/compute_TIS_analysis.py
        python analysis/analysis_nocodazole/compute_TIS_periph_analysis.py
        python analysis/analysis_nocodazole/plot_figures_MTOC.py
```
###### CytoD analysis
Adaptations of the previous methods to the case of nocodazole treated cells.
```
        python analysis/analysis_cytoD/search_enriched_quad.py
        python analysis/analysis_cytoD/main_cyt_spread.py
        python analysis/analysis_cytoD/main_cyt_total.py
        python analysis/analysis_cytoD/main_periph_frac.py
        python analysis/analysis_cytoD/plot_figures_MTOC.py
```

###### Muscle analysis
Adaptations of the previous methods to the case of muscle cells.
```
        python analysis/analysis_muscle_data/search_enriched_quad.py
        python analysis/analysis_muscle_data/main_cyt_spread.py
        python analysis/analysis_muscle_data/main_cyt_total.py
        python analysis/analysis_muscle_data/main_periph_frac.py
        python analysis/analysis_muscle_data/plot_figures_MTOC.py
```

## Outputs:
The output files produced by DypFish will be stored in the corresponding analysis figures and dataframe folders. 
As an example here are the outputs produced by the cytoplasmic spread analysis :
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

# LICENSE

    Copyright (c) 2018 Macha Nikolski (1,2) (macha.nikolski@u-bordeaux.fr) 
    Benjamin Dartigues (1) (benjamin.dartigues@u-bordeaux.fr) 
    Emmanuel Bouilhol (1,2) (emmanuel.bouilhol@u-bordeaux.fr)
                
    
        (1) CBiB - University of Bordeaux,
        146, rue Leo Saignat, 33076 Bordeaux, France

        (2) CNRS, LaBRI - University of Bordeaux,
        351 cours de la Libération, 33400 Talence, France

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
