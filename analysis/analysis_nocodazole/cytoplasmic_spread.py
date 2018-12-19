#!/usr/bin/python
# encoding: UTF-8

import logging
import sys
import h5py
import src.acquisition_descriptors as adsc
import src.path as path
import src.statistical_analysis as stan
import src.helpers as helps
from src.utils import check_dir

logger = logging.getLogger('DYPFISH_HELPERS')
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(filename)s - %(message)s', "%Y-%m-%d %H:%M:%S")
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.info("Running %s", sys.argv[0])

def cytoplasmic_spread(file_handler, molecule_type,genes,path_data,colors,y_lim):
    cyt_spreads = []
    figname = check_dir(path.analysis_dir + '/analysis_nocodazole/figures/cyt_spread/') + molecule_type[0].split('/')[1] + '_cytoplasmic_spread.png'
    for gene in genes:
        image_list = helps.preprocess_image_list2(file_handler, molecule_type[0], gene)
        cyt_spreads.append(adsc.compute_cytoplasmic_spread(image_list, file_handler, path_data))
    stan.plot_bar_profile(cyt_spreads, genes, y_lim, 'Cytoplasmic spread', figname, colors)

if __name__ == "__main__":
    # Required descriptors: spots, IF, zero level, cell mask, nucleus_centroid and height_map
    # Import basics descriptors in H5 Format using 'import_h5.sh' or use own local file
    # This import script takes username and password arguments to connect to remote server bb8
    ''' 
    1-You need to create a password.txt file before running to connect via ssh
    '''
    basic_file_basename = 'basic'
    sec_file_basename = 'secondary'
    path_data = path.raw_data_dir
    basic_file_path = path.analysis_data_dir + basic_file_basename + '.h5'
    secondary_file_path = path.analysis_data_dir + sec_file_basename + '.h5'
    colors = ['#F16c1b', '#f1bc1b']
    genes = ["pard3", "pard3_nocodazole"]
    timepoints = [ "3h", "5h"]
    timepoints_protein = ["3h",  "5h"]
    cell_type = 'micropatterned/'

    with h5py.File(basic_file_path, "r") as file_handler:
        molecule_type = ['/mrna']
        cytoplasmic_spread(file_handler,molecule_type,genes,path_data,colors,0.4)

        molecule_type = ['/protein']
        cytoplasmic_spread(file_handler,molecule_type,genes,path_data,colors,0.9)