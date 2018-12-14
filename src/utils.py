import logging
import sys

import os

plot_colors = ['#0A3950', '#1E95BB', '#A1BA6D', '#F16C1B', '#C02A18', '#E9CB45']
cell_type_micropatterned = 'micropatterned'



# TODO fonction d'extractions des timepoints depuis les h5

def enable_logger():
    logger = logging.getLogger('DYPFISH_HELPERS')
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(filename)s - %(message)s', "%Y-%m-%d %H:%M:%S")
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    # if "log" not in globals():
    # logger = Logger.init_logger('REFORMAT_%s'%(cfg.language_code), load_config())
    logger.info("Running %s", sys.argv[0])


def check_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

    return path
