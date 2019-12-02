import logging
import sys
import os
import path
import json

plot_colors = ['#0A3950', '#1E95BB', '#A1BA6D', '#F16C1B', '#C02A18', '#E9CB45']
plot_colors_chx = ['#0A3950', '#1E95BB', '#A1BA6D', '#d3deba']
plot_colors_cytoD = ['#1E95bb', '#1ec5d4']
my_pal = {"MTOC": "#66b2ff", "Non MTOC": "#003366", "MTOC leading edge": "#0080ff"}

cell_type_micropatterned = 'micropatterned'

def enable_logger():
    logger = logging.getLogger('DYPFISH_HELPERS')
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(filename)s - %(message)s', "%Y-%m-%d %H:%M:%S")
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.info("Running %s", sys.argv[0])
    return logger

def check_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)
    return path


def loadconfig(data_label):
    data={}
    print(path.analysis_data_dir)
    with open(path.config_data_dir+'/config_'+data_label+'.json', 'r') as fichier:
        data = json.load(fichier)
    return data
