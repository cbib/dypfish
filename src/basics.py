# -*- coding: utf-8 -*-

import os
import socket
from configobj import ConfigObj
from path import config_dir

def get_hostname():
    """
    Get host name
    :return: host name in lower case
    """
    return socket.gethostname().lower()

def load_config(filepath=None):
    """
    Load config object from a file
    :param filepath: path of config file
    :return: a config object
    """
    if not filepath:
        filename = '%s.ini' % get_hostname()
        filepath = os.path.join(config_dir, filename)
        if not os.path.exists(filepath):
            filepath = os.path.join(config_dir, 'default.ini')
    return ConfigObj(filepath)
