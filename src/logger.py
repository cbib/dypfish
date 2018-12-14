# -*- coding: utf-8 -*-

import os
import logging
from logging.handlers import RotatingFileHandler
from path import log_dir


class Logger(object):
    @classmethod
    def init_class_logger(cls, obj):
        """
        Init logger for a class
        :param obj: object instance
        :return: logger
        """
        logger = logging.getLogger('%s.%s' % (obj.__class__.__module__, obj.__class__.__name__))
        Logger.setup_logger(logger, obj.config)
        return logger

    @staticmethod
    def init_logger(name, config):
        print "init logger %s" % name
        logger = logging.getLogger('%s' % name)
        Logger.setup_logger(logger, config)
        return logger

    @staticmethod
    def setup_logger(logger, config):
        """
        Setup a logger from a config
        :param logger: logger to setup
        :param config: config to use
        :return:
        """
        try:
            level = config['LOGGER']['LEVEL']
        except KeyError:
            level = 'debug'
        logger.setLevel(getattr(logging, level.upper()))

        try:
            handlers = config['LOGGER']['HANDLERS']
            if isinstance(handlers, str):
                handlers = [handlers]
        except KeyError:
            handlers = ['console']

        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        for h in handlers:
            hdlr = h.lower()

            if hdlr == 'console':
                console_handler = logging.StreamHandler()
                console_handler.setFormatter(formatter)
                logger.addHandler(console_handler)

            elif hdlr == 'file':
                try:
                    tag = '_%s' % config['LOGGER']['TAG']
                except KeyError:
                    tag = ''
                logfile = os.path.join(log_dir, 'skills%s.log' % tag)
                try:
                    max_bytes = config['LOGGER']['ROTATE_MAX_BYTES']
                except KeyError:
                    max_bytes = 100000
                try:
                    count = config['LOGGER']['ROTATE_COUNT']
                except KeyError:
                    count = 3
                handler = RotatingFileHandler(logfile, maxBytes=max_bytes, backupCount=count)
                handler.setFormatter(formatter)
                logger.addHandler(handler)
