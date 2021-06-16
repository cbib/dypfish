#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

from image3d import Image3dMultiNucleus
from image3dWithIntensities import Image3dWithIntensitiesAndMTOC
from image3dWithSpots import Image3dWithSpotsAndMTOC, Image3dWithSpots
from imageWithSpots import ImageWithSpots, ImageWithSpotsAndMTOC
from imageWithIntensities import ImageWithIntensities, ImageWithIntensitiesAndMTOC
from repository import Repository


class ImageWithSpotsAndIntensities(ImageWithSpots, ImageWithIntensities):
    @staticmethod
    def is_a(repo: Repository, path: str):
        return ImageWithSpots.is_a(repo, path) and ImageWithIntensities.is_a(repo, path)


class ImageWithSpotsAndIntensitiesAndMTOC(ImageWithSpotsAndMTOC, ImageWithIntensitiesAndMTOC):
    @staticmethod
    def is_a(repo: Repository, path: str):
        return ImageWithSpotsAndMTOC.is_a(repo, path) and ImageWithIntensitiesAndMTOC.is_a(repo, path)


class Image3dWithSpotsAndIntensitiesAndMTOC(Image3dWithSpotsAndMTOC, Image3dWithIntensitiesAndMTOC):
    @staticmethod
    def is_a(repo: Repository, path: str):
        return Image3dWithSpotsAndMTOC.is_a(repo, path) and Image3dWithIntensitiesAndMTOC.is_a(repo, path)


class Image3dMultiNucleusWithSpots(Image3dMultiNucleus, Image3dWithSpots):
    """
    Represents an image with identified spots (e.g. from FISH), has to have spots descriptor
    """

    @staticmethod
    def is_a(repo: Repository, path: str):
        return Image3dMultiNucleus.is_a(repo, path) and Image3dWithSpots.is_a(repo, path)
