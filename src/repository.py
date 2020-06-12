#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

import os
import pathlib
from typing import List, Union, Optional
import re

import h5py
import numpy as np
from loguru import logger


class Repository:
    def __init__(self, repo_path: pathlib.Path, ):
        self.repo_path = repo_path

    def is_present(self, image_path) -> bool:
        raise NotImplemented("ABC class")

    def is_include(self, image_path, path) -> bool:
        raise NotImplemented("ABC class")

    def is_multiple(self, image_path, path) -> bool:
        raise NotImplemented("ABC class")

    def has_attribute(self, image_path, attribute_name) -> bool:
        raise NotImplemented("ABC class")

    def get_attribute(self, path: str, attribute_name: str):
        raise NotImplemented("ABC class")

    def get(self, path):
        raise NotImplemented("ABC class")

    def get_by_regex(self, image_path, path):
        raise NotImplemented("ABC class")

    def get_multiple(self, image_path, path):
        raise NotImplemented("ABC class")

    def get_image_path_list(self, path_list: List[str]) -> List[str]:
        raise NotImplemented("ABC class")

    def save(self, image):
        raise NotImplemented("ABC class")

    def save_descriptor(self, descriptor_path: str, value: Union[np.ndarray, float, int], dtype: type):
        raise NotImplemented("ABC class")

    def clear(self):
        raise NotImplemented("ABC class")


class H5RepositoryWithCheckpoint(Repository):
    def __init__(self, repo_path: pathlib.Path, secondary_repo_path: Optional[pathlib.Path] = None):
        super(H5RepositoryWithCheckpoint, self).__init__(repo_path=repo_path)
        self.primary_repo_path = repo_path
        self.primary_repo = H5Repository(repo_path)
        if secondary_repo_path is None:
            secondary_repo_path = pathlib.Path(repo_path.parent, repo_path.stem + "_secondary.h5")
            logger.info("Reverted to default location for secondary repo path : {}", secondary_repo_path)
        self.secondary_repo_path = secondary_repo_path

        self.secondary_repo = H5Repository(pathlib.Path(self.secondary_repo_path), "a")

    def is_present(self, image_path) -> bool:
        return self.primary_repo.is_present(image_path) or self.secondary_repo.is_present(image_path)

    def is_include(self, image_path, path) -> bool:
        return self.primary_repo.is_include(image_path, path)

    def is_multiple(self, image_path, path) -> bool:
        return self.primary_repo.is_multiple(image_path, path)

    def has_attribute(self, image_path, attribute_name) -> bool:
        return self.primary_repo.has_attribute(image_path=image_path, attribute_name=attribute_name) \
               or self.secondary_repo.has_attribute(image_path=image_path, attribute_name=attribute_name)

    def get(self, path):
        if self.primary_repo.is_present(path):
            return self.primary_repo.get(path)
        elif self.secondary_repo.is_present(path):
            return self.secondary_repo.get(path)
        else:
            raise LookupError

    def get_by_regex(self, image_path, path):
        if self.primary_repo.is_include(image_path, path):
            return self.primary_repo.get_by_regex(image_path, path)
        elif self.secondary_repo.is_include(image_path, path):
            return self.secondary_repo.get_by_regex(image_path, path)
        else:
            raise LookupError
    def get_multiple(self, image_path, path):
        if self.primary_repo.is_include(image_path, path):
            return self.primary_repo.get_multiple(image_path, path)
        else:
            raise LookupError

    def get_attribute(self, path: str, attribute_name: str):
        if self.primary_repo.has_attribute(path, attribute_name):
            return self.primary_repo.get_attribute(path, attribute_name)
        elif self.secondary_repo.has_attribute(path, attribute_name):
            return self.primary_repo.get_attribute(path, attribute_name)

    def save_descriptor(self, descriptor_path: str, value: Union[np.ndarray, float, int], dtype: type):
        self.secondary_repo.save_descriptor(descriptor_path, value=value, dtype=dtype)

    def clear(self):
        self.secondary_repo.clear()

    def get_image_path_list(self, path_list: List[str]) -> List[str]:
        return self.primary_repo.get_image_path_list(path_list)


class H5Repository(Repository):
    repo: h5py

    def __init__(self, repo_path: pathlib.Path, rights: str = "r"):
        super(H5Repository, self).__init__(repo_path=repo_path)
        self.repo_path = repo_path
        self.repo = h5py.File(repo_path, rights)
        self.rights = rights
        logger.info("Opened repo at {}", self.repo_path)

    def __eq__(self, r2: Repository) -> bool:
        return self.repo_path == r2.repo_path

    def is_present(self, path):
        # print("is present", path)
        # print("path", path)
        # print(list(self.repo.keys()))
        if path not in self.repo.keys():
            return False
        return True

    def is_include(self, image_path, path):
        r = re.compile(path.split('/')[::-1][0]+".*")
        match_list = list(filter(r.match, self.repo[image_path]))
        if len(match_list) == 0:
            return False
        return True

    def is_multiple(self, image_path, path):
        r = re.compile(path.split('/')[::-1][0] + ".*")
        match_list = list(filter(r.match, self.repo[image_path]))
        if len(match_list) == 1:
            return False
        return True

    def has_attribute(self, image_path: str, attribute_name: str) -> bool:
        attributes = self.repo[image_path].attrs.keys()
        return attribute_name in attributes

    def get(self, path):
        return self.repo[path]

    def get_by_regex(self,image_path,  path):
        r = re.compile(path.split('/')[::-1][0] + ".*")
        match_list = list(filter(r.match, self.repo[image_path]))
        return self.repo[image_path+ '/' + match_list[0]]

    def get_multiple(self, image_path, path):
        r = re.compile(path.split('/')[::-1][0] + ".*")
        match_list = list(filter(r.match, self.repo[image_path]))
        #print(len([self.repo[image_path+ '/' + x] for x in match_list]))
        return [self.repo[image_path+ '/' + x] for x in match_list]

    def get_attribute(self, path: str, attribute_name: str):
        return self.repo[path].attrs[attribute_name]

    def get_image_path_list(self, path_list: List[str]) -> List[str]:
        all_path_keys = []
        for path in path_list:
            if not path.endswith("/"):
                raise ValueError("All paths must end with a /")
            path_keys = []
            self.repo[path].visit(path_keys.append)
            path_keys = [path + k for k in path_keys]
            images_in_path = set(
                [self.repo[p].parent.name for p in path_keys if isinstance(self.repo[p], h5py.Dataset)])
            all_path_keys.extend(list(images_in_path))
        return all_path_keys

    def save(self, image):
        raise PermissionError("Repository %s can't write to a secondary checkpoint HDF5 file" % self.repo_path)

    def save_descriptor(self, descriptor_path: str, value: Union[np.ndarray, float, int], dtype: type):
        if self.rights == 'r':
            raise PermissionError("Repository %s can't write to a secondary checkpoint HDF5 file" % self.repo_path)
        else:
            self.repo.create_dataset(descriptor_path, data=value, dtype=dtype)

    def clear(self):
        if self.rights == 'r':
            raise PermissionError("Repository %s is read only" % self.repo_path)
        else:
            os.remove(str(self.repo_path))
            logger.info("Deleted repo at {}", self.repo_path)
