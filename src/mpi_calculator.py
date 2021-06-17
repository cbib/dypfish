#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Credits: Benjamin Dartigues, Emmanuel Bouilhol, Hayssam Soueidan, Macha Nikolski

from dataclasses import dataclass
from typing import Optional, List, Any

import numpy as np
import pandas as pd
from loguru import logger
from scipy.stats import mannwhitneyu


@dataclass
class MTOCPolarityIndex(object):
    index: float
    pvalue: float
    errs: Optional[List[Any]] = None
    is_random: bool = False

    def compute_error(self):
        return np.std([x.index for x in self.errs])

    def upper_envelope(self):
        return self.index + self.compute_error()

    def lower_envelope(self):
        return self.index - self.compute_error()


@dataclass
class DensityStats(object):
    df: pd.DataFrame
    group_key: List[str]
    mpi_sample_size: int
    quadrant_labels: List[str]
    bootstrap_num: int = 100
    mtoc_quadrant_label: str = "MTOC"

    def make_labels(self):
        labels = []
        grouped_df = self.df.groupby(self.group_key, as_index=False)
        for label in grouped_df.groups.keys():
            if len(self.group_key) > 1:
                label = '-'.join(label)
            labels.append(label)
        return labels

    def update_group_key(self, group_key: List[str]):
        self.group_key = group_key

    def subset_stats(self, column_label, column_value):
        return DensityStats(
            df=self.df[self.df[column_label] == column_value],
            group_key=self.group_key,
            mpi_sample_size=self.mpi_sample_size,
            quadrant_labels=self.quadrant_labels,
            bootstrap_num=self.bootstrap_num,
            mtoc_quadrant_label=self.mtoc_quadrant_label
        )

    def mpi(self, use_mean=False):
        mpis, errs = [], []
        grouped_df = self.df.groupby(self.group_key, as_index=False)
        for gene in sorted(grouped_df.groups.keys()):
            logger.info("Running MPI calculation with bootstrap for {}", gene)
            k_grouped_df = grouped_df.get_group(gene).copy()
            mpi = compute_mpis(
                k_grouped_df,
                self.mpi_sample_size,
                self.quadrant_labels, self.mtoc_quadrant_label, use_mean=use_mean)
            mpis.append(mpi.index)
            errs.append(mpi.compute_error())
        return mpis, errs

    def ratios(self):
        return self.df[self.mtoc_quadrant_label] / self.df[self.quadrant_labels].mean(axis=1)


def calculate_mpi(mtoc: list, quadrants: list, use_mean=False) -> MTOCPolarityIndex:
    """
    Given two lists of floats, calculate the MTOC Polarity Index
    It represents the number of times the mtoc list values are greater than those in the other quadrants
    """
    if use_mean:
        adjusted_mtoc = mtoc - np.mean(quadrants)  # TODO : why nanmedian and not just median?
    else:
        adjusted_mtoc = mtoc - np.median(quadrants)  # TODO : why nanmedian and not just median?
    npos = sum(x > 0 for x in adjusted_mtoc)
    mpi = ((float(npos) / len(adjusted_mtoc)) * 2) - 1
    # TODO it's just for displaying
    if mpi == 0.0:
        mpi = 0.01
    statistic, p = mannwhitneyu(mtoc, quadrants, alternative='two-sided')
    return MTOCPolarityIndex(index=mpi, pvalue=p)


def calculate_random_mpi(mtoc, quadrants, mpi_sub_sample_size) -> MTOCPolarityIndex:
    adjusted_mtoc = []
    if mpi_sub_sample_size >= len(quadrants):
        mpi_sub_sample_size = int(2 * len(quadrants) / 3)
    for value in mtoc:
        tmp_list = np.random.choice(np.array(quadrants), mpi_sub_sample_size)
        med = np.median(tmp_list)
        adjusted_mtoc.append(value - med)
    assert (len(mtoc) == len(adjusted_mtoc))

    npos = sum(x > 0 for x in adjusted_mtoc)
    mpi = (float(npos) / len(adjusted_mtoc)) * 2 - 1
    statistic, p = mannwhitneyu(mtoc, quadrants, alternative='two-sided')
    return MTOCPolarityIndex(index=mpi, pvalue=p, is_random=True)


def compute_mpis(df, mpi_sub_sample_size, quadrant_labels, mtoc_quadrant_label='MTOC',
                 bootstrap_num=100, use_mean=False) -> MTOCPolarityIndex:
    """
    Compute both the bootstrapped (random) and the non-random MTOC Polarity Indices
    """
    mtoc = df[mtoc_quadrant_label].values
    nonmtoc = list(df[quadrant_labels].values.flatten())
    mpi = calculate_mpi(mtoc, nonmtoc, use_mean=use_mean)
    mpi.errs = [calculate_random_mpi(mtoc, nonmtoc, mpi_sub_sample_size) for _ in range(bootstrap_num)]

    return mpi
