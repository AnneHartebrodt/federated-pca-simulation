import numpy as np
import pandas as pd
import scipy as sc
import random as r
import math as math
import scipy.linalg as la
import scipy.sparse.linalg as lsa
import random as rnd
import os as os
import logging
import copy
import scipy.spatial.distance as d
import argparse as ap
from statsmodels import robust

class OutlierRemoval():

    def __init__(self, logger = None, file='/tmp/'):
        if logger is None:
            self.logger = logging.getLogger('default')
        else:
            self.logger = logger
        self.file = file

    # https://stackoverflow.com/questions/8930370/where-can-i-find-mad-mean-absolute-deviation-in-scipy
    def outlier_removal_mad(self, data, mulitplier=6, nr_pcs=6):
        flagged_pc = list()
        for column in range(min(nr_pcs,data.shape[1])):
            med = sc.median(data[:, column])
            mad = robust.mad(data[:, column])
            for row in range(data.shape[0]):
                if abs(data[row, column]-med)> (mulitplier*mad):
                    flagged_pc.append(row)
        flagged_pc = np.unique(flagged_pc)
        return flagged_pc
