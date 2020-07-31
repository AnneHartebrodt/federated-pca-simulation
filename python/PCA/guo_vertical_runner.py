import scipy.sparse.linalg as lsa
import matplotlib.pyplot as plt
import comparison as co
import os.path as path

import math as math
import numpy as np
import scipy.linalg as la
import scipy.linalg as lsa

import fast_power as powerit
import import_export.easy_import as easy
import import_export.import_data as imp
import guo_vertical as gv


file = '/home/anne/Documents/data/seeds_dataset.txt'
outfile = '/home/anne/Documents/featurecloud/gwas/chr10'
header = 0
rownames = None
center = False
scale = False
scale_var = False
scale01 = False
scale_unit = False
p = 23
transpose = False
sep = ','
#
import pandas as pd

data = easy.easy_import(file, header=header, rownames=rownames, center=center, scale_var=scale_var,
                        scale01=scale01, scale_unit=scale_unit,
                        outfile=outfile, transpose=transpose, sep=sep)
data = imp.CustomDataImporter().scale_center_data(data, center=True)
data = data[:,0:(data.shape[1]-1)]

g = gv.standalone(data)

s,v,d = np.linalg.svd(data)