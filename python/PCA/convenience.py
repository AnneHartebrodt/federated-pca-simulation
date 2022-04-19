import os as os
import os.path as path
import shutil as sh
import time as time
import pandas as pd
import numpy as np

"""
    A few convenience functions for data writing
"""

def collapse_array_to_string(a, study_id):
    res = study_id + '\t'
    for i in range(len(a)-1):
        e = a[i]
        res = res + str(e) + '\t'
    res = res + str(a[len(a)-1]) +'\n'
    return res

def collapse_array_to_string_nrdropped(a, study_id, nr_dropped):
    res = study_id+'\t'+str(nr_dropped)+'\t'
    for e in a:
        res =  res + str(e)+'\t'
    res = res+'\n'
    return res


def make_eigenvector_path(inputfile, foldername):
    """
    creates a folder called eigenvectors in the input directory
    if the eigenvector folder already exists a folder named
    'eigenvectors<currentTimestamp> will be created
    Args:
        inputfile: name of the inputfile
        foldername:

    Returns: a pathname according to the stated rules

    """
    print(path.dirname(inputfile))
    if not os.path.exists(path.join(path.dirname(inputfile), foldername)):
        pn = path.join(path.dirname(inputfile), foldername)
        os.makedirs(pn)
    else:
        print('Path exists')
        pn = path.join(path.dirname(inputfile), foldername + str(time.process_time()))
        os.makedirs(pn)
    return pn


def delete(eigenvector_path):
    '''
    CAREFUL with this.
    To save disk space, the eigenvectors will be remove after the angles have been computed
    Args:
        eigenvector_path: directory which will be removed including all subdirectories.

    Returns:

    '''
    sh.rmtree(eigenvector_path)


def parse_array(value_str):
    values_str = value_str.split(',')
    values_int = []
    for v in values_str:
        values_int.append(float(v))
    return values_int


def write_summary(res, header, outfile):
    try:
        os.makedirs(path.dirname(outfile))
    except OSError:
        print(path.dirname(outfile) + ' could not be created')
    else:
        print(path.dirname(outfile) + ' was created')

    if not path.exists(outfile):
        with open(outfile, 'w') as handle:
            handle.write(header + '\n')

    with open(outfile, 'a+') as handle:
        handle.write(res + '\n')

