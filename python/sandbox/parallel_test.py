import numpy as np
from Pyfhel import Pyfhel, PyPtxt, PyCtxt
import tempfile
from pathlib import Path
import scipy.linalg as la
import scipy.linalg as la
import python.PCA.shared_functions as sh
import time
import python.PCA.convenience as cv
import python.PCA.comparison as co
import os.path as path
import sys

#import ray


@ray.remote
def encrypt(x):
    return HE.encryptFrac(x)

if __name__ == '__main__':

    setup_time = time.monotonic()
    secure_channel = tempfile.TemporaryDirectory()
    sec_con = Path(secure_channel.name)
    pk_file = sec_con / "mypk.pk"
    contx_file = sec_con / "mycontx.con"

    ##### CLIENT
    # HE Object Creation, including the public and private keys
    HE = Pyfhel()
    HE.contextGen(p=65537, m=2 ** 12)
    HE.keyGen()  # Generates both a public and a private key

    # Saving only the public key and the context
    HE.savepublicKey(pk_file)
    HE.saveContext(contx_file)
    setup_time = time.monotonic() - setup_time

    mat = sh.generate_random_gaussian(300,1)
    mat = mat.tolist()

    ray.init()

    encr = [encrypt.remote(v) for v in mat]