import numpy as np
import python.PCA.vertical_pca_library as gv
import python.PCA.shared_functions as sh
import scipy.linalg as la
#import import_export.easy_import as easy
import scipy.sparse.linalg as lsa
import argparse as ap
import pandas as pd
import os.path as path
import python.import_export.mnist_import as imnist
import python.PCA.comparison  as co
from scipy.sparse import coo_matrix
import os
from Pyfhel import Pyfhel, PyPtxt, PyCtxt
import tempfile
from pathlib import Path


def simulate_federated_qr(local_data,  encrypt):

    # Using a temporary dir as a "secure channel"
    # This can be changed into real communication using other python libraries.
    if encrypt:
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
    # vector 2 norm
    alist = []
    ortho = []
    for d in local_data:
        #a = d[:,0]/sum
        a = d[:,0]
        alist.append(a)
    ortho.append(alist)

    norms = []
    for i in range(1,local_data[0].shape[1]-1):
        sums = []
        if encrypt:
            norm = HE.encryptFrac(0)
        else:
            norm = 0
        for o in ortho[i-1]:
            if encrypt:
                norm = norm+ HE.encryptFrac(np.dot(o,o))
            else:
                norm = norm + np.dot(o, o)
            print(norm)
        aplist = []
        norms.append(norm)
        print(len(ortho))

        for (o, nn) in zip(ortho, norms):
            if encrypt:
                n = HE.decryptFrac(nn)
                sum = HE.encryptFrac(0)
            else:
                n = nn
                sum = 0
            for d, o1 in zip(local_data, o):
                if encrypt:
                    sum = sum + HE.encryptFrac(np.dot(o1, d[:, i])/n)
                else:
                    sum = sum + np.dot(o1, d[:, i]) / n
            sums.append(sum)
        for d in range(len(local_data)):
            # vi
            ap = local_data[d][:, i]
            for j in range(len(sums)):
                if encrypt:
                    ap = ap - HE.decryptFrac(sums[j]) * ortho[j][d]
                else:
                    ap = ap - sums[j] * ortho[j][d]
            aplist.append(ap)
        ortho.append(aplist)
        #aplist = aplist/np.linalg.norm(aplist)

    oo = []
    for o in ortho:
        a = np.concatenate(o)
        a = a/np.linalg.norm(a)
        oo.append(a)
    ortho = np.stack(oo, axis=1)
    return ortho


def simulate_federated_qr_stabilised(local_data, ndim):
    # vector 2 norm
    alist = []
    ortho = []
    sum = 0
    # orthonnormalised first vector
    for d in local_data:
        sum = sum + np.sum(np.square(d[:,0]))
    sum = np.sqrt(sum)
    print(sum)
    # u1 = v
    for d in local_data:
        #a = d[:,0]/sum
        a = d[:,0]
        alist.append(a)
    ortho.append(alist)

    # orthonormalise 2nd vecctor
    for i in range(1, ndim):
        norm = 0
        ortho_loc = []
        for o in ortho[i-1]:
            norm = norm + np.sum(np.square(o))
        norm = np.sqrt(norm)
        print(norm)
        sum = 0
        for o1, o in zip(ortho[i-1], ortho[i]):
            sum = sum + np.dot(o1, o)/norm
        # no  = 0
        for o1, o in zip(ortho[i-1], ortho[i]):
            orth = o - sum * o[i]
            ortho_loc.append((orth))
        #     no = no + np.sum(np.square(orth))
        # no = np.sqrt(no)
        # for d in range(len(ortho_loc)):
        #     ortho_loc[d] = ortho_loc[d]/no
        ortho.append(ortho_loc)

    oo = []
    for o in ortho:
        a = np.concatenate(o)
        oo.append(a)
    ortho = np.stack(oo, axis=1)
    return ortho


if __name__ == '__main__':
    import pandas as pd
    import numpy.linalg as la
    import python.PCA.comparison
    import python.import_export.mnist_import as imnist
    #data = pd.read_table('/home/anne/Documents/data/seeds/seeds_dataset.txt', sep=',')
    data , labels = imnist.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw/', 'train')
    data = coo_matrix.asfptype(data)
    #u, v, w = lsa.svds(data, k=20)
    #u = np.flip(u, 1)
    import python.PCA.shared_functions as sh
    data_list = sh.partition_data_horizontally(data, 2)
    q,r = la.qr(data)

    ee = np.dot(data[:,0], data[:,0])

    import  time
    start = time.monotonic()
    ortho1 = simulate_federated_qr(data_list, encrypt=False)
    print(time.monotonic()-start)
    print(co.compute_angles(q, ortho1))
    start = time.monotonic()
    ortho = simulate_federated_qr(data_list, encrypt=True)
    print(time.monotonic()-start)
    print(co.compute_angles(q, ortho))
    import scipy.sparse.linalg
    print(ortho -ortho1)
