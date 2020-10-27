import pandas as pd
import numpy as np
import scipy.linalg as la
import pandas_plink as pp
import scipy.sparse.linalg as lsa

plink = pd.read_table('/home/anne/Documents/featurecloud/pca/vertical-pca/results/1000g/chr1/plink/chr1.thin.eigenvec.values')

plink = plink.iloc[:, 2:].values
approx = np.dot(plink, plink.T)

pd.DataFrame(approx).to_csv('/home/anne/Documents/featurecloud/pca/vertical-pca/results/1000g/chr1/plink/rel_approx.tsv')

data = pd.read_table('/home/anne/Documents/featurecloud/pca/vertical-pca/data/1000g/raw/chr1/chr1.thin.rel', sep='\t', header=None)


chr1 = pp.read_plink1_bin('/home/anne/Documents/featurecloud/pca/vertical-pca/data/1000g/raw/chr1/chr1.thin.bed', '/home/anne/Documents/featurecloud/pca/vertical-pca/data/1000g/raw/chr1/chr1.thin.bim','/home/anne/Documents/featurecloud/pca/vertical-pca/data/1000g/raw/chr1/chr1.thin.fam')


chr1val = np.abs(chr1.values-2)
chr1 = None

chr1valsub = chr1val[0:1250, :]
chr1valsub[np.isnan(chr1valsub)] = 0
chr1valsub = chr1valsub.T
chr1valsubinv = la.inv(chr1valsub)
linv_A = np.linalg.solve(chr1valsub.T.dot(chr1valsub), chr1valsub.T)


val = approx.flatten()



val = np.histogram(approx.flatten())


seed = pd.read_table('/home/anne/Documents/data/seeds/seeds_dataset.txt', sep=',')
seed = seed.iloc[:, 0:7].values
seedcov = np.dot(seed.T, seed)/210
covAA = seedcov[0:3, 0:3]
covAB = seedcov[0:3, 3:7]
invA= np.linalg.pinv(seed[:,0:3].T)
Bprime= np.dot(invA, covAB*210)
sum(np.abs(Bprime-seed[:,3:7]))


u,s,v = la.svd(seed.T)

seed2 = np.dot(np.dot(u[:,0:7], np.diag(s[0:7])), v[:,0:7].T)
seed2 = seed2.T
sum(np.abs(seed2-seed))

covAA = seedcov2[0:3, 0:3]
covAB = seedcov2[0:3, 3:7]
invA= np.linalg.pinv(seed[:,0:3].T)
Bprime= np.dot(invA, covAB)
sum(np.abs(Bprime-seed[:,3:7]))




seed = pd.read_table('/home/anne/Documents/data/mnist/flat.csv', sep=',')
seed = seed.values
seedcov = np.dot(seed.T, seed)/10000
covAA = seedcov[0:300, 0:300]
covAB = seedcov[0:300, 300:785]
invA= np.linalg.pinv(seed[:,0:300].T)
Bprime= np.dot(invA, covAB)
ex = np.sum(np.abs(Bprime-seed[:,300:785]))


u,s,v = lsa.svds(seed.T, k=100)

#seed = seed.values
u = np.flip(u, axis=1)
v = np.flip(v, axis=0)
s = np.flip(s)
seed2 = np.dot(np.dot(u, np.diag(s)), v)
seed2 = seed2.T
seedcov = np.dot(u[:,0:1], u[:,0:1].T)
seedcov2 = np.dot(u, u.T)
covAA = seedcov[0:300, 0:300]
covAB = seedcov[0:300, 300:785]
invA= np.linalg.pinv(seed[:,0:300].T)
Bprime= np.dot(invA, covAB)
app = np.sum(np.abs(Bprime-seed[:,300:785]))

from phe import paillier
import pandas as pd
import time
seed = pd.read_table('/home/anne/Documents/data/mnist/flat.csv', sep=',')
seed = seed.values
public_key, private_key = paillier.generate_paillier_keypair()
secret_number_list = [3.141592653, 300, -4.6e-12]
encrypted_number_list = [private_key.encrypt(x) for x in secret_number_list]

start = time.monotonic()
seed2_enc = [public_key.encrypt(x) for x in seed[:, 0].tolist()[1:1000]]
seed2_sec = [private_key.decrypt(x) for x in seed2_enc]
end = time.monotonic()
end-start

start = time.monotonic()
seed2_doub = [x+x for x in seed2_enc]
seed2_secd = [private_key.decrypt(x) for x in seed2_doub]
end = time.monotonic()
end-start

from Pyfhel import Pyfhel, PyPtxt, PyCtxt
import tempfile
from pathlib import Path

# Using a temporary dir as a "secure channel"
# This can be changed into real communication using other python libraries.
secure_channel = tempfile.TemporaryDirectory()
sec_con = Path(secure_channel.name)
pk_file = sec_con / "mypk.pk"
contx_file = sec_con / "mycontx.con"


##### CLIENT
#HE Object Creation, including the public and private keys
HE = Pyfhel()
HE.contextGen(p=65537, m=2**12)
HE.keyGen() # Generates both a public and a private key

# Saving only the public key and the context
HE.savepublicKey(pk_file)
HE.saveContext(contx_file)

# Serializing two float values
a = 1.5
b = 2.5
ca = HE.encryptFrac(a)
cb = HE.encryptFrac(b)

start = time.monotonic()
seed2_enc = [HE.encryptFrac(x) for x in seed[:, 0].tolist()[1:10000]]
seed2_sec = [HE.decrypt(x) for x in seed2_enc]
end = time.monotonic()
end-start

start = time.monotonic()
seed2_doub = [x+x for x in seed2_enc]
seed2_secd = [HE.decryptFrac(x) for x in seed2_doub]
end = time.monotonic()
end-start

start = time.monotonic()
seed2_doub = [x*x for x in seed2_enc]
seed2_secd = [HE.decryptFrac(x) for x in seed2_doub]
end = time.monotonic()
end-start