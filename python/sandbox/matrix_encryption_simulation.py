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


    with open('/home/anne/Documents/featurecloud/pca/vertical-pca/results/matrix_encryption/encryption_time.tsv',
              'a') as handle:
        with open('/home/anne/Documents/featurecloud/pca/vertical-pca/results/matrix_encryption/encryption_accuracy.tsv', 'a') as handle2:
            for n in [5000, 1000]:
                print(n)
                prev_enc = None
                for it in range(10):
                    k=10
                    data = sh.generate_random_gaussian(n, k)
                    data_encrypted = []
                    start = time.monotonic()
                    for i in range(data.shape[0]):
                        for j in range(data.shape[1]):
                            data_encrypted.append(HE.encryptFrac(data[i, j]))
                    end = time.monotonic()
                    handle.write(str(it)+'\t'+str(n)+'\t'+str(k)+'\t'+'encryption'+'\t'+str(end-start)+'\t'+str(sys.getsizeof(data))+'\t'+str(sys.getsizeof(data_encrypted))+'\n')

                    data_added = []
                    start = time.monotonic()
                    if prev_enc is not None:
                        for i in range(len(data_encrypted)):
                            data_added.append (data_encrypted[i]+prev_enc[i])
                    prev_enc = data_encrypted
                    prev = data
                    end = time.monotonic()
                    handle.write(
                        str(it) + '\t' + str(n) + '\t' + str(k) + '\t' + 'addition' + '\t' + str(end - start) + '\t' + str(
                            sys.getsizeof(data)) + '\t' + str(sys.getsizeof(data_encrypted)) + '\n')

                    data_decrypted = []
                    start = time.monotonic()
                    for i in range(len(data_encrypted)):
                        data_decrypted.append(HE.decryptFrac(data_encrypted[i]))
                    end = time.monotonic()
                    handle.write(
                        str(it) + '\t' + str(n) + '\t' + str(k) + '\t' + 'decryption' + '\t' + str(end - start) + '\t' + str(
                            sys.getsizeof(data)) + '\t' + str(sys.getsizeof(data_encrypted)) + '\n')

                    start = time.monotonic()
                    data_added_decrypted = []
                    if it > 0:
                        for i in range(len(data_encrypted)):
                            data_added_decrypted.append(HE.decryptFrac(data_added[i]))
                    end = time.monotonic()
                    handle.write(
                        str(it) + '\t' + str(n) + '\t' + str(k) + '\t' + 'decryption' + '\t' + str(end - start) + '\t' + str(
                            sys.getsizeof(data)) + '\t' + str(sys.getsizeof(data_encrypted)) + '\n')

                    acc1 = np.sum(np.abs(data.flatten()-data_decrypted))
                    handle2.write(str(it) + '\t' + str(n) + '\t' + str(k) + '\t' + 'enc_dec_loss' + '\t' + str(acc1)+'\n')
                    if it > 0 :
                        data_added_orig = data+prev
                        acc = np.sum(np.abs(data_added_decrypted-data_added_orig.flatten()))
                        handle2.write(str(it) + '\t' + str(n) + '\t' + str(k) + '\t'+ 'addition_loss'+'\t'+str(acc)+'\n')

            import struct

            value = 5.1

            ba = bytearray(struct.pack("f", value))

            sum = 0
            for e in data_encrypted:
                sum = sum+sys.getsizeof(e)

            sum2 = 0
            sum3 = 0
            for s in data_decrypted:
                sum3 = sum3+sys.getsizeof(float)
                sum2 = sum2 + sys.getsizeof(bytearray(struct.pack("f", s)))