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
#from python.PCA.vertical_pca_benchmark import log_transmission


def log_transmission(logfile, log_entry_string, iterations, counter, element, eigenvector=10):
    with open(logfile + '.transmission', 'a') as handle:
        if type(element)=='numpy.ndarray':
            try:
                handle.write(log_entry_string + '\t' + str(iterations)
                     + '\t' + str(counter) +'\t' + str(eigenvector)+'\t' + str(sys.getsizeof(element.tobytes()))+'\n')
            except AttributeError:
                handle.write(log_entry_string + '\t' + str(iterations)
                             + '\t' + str(counter) + '\t' + str(eigenvector) + '\t' + str(
                    sys.getsizeof(element)) + '\n')
        else:
            handle.write(log_entry_string + '\t' + str(iterations)
                         + '\t' + str(counter) + '\t' + str(eigenvector) + '\t' + str(
                sys.getsizeof(element)) + '\n')

def log_costs(filename, action, duration, split, repeat):
    with open(filename, 'a') as handle:
            handle.write(action+'\t'+str(split)+'\t'+str(repeat)+'\t'+str(duration)+'\n')

def simulate_federated_qr(local_data,  encrypt, filename=None, split=None, repeat=None):
    if filename is not None:
        log = True
    else:
        log = False
    start = time.monotonic()
    # Using a temporary dir as a "secure channel"
    # This can be changed into real communication using other python libraries.
    if encrypt:
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
        setup_time = time.monotonic()-setup_time

    encrypted_values = 0
    decrypted_values = 0
    encrypted_additions = 0

    encryption_time = 0
    decryption_time = 0
    addition_time = 0

    # vector 2 norm
    alist = []
    ortho = []
    if encrypt:
        s = time.monotonic()
        sum = HE.encryptFrac(0)
        encryption_time = time.monotonic()-s
        encrypted_values +=1

    else:
        sum = 0

    for d in range(len(local_data)):
        if encrypt:
            s = time.monotonic()
            se = HE.encryptFrac(np.dot(local_data[d][:, 0], local_data[d][:,0]))
            encryption_time += time.monotonic()-s
            encrypted_values+=1

            s = time.monotonic()
            sum = sum + se
            addition_time += time.monotonic() - s
            encrypted_additions+=1
        else:
            se = np.dot(local_data[d][:, 0], local_data[d][:,0])
            sum = sum+se
        alist.append(local_data[d][:,0])
        if log:
            log_transmission(filename, "qr_local_norm=CS", repeat, d, se)
    if log:
        log_transmission(filename, "qr_global_norm=SC", repeat, 1, sum)
    ortho.append(alist)


    # ortho [eigenvector rank] [data set]
    # list of lists containing the already
    #  orthogonal eigenvectors
    norms = [sum]
    # iterate over the eigenvectors
    for i in range(1,local_data[0].shape[1]):
        sums = []
        aplist = []

        if encrypt:
            s = time.monotonic()
            norm = HE.encryptFrac(0)
            encryption_time += time.monotonic() - s
            encrypted_values += 1
        else:
            norm = 0

        # iterate over the all already orthonormal eigenvector snippet
        # lists and the associated norms
        # decrypt norms
        for di in range(len(ortho)):
            o = ortho[di]
            nn = norms[di]
            if encrypt:
                s = time.monotonic()
                n = HE.decryptFrac(nn)
                decryption_time += time.monotonic()-s
                decrypted_values += 1
                s = time.monotonic()
                sum = HE.encryptFrac(0)
                encryption_time += time.monotonic() - s
                encrypted_values += 1
            else:
                n = nn
                sum = 0
            # iterate over the local data sets
            # combined with the local eigenvector snippets
            # o takes all othonormal ranks
            # i is the currently to be orthonomalised rank
            for ik in range(len(local_data)):
                d = local_data[ik]
                o1 = o[ik]
                if encrypt:
                    s = time.monotonic()
                    # this is the projection operation
                    se = HE.encryptFrac(np.dot(o1, d[:, i])/n)
                    encryption_time += time.monotonic() - s
                    encrypted_values += 1
                    s = time.monotonic()
                    # this is the aggregation at the server
                    sum = sum + se
                    addition_time += time.monotonic()-s
                    encrypted_additions += 1
                else:
                    se = np.dot(o1, d[:, i]) / n
                    sum = sum + se

                if log:
                    log_transmission(filename, "qr_local_dp=CS", repeat, ik, se)
            if log:
                log_transmission(filename, "qr_global_dp=SC", repeat, 1, sum)
            sums.append(sum)

        for d in range(len(local_data)):
            ap = local_data[d][:, i]
            for j in range(len(sums)):
                if encrypt:
                    s = time.monotonic()
                    # this is the local normalisation step
                    ap = ap - HE.decryptFrac(sums[j]) * ortho[j][d]
                    decryption_time += time.monotonic() - s
                    decrypted_values += 1
                else:
                    ap = ap - sums[j] * ortho[j][d]

            # compute the local norm of the freshly orthogonalised
            # eigenvector
            if encrypt:
                s = time.monotonic()
                se = HE.encryptFrac(np.dot(ap, ap))
                encryption_time += time.monotonic() - s

                s = time.monotonic()
                norm = norm + se
                addition_time += time.monotonic() - s
                encrypted_additions += 1
                encrypted_values += 1
            else:
                se = np.dot(ap, ap)
                norm = norm + se
            if log:
                log_transmission(filename, "qr_local_norm=CS", repeat, d, se)
            aplist.append(ap)
        if log:
            log_transmission(filename, "qr_global_norm=SC", repeat, 1, norm)
        norms.append(norm)
        ortho.append(aplist)

    G_list = []
    # normalise the vector norms to unit norm.
    for d in range(len(local_data)):
        oo = []
        for i in range(len(ortho)):
            if encrypt:
                s = time.monotonic()
                no = HE.decryptFrac(norms[i])
                decryption_time += time.monotonic() - s
                decrypted_values += 1
                oo.append(ortho[i][d]/np.sqrt(no))
            else:
                oo.append(ortho[i][d] / np.sqrt(norms[i]))
        oo = np.stack(oo, axis = 1)
        G_list.append(oo)

    end = time.monotonic()
    # log time for run, and time for encryption
    if filename is not None:
        if encrypt:
            logentry = 'total_encrypted'
            log_costs(filename, logentry, end - start, split=split, repeat=repeat)
            log_costs(filename, 'setup', setup_time, split=split, repeat=repeat)
            log_costs(filename, 'encryption', encryption_time, split=split, repeat=repeat)
            log_costs(filename, 'decryption', decryption_time, split=split, repeat=repeat)
            log_costs(filename, 'addition', addition_time, split=split, repeat=repeat)

        else:
            logentry = 'total_unencrypted'
            log_costs(filename, logentry, end-start, split=split, repeat=repeat)

    ortho = np.concatenate(G_list, axis=0)
    return ortho, G_list




def simulate_federated_qr_stabilised(local_data,  encrypt):

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
    sum  = 0
    # first verctor simply scaled to unit norm
    for d in local_data:
        sum  = sum+ np.sum(np.square(d[:,0]))

    sum = np.sqrt(sum)
    for d in local_data:
        a = d[:,0]/sum
        alist.append(a)

    ortho.append(alist)

    norms = []
    prev = alist
    for i in range(1,local_data[0].shape[1]):
        sums = []
        norm = 0
        for o in ortho[i-1]:
            norm = norm + np.dot(o, o)
            #print(norm)
        aplist = []
        norms.append(norm)
        #print(len(ortho))


        sum = 0
        no = 0
        for k in range(len(local_data)):
            sum = sum + np.dot(prev[k], local_data[k][:, i]) /norm
        for d in range(len(local_data)):
            # vi
            ap = prev[d]
            ap = ap - sum * prev[d]
            no = no + np.sum(np.sum(ap))
            aplist.append(ap)
        no = np.sqrt(no)
        for a in range(len(aplist)):
            aplist[a] = aplist[a]/no

        prev=aplist
        ortho.append(aplist)
        #aplist = aplist/np.linalg.norm(aplist)

    oo = []
    for o in ortho:
        a = np.concatenate(o)
        a = a/np.linalg.norm(a)
        oo.append(a)
    ortho = np.stack(oo, axis=1)
    return ortho

def all_angles_against_all(ortho1):
    oo1 = []
    for i in range(ortho1.shape[1]-2):
        for j in range(i + 1, ortho1.shape[1]-1):
            oo1.append(co.angle(ortho1[:, i], ortho1[:, j].T))
    return oo1

def benchmark_encryption(repeats, splits, n, m, filename):
    fn = filename+'_encryption'
    filename1 = filename + '_angles_11.tsv'
    filename2 = filename + '_angles_22.tsv'
    filename12 = filename + '_angles_12.tsv'
    filenameq1 = filename + '_angles_scipy_fed.tsv'

    for i in range(repeats):

        for s in splits:
            print(s)
            # generate an equal number of samples per site for all rounds
            data = sh.generate_random_gaussian(n*s, m)
            q, r = la.qr(data, mode='economic')
            data_list, choice = sh.partition_data_horizontally(data, s)
            ortho1, G_list = simulate_federated_qr(data_list, encrypt=False, filename=fn, split=s, repeat=i)
            ortho2, G_list = simulate_federated_qr(data_list, encrypt=True, filename=fn+'_encrypted', split=s, repeat=i)

            # Compute angles between two consecutive columns,
            # they should be 90
            oo1 = all_angles_against_all(ortho1)
            log_angles(filename1, oo1, 'unencrypted', i, s)


            oo2 = all_angles_against_all(ortho2)
            log_angles(filename2, oo2, 'encrypted', i, s)

            # compute angles between scipy linalg and fed
            ooq1 = co.compute_angles(q, ortho1)
            log_angles(filenameq1, ooq1, 'fed_vs_scipy', i, s)

            # compute angles between encrypted run and decrypted run
            oo12 = co.compute_angles(ortho2, ortho1)
            log_angles(filename12, oo12, 'encrypted_vs_uncrypted', i, s)


def log_angles(filename, angles, encrypted, repeat, splits):
    with open(filename, 'a') as handle:
        id = str(repeat)+'\t'+str(splits)+'\t'+ encrypted
        handle.write(cv.collapse_array_to_string(angles, id))

if __name__ == '__main__':

    dirname = '/home/anne/Documents/featurecloud/pca/vertical-pca/results/qr_encryption'
    filename = 'benchmark_encryption_qr'
    filename = path.join(dirname, filename)

    # data = sh.generate_random_gaussian(5000, 20)
    # q, r = la.qr(data, mode='economic')
    #
    # data_list, choice = sh.partition_data_horizontally(data, 2, randomize=False)
    # ortho1, G_list = simulate_federated_qr(data_list, encrypt=False)
    # co.compute_angles(ortho1, ortho1[:, 1:])
    # co.compute_angles(ortho1, q)
    benchmark_encryption(100, [2,3,5,10], 10000, 20, filename)