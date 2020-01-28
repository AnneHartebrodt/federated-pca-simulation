import pandas as pd
import math as math
import os as os
import shutil as sh
import time as time
import scipy.sparse.linalg as lsa
import scipy.linalg as la
from master import Distributed_DP_PCA
from PCA_runner import SimulationRunner
from outlier_removal import OutlierRemoval
from import_export.import_data import CustomDataImporter
import os.path as path
import numpy as np
import convenience as cv
import itertools as it

class Dropout():
    def __init__(self, file = '/tmp/'):
        self.ddppca = Distributed_DP_PCA(file = file)
        self.importer = CustomDataImporter()
        self.simulation = SimulationRunner()
        self.outlier = OutlierRemoval()


    def extract_eigenvals(self, E):
        '''
        Eigendecomposition from scipy.linalg.sparse returns eigenvalues ordered in
        increasing order, followed by eigenvalues which are 0.
        Eigenvalues are returned in decreasing order ommiting th 0s alltogether
        Args:
            E: Eigenvalues from a sparse singular value decomposition.

        Returns: Eigenvalue vector in decreasing order, without 0s.

        '''
        indz = np.where(E == 0)
        E = np.flip(E)
        E = E[E != 0]
        return E, indz



    def collapse_array_to_string(self, a, study_id, nr_dropped):
        res = study_id+'\t'+str(nr_dropped)+'\t'
        for e in a:
            res =  res + str(e)+'\t'
        res = res+'\n'
        return res


    def angle(self, v1, v2):
        """
        Calculates the angle between to n-dimensional vectors
        and returns the result in degree
        Args:
            v1: first vector
            v2: second vector

        Returns: angle in degree or NaN

        """
        dp = np.dot(v1,v2)
        norm_x = la.norm(v1)
        norm_y = la.norm(v2)
        theta = np.arccos(dp/(norm_x * norm_y))
        angle = math.degrees(theta)
        # angle can be Nan
        #print(angle)
        if math.isnan(angle):
            return angle
        # return the canonical angle
        if angle>90:
            return np.abs(angle-180)
        else:
            return angle

    def superset(self, set):
        la = []
        if len(set)<5:
            ll = list(it.product([0,1], repeat=len(set)))
        else:
            ll = []
            for i in range(20):
                ll.append(np.random.choice([0, 1], size=len(set)))
        for l in range(len(ll)):
            tmp = []
            for i in range(len(ll[l])):
                if ll[l][i] == 1:
                    tmp.append(set[i])
            if len(tmp)>1 and len(tmp)<len(set):
                la.append(tmp)
        for s in set:
            la.append([s])
        return la

    def run_dropout(self, data_scaled,  outdir, tempdir, nr_dropped, study_id, outliers=[], mode='outlier_free', ndims = 100,header=None, rownames=None, center=True, scale_var=True, scale01=False, scale_unit=False, transpose= False, sep='\t'):
        '''
        This function runs a 'leave-k'out' simulation. An eigenvalue decomposition
        is calculated n times on a datamatrix, leaving out nr_dropped samples at each time.
        The indices of the dropped samples, the eigenvalues and the eigenvectors are saved
        (up to 20).
        Args:
            data_scaled: The scaled data matrix
            nr_dropped: The number of samples to drop at each pass

        Returns: None, the according results are saved in files

        '''

        n = data_scaled.shape[0]

        # Outlier free mode means, we can estimate the influence of a sample
        # under ideal conditions
        if mode == 'outlier_free':
            if len(outliers) != 0:
                data_scaled = np.delete(data_scaled, outliers, 0)


        data_scaled_copy = np.copy(data_scaled)
        # the number of draws is 20% of the number of samples,
        # but a minimum of 10 samples, unless the dataframe is smaller
        # in which case very sample is dropped once.

        draws = max(int(np.ceil(n / 5)), min(10, n))

        if mode == 'cheat':
            # gene
            droprows = self.superset(outliers)
            nr_dropped = 1
        else:
            droprows = np.random.choice(n, size=draws, replace=False)
        index = 0

        #
        for row in range(0, len(droprows), nr_dropped):
            if mode == 'cheat':
                data_scaled = np.delete(data_scaled, droprows[row], 0)
            else:
                # remove the samples rows= patients from the data
                e = min(row+nr_dropped, data_scaled.shape[0])
                data_scaled = np.delete(data_scaled, droprows[row:e], 0)
            #calculate covariance matrix and svd
            #print('Computing Singular Value Decomposition')
            # It should be sufficient to calculate at most n PCs, because the others will contain noise
            nd = min(min(n - nr_dropped, data_scaled.shape[0]-1), ndims)
            pca, U, E = self.simulation.run_standalone(data_scaled, outdir, dims=ndims, header=header, rownames=rownames,
                                                         center=center, scale_var=scale_var, scale01=scale01,
                                                         scale_unit=scale_unit, transpose=transpose, sep=sep,
                                                         filename='/pca.after_outlier_removal', log=True,
                                                         exp_var=1.0)
            E= E[0:20]
            U = U[:, 0:20]

            # print the 20 first eigenvectors
            res_eigen = self.collapse_array_to_string(E, study_id, nr_dropped)
            # save the dropped patients, to trace back big changes later
            res_dropped = self.collapse_array_to_string(droprows[row:(row+nr_dropped)],study_id, nr_dropped)

            # write eigenvalues, dropped samples and eigenvectors to file
            with open(outdir+'/eigenvalues.tsv', 'a+') as handle:
                handle.write(res_eigen)
            with open(outdir+'/dropped_samples.tsv', 'a+') as handle:
                handle.write(res_dropped)
            # save the eigenvectors in a separate directory. They likely will be deleted after
            # calculating the angles due to missing diskspace
            pd.DataFrame(U).to_csv(tempdir+'/eigenvectors'+str(index)+'.tsv', sep='\t',header=None, index=None)
            # increment the index
            index = index+1
            # get the correct data back.
            data_scaled = np.copy(data_scaled_copy)



    def run_standalone(self, data, outdir=None, dims=100, header=None, rownames=None, center=True, scale_var=True, scale01=False, scale_unit=False,transpose=False, sep = '\t', reported_angles = 20, exp_var = 1, study_id='ID1'):


        pca, W1, E1 = self.simulation.run_standalone(data, outdir, dims=dims, header=header, rownames=rownames,
                                                     center=center, scale_var=scale_var, scale01=scale01,
                                                     scale_unit=scale_unit, transpose=transpose, sep=sep,
                                                     filename='/pca.before_outlier_removal', log=True, exp_var=exp_var)

        res_eigen = self.collapse_array_to_string(E1[0:reported_angles],study_id, '0')
        with open(outdir + '/eigenvalues_with_outliers.tsv', 'a+') as handle:
            handle.write(res_eigen)
        #pd.DataFrame(W1[:,0:20]).to_csv(ev_path + '/eigenvectors_with_outlier.tsv', sep='\t', header=None, index=None)


        print('Logging outliers')
        outliers = self.outlier.outlier_removal_mad(pca, 6, 3)
        print(outliers)
        with open(outdir + '/removed_outliers.tsv', 'a+') as handle:
            handle.write(cv.collapse_array_to_string(outliers, study_id))

        # delete outliers manually
        if len(outliers) != 0:
            data = np.delete(data, outliers, 0)

        print('Standalone PCA after outlier removal')
        pca, W1, E1 = self.simulation.run_standalone(data, outdir, dims=dims, header=header, rownames=rownames,
                                                     center=center, scale_var=scale_var, scale01=scale01,
                                                     scale_unit=scale_unit, transpose=transpose, sep=sep,
                                                     filename='/pca.after_outlier_removal', log=True,
                                                     exp_var=exp_var)

        res_eigen = self.collapse_array_to_string(E1[0:reported_angles],study_id, '0')

        # write eigenvalues, dropped samples and eigenvectors to file
        with open(outdir + '/eigenvalues_reference.tsv', 'a+') as handle:
            handle.write(res_eigen)
        pd.DataFrame(W1[:,0:reported_angles]).to_csv(outdir + '/eigenvectors_reference.tsv', sep='\t', header=None, index=None)
        return outliers

    def run_study(self, datafile, outdir, tempdir, header = None, rownames = None, sep = '\t', dims = 100, mode = 'outlier_free', nr_dropped=1):
        # ID is the folder name
        study_id = path.basename(path.dirname(datafile))

        data, varnames, sampleids = self.importer.data_import(datafile, header=header, rownames=rownames, outfile=outdir,transpose=False, sep=sep)
        dims = min(dims, data.shape[0])
        outliers = self.run_standalone(data, outdir = outdir, study_id=study_id)
        self.run_dropout(data, outdir, tempdir, nr_dropped=nr_dropped, study_id=study_id, outliers=outliers, mode=mode, ndims=dims)
        self.compute_angles(tempdir, outdir, nr_dropped)



    def compute_angles(self, tempdir, outdir, nr_dropped):
        '''
        computes the angle in degrees for all the input files
        in the given directory vs the canoncial pca. Angles are between the matching vectors
        (columnwise)

        Args:
            eigenvector_path: folder containing files, which contain
            columnwise eigenvectors in decreaseing order according to the
            eigenvector
            other_folder: a different folder, where the outputs are saved.
            Make sure it is a different folder in case the eigenvector folder
            is deleted.

        Returns:

        '''

        # read reference
        reference = pd.read_csv(filepath_or_buffer=outdir + '/eigenvectors_reference.tsv', sep='\t', header=None, engine='python')
        angles = []
        i = 1
        for file in os.listdir(tempdir):
            #read one of the dropped out eigenvalue files
            current = pd.read_csv(filepath_or_buffer=tempdir+'/'+file, sep='\t', header = None, engine='python')
            # calculate angles against reference
            for k in range(max(current.shape[1], reference.shape[1])):
                ll = [i, 0, k, self.angle(current.iloc[:,k], reference.iloc[:,k])]
                angles.append(ll)
            i = i+1
        pd.DataFrame(angles).to_csv(outdir+'/angles_dropout'+str(nr_dropped)+'.tsv', sep = '\t', header = None, index=None)

    def make_eigenvector_path(self, inputfile, foldername):
        """
        creates a folder called eigenvectors in the input directory
        if the eigenvector folder already exists a folder named
        'eigenvectors<currentTimestamp> will be created
        Args:
            inputfile: name of the inputfile

        Returns: a pathname according to the stated rules

        """
        print(path.dirname(inputfile))
        if not os.path.exists(path.dirname(inputfile)+'/'+foldername):
            pn = path.dirname(inputfile)+'/'+foldername
            os.makedirs(pn)
        else:
            print('Path exists')
            pn = path.dirname(inputfile)+'/'+foldername+str(time.clock())
            os.makedirs(pn)
        return pn

    def delete(self, eigenvector_path):
        '''
        CAREFUL with this.
        To save disk space, the eigenvectors will be remove after the angles have been computed
        Args:
            eigenvector_path: directory which will be removed including all subdirectories.

        Returns:

        '''
        sh.rmtree(eigenvector_path)


if __name__=="__main__":
    print('run eigengap script')

    parser = ap.ArgumentParser(description='Eigengap calculation')
    parser.add_argument('-f', metavar='file', type=str, help='filename of data file; file should be tab separated')
    parser.add_argument('-o', metavar='outfile', type=str, help='output file')
    parser.add_argument('-d', metavar='dims', type=int, help='field delimiter', default=100)
    parser.add_argument('-s', metavar='sep', type=str, help='field delimiter')
    parser.add_argument('-n', metavar='nr_dropped', type=int, help='number of rows dropped')
    parser.add_argument('-m', metavar='mode', type=str, help='cheat or outflier_free')
    args = parser.parse_args()

    inputfile = args.f
    outfile = args.o
    dims = args.d
    sep = args.s
    nr_dropped = args.n
    mode = args.m

    #inputfile ='/home/anne/Documents/featurecloud/data/tcga/data_clean/BEATAML1/coding_trunc.tsv'
    #outfile = '/home/anne/Documents/featurecloud/results/gexp_stats/leaveoone1out/'
    #dims = 100
    #sep = ','
    #nr_dropped = 1
    #mode = 'outlier_free'

    d = Dropout()

    dname = path.basename(path.dirname(inputfile))+'/'+mode+'_'+str(nr_dropped)+'/'
    summaryfile = d.make_eigenvector_path(outfile, dname)
    ev_path  = d.make_eigenvector_path(summaryfile, 'eigenvectors')
    d.run_study(inputfile, summaryfile, ev_path, sep=sep, header = 0, dims=dims, mode=mode)





