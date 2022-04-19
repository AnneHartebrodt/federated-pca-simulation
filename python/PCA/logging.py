import python.PCA.convenience as cv
import python.PCA.comparison as co
import numpy as np
import os.path as op
import time
import codetiming as ct
import pandas as pd


# Extending Timer to keep track of the total time
class TimerCounter(ct.Timer):
    def __init__(self):
        self._total_timer = 0

    def total(self) -> float:
        return self._total_timer

    def stop(self):
        """Stop the timer, and report the elapsed time"""
        if self._start_time is None:
            raise ct.TimerError(f"Timer is not running. Use .start() to start it")

        elapsed_time = time.perf_counter() - self._start_time
        self.last = elapsed_time
        self._start_time = None
        self._total_timer = self._total_timer + self.last
        return elapsed_time

class FileNotOpenError(Exception):
    """ File Not open"""

class TransmissionLogger:
    def __init__(self):
        self.logs = []

    def open(self, filename):
        self.is_open = True
        self.filename = filename

    def close(self):
        if not self.is_open:
            raise FileNotOpenError
        logs = pd.DataFrame(self.logs)
        if self.filename is not None:
            logs.to_csv(self.filename + '.transmission', header=False, sep='\t', mode='a')

    def log_transmission(self, key, iterations, site, element, eigenvector=1):
        size = self.determine_size(element)
        self.logs.append([key, iterations, site, size, eigenvector])

    @staticmethod
    def determine_size(element):
        if isinstance(element, np.ndarray):
            size = len(element.flatten())
        elif isinstance(element, float) or isinstance(element, int):
            size = 1
        elif isinstance(element, list):
            # cast to numpy array and flatten
            # in case nested list.
            size = len(np.asanyarray(element).flatten())
        else:
            # in case something goes wrong
            size = -1
        return size

####### LOGGING FUNCTIONS #######
def create_filename(outdir, dataset_name, splits, counter, k, maxit, time):
    """
    Make a file name
    Args:
        dataset_name: Name of the dataset currently running
        splits: The current split
        counter: Current repeat experiment
        k: Number of targeted dimensions
        maxit: maximal iterations

    Returns: A concatenated filename with filestamp

    """
    fn = dataset_name + '_' + str(splits) + '_' + str(counter) + '_' + str(k) + '_' + str(maxit) + '_' + str(time)
    fn = op.join(outdir, fn)
    return fn

class AccuracyLogger:
    def __init__(self):
        self.angles_u = []
        self.mev_u = []
        self.cor_u =[]
        self.angles_v = []
        self.mev_v = []
        self.cor_v = []
        self.eigenvals_list = []
        self.conv_list = []
        self.prec_angles = []

    def open(self, filename):
        self.filename = filename

    def info(self, current_ev, current_iteration):
        if current_ev is not None:
            info = [current_ev, current_iteration]
        else:
            info = [current_iteration]
        return info

    def log_u(self, u,G_i, choices, info):
        if not u is None and not G_i is None:
            self.angles_u.append(info + co.compute_angles(u[choices, :], G_i))
            self.mev_u.append(info + [co.mev(u[choices, :], G_i)])
            self.cor_u.append(info +co.compute_correlations(u[choices, :], G_i))

    def log_v(self, v, H_i, info):
        if not v is None and not H_i is None:
            self.angles_v.append(info + co.compute_angles(v, H_i))
            self.mev_v.append(info + [co.mev(v, H_i)])
            self.cor_v.append(info +co.compute_correlations(v, H_i))


    def log_current_accuracy(self, u, G_i, current_iteration, choices=None, precomputed_pca=None,
                             current_ev=None,  v=None, H_i=None, eigenvals=None, conv=None,):
        """
        Log the current iterations angle to the canonical
        Args:
            u: column vector matrix with canonical eigenvectors
            G_i: column vector based matrix with current eigenvector estimation
            current_iteration: iteration index
            filename: output filename prefix > out will be saved to out.angles, and out.cor

        Returns: None

        """

        info = self.info(current_ev, current_iteration)
        self.log_u(u=u,G_i=G_i, choices=choices, info=info )
        self.log_v(v=v, H_i=H_i, info=info)
        self.log_conv(conv, info)
        self.log_eigenvals(eigenvals, info)
        self.log_precomputed(precomputed_pca, G_i, info, choices)


    def log_eigenvals(self, eigenvals, info):
        if eigenvals is not None:
            self.eigenvals_list.append(info+eigenvals)

    def log_conv(self, conv, info):
        if conv is not None:
            self.conv_list.append(info+conv)


    def log_precomputed(self, precomputed_pca, G_i, info, choices):
        if precomputed_pca is not None:
            self.prec_angles.append(info+co.compute_angles(precomputed_pca[choices, :], G_i))

    def write_to_pandas(self, log_list, suffix):
        if self.filename is not None:
            if len(log_list) > 0:
                pd.DataFrame(log_list).to_csv(self.filename + suffix, mode='a', sep='\t', header=False, index=False)


    def close(self):
        self.write_to_pandas(self.angles_u, '.angles.u')
        self.write_to_pandas(self.angles_v, '.angles.v')
        self.write_to_pandas(self.mev_u, '.mev.u')
        self.write_to_pandas(self.mev_v, '.mev.v')
        self.write_to_pandas(self.eigenvals_list, '.eigenvals')
        self.write_to_pandas(self.conv_list, '.conv')
        self.write_to_pandas(self.prec_angles, '.prec_angles.u')




def log_choices(logfile, filename, choices):
    """
    Log the permutation of the data sets.
    Args:
        logfile: Name of the log file
        filename: Filename of the result file, this permutation belongs to
        choices: the actual choice array

    Returns: None

    """
    with open(logfile, 'a+') as handle:
        handle.write(cv.collapse_array_to_string(choices, filename))








def log_time(logfile, algorithm, time, split, repeat):
    """
    Log the permutation of the data sets.
    Args:
        logfile: Name of the log file
        filename: Filename of the result file, this permutation belongs to
        choices: the actual choice array

    Returns: None

    """
    with open(logfile, 'a+') as handle:
        handle.write(algorithm + '\t' + str(split) + '\t' + str(repeat) + '\t' + str(time) + '\n')

def log_time_keywords(logfile, keyword, time ):
    """
    Log the permutation of the data sets.
    Args:
        logfile: Name of the log file
        filename: Filename of the result file, this permutation belongs to
        choices: the actual choice array

    Returns: None

    """
    with open(logfile + '.eo', 'a+') as handle:
        handle.write(keyword + '\t' + str(time) + '\n')
