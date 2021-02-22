import python.PCA.convenience as cv
import python.PCA.comparison as co
import numpy as np

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


def log_current_accuracy(u, G_i, eigenvals, conv, current_iteration, filename, choices, precomputed_pca=None,
                         current_ev=None, gi_delta_obj=None, v=None, H_i=None):
    """
    Log the current iterations angle to the canonical
    Args:
        u: column vector matrix with canonical eigenvectors
        G_i: column vector based matrix with current eigenvector estimation
        current_iteration: iteration index
        filename: output filename prefix > out will be saved to out.angles, and out.cor

    Returns: None

    """
    if current_ev is not None:
        info_string = str(current_ev) + '\t' + str(current_iteration)
    else:
        info_string = str(current_iteration)

    if not u is None:
        with open(filename + '.angles.u', 'a+') as handle:
            angles = co.compute_angles(u[choices, :], G_i)
            if angles is not None and len(angles) > 0:
                info = cv.collapse_array_to_string(angles, info_string)
                handle.write(info)

        with open(filename + '.cor', 'a+') as handle:
            correlations = co.compute_correlations(u[choices, :], G_i)
            if correlations is not None and len(correlations) > 0:
                info = cv.collapse_array_to_string(correlations, info_string)
                handle.write(info)

    if not v is None:
        with open(filename + '.angles.v', 'a+') as handle:
            angles = co.compute_angles(v, H_i)
            if angles is not None and len(angles) > 0:
                info = cv.collapse_array_to_string(angles, info_string)
                handle.write(info)

    with open(filename + '.eigenval', 'a+') as handle:
        info = cv.collapse_array_to_string(eigenvals, info_string)
        handle.write(info)

    if conv is not None:
        with open(filename + '.conv', 'a+') as handle:
            conv = cv.collapse_array_to_string(conv, info_string)
            handle.write(conv)

    if precomputed_pca is not None:
        with open(filename + '.angles_precomp', 'a+') as handle:
            angles = co.compute_angles(precomputed_pca[choices, :], G_i)
            info = cv.collapse_array_to_string(angles, info_string)
            handle.write(info)

        with open(filename + '.cor_precomp', 'a+') as handle:
            correlations = co.compute_correlations(precomputed_pca[choices, :], G_i)
            info = cv.collapse_array_to_string(correlations, info_string)
            handle.write(info)
    if gi_delta_obj is not None:
        with open(filename + '.eigenvector_convergence', 'a+') as handle:
            conv = cv.collapse_array_to_string(gi_delta_obj, info_string)
            handle.write(conv)

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


def log_transmission(logfile, log_entry_string, iterations, counter, element, eigenvector=10):

    """
    Dumps object to json to estimate the size.
    Args:
        logfile:
        log_entry_string:
        iterations:
        counter:
        element:
        eigenvector: Dummy value 10 if not put

    Returns:

    """
    with open(logfile + '.transmission', 'a+') as handle:

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
        handle.write(log_entry_string + '\t' + str(iterations)+ '\t' + str(counter) + '\t' + str(eigenvector) + '\t' + str(size) + '\n')


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
