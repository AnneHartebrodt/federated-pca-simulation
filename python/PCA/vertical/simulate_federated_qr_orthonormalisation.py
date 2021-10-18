"""
    Copyright (C) 2020 Anne Hartebrodt

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    Authors: Anne Hartebrodt

"""

##########################################################################################
#
#    If you just want some pseudo code for federated QR orthonormalisation,
#    this here is your script. No unnecessary logging or other stuff.
#
##########################################################################################


import numpy as np
import scipy.linalg as la
import python.PCA.shared_functions as sh
import python.PCA.comparison as co



def simulate_federated_qr(local_data):
    """
    Simulate federated orthonormalisation of a horizontally split
    data matrix.
    Args:
        local_data: Local datasets formatted as numpy arrays.

    Returns: The orthogonal vectors, as concatenated matrix and individual snippets.

    """

    alist = []
    ortho = []
    sum = 1e-16

    # Compute first squared eigenvector norm
    for d in range(len(local_data)):
        se = np.dot(local_data[d][:, 0], local_data[d][:,0])
        sum = sum+se
        alist.append(local_data[d][:,0])
    ortho.append(alist)


    # ortho [eigenvector rank] [data set]
    # list of lists containing the already
    #  orthogonal eigenvectors
    norms = [sum]
    # iterate over the eigenvectors
    for i in range(1,local_data[0].shape[1]):
        # conorms we want to calculate
        sums = []
        aplist = []
        # iterate over the all already orthonormal eigenvector snippet
        # lists and the associated norms
        # decrypt norms
        for di in range(len(ortho)):
            # local eigenvector snippets
            o = ortho[di]
            # eigenvector norms
            nn = norms[di]
            n = nn + 1e-16
            sum = 0
            # iterate over the local data sets
            # combined with the local eigenvector snippets
            # o takes all othonormal ranks
            # i is the currently to be orthonomalised rank
            for ik in range(len(local_data)):
                d = local_data[ik]
                o1 = o[ik]
                # Compute conorm
                se = np.dot(o1, d[:, i]) / n
                sum = sum + se
            sums.append(sum)

        # init norm of the current vector
        norm = 0
        for d in range(len(local_data)):
            # ap = newly reorthogonalised eigenvector snippet
            ap = local_data[d][:, i]
            for j in range(len(sums)):
                # reorthonogonalise
                ap = ap - sums[j] * ortho[j][d]

            # compute the local norm of the freshly orthogonalised
            # eigenvector snippet
            se = np.dot(ap, ap)
            norm = norm + se
            aplist.append(ap)
        norms.append(norm+1e-15)
        ortho.append(aplist)

    G_list = []

    # normalise the vector norms to unit norm.
    for d in range(len(local_data)):
        oo = []
        for i in range(len(ortho)):
            # norms are still squared
            print(i)
            oo.append(ortho[i][d] / np.sqrt(norms[i]))
        oo = np.stack(oo, axis = 1)
        G_list.append(oo)

    # just for convenience stack the data
    ortho = np.concatenate(G_list, axis=0)
    return ortho, G_list


if __name__ == '__main__':


    data = sh.generate_random_gaussian(50000, 10)
    q, r = la.qr(data, mode='economic')
    data_list, choice = sh.partition_data_horizontally(data, 3)
    ortho, G_list = simulate_federated_qr(data_list)
    angles = co.compute_angles(q, ortho)
    print(angles)