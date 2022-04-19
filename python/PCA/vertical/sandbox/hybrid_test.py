
####### HYBRID POWER ITERATION SCHEME #######
def hybrid_scheme(data_list, k, maxit, filename=None, filename2=None, scipy=None, choices=None, precomputed_pca=None, federated_qr=False):
    ug, eigenvalues, counter = better_hybrid_scheme(local_data=data_list, k=k, maxit=maxit, filename=filename, choices=choices, scipy=scipy, precomputed_pca=precomputed_pca, federated_qr=federated_qr)

    #ug, ev, eigenvalues = simulate_guo_benchmark(local_data=data_list, k=k, maxit=maxit, filename=filename, choices=choices, scipy=scipy, precomputed_pca=precomputed_pca, fractev=0.75)
    ug1 = ug.copy()
    #restart = assure_consecutive(eigenvalues)
    restart = counter
    print(restart)
    for i in range(restart, k - 1):
        next_eigenvec = simulate_guo(data_list, maxit=maxit, V_k=ug1[:, 0:i],
                                     filename=filename2, scipy=scipy, choices=choices, precomputed_pca=None, federated_qr = False, starting_vector= ug1[:, i])
        ug = np.concatenate((ug1[:, 0:i], next_eigenvec), axis=1)
    return ug





def better_hybrid_scheme(local_data, k, maxit, filename, scipy, choices, precomputed_pca=None, federated_qr=False):
    '''
       Simulate a federated run of principal component analysis using Guo et als algorithm in a modified version.

       Args:
           local_data: List of numpy arrays containing the data. The data has to be scaled already.
           k: The number of dimensions to retrieve
           maxit: Maximal number of iterations

       Returns: A column vector array containing the global eigenvectors

       '''
    G_list = []
    iterations = 0
    converged = False
    total_len = 0
    # generate an intitial  orthogonal noise matrix
    for d in local_data:
        total_len = total_len + d.shape[1]
    start = 0
    G_i = sh.generate_random_gaussian(total_len, k)
    G_i, R = la.qr(G_i, mode='economic')

    # send parts to local sites
    for i in range(len(local_data)):
        G_list.append(G_i[start:start + local_data[i].shape[1], :])
        log_transmission(filename, "G_i=SC", iterations, i, G_list[i])

        start = start + local_data[i].shape[1]

    H_i_prev = sh.generate_random_gaussian(local_data[0].shape[0], k)  # dummy init
    Gi_prev = G_i

    G_conv = None  # converged eigenvector matrix
    G_conv_list = []
    converged_counter = 0
    annormaly = [0] * k  # init counter for small eigengaps
    global_eigenvalues = []
    global_deltas = []
    while not converged and iterations < maxit and annormaly[converged_counter] < 10:
        iterations = iterations + 1
        H_i = np.zeros((local_data[0].shape[0], k - converged_counter))

        for i in range(len(local_data)):
            H_local = np.dot(local_data[i], G_list[i])
            log_transmission(filename, "H_local=CS", iterations, i, H_local)
            H_i = H_i + H_local

        log_transmission(filename, "H_global=SC", iterations, 1, H_i)

        converged, sum_of_delta, converged_eigenvals, delta = gv.convergence_checker(H_i, H_i_prev, return_converged=True)
        #print(delta)
        for i in range(len(G_list)):
            G_list[i] = np.dot(local_data[i].T, H_i) + G_list[i]

        G_i = np.concatenate(G_list, axis=0)
        gi_delta_obj = sh.eigenvector_convergence_checker(G_i, Gi_prev)

        # get the converged eigenvalues
        eigenvals = global_eigenvalues.copy()
        delta = np.append(global_deltas, delta)
        diff = []
        # iterate over non converged eigenvectors
        for col in range(G_i.shape[1]):
            eigenvals.append(np.sqrt(np.linalg.norm(G_i[:, col])))
        # compute eigengaps for all eigenvalues
        for e in range(1, len(eigenvals)):
            diff.append(np.log(np.abs(eigenvals[e - 1] - eigenvals[e])))
        # compute mean and standard deviation over all eigengaps
        mm = np.mean(diff)
        ssd = np.sqrt(np.var(diff))

        # flag eigenvectors that exhibit a very small eigengap
        for d in range(len(diff)):
            if diff[d] < mm - ssd:
                annormaly[d] = annormaly[d] + 1
        #print(annormaly)



        if federated_qr:
            temp = []
            if len(G_conv_list)!=0:
                for i in range(len(G_list)):
                    temp.append(np.concatenate([G_conv_list[i], G_list[i]], axis=1))
                G_i, G_list = qr.simulate_federated_qr(temp, encrypt=False)
            else:
                G_i, G_list = qr.simulate_federated_qr(G_list, encrypt=False)
        else:
            # Orthonormalise based on all eigenvectors
            if G_conv is not None and G_conv.shape[1] > 0:
                G_i = np.concatenate([G_conv, G_i], axis=1)
            G_i, R = la.qr(G_i, mode='economic')

        # current iteration led to some eigenvectors converging
        if len(converged_eigenvals) > 0:
            converged_counter = converged_counter + len(converged_eigenvals)
            # slice H_i and reset annomaly counter
            H_i_prev = H_i[:, len(converged_eigenvals):]
            annormaly = k * [0]
            # update converged eigenvector array
            global_eigenvalues = eigenvals[0:converged_counter]
            global_deltas = delta[0:converged_counter]
        else:
            H_i_prev = H_i

        log_current_accuracy(scipy=scipy, G_i=G_i, eigenvals=eigenvals, conv=delta, current_iteration=iterations, filename=filename, choices=choices, precomputed_pca=precomputed_pca, gi_delta_obj=gi_delta_obj)
        # in any case reslice G_i into converged and
        # unconverged
        #print(converged_counter)

        if federated_qr:
            G_conv_list = []
            for i in range(len(G_list)):
                G_conv_list.append(G_list[i][:, 0:converged_counter])
                G_list[i] = G_list[i][:, converged_counter:]

        else:
            G_conv = G_i[:, 0:converged_counter]
            G_i = G_i[:, converged_counter:]
            # redistribute the eigenvector parts
            start = 0
            for i in range(len(local_data)):
                G_list[i] = G_i[start:start + local_data[i].shape[1], :]
                log_transmission(filename, "G_i=SC", iterations, i, G_list[i])
                start = start + local_data[i].shape[1]


    #print(iterations)
    if federated_qr:
        temp = []
        if len(G_conv_list) != 0:
            for i in range(len(G_list)):
                temp.append(np.concatenate([G_conv_list[i], G_list[i]], axis=1))
        G_conv = np.concatenate(temp, axis=0)
    else:
        if G_i.shape[1] > 0:
            G_conv = np.concatenate([G_conv, G_i], axis= 1)
    return G_conv, global_eigenvalues, converged_counter


def assure_consecutive(arr):
    if len(arr) == 0:
        return -1
    i = 0
    while i < (len(arr) - 1) and arr[i] + 1 == arr[i + 1]:
        i = i + 1
    return i


####### END HYBRID POWER ITERATION SCHEME #######
