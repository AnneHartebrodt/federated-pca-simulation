import pandas as pd
import scipy.sparse.linalg as lsa
from python.PCA.vertical.vertical_pca_benchmark import *



if __name__ == '__main__':
    # this is here to avoid circular import
    from python.PCA.horizontal.horizontal_pca_benchmark import read_presplit_data_folders, compute_canonical, scale_datasets


    # MNIST for reference
    data, test_lables = mi.load_mnist('/home/anne/Documents/featurecloud/pca/vertical-pca/data/mnist/raw', 'train')
    # data, test_labels = mi.load_mnist(input_dir, 'train')
    data = coo_matrix.asfptype(data)

    #ds = si.scale_center_data_columnwise(data, center=True, scale_variance=True)

    data_list, choice = sh.partition_data_vertically(data, splits=2, randomize=False)
    #data_list = scale_datasets(data_list)
    u,s,v = lsa.svds(np.concatenate(data_list,axis=1))
    u = np.flip(u)
    v = np.flip(v.T, axis=1)
    data_list = [d for d in data_list]

    print('randomized')
    uu = run_randomized(data_list, 10,10, 500, use_approximate=False)
    print(co.compute_angles(uu, v))
    uu, r = la.qr(uu)
    print(co.mev(uu,v))

    print('approx+randomized')
    uuuu = run_randomized_2(data_list, 10,10, 500, use_approximate=True)
    print(co.compute_angles(uuuu, v))
    uuuu, r = la.qr(uu)
    print(co.mev(uuuu, v))

    print('randomized+approx')
    uuu = run_randomized(data_list, 10,10, 500, use_approximate=True)
    print(co.compute_angles(uuu, v))
    uuu, r = la.qr(uuu)
    print(co.mev(uuu,v))


    print('exact')
    a,b,c,d,e,f = simulate_subspace_iteration(data_list,10, 10000 )
    print(co.compute_angles(v,a))
    a, r = la.qr(a,mode='economic')
    print(co.mev(v,a))


    # G_i, eigenvals, H_i, iterations, H_stack, G_list = simulate_subspace_iteration(data_list, 2 * 10, 10)
    # H_stack = np.concatenate(H_stack, axis=1)
    # H, S, G = lsa.svds(H_stack, k=199)
    # H = np.flip(H)
    # pp = np.dot(H.T, data)
    # up, sp, vp = lsa.svds(pp)
    # co.compute_angles(np.flip(vp.T, axis=1), v)
    #
    # co.mev(np.flip(vp.T, axis=1), v)
    # a,b,c,d,e,f = simulate_subspace_iteration(data_list, 10, 100000)
    # co.compute_angles(u, a)
    # ff = simulate_tt_pca(data_list)
    # # sample_count = [d.shape[0] for d in data_list]
    # # total_samples = sum(sample_count)
    # # weights = [sc/total_samples for sc in sample_count]
    # # u,s,v = compute_canonical(data_list, k=20)
    # uu, ss, vv = lsa.svds(np.concatenate(data_list, axis=0))
    # vv = np.flip(vv.T, axis=1)
    # x, e = simulate_federated_horizontal_pca(data_list, k=10, factor_k=3)
    # co.compute_angles(vv, x)
    # ff = simulate_tt_pca(data_list)
    # co.compute_angles(ff, x)
    # x2, e = simulate_federated_horizontal_pca(data_list, k=20, factor_k=2)
    # v1 = v
    # v1[0,0]=-8.36060474094277e-05
    # v1[5,0] =-0.00056577135205741441
    # co.compute_angles(x, v)[12]
    # np.sum(x[:,12]-v[:,12])/717
    # co.compute_angles(x, v)
    # print(co.subspace_reconstruction_error(np.concatenate(data_list, axis=0), v))
    # #print(co.subspace_reconstruction_error(np.concatenate(data_list, axis=0), vv))
    # print(co.subspace_reconstruction_error(np.concatenate(data_list, axis=0), x))
    #

