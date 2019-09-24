import PCA.master
import scipy as sc
import pandas as pd
import PCA.simulation_runner as simul

ddpca = PCA.master.Distributed_DP_PCA()
############IRIS#############
original, labels = ddpca.prepare_data()
orig = original

# standalone PCA
pca, W, s, d = simul.perform_standalone_pca(orig, ndims=4, dppca=ddpca, center=True, scale=True, scale01=True)
# #pca = None
pca = pd.DataFrame(pca)
pca.to_csv('/home/anne/Documents/featurecloud/featurecloud-test-data/PCA01/IRIS_standalonePCA.tsv', sep='\t', header=None,
           index=False)

o = sc.copy(orig)
scaled = ddpca.scale_data(orig, center=True, scale_variance=True, scale01=True)
# run a 5 split without noise
for i in range(1, 6):
    sim = ddpca.simulate_multisite_PCA(orig, i, noise=False, ndims=4)
    print('A' + str(i))
    print(sim)
    projection = sc.dot(scaled, sim[:, 0:2])
    print(scaled[1])
    print(projection[1])
    projection = pd.DataFrame(projection)
    projection.to_csv(
        '/home/anne/Documents/featurecloud/featurecloud-test-data/PCA01/IRIS_multi_no_noise_' + str(
            i) + '.tsv', sep='\t', header=None, index=False)
    sim = pd.DataFrame(sim)
    sim.to_csv(
        '/home/anne/Documents/featurecloud/featurecloud-test-data/PCA01/IRIS_multi_no_noise_' + str(i) + 'eigen.tsv',
        sep='\t', header=None, index=False)

original, labels = ddpca.prepare_data()
orig = original
scaled = ddpca.scale_data(orig, center=True, scale_variance=True, scale01=True)
for i in range(1, 6):
    sim2 = ddpca.simulate_multisite_PCA(orig, i, noise=False, ndims=4)
    print('B' + str(i))
    print(sim2)
    projection = sc.dot(scaled, sim2[:, 0:2])
    print(scaled[1])
    print(projection[1])

original, labels = ddpca.prepare_data()
orig = original
scaled = ddpca.scale_data(orig, center=True, scale_variance=True, scale01=True)
for i in range(1, 6):
    print('C' + str(i))
    sim3 = ddpca.simulate_multisite_PCA(orig, i, noise=False, ndims=4)
    print(sim3)
    projection = sc.dot(scaled, sim3[:, 0:2])
    print(scaled[1])
    print(projection[1])

# multiple simulations with noise
print('middle')
orig = original

dm = 4
results = sc.empty(shape=(0, dm + 4))
scaled1 = ddpca.scale_data(orig, center=True, scale_variance=True, scale01=True)

res = ddpca.simulate_multisite_PCA(orig, 1, noise=False, ndims=4)
res2 = ddpca.simulate_multisite_PCA(orig, 1, noise=False, ndims=4)
for split in range(1, 6):
    res = ddpca.simulate_multisite_PCA(orig, split, noise=False, ndims=4)
    print(res)
    proj = sc.dot(scaled1, res[:, 0:dm])
scaled1 = pd.DataFrame(scaled1)
scaled1.to_csv('/home/anne/Documents/featurecloud/featurecloud-test-data/PCA01/IRIS_Scaled.tsv',
               sep='\t', header=None, index=False)

original, labels = ddpca.prepare_data()
orig = original
epsilons = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 100]
deltas = [0.01]
scaled = ddpca.scale_data(orig, center=True, scale_variance=True, scale01=True)
res = simul.run_multiple_simulations(orig, dims=4, dppca=ddpca, noise=True, nrSamples=10, splits=5, epsilons=epsilons,
                               deltas=deltas)
res = pd.DataFrame(res)
res.to_csv('/home/anne/Documents/featurecloud/featurecloud-test-data/PCA01/IRIS_simulation_multiple.tsv',
           sep='\t', header=None, index=False)
labels = pd.DataFrame(labels)
labels.to_csv('/home/anne/Documents/featurecloud/featurecloud-test-data/PCA01/IRIS_labels.tsv',
           sep='\t', header=None, index=False)