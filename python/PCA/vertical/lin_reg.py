from sklearn import datasets
import pandas as pd
import python.PCA.vertical.simulate_federated_qr_orthonormalisation as fed
import numpy as np
import scipy.linalg as la
import python.PCA.shared_functions as sh
import python.PCA.comparison as co
import python.PCA.convenience as cv
import scipy.stats as stats


def simulate_linear_regression(data_list, label_list):
    # Step 1: Compute QR factorization
    ortho, G_list, r2, rl = fed.simulate_federated_qr(data_list)

    # Step 2 Compute y (Client side)
    ylist = [g.T @ l for g,l in zip(G_list, label_list)]

    # Aggregate y (Aggregator side)
    y = np.sum(ylist, axis=0)
    # Solve x (Aggregator side)
    x = la.solve_triangular(r2, y)

    # Step 3:
    # Compute squared sum of predicted labels (client side)
    mss = [np.sum((d @ x)**2) for d, y in zip(data_list, label_list)]
    # Aggregate
    mss = np.sum(mss, axis=0)

    # Compute sum of squared residuals (client)
    residuals = [(d @ x)-y for d, y in zip(data_list, label_list)]
    sses = [np.sum(np.square(r)) for r in residuals]
    # Aggregate
    sses = np.sum(sses)

    # Compute R squared
    # df = degree of freedom = # samples - # features -1
    df = np.sum([d.shape[0] for d in data_list]) - r2.shape[0] +1
    print(df)
    rsq = np.sum(mss)/(sses+mss)
    rsq

    # invert covariance matric
    rinv  =np.linalg.inv(r2.T @ r2)
    # Compute T statistic
    sigma = sses/(df)*rinv
    se= np.sqrt(np.diag(sigma))
    T = x/se

    # find p-value
    pval = (stats.t.sf(np.abs(T), df))*2

    return x, rsq, pval



if __name__ == '__main__':

    # Diabetes
    diabetes_X, diabetes_y = datasets.load_diabetes(return_X_y=True)

    pd.DataFrame(diabetes_X).to_csv('/home/anne/Documents/featurecloud/pca/qr/data/diabetes/diabetes.tsv', header=False, index=False)
    pd.DataFrame(diabetes_y).to_csv('/home/anne/Documents/featurecloud/pca/qr/data/diabetes/diabetes_lab.tsv', header=False, index=False)

    # Use only one feature
    diabetes_X = diabetes_X

    # Split the data into training/testing sets
    diabetes_X_train = diabetes_X[:-20]
    diabetes_X_test = diabetes_X[-20:]
    # Split the targets into training/testing sets
    diabetes_y_train = diabetes_y[:-20]
    diabetes_y_test = diabetes_y[-20:]

    data_list, choice = sh.partition_data_horizontally(diabetes_X_train, 3)
    q, r = la.qr(diabetes_X_train, mode='economic')
    ortho, G_list, r2, rl = fed.simulate_federated_qr(data_list)
    print(np.linalg.norm(np.abs(q) - np.abs(ortho)))
    print(np.linalg.norm(np.abs(r) - np.abs(r2)))

    interval = 0
    label_list = []
    for d in range(len(data_list)):
        label_list.append(diabetes_y_train[interval:interval+data_list[d].shape[0]])
        interval = interval+data_list[d].shape[0]
    x, rsquared, pval = simulate_linear_regression(data_list, label_list)


    with open('/home/anne/Documents/featurecloud/pca/qr/results/diabetes.txt', 'w') as handle:
        handle.write(cv.collapse_array_to_string(x,'Coefficients:'))
        handle.write('r-squared:\t'+ str(rsquared)+'\n')
        handle.write(cv.collapse_array_to_string(pval,'p-values:'))

    # Fish market
    fish = pd.read_csv('/home/anne/Documents/featurecloud/pca/qr/data/fish/Fish.csv')
    fish = fish.values
    weight = fish[:,1]
    weight = weight.astype(float)
    fish = fish[:,2:]
    fish = fish.astype(float)
    data_list, choice = sh.partition_data_horizontally(fish, 3)
    q, r = la.qr(fish, mode='economic')
    ortho, G_list, r2, rl = fed.simulate_federated_qr(data_list)
    print(np.linalg.norm(np.abs(q) - np.abs(ortho)))
    print(np.linalg.norm(np.abs(r) - np.abs(r2)))

    interval = 0
    label_list = []
    for d in range(len(data_list)):
        label_list.append(weight[interval:interval + data_list[d].shape[0]])
        interval = interval + data_list[d].shape[0]
    x, rsquared, pval = simulate_linear_regression(data_list, label_list)
    print(x)
    print(rsquared)
    print(pval)
    with open('/home/anne/Documents/featurecloud/pca/qr/results/fish.txt', 'w') as handle:
        handle.write(cv.collapse_array_to_string(x,'Coefficients:'))
        handle.write('r-squared:\t'+ str(rsquared)+'\n')
        handle.write(cv.collapse_array_to_string(pval,'p-values:'))

    ### WHO
    print('WHO')
    who = pd.read_csv('/home/anne/Documents/featurecloud/pca/qr/data/who/who.csv')
    who = who.interpolate()
    who.to_csv('/home/anne/Documents/featurecloud/pca/qr/data/who/who2.csv', index=False)
    who = who[(who['Year']==2012) & (who['Status']=='Developing')]
    outcome = who['Life expectancy ']
    who = who.iloc[:, 4:]

    who = who.values
    who = who.astype(float)

    # mean impute

    outcome =outcome.values
    outcome = outcome.astype(float)

    data_list, choice = sh.partition_data_horizontally(who, 3)
    q, r = la.qr(who, mode='economic')
    ortho, G_list, r2, rl = fed.simulate_federated_qr(data_list)
    print(np.linalg.norm(np.abs(q) - np.abs(ortho)))
    print(np.linalg.norm(np.abs(r) - np.abs(r2)))

    interval = 0
    label_list = []
    for d in range(len(data_list)):
        label_list.append(outcome[interval:interval + data_list[d].shape[0]])
        interval = interval + data_list[d].shape[0]
    x, rsquared, pval = simulate_linear_regression(data_list, label_list)
    print(x)
    print(rsquared)
    print(pval)

    with open('/home/anne/Documents/featurecloud/pca/qr/results/who.txt', 'w') as handle:
        handle.write(cv.collapse_array_to_string(x,'Coefficients:' ))
        handle.write('r-squared:\t'+ str(rsquared)+'\n')
        handle.write(cv.collapse_array_to_string(pval,'p-values:'))



