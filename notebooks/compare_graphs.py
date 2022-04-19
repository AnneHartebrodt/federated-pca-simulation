
#%%

import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc
import numpy as np
import sklearn as sk
import sklearn.metrics.cluster as m
#%%

def find_cluster_matching(true_labels, predictions):
    """
    Compute a cluster matching between the two labellings.
    This is done by computing a contingency matrix (computing the co-occurence of label
    pairs) and choosing the argmax for each label/
    :param true_labels: The gold standard labels
    :param predictions: The predictions
    :return: A matching beween the true labels and the predicted labels in form of a
    dictionary {true_label: predicted_label, ...}
    """
    # true_labels = [1,1,2,2,3,4,5,5]
    # predictions = [2,2,1,1,4,3,6,6]

    rownames = np.unique(true_labels)
    colnames = np.unique(predictions)
    cont = m.contingency_matrix(true_labels, predictions)
    # print(cont)
    matching = []
    for row in range(cont.shape[0]):
        matching.append(np.argmax(cont[row, :]))

    label_match = {}
    for a, b in zip(rownames, colnames[matching]):
        label_match[b] = a

    # Fallback if there is not an optimal cluster correspondence
    for c in colnames:
        if c not in label_match.keys():
            label_match[c] = -1

    # print(colnames)
    # print(colnames[matching])
    #print(label_match)
    new_prediction_labels = [label_match[b] for b in predictions]

    return label_match, new_prediction_labels


def confusion_matrix(true_labels, predictions):
    """
    Compute the confusion matrix for each class
    :param true_labels: Gold standard lables
    :param predictions: predicted lables
    :param labels: matching between the class names
    :return: true positives, false positives, true negtiaves, false negatives in that order
    """
    tp = {}
    fp = {}
    fn = {}
    tn = {}

    for cl in np.unique(true_labels):
        tp[cl] = 0
        fp[cl] = 0
        fn[cl] = 0
        tn[cl] = 0

    for tl, pred in zip(true_labels, predictions):
        for cl in np.unique(true_labels):
            if tl == pred and tl == cl:
                tp[cl] = tp[cl] + 1
            elif tl == cl and pred != cl:
                fn[cl] = fn[cl] + 1
            elif tl != cl and pred != cl:
                tn[cl] = tn[cl] + 1
            else:
                fp[cl] = fp[cl] + 1

    return tp, fp, tn, fn


def macro_f1(tp, fp, tn, fn, true_labels):
    """
    Compute macro f1 score and additional evaluation metrics for the clustering
    :param tp: True positives
    :param fp: False positives
    :param tn: True negatives
    :param fn: False negatives
    :param true_labels: The true labels
    :return: Macro F1 score, dictionaries of class-wise f1 scores, precision and recall
    """
    precision = {}
    recall = {}
    f1 = {}
    macro_f1 = 0
    for cl in np.unique(true_labels):
        if tp[cl] > 0 or fp[cl] > 0:
            precision[cl] = tp[cl] / (tp[cl] + fp[cl])
        else:
            precision[cl] = 0.0

        if tp[cl] > 0 or fn[cl] > 0:
            recall[cl] = tp[cl] / (tp[cl] + fn[cl])
        else:
            recall[cl] = 0.0
        if precision[cl] or recall[cl] > 0:
            f1[cl] = 2 * precision[cl] * recall[cl] / (precision[cl] + recall[cl])
        else:
            f1[cl] = 0.0
        macro_f1 = macro_f1 + f1[cl]

    macro_f1 = macro_f1 / len(np.unique(true_labels))

    return macro_f1, f1, precision, recall


def evaluate(true_labels, predictions):
    labels, matching = find_cluster_matching(true_labels, predictions)
    tp, fp, tn, fn = confusion_matrix(true_labels, predictions)
    ma_f1, f1, precision, recall = macro_f1(tp, fp, tn, fn, true_labels)
    #print(ma_f1, f1, precision, recall)
    return ma_f1, f1, precision, recall


def simplified_silhouette_scores(data, label, centroids):
    # cannot use sklearn silhouette coefficient here,
    # because it relies solely on the labels, which means, the centroids
    # will be computed in a wrong manner.

    # need to compute simpliefied silhouette coefficient based on the
    # global centroids

    if len(np.unique(label))<=1:
        return [0] * data.shape[0]

    nearest_centroid = []
    nearest_distance = []

    second_best_centroid = []
    second_distance = []


    for i in range(data.shape[0]):
        # compute second best centroid for each data point
        current_nearest = 0
        current_nearest_distance = np.inf

        current_second = 1
        current_second_distance = np.inf


        for c in range(centroids.shape[0]):
            # compute distance
            cd = np.linalg.norm(data[i, :]-centroids[c, :])
            if cd < current_nearest_distance:
                if current_nearest_distance < current_second_distance:
                    current_second_distance = current_nearest_distance
                    current_second = current_nearest
                current_nearest_distance = cd
                current_nearest = c
            elif cd < current_second_distance:
                current_second_distance = cd
                current_second = c

        nearest_centroid.append(current_nearest)
        nearest_distance.append(current_nearest_distance)
        second_best_centroid.append(current_second)
        second_distance.append(current_second_distance)

    silhouette_value = []
    for i in range(data.shape[0]):
        try:
            sv = (second_distance[i]-nearest_distance[i])/max(second_distance[i], nearest_distance[i])
        except:
            sv = -1
        silhouette_value.append(sv)

    # set silhouette for singletons to 0
    bc = np.bincount(nearest_centroid)
    for i in np.where(bc==1)[0]:
        silhouette_value[np.where(nearest_centroid==i)[0][0]]=0

    return silhouette_value

def simplified_silhouette_coefficient(data, label, centroids):
    silhouette_values = simplified_silhouette_scores(data, label, centroids)
    return np.mean(silhouette_values)



#%%

#sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
#sc.logging.print_versions()
#results_file = './write/paul15.h5ad'
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor='white')
sc.settings.figdir='/home/anne/Documents/featurecloud/pca/horizontal-pca/figures/single-cell'# low dpi (dots per inch) yields small inline figures


#%%
orig  = sc.read_h5ad('/home/anne/Documents/featurecloud/pca/horizontal-pca/figures/orig.h5ad')
approx = sc.read_h5ad('/home/anne/Documents/featurecloud/pca/horizontal-pca/figures/approx.h5ad')

#%%
approx.obs
#%%

orig.obs

#%%

#%%
match, approx_labels = find_cluster_matching(orig.obs['leiden'].values, approx.obs['leiden'].values)
#%%
tp, fp, tn, fn = confusion_matrix(orig.obs['leiden'].values, approx_labels)
#%%
mf1, f1, precision, recall  = macro_f1(tp, fp, tn, fn, orig.obs['leiden'].values)

#%%
df = []
precision =  precision
for k,p in zip(precision.keys(), precision.values()):
    df.append([int(k), np.round(float(p),2), np.round(float(recall[k]),2)])

#%%
df = np.stack(df, axis=1).T
#%%
df = pd.DataFrame(df)
df=df.iloc[np.argsort(df.iloc[:, 0])]

df.to_latex('/home/anne/Documents/manuscripts/horizontal-pca-bioinv-adv-clean/supplementary_figures/prec.tex', index=False)
