# This file reads the 3d co-ordinates from the R output and then plots a 3d scatter plot of the entire data.

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.patches as mpatches
from sklearn.cluster import KMeans as km
from sklearn.metrics import silhouette_score
import sys
from scipy.spatial.distance import cdist
import numpy as np


df = pd.read_csv("data_dist.csv.virus")
xs = list(df['x'])
ys = list(df['y'])
xs = xs - min(xs)
ys = ys - min(ys)
X = np.matrix(zip(xs, ys))
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
X = scaler.fit_transform(X)
stat = open('kstatistics.csv', 'w')
for nc in range(2, 13):
    kmeans = km(n_clusters=nc, random_state=10)
    cluster_labels = kmeans.fit_predict(X)
    silhouette_avg = silhouette_score(X, cluster_labels)
    #print("For n_clusters =", no_cluster, "The average silhouette_score is :", silhouette_avg)
    labels = kmeans.labels_
    dist = kmeans.transform(X)
    distortion = 0
    distortion = (
        sum(np.min(cdist(X, kmeans.cluster_centers_, 'euclidean'), axis=1)) / X.shape[0])
    # print kmeans.cluster_centers_, sum(np.min(cdist(X, kmeans.cluster_centers_, 'euclidean'), axis=1))/13

    print nc, distortion, silhouette_avg
    print >>stat, nc, distortion, silhouette_avg