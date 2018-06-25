import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans as km
import sys

df = pd.read_csv("data_dist.csv.virus")
no_cluster = int(sys.argv[1])

xs = list(df['x'])
ys = list(df['y'])
xs = xs - min(xs)
ys = ys - min(ys)
X = np.matrix(zip(xs, ys))
kmeans = km(n_clusters=no_cluster).fit(X)
cluster_labels = kmeans.fit_predict(X)

colour = {0: 'blue', 1: 'red', 2: 'green',
          3: 'black', 4: 'magenta', 5: 'orange'}
# print kmeans.predict(X)
labels = kmeans.labels_
dist = kmeans.transform(X)
rep = {}
for i in range(no_cluster):
    rep[i] = []
for i in range(len(labels)):
    rep[labels[i]].append(i)
print rep
choice = []
distortion = 0
for i in rep.keys():
    mindist = 100
    minpoint = rep[i][0]
    print "Cluster " + str(i) + ": ",
    for j in rep[i]:
        if j != rep[i][-1]:
            print str(list(df.index)[j]) + ",",
        else:
            print str(list(df.index)[j])
        if dist[j][i] < mindist:
            mindist = dist[j][i]
            minpoint = j
    choice.append(minpoint)
    print "Representative: ", str(list(df.index)[minpoint])
    print "\n"

fig = plt.figure()
ax = fig.add_subplot(111)
centroid = np.array(kmeans.cluster_centers_)
for i in range(len(xs)):
    # print list(df.index)
    if i in choice:
        ax.scatter(xs[i], ys[i], s=50, c=colour[kmeans.labels_[i]],
                   marker='^', linewidths=0.5, edgecolors='face')
        plt.text(xs[i], ys[i] + .8, str(list(df.index)[i]), ha='center',
                 va='center', color=colour[kmeans.labels_[i]], fontsize=15)
    else:
        ax.scatter(xs[i], ys[i], s=20, c=colour[kmeans.labels_[i]],
                   marker='^', linewidths=None, edgecolors='face')
        plt.text(xs[i], ys[i] + .8, str(list(df.index)[i]), ha='center',
                 va='center', color=colour[kmeans.labels_[i]], fontsize=12)

# ax.scatter(centroid[:,0],centroid[:,1],s=20,marker='o',linewidths=None,edgecolors='face')
plt.savefig('map_rep' + str(no_cluster) + '.eps', format='eps', dpi=800)
plt.show()
