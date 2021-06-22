from __future__ import print_function
import os.path
import math
import copy
import numpy as np
import heapq
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D
from sklearn.metrics import silhouette_score
from scipy.spatial.distance import cdist
from sklearn.cluster import KMeans as km
import gap
import logging
import subprocess
import os
import sys


def runmap(name, dim, loo):
    newname = "../processing/" + name
    wm = newname + "_wm.csv"
    inputfile = newname + "_dist.csv"
    distancefile = newname + "-distances.txt"
    print(newname, name, wm, inputfile, distancefile)
    os.system(
        "Rscript nmds.R " +
        inputfile +
        " " +
        distancefile +
        " " +
        wm +
        " " +
        dim +
        " " +
        loo)


def readdistance(name):
    data = open(name + "-distances.txt", 'r')
    dist = {}
    for l in data.readlines():
        line = l.strip().split()
        dist[(line[0].upper(), line[1].upper())] = [
            float(line[2]), float(line[3])]
        dist[(line[1].upper(), line[0].upper())] = [
            float(line[2]), float(line[3])]
    data.close()
    return dist


def distancefromlayout(name):  # read distance from layout file of acmacs
    data = open(name + ".layout.txt", 'r')
    sera = {}
    strain = {}
    for l in data.readlines():
        line = l.strip().split()
        if line[0] == 'AG':
            sera[line[1]] = (float(line[2]), float(line[3]))
        if line[0] == 'SR':
            strain[line[1]] = (float(line[2]), float(line[3]))
    data.close()
    dist = {}
    sortedstrain = sorted(strain.keys())
    for i in sorted(sera.keys()):
        for j in sortedstrain:
            dist[(i.upper(), j.upper())] = math.hypot(
                strain[j][0] - sera[i][0], strain[j][1] - sera[i][1])
            dist[(j.upper(), i.upper())] = dist[(i, j)]
    return dist


def readmatrix(name):  # read the original matrix
    sera = []
    strain = []
    hi = []
    table = open("../input/" + name + ".csv", 'r')
    counter = 0
    for line in table:
        if counter == 0:
            # print string.split(line)
            strain = line.strip().split(",")
        else:
            tmp = line.split(",")
            sera.append(tmp[0].upper())
            hi.append(tmp[1:])
        counter += 1
    return sera, strain, hi


def readdistmatrix(name, strain, sera):
    strain = strain
    table = open(name + "-distances.txt", 'r')
    dist = {}
    lines = table.readlines()
    for i in range(1, len(lines)):
        data = lines[i].split()
        for j in range(1, len(data)):
            dist[(sera[i - 1].strip(), strain[j - 1].strip().upper())
                 ] = float(data[i])
            dist[(strain[j - 1].strip().upper(),
                  sera[i - 1].strip())] = float(data[j])
    return dist


def findmin(hi):  # find minimum across each column
    mi = [100] * len(hi[0])
    for i in range(len(hi[0])):
        for j in range(len(hi)):
            # print (hi[j][i])
            if (str(hi[j][i]).strip() == "*" or
                    str(hi[j][i]).strip() == "<" or
                    str(hi[j][i]).strip() == ">"):
                continue
            else:
                # print max
                if float(hi[j][i]) <= mi[i]:
                    mi[i] = float(hi[j][i])
    print("minimum column values: ", mi)
    return mi


def findmax(hi):  # find maximum across each column
    mx = [0] * len(hi[0])
    for i in range(len(hi[0])):
        for j in range(len(hi)):
            if float(hi[j][i]) >= mx[i]:
                mx[i] = float(hi[j][i])
    print("maximum column values: ", mx)
    return mx


def addoptions():
    parser = argparse.ArgumentParser(
        description=("Take neutralization matrix as input and"
                     "give the final antigenic plot"),
        formatter_class=argparse.RawTextHelpFormatter,
        usage=("%(prog)s inputfile [-h]"
               "[-norm NORMALIZE] [-w WEIGHT] [-dim DIM]"))
    parser.add_argument("i", help=" Input file with neutralization matrix")
    parser.add_argument("-norm", "--normalize",
                        help=(" Normalization method \n"
                              " 0: No Normalization (Default) \n"
                              " 1: take column minimum \n"
                              " 2: take global minimum \n"
                              ":>=3: minimum column basis \n"))
    parser.add_argument("-w", "--weight",
                        help=(" Weighing method \n"
                              " 0: No weights (Default)\n"
                              " 1: 1/dij as weights \n"))
    parser.add_argument("-dim", "--dim",
                        help=(" Number of dimensions \n"
                              " 2: Two dimensions (Default) \n"
                              " 3: Three dimensions"))
    parser.add_argument("-loo", "--loo",
                        help=(" Leave-one-out tests \n"
                              " 0: No loo test (Default)) \n"
                              " 1: Compute loo error"))
    args = parser.parse_args()
    return (args.i, args.normalize, args.weight, args.dim, args.loo)


def rnormalize(hi, basis):  # normalize the row neutralization
                            # matrix into distance matrix
    mi = findmin(hi)
    print("basis:", basis)
    # print (min)
    hinew = copy.deepcopy(hi)
    if basis == 0:
        for i in range(len(hi)):
            for j in range(len(mi)):
                hinew[i][j] = float(hinew[i][j])
        print("No normalization")
    if basis == 1:
        for i in range(len(hi)):
            for j in range(len(mi)):
                # hinew[i][j]=(float(hinew[i][j])-min(mi))/10 #enforcing
                # minimum column basis\
                hinew[i][j] = float(hinew[i][j]) - min(mi)
                # print("minimum basis:",basis)
    if basis == 2:
        for i in range(len(hi)):
            for j in range(len(mi)):
                hinew[i][j] = float(hinew[i][j]) - mi[j]
    if basis >= 3:
        for i in range(len(hi)):
            for j in range(len(mi)):
                if mi[j] <= basis:
                    hinew[i][j] = (float(hinew[i][j]) - (mi[j]))
                else:
                    hinew[i][j] = (float(hinew[i][j]) - basis)

    return hinew


def filterbad(hi, bad):
    hinew = copy.deepcopy(hi)
    if bad > 0:
        for i in range(len(hi)):
            for j in range(len(hi[0])):
                if float(hinew[i][j]) >= bad:
                    hinew[i][j] = float("nan")
                else:
                    continue
        else:
            print("All values considered")
    return hinew


''' def addoptions():
    sys.option[0]("-i"+file+"-h Please enter the input file.")
    sys.option[1]("-b") '''


# write distance matrix for the original data
def writeoriginalmatrix(name, sera, strain, hi, weight):
    distfile = open("../processing/" + name + "_dist.csv", "w")
    print(*strain, sep=',', file=distfile)
    tmp = copy.deepcopy(hi)
    # tmp[icheck][jcheck]="Nan"
    for i in range(len(sera)):
        print(sera[i], *tmp[i], sep=',', file=distfile)
    distfile.close()

    # tmp1=[[]*5]
    tmp1 = [[1 for i in range(len(hi[0]))] for j in range(len(hi))]
    # print (tmp1)
    for i in range(len(tmp)):
        for j in range(len(tmp[0])):
            if weight == 1:
                if str(
                    tmp[i][j]) == str("nan") or math.isnan(
                    float(
                        tmp[i][j])):
                    tmp1[i][j] = 1
                    # print("testthis works")
                if str(
                    tmp[i][j]) != str("nan") or math.isnan(
                    float(
                        tmp[i][j])) is False:
                    try:
                        tmp1[i][j] = 1 / tmp[i][j]
                    except BaseException:
                        tmp1[i][j] = 1
            else:
                tmp1[i][j] = 1

    wmfile = open("../processing/" + name + "_wm.csv", "w")
    print(*strain, sep=',', file=wmfile)
    # tmp[icheck][jcheck]="Nan"
    for i in range(len(sera)):
        print(sera[i], *tmp1[i], sep=',', file=wmfile)


# write new matrix with leaving one out as NA
def writematrix(name, sera, strain, hi, icheck, jcheck):
    outname = name + str(icheck) + str(jcheck)
    distfile = open(outname + "_dist.csv", "w")
    print(*strain, sep=',', file=distfile)
    tmp = copy.deepcopy(hi)
    tmp[icheck][jcheck] = "Nan"
    for i in range(len(sera)):
        print(sera[i], *tmp[i], sep=',', file=distfile)
    return outname


def plot2d(name):  # plot 2dimensional neutralization plot
    newname = "../processing/" + name
    dfs = pd.read_csv(newname + "_dist.csv.serum")
    df = pd.read_csv(newname + "_dist.csv.virus")
    fig = plt.figure()
    ax = fig.add_subplot(111)
    xs = list(df['x'])
    ys = list(df['y'])
    xss = list(dfs['x'])
    yss = list(dfs['y'])
    joinedx = xs + xss
    joinedy = ys + yss
    print (xs, min(joinedx))
    xs = [x - min(joinedx) for x in xs]
    ys = [x - min(joinedy) for x in ys]
    xss = [x - min(joinedx) for x in xss]
    yss = [x - min(joinedy) for x in yss]
    print (xs, min(joinedx))
    # print xs

    major_ticks = np.arange(-200, 200, 10)
    minor_ticks = np.arange(-200, 200, 5)

    ax.set_xticks(major_ticks)
    ax.set_xticks(minor_ticks, minor=True)
    ax.set_yticks(major_ticks)

    ax.set_yticks(minor_ticks, minor=True)

    colour = {
        1: 'blue',
        2: 'orange',
        3: 'green',
        4: 'purple',
        5: 'magenta',
        6: 'red'}

    # For each set of style and range settings, plot n random points in the box
    # defined by x in [23, 32], y in [0, 100], z in [zlow, zhigh].
    serum_color = open("serum_gt.csv", 'r').readlines()
    sc = {}
    for i in serum_color:
        sc[int(i.split(',')[0])] = int(i.split(',')[1])

    cols = []
    for i in list(dfs.index):
        cols.append(sc[i])

    for i in range(len(xss)):
        ax.scatter(xss[i], yss[i], s=2, c=colour[cols[i]],
                   marker='.', linewidths=None, edgecolors='face')

        # ax.text(xs[i]-0.055,ys[i]-.16,str(list(df.index)[i]),size=6.75,color=colour[col[i]])
        texts = plt.text(xss[i],
                         yss[i] - .11,
                         str(list(dfs.index)[i]),
                         ha='center',
                         va='center',
                         color=colour[cols[i]],
                         fontsize=8)

    stressdata = open(newname + "_dist.csv.stress", 'r')
    stress = stressdata.readlines()[0]

    col_assign = {
        '2a.3': 2,
        '2b.J8': 2,
        '2b.4': 2,
        '2b.5': 2,
        '2k': 2,
        '2r': 2,
        'Con1': 1,
        '3a': 3,
        '4a': 4,
        '5a': 5,
        'H77': 1,
        '1b.J4': 1,
        '2a.Jc1': 2}
    col = []
    for i in list(df.index):
        col.append(col_assign[i])
    for i in range(len(xs)):
        # print list(df.index)
        ax.scatter(xs[i], ys[i], s=20, c=colour[col[i]],
                   marker='^', linewidths=None, edgecolors='face')

        texts = plt.text(xs[i],
                         ys[i] + .11,
                         str(list(df.index)[i]),
                         ha='center',
                         va='center',
                         color=colour[col[i]],
                         fontsize=9)
    ax.set_xlabel('residual infectivity')
    ax.set_ylabel('residual infectivity')
    ax.set_title("Neutralization map", fontsize=15)
    # ax.tick_params(grid_alpha=0.5)
    ax.grid(
        which='major',
        color='black',
        linestyle='--',
        linewidth=0.5,
        alpha=0.5)
    ax.grid(
        which='minor',
        color='blue',
        linestyle=':',
        linewidth=0.5,
        alpha=0.5)
    plt.axes().set_aspect('equal', 'datalim')
    gt1 = mpatches.Patch(color='blue', label='GT1')
    gt2 = mpatches.Patch(color='orange', label='GT2')
    gt3 = mpatches.Patch(color='green', label='GT3')
    gt4 = mpatches.Patch(color='purple', label='GT4')
    gt5 = mpatches.Patch(color='magenta', label='GT5')
    gt6 = mpatches.Patch(color='red', label='GT6')
    ax.text(
        0.9,
        0.05,
        'Stress: ' +
        stress,
        ha='center',
        va='top',
        transform=ax.transAxes,
        fontsize=8)
    plt.legend(handles=[gt1, gt2, gt3, gt4, gt5, gt6], fontsize='x-small')

    plt.savefig('../results/' + name + '_nmds_plot2d.eps',
                format='eps', dpi=800)
    plt.savefig('../results/' + name + '_nmds_plot2d.png',
                format='png', dpi=800)
    plt.show(block=False)


def plot3d(name):  # plot 3-dimesional neutralization plot
    dfs = pd.read_csv(name + "_dist.csv.serum")
    df = pd.read_csv(name + "_dist.csv.virus")
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    xs = list(df['x'])
    ys = list(df['y'])
    zs = list(df['z'])
    xs = list(df['x'])
    ys = list(df['y'])
    zs = list(df['z'])
    xss = list(dfs['x'])
    yss = list(dfs['y'])
    zss = list(dfs['z'])
    joinedx = xs + xss
    joinedy = ys + yss
    joinedz = zs + zss
    xs = xs - min(joinedx)
    ys = ys - min(joinedy)
    zs = zs - min(joinedz)
    xss = xss - min(joinedx)
    yss = yss - min(joinedy)
    zss = zss - min(joinedz)
    # print xs,ys,zs

    colour = {
        1: 'blue',
        2: 'orange',
        3: 'green',
        4: 'purple',
        5: 'magenta',
        6: 'red'}

    # For each set of style and range settings, plot n random points in the box
    # defined by x in [23, 32], y in [0, 100], z in [zlow, zhigh].
    serum_color = open("serum_gt.csv", 'r').readlines()
    sc = {}
    for i in serum_color:
        sc[int(i.split(',')[0])] = int(i.split(',')[1])

        cols = []
    for i in list(dfs.index):
        cols.append(sc[i])
    for i in range(len(xss)):
        ax.text(xss[i], yss[i], zss[i], str(list(dfs.index)[i]),
                color=colour[cols[i]], fontsize=8)

    stressdata = open("data_dist.csv.stress", 'r')
    stress = stressdata.readlines()[0]
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_zlim(0, 100)

    col_assign = {
        '2a.3': 2,
        '2b.J8': 2,
        '2b.4': 2,
        '2b.5': 2,
        '2k': 2,
        '2r': 2,
        'Con1': 1,
        '3a': 3,
        '4a': 4,
        '5a': 5,
        'H77': 1,
        '1b.J4': 1,
        '2a.Jc1': 2}
    col = []
    for i in list(df.index):
        col.append(col_assign[i])
    for i in range(len(xs)):
        # print list(df.index)
        ax.scatter(xs[i], ys[i], zs[i], s=20, c=colour[col[i]],
                   marker='^', linewidths=None, edgecolors='face')

    ax.text(xs[i], ys[i], zs[i] + .21, str(list(df.index)[i]),
            ha='center', va='center', color=colour[col[i]], fontsize=12)
    ax.set_xlabel('residual infectivity')
    ax.set_ylabel('residual infectivity')
    ax.set_zlabel('residual infectivity')
    ax.set_title("Neutralization map", fontsize=15)
    gt1 = mpatches.Patch(color='blue', label='GT1')
    gt2 = mpatches.Patch(color='orange', label='GT2')
    gt3 = mpatches.Patch(color='green', label='GT3')
    gt4 = mpatches.Patch(color='purple', label='GT4')
    gt5 = mpatches.Patch(color='magenta', label='GT5')
    gt6 = mpatches.Patch(color='red', label='GT6')
    ax.text(
        10,
        100,
        10,
        'Stress: ' +
        stress,
        ha='center',
        va='top',
        transform=ax.transAxes,
        fontsize=12)
    plt.legend(handles=[gt1, gt2, gt3, gt4, gt5, gt6], fontsize='x-small')
    plt.savefig('../results/' + name + '_nmds_plot3d.eps',
                format='eps', dpi=800)
    plt.savefig('../results/' + name + '_nmds_plot3d.png',
                format='png', dpi=800)
    plt.show(block=False)


def silhoutte(name):
    df = pd.read_csv("../processing/" + name + "_dist.csv.virus")
    xs = list(df['x'])
    ys = list(df['y'])
    xs = [x - min(xs) for x in xs]
    ys = [x - min(ys) for x in ys]
    X = np.matrix(zip(xs, ys))
    stat = open('kstatistics.csv', 'w')
    ncluster = []
    distortion_set = []
    silh = []
    for nc in range(2, 13):
        kmeans = km(n_clusters=nc).fit(X)
        cluster_labels = kmeans.fit_predict(X)
        silhouette_avg = silhouette_score(X, cluster_labels)
        distortion = 0
        distortion = (sum(np.min(cdist(X, kmeans.cluster_centers_,
                      'euclidean'), axis=1)) / X.shape[0])
        ncluster.append(nc)
        distortion_set.append(distortion)
        silh.append(silhouette_avg)

        # print kmeans.cluster_centers_, sum(np.min(cdist(X,
        # kmeans.cluster_centers_, 'euclidean'), axis=1))/13

        print(nc, distortion, silhouette_avg)
        print(nc, distortion, silhouette_avg, file=stat)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(ncluster, distortion_set, '-o')
    ax.set_xlabel('No.of clusters')
    ax.set_ylabel('Distortion')
    ax.set_title("Selecting K with Elbow method", fontsize=10)
    plt.show(block=False)

    ''' fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    plt.plot(ncluster, silh, '-o')
    ax1.set_xlabel('No.of clusters')
    ax1.set_ylabel('Silhoutte error')
    ax1.set_title("Selecting K with Silhoutte Error",fontsize=10)
    plt.show(block = False) '''

    # fig2 = plt.figure()


def gapsstatistics(name):
    df = pd.read_csv("../processing/" + name + "_dist.csv.virus")
    xs = list(df['x'])
    ys = list(df['y'])
    X = np.array(zip(xs, ys))
    # print(X)
    # n_clusters = optimalK(df, cluster_array=np.arange(1, 5))
    # print (optimalK.gap_df.head())
    avg_gap = []
    for i in range(1, 5):
        gaps, s_k, K = gap.gap_statistic(
            X, refs=None, B=10, K=range(
                1, 13), N_init=10)
        avg_gap.append(gaps)

    avg_gap = np.array(avg_gap)
    nc = []
    avg_gaperror = []
    max_gap = -100
    max_gap_nc = 1
    for i in range(12):
        gerror = np.mean(avg_gap[:, i])
        nc.append(i + 1)
        avg_gaperror.append(gerror)
        if gerror >= max_gap:
            max_gap = gerror
            max_gap_nc = i + 1
        print(i + 1, np.mean(gerror))
    fig = plt.figure()
    # min_gap = min(avg_gaperror)
    ax = fig.add_subplot(111)
    plt.plot(nc, avg_gaperror, '-o')
    ax.set_xlabel('No.of clusters')
    ax.set_ylabel('Gap error')
    ax.set_title("Selecting K with Gap statistics", fontsize=10)
    ax.arrow(max_gap_nc,
             plt.ylim()[0],
             0,
             max_gap + abs(plt.ylim()[0]),
             width=0.02,
             color='red',
             head_length=0.0,
             head_width=0.0)
    ax.arrow(
        0,
        max_gap,
        max_gap_nc,
        0,
        width=0.0005,
        color='red',
        head_length=0.0,
        head_width=0.0)
    print("Best k:", max_gap_nc)
    plt.show(block=False)
    return max_gap_nc


def cluster(name, k, dim):   # clutser viruses based on k-means clustering
    if dim >= 4:
        print("Dimension is greater than 3. No clustering")
        return None
    df = pd.read_csv("../processing/" + name + "_dist.csv.virus")
    no_cluster = k
    X = []
    if dim == 2:
        xs = list(df['x'])
        ys = list(df['y'])
        X = np.matrix(zip(xs, ys))

    if dim == 3:
        xs = list(df['x'])
        ys = list(df['y'])
        zs = list(df['z'])
        X = np.matrix(zip(xs, ys, zs))
    kmeans = km(n_clusters=no_cluster).fit(X)
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
    print(rep)
    choice = []
    for i in rep.keys():
        mindist = 100
        minpoint = rep[i][0]
        print("Cluster " + str(i) + ": ", end="")
        for j in rep[i]:
            if j != rep[i][-1]:
                print(str(list(df.index)[j]) + ", ", end="")
            else:
                print(str(list(df.index)[j]))
            if dist[j][i] < mindist:
                mindist = dist[j][i]
                minpoint = j
        choice.append(minpoint)
        print("Representative: ", str(list(df.index)[minpoint]))
        print("\n")

    fig = plt.figure()
    ax = fig.add_subplot(111)
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
    plt.savefig('../results/' + name + '_cluster.png', format='png', dpi=800)
    plt.show(block=False)


def main():
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
    tee = subprocess.Popen(["tee", "../log/log.txt"], stdin=subprocess.PIPE)
    os.dup2(tee.stdin.fileno(), sys.stdout.fileno())
    os.dup2(tee.stdin.fileno(), sys.stderr.fileno())
    norm = 0
    weight = 0
    dim = 2
    loo = 0
    # print ("ndwndwi")

    # Read options from commandline
    original_matrixfile, norm, weight, dim, loo = addoptions()

    if norm is not None:
        norm = float(norm)
    if weight is not None:
        weight = float(weight)
    if dim is not None:
        dim = int(dim)
    if loo is not None:
        loo = int(loo)
    if norm is None:
        norm = 0
    if weight is None:
        weight = 0
    if dim is None:
        dim = 2
    if loo is None:
        loo = 0

    print(norm, weight, dim)
    # name of the file with the original neutralization matrix
    original_name = os.path.splitext(original_matrixfile)[0]
    sera, strain, original_hi = readmatrix(original_name)
    # remove unreliable values i.e. values > 100 or 90 etc.
    original_hi = filterbad(original_hi, 101)

    # normlaize the neutralization matrix into a distance matrix
    original_hidist = rnormalize(original_hi, norm)
    # print (original_hi,original_hidist)

    # write the distance matrix for the original neutralization data
    writeoriginalmatrix(original_name, sera, strain, original_hidist, weight)

    runmap(original_name,
           str(dim),
           str(loo))  # run mds for the original data
    if dim == 2:
        plot2d(original_name)
    if dim == 3:
        plot3d(original_name)

    print("------------------------------Results"
          "-----------------------------------------\n")

    stresserror_file = open("../processing/" + original_name +
                            "_dist.csv.stress", 'r')
    stresserror = stresserror_file.readlines()[0].strip()
    stresserror_file.close()

    kerror_file = open("../processing/" + original_name +
                       "_dist.csv.kerror", 'r')
    kerror = kerror_file.readlines()[0].strip()
    kerror_file.close()

    print("Stress:", stresserror)
    print("Kruskal's error in %:", kerror, "\n")

    if loo == 1:

        looerror_file = open("../processing/" +original_name + "_dist.csv.loo", 'r')
        looerror = looerror_file.readlines()[0].strip()
        looerror_file.close()
        print("Leave-one-out error:", looerror, "\n")

    print("-------------------------Cluster Analysis"
          "-------------------------------------\n")
    silhoutte(original_name)
    print("\n")
    k = gapsstatistics(original_name)
    print("Best k from gap statistics:", k)
    print(" Please enter the choice of k:")
    userk = int(raw_input())

    cluster(original_name, userk, dim)

    # print("Pausing works")
    print("Do you want to see silhoutte plots? Enter yes/no:")
    choice = raw_input()
    if choice == "yes" or choice == "Yes" or choice == "YES":
        os.system("python2 silhoutte.py")
    else:
        print("wrong choice")
    plt.show()


if __name__ == "__main__":
    main()
