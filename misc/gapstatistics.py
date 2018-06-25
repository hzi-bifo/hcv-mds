import numpy as np
import pandas as pd
import gap

# optimalK = OptimalK(parallel_backend='joblib')
df = pd.read_csv("data_dist.csv.virus")
# sno_cluster = int(sys.argv[1])

xs = list(df['x'])
ys = list(df['y'])
X = np.array(zip(xs, ys))
# print(X)
# n_clusters = optimalK(df, cluster_array=np.arange(1, 5))
# print (optimalK.gap_df.head())
avg_gap = []
for i in range(1, 2):
    gaps, s_k, K = gap.gap_statistic(
        X, refs=None, B=10, K=range(1, 13), N_init=10)
    avg_gap.append(gaps)

avg_gap = np.array(avg_gap)
for i in range(12):
    print i + 1, np.mean(avg_gap[:, i])