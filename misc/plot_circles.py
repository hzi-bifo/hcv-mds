# This file reads the 3d co-ordinates from the R output and then plots a 3d scatter plot of the entire data.

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.patches as mpatches
import sys

dfs=pd.read_csv("data_dist.csv.serum")
fig = plt.figure()
ax = fig.add_subplot(111)


major_ticks = np.arange(-10, 10, 2)
minor_ticks = np.arange(-10, 10, 1)

ax.set_xticks(major_ticks)
ax.set_xticks(minor_ticks, minor=True)
ax.set_yticks(major_ticks)

ax.set_yticks(minor_ticks, minor=True)

colour={1:'blue',2:'orange',3:'green',4:'purple',5:'magenta',6:'red'}

# For each set of style and range settings, plot n random points in the box
# defined by x in [23, 32], y in [0, 100], z in [zlow, zhigh].
serum_color=open("serum_gt.csv",'r').readlines()
sc={}
for i in serum_color:
	sc[int(i.split(',')[0])]=int(i.split(',')[1])

xss = list(dfs['x'])
yss = list(dfs['y'])
cols= []
for i in list(dfs.index):
	cols.append(sc[i])

for i in range(len(xss)):
	ax.scatter(xss[i], yss[i],s=2,c=colour[cols[i]],marker='.',linewidths=None,edgecolors='face')

	#ax.text(xs[i]-0.055,ys[i]-.16,str(list(df.index)[i]),size=6.75,color=colour[col[i]])
	texts=plt.text(xss[i],yss[i]-.11,str(list(dfs.index)[i]), ha='center', va='center',color=colour[cols[i]],fontsize=8)

df=pd.read_csv("data_dist.csv.virus")
stressdata=open("data_dist.csv.stress",'r')
stress=stressdata.readlines()[0]

xs = list(df['x'])
ys = list(df['y'])
col_assign={'2a.3': 2, '2b.J8':2,'2b.4':2, '2b.5':2,'2k':2,'2r':2,'Con1':1, '3a':3, '4a':4, '5a':5,'H77':1,'1b.J4':1,'2a.Jc1':2}
col=[]
for i in list(df.index):
	col.append(col_assign[i])
for i in range(len(xs)):
	#print list(df.index)
	ax.scatter(xs[i], ys[i],s=20,c=colour[col[i]],marker='^',linewidths=None,edgecolors='face')

	texts=plt.text(xs[i],ys[i]+.11,str(list(df.index)[i]), ha='center', va='center',color=colour[col[i]],fontsize=9)
ax.set_xlabel('Antigenic distance')
ax.set_ylabel('Antigenic distance')
ax.set_title("Neutralization map",fontsize=15)
#ax.tick_params(grid_alpha=0.5)
ax.grid(which='major',color='black', linestyle='--', linewidth=0.5,alpha=0.5)
ax.grid(which='minor',color='blue', linestyle=':', linewidth=0.5,alpha=0.5)
plt.axes().set_aspect('equal', 'datalim')
gt1 = mpatches.Patch(color='blue', label='GT1')
gt2 = mpatches.Patch(color='orange', label='GT2')
gt3 = mpatches.Patch(color='green', label='GT3')
gt4 = mpatches.Patch(color='purple', label='GT4')
gt5 = mpatches.Patch(color='magenta', label='GT5')
gt6 = mpatches.Patch(color='red', label='GT6')
ax.text(0.9,0.05, 'Stress: '+stress, ha = 'center',va = 'top', transform = ax.transAxes, fontsize=8)
plt.legend(handles=[gt1, gt2, gt3, gt4, gt5, gt6],fontsize = 'x-small')

vn=int(sys.argv[1])
circle1=plt.Circle((xs[vn],ys[vn]),50,color=colour[col_assign[str(list(df.index)[i])]],fill=False)
ax.add_artist(circle1)
plt.savefig('map_col_'+str(vn)+'.eps', format='eps', dpi=800)
plt.savefig('map_col_'+str(vn)+'.svg', format='svg',figsize=(10,10),dpi=800)
plt.show()
