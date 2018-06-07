from __future__ import print_function
import sys
import time
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

def runmap(name,wm,dim,loo):
	inputfile=name+"_dist.csv"
	distancefile=name+"-distances.txt"
	os.system("Rscript nmds.R "+inputfile+" "+distancefile+" "+wm+" "+dim+" "+loo)

def readdistance(name):
	data=open(name+"-distances.txt",'r')
	dist={}
	for l in data.readlines():
		line=l.strip().split()
		dist[(line[0].upper(),line[1].upper())]=[float(line[2]),float(line[3])]
		dist[(line[1].upper(),line[0].upper())]=[float(line[2]),float(line[3])]
	data.close()
	return dist

def distancefromlayout(name): #read distance from layout file of acmacs output
	data=open(name+".layout.txt",'r')
	sera={}
	strain={}
	for l in data.readlines():
		line=l.strip().split()
		if line[0]=='AG':
			sera[line[1]]=(float(line[2]),float(line[3]))
		if line[0]=='SR':
			strain[line[1]]=(float(line[2]),float(line[3]))
	data.close()
	dist={}
	sortedstrain=sorted(strain.keys())
	for i in sorted(sera.keys()):
		for j in sortedstrain:
			dist[(i.upper(),j.upper())]=math.hypot(strain[j][0]-sera[i][0],strain[j][1]-sera[i][1])
			dist[(j.upper(),i.upper())]=dist[(i,j)]
	return dist

def readmatrix(name): #read the original matrix
	sera = []
	strain = []
	hi = []
	table = open(name+".csv", 'r')
	counter = 0
	for line in table:
		if counter == 0:
			#print string.split(line)
			strain = line.strip().split(",")
		else:
			tmp = line.split(",")
			sera.append(tmp[0].upper())
			hi.append(tmp[1:])
		counter += 1
	return sera,strain,hi



def readdistmatrix(name, strain, sera):
	strain = strain
	hidist = []
	table = open(name+"-distances.txt", 'r')
	dist = {}
	lines = table.readlines()
	for i in range (1,len(lines)):
		data = lines[i].split()
		for j in range(1, len(data)):
			#print (data[0].strip().upper(), strain[j-1].upper(), float(data[j]))
			dist[(sera[i-1].strip(), strain[j-1].strip().upper())] = float(data[i])
			dist[(strain[j-1].strip().upper(), sera[i-1].strip())] = float(data[j])
	return dist





def findmin(hi): #find minimum across each column
	mi=[100]*len(hi[0])
	for i in range(len(hi[0])):
		for j in range(len(hi)):
			#print (hi[j][i])
			if str(hi[j][i]).strip()=="*" or str(hi[j][i]).strip()=="<" or str(hi[j][i]).strip()==">" :
				continue
			else:
				#print max
				if float(hi[j][i])<=mi[i]:
					mi[i]=float(hi[j][i])
	print ("minimum column values: ",mi)
	return mi

def findmax(hi): #find maximum across each column
	mx = [0]*len(hi[0])
	for i in range(len(hi[0])):
		for j in range(len(hi)):
			if float(hi[j][i]) >= mx[i]:
				mx[i] = float(hi[j][i])
	print ("maximum column values: ", mx)
	return mx

def addoptions():
	parser = argparse.ArgumentParser(description = "Take neutralization matrix as input and give the final antigenic plot",formatter_class = argparse.RawTextHelpFormatter, usage = "%(prog)s inputfile [-h] [-norm NORMALIZE] [-w WEIGHT] [-dim DIM]")
	parser.add_argument("i", help = " Input file with neutralization matrix")
	parser.add_argument("-norm", "--normalize", help = " Normalization method \n   0: No Normalization (Default) \n   1: take column minimum \n   2: take global minimum \n >=3: minimum column basis \n")
	parser.add_argument("-w", "--weight", help = " Weighing method \n   0: No weights (Default) 	\n   1: 1/dij as weights \n")
	parser.add_argument("-dim", "--dim", help = " Number of dimensions \n   2: Two dimensions (Default) \n   3: Three dimensions")
	parser.add_argument("-loo", "--loo", help = " Leave-one-out tests \n   0: No loo test (Default)) \n   1: Compute loo error")
	args = parser.parse_args()
	return (args.i, args.normalize, args.weight, args.dim, args.loo)


def rnormalize(hi, basis): #normalize the row neutralization matrix into distance matrix
	mi = findmin(hi)
	print ("basis:",basis)
	#print (min)
	hinew=copy.deepcopy(hi)
	if basis == 0:
		for i in range(len(hi)):
			for j in range(len(mi)):
				hinew[i][j] = float(hinew[i][j])
		print ("No normalization")
	if basis == 1:
		for i in range(len(hi)):
			for j in range(len(mi)):
				#hinew[i][j]=(float(hinew[i][j])-min(mi))/10 #enforcing minimum column basis\
				hinew[i][j] = float(hinew[i][j]) - min(mi)
				#print("minimum basis:",basis)
	if basis == 2:
		for i in range(len(hi)):
			for j in range(len(mi)):
				hinew[i][j] = float(hinew[i][j]) - mi[j]
	if basis >= 3:
		for i in range(len(hi)):
			for j in range(len(mi)):
				if mi[j]<=basis:
					hinew[i][j] = (float(hinew[i][j]) - (mi[j]))
				else:
					hinew[i][j] = (float(hinew[i][j]) - basis)


	return hinew



def filterbad(hi,bad):
	hinew = copy.deepcopy(hi)
	if bad > 0:
		for i in range(len(hi)):
			for j in range(len(hi[0])):
				if float(hinew[i][j]) >= bad:
					hinew[i][j] = float("nan")
				else:
					continue
		else:
			print ("All values considered")
	return hinew

'''def addoptions():
	sys.option[0]("-i"+file+"-h Please enter the input file.")
	sys.option[1]("-b")'''

def writeoriginalmatrix(name, sera, strain, hi, weight): #write distance matrix for the original data
	distfile = open(name + "_dist.csv", "w")
	print (*strain, sep = ',', file = distfile)
	tmp = copy.deepcopy(hi)
	#tmp[icheck][jcheck]="Nan"
	for i in range(len(sera)):
		print (sera[i], *tmp[i], sep = ',', file = distfile)
	distfile.close()

	#tmp1=[[]*5]
	tmp1 = [[1 for i in range(len(hi[0]))] for j in range(len(hi))]
	#print (tmp1)
	mx = max(findmax(tmp))
	mn = min(findmin(tmp))
	for i in range(len(tmp)):
		for j in range(len(tmp[0])):
			if weight==1:
				if str(tmp[i][j]) == str("nan") or math.isnan(float(tmp[i][j])):
					tmp1[i][j] = 1
					#print("testthis works")
				if str(tmp[i][j]) != str("nan") or math.isnan(float(tmp[i][j])) == False:
					try:
						tmp1[i][j] = 1/tmp[i][j]
					except:
						tmp1[i][j] = 1
			else:
				tmp1[i][j] = 1
				

	wmfile=open(name+"_wm.csv","w")
	print (*strain,sep=',',file=wmfile)
	#tmp[icheck][jcheck]="Nan"
	for i in range(len(sera)):
		print (sera[i],*tmp1[i],sep=',',file=wmfile)




def writematrix(name, sera, strain, hi, icheck, jcheck): #write new matrix with leaving one out as NA
	outname = name + str(icheck) + str(jcheck)
	distfile = open(outname + "_dist.csv", "w")
	print (*strain, sep = ',', file = distfile)
	tmp = copy.deepcopy(hi)
	tmp[icheck][jcheck] = "Nan"
	for i in range(len(sera)):
		print (sera[i], *tmp[i], sep = ',', file = distfile)
	return outname

def plot2d(name):


	dfs = pd.read_csv(name+"_dist.csv.serum")
	df = pd.read_csv(name+"_dist.csv.virus")
	fig = plt.figure()
	ax = fig.add_subplot(111)
	xs = list(df['x'])
	ys = list(df['y'])
	xss = list(dfs['x'])
	yss = list(dfs['y'])
	joinedx = xs + xss
	joinedy = ys + yss
	#print xs
	xs = xs - min(joinedx)
	ys = ys - min(joinedy)
	xss = xss - min(joinedx)
	yss = yss - min(joinedy)
	#print xs


	major_ticks = np.arange(-200, 200, 10)
	minor_ticks = np.arange(-200, 200, 5)

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


	cols= []
	for i in list(dfs.index):
		cols.append(sc[i])

	for i in range(len(xss)):
		ax.scatter(xss[i], yss[i],s=2,c=colour[cols[i]],marker='.',linewidths=None,edgecolors='face')
		
		#ax.text(xs[i]-0.055,ys[i]-.16,str(list(df.index)[i]),size=6.75,color=colour[col[i]])
		texts=plt.text(xss[i],yss[i]-.11,str(list(dfs.index)[i]), ha='center', va='center',color=colour[cols[i]],fontsize=8)


	stressdata=open("data_dist.csv.stress",'r')
	stress=stressdata.readlines()[0]


	col_assign={'2a.3': 2, '2b.J8':2,'2b.4':2, '2b.5':2,'2k':2,'2r':2,'Con1':1, '3a':3, '4a':4, '5a':5,'H77':1,'1b.J4':1,'2a.Jc1':2}
	col=[]
	for i in list(df.index):
		col.append(col_assign[i])
	for i in range(len(xs)):
		#print list(df.index)
		ax.scatter(xs[i], ys[i],s=20,c=colour[col[i]],marker='^',linewidths=None,edgecolors='face')
		
		texts=plt.text(xs[i],ys[i]+.11,str(list(df.index)[i]), ha='center', va='center',color=colour[col[i]],fontsize=9)
	ax.set_xlabel('residual infectivity')
	ax.set_ylabel('residual infectivity')
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


	plt.savefig('map.eps', format='eps', dpi=800)
	plt.savefig('myimage.svg', format='svg',figsize=(10,10),dpi=800)
	plt.show(block = False)


def plot3d(name):
	dfs = pd.read_csv(name+"_dist.csv.serum")
	df = pd.read_csv(name+"_dist.csv.virus")
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
	joinedz = zs +  zss
	xs = xs - min(joinedx)
	ys = ys - min(joinedy)
	zs = zs - min(joinedz)
	xss = xss - min(joinedx)
	yss = yss - min(joinedy)
	zss = zss - min(joinedz)
	#print xs,ys,zs

	colour = {1:'blue', 2:'orange', 3:'green', 4:'purple', 5:'magenta', 6:'red'}

	# For each set of style and range settings, plot n random points in the box
	# defined by x in [23, 32], y in [0, 100], z in [zlow, zhigh].
	serum_color = open("serum_gt.csv", 'r').readlines()
	sc = {}
	for i in serum_color:
		sc[int(i.split(',')[0])] = int(i.split(',')[1])


		cols = []
	for i in list(dfs.index):
		cols.append(sc[i])
		#ax.scatter(xss, yss,zss,s=2,c=colour[cols[i]],marker='.',linewidths=None,edgecolors='face')
	for i in range(len(xss)):
		ax.text(xss[i], yss[i], zss[i], str(list(dfs.index)[i]), color = colour[cols[i]], fontsize = 8)

	stressdata = open("data_dist.csv.stress", 'r')
	stress = stressdata.readlines()[0]
	ax.set_xlim(0, 100)
	ax.set_ylim(0, 100)
	ax.set_zlim(0, 100)

	col_assign={'2a.3': 2, '2b.J8':2,'2b.4':2, '2b.5':2,'2k':2,'2r':2,'Con1':1, '3a':3, '4a':4, '5a':5,'H77':1,'1b.J4':1,'2a.Jc1':2}
	col=[]
	for i in list(df.index):
		col.append(col_assign[i])
	for i in range(len(xs)):
		#print list(df.index)
		ax.scatter(xs[i], ys[i], zs[i], s=20,c=colour[col[i]],marker='^',linewidths=None,edgecolors='face')

	ax.text(xs[i],ys[i], zs[i]+.21, str(list(df.index)[i]), ha='center', va='center',color=colour[col[i]],fontsize=12)
	ax.set_xlabel('residual infectivity')
	ax.set_ylabel('residual infectivity')
	ax.set_zlabel('residual infectivity')
	ax.set_title("Neutralization map",fontsize=15)
	gt1 = mpatches.Patch(color='blue', label='GT1')
	gt2 = mpatches.Patch(color='orange', label='GT2')
	gt3 = mpatches.Patch(color='green', label='GT3')
	gt4 = mpatches.Patch(color='purple', label='GT4')
	gt5 = mpatches.Patch(color='magenta', label='GT5')
	gt6 = mpatches.Patch(color='red', label='GT6')
	ax.text(10,100,10, 'Stress: '+stress, ha = 'center',va = 'top', transform = ax.transAxes, fontsize=12)
	plt.legend(handles=[gt1, gt2, gt3, gt4, gt5, gt6],fontsize = 'x-small')
	plt.savefig('map3d.eps', format='eps', dpi=800)
	plt.savefig('myimage3d.png', format='png',figsize=(10,10),dpi=800)
	plt.show(block = False)


def silhoutte(name):
	df = pd.read_csv(name + "_dist.csv.virus")
	xs = list(df['x'])
	ys = list(df['y'])
	xs = xs - min(xs)
	ys = ys - min(ys)
	X = np.matrix(zip(xs,ys))
	stat=open('kstatistics.csv','w')
	ncluster=[]
	distortion_set=[]
	silh=[]
	for nc in range(2,13):
		kmeans = km(n_clusters=nc).fit(X)
		cluster_labels = kmeans.fit_predict(X)
		silhouette_avg = silhouette_score(X, cluster_labels)
		#print("For n_clusters =", no_cluster, "The average silhouette_score is :", silhouette_avg)
		labels = kmeans.labels_
		dist = kmeans.transform(X)
		distortion=0
		distortion= (sum(np.min(cdist(X, kmeans.cluster_centers_, 'euclidean'), axis=1)) / X.shape[0])
		ncluster.append(nc)
		distortion_set.append(distortion)
		silh.append(silhouette_avg)

		#print kmeans.cluster_centers_, sum(np.min(cdist(X, kmeans.cluster_centers_, 'euclidean'), axis=1))/13

		print (nc, distortion, silhouette_avg)
		print (nc, distortion, silhouette_avg, file = stat)
	fig = plt.figure()
	ax = fig.add_subplot(111)
	plt.plot(ncluster, distortion_set, '-o')
	ax.set_xlabel('No.of clusters')
	ax.set_ylabel('Distortion')
	ax.set_title("Selecting K with Elbow method",fontsize=10)
	plt.show(block = False)

	'''fig1 = plt.figure()
	ax1 = fig1.add_subplot(111)
	plt.plot(ncluster, silh, '-o')
	ax1.set_xlabel('No.of clusters')
	ax1.set_ylabel('Silhoutte error')
	ax1.set_title("Selecting K with Silhoutte Error",fontsize=10)
	plt.show(block = False)'''

	fig2 = plt.figure()


def gapsstatistics(name):
	df = pd.read_csv(name + "_dist.csv.virus")
	xs = list(df['x'])
	ys = list(df['y'])
	X = np.array(zip(xs,ys))
	#print(X)
	#n_clusters = optimalK(df, cluster_array=np.arange(1, 5))
	#print (optimalK.gap_df.head())
	avg_gap = []
	for i in range(1, 5):
		gaps, s_k, K = gap.gap_statistic(X, refs = None, B = 10, K = range(1,13), N_init = 10)
		avg_gap.append(gaps)

	avg_gap = np.array(avg_gap)
	nc=[]
	avg_gaperror=[]
	max_gap=0
	max_gap_nc=1
	for i in range(12):
		gerror=np.mean(avg_gap[:,i])
		nc.append(i + 1)
		avg_gaperror.append(gerror)
		if gerror >= max_gap:
			max_gap = gerror
			max_gap_nc = i+1
		print (i+1, np.mean(gerror))
	fig = plt.figure()
	min_gap=min(avg_gaperror)
	ax = fig.add_subplot(111)
	plt.plot(nc, avg_gaperror, '-o')
	ax.set_xlabel('No.of clusters')
	ax.set_ylabel('Gap error')
	ax.set_title("Selecting K with Gap statistics",fontsize=10)
	ax.arrow(max_gap_nc,plt.ylim()[0],0,max_gap+abs(plt.ylim()[0]),width=0.02,color='red',head_length=0.0,head_width=0.0)
	ax.arrow(0,max_gap,max_gap_nc,0,width=0.0005,color='red',head_length=0.0,head_width=0.0)
	print ("Best k:",max_gap_nc)
	plt.show(block = False)





def main():
	norm = 0
	weight = 0
	dim = 2
	loo = 0
	#print ("ndwndwi")


	#Read options from commandline
	original_matrixfile, norm, weight, dim, loo = addoptions()
	
	if norm != None:
		norm = float(norm)
	if weight != None:
		weight = float(weight)
	if dim != None:
		dim = int(dim)
	if loo != None:
		loo =int(loo)
	if norm is None:
		norm = 0
	if weight is None:
		weight = 0
	if dim is None:
		dim = 2
	if loo is None:
		loo = 0
	
	print (norm, weight, dim)
	original_name = os.path.splitext(original_matrixfile)[0] #name of the file with the original neutralization matrix
	sera,strain,original_hi = readmatrix(original_name)
	original_hi = filterbad(original_hi,101) #remove unreliable values i.e. values > 100 or 90 etc.

	original_hidist = rnormalize(original_hi,norm) #normlaize the neutralization matrix into a distance matrix
	#print (original_hi,original_hidist)
	
	writeoriginalmatrix(original_name,sera,strain,original_hidist,weight) #write the distance matrix for the original neutralization data
	
	runmap(original_name,original_name + "_wm.csv",str(dim),str(loo)) #run mds for the original data
	if dim == 2:
		plot2d(original_name)
	if dim == 3:
		plot3d(original_name)

	print("------------------------------Results-----------------------------------------\n")

	stresserror_file = open(original_name + "_dist.csv.stress", 'r')
	stresserror = stresserror_file.readlines()[0].strip()
	stresserror_file.close()

	kerror_file = open(original_name + "_dist.csv.kerror", 'r')
	kerror = kerror_file.readlines()[0].strip()
	kerror_file.close()

	print("Stress:", stresserror)
	print("Kruskal's error in %:", kerror,"\n")
	
	if loo==1:

		looerror_file = open(original_name + "_dist.csv.loo", 'r')
		looerror = looerror_file.readlines()[0].strip()
		looerror_file.close()
		print("Leave-one-out error:", looerror,"\n")

	print("-------------------------Cluster Analysis-------------------------------------\n")
	silhoutte(original_name)
	print("\n")
	gapsstatistics(original_name)

	
	#print("Pausing works")
	print ("Do you want to see silhoutte plots? Enter yes/no:")
	choice = raw_input()
	if choice == "yes" or choice == "Yes" or choice == "YES":
		os.system("python2 silhoutte.py")
	else:
		print("wrong choice")
	plt.show()


if __name__ == "__main__":
	main()
