name="Reduce-pocn[2,72]-size[1,13]-K-ary.out"
f=open(name)
TEXT=f.readlines()
starti=0;
KMap={}
startNode=50
EndNode=950
NodeStep=50
loopn=5
sizeN=13
for k in range(startNode,EndNode,NodeStep):
	KMap[k]=[]
for noden in range(startNode,EndNode,NodeStep):
	for k in range(2,17,2):
		mpivec=[0.0]*sizeN
		glexvec=[0.0]*sizeN
		for loop in range(0,loopn):
			starti+=1
			for size in range(0,sizeN):
				strs=TEXT[starti].split()
				for str1 in strs:
					str1=str1.strip()
				# print(strs)
				mpivec[size]+=float(strs[1])/loopn
				glexvec[size]+=float(strs[2])/loopn
				starti+=1
		improvementvec=[0.0]*sizeN
		for i in range(0,sizeN):
			improvementvec[i] = (mpivec[i]-glexvec[i])/glexvec[i]
		# print(mpivec,glexvec)
		# print(improvementvec)
		KMap[k].append(improvementvec)
		avgVal = sum(improvementvec)

#接下来找到最好的K值

maxK=-1
maxval=-10000.0
for k in range(2,17,2):
	sumall = 0.0
	for vec in KMap[k]:
		sumall+=sum(vec)
	if sumall > maxval:
		maxval = sumall
		maxK = k
#接下来对每一个K值的不同节点数和消息大小绘制一张三维图
f.close()

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import numpy as np
from matplotlib import cm
print(maxK)
for k in range(maxK,maxK+1,2):
	fig=plt.figure()
	ax=fig.add_subplot(111,projection='3d')
	noden=np.arange(0,71)
	size=np.arange(0,13)
	x,y=np.meshgrid(size,noden)
	z=np.array(KMap[k])
	improve=x*y
	# print(improve)
# 
# z=x**2+y**2
	ax.plot_wireframe(x,y,z, rstride=1, cstride=1)
	#ax.plot_surface(x,y,z,rstride=1,cstride=1,cmap=plt.get_cmap('Blues_r'))
	ax.set_xlabel('消息大小2^x')
	ax.set_ylabel('节点数量')
	ax.set_zlabel('对比MPI性能提升')
	tile=str(maxK)
	ax.set_title(tile+'-ary RDMA Reduce树对比MPI Reduce')
	plt.rcParams[
     'font.sans-serif']=[
     'SimHei'] 
	# plt.show()
	# exit(0)

# for k in range(2,17,2):
# 	for noden in range(0,71):
# 		plt.clf()
# 		print("-------------------------k="+str(k)+" noden="+str(noden+2)+"-------------------------")
# 		title="k="+str(k)+" noden="+str(noden+2)+"-msgsize-improvment"
# 		vec = KMap[k][noden]
# 		# for v in vec:
# 		# 	print(v,end=' ')
# 		# print("")
# 		x_data = range(1,14,1)
# 		y_data = KMap[k][noden]
# 		plt.plot(x_data,y_data,color='red',linewidth=2.0,linestyle='--')
# 		ax.set_xlabel('msgsize2^x')
# 		ax.set_ylabel('Improvement compared to MPI')
# 		ax.set_title(title)
# 	# plt.show()
#       	#plt.plot(x_data,y_data2,color='blue',linewidth=3.0,linestyle='-.')
# 		plt.savefig("./消息大小-性能提升图/"+title+".pdf")

for k in range(2,17,2):
	for msgsize in range(0,13):
		title="k="+str(k)+" msgsize=2^"+str(msgsize+1)+"-msgsize-improvment"
		x_data = range(2,73)
		y_data = [x[msgsize] for x in KMap[k]]
		# print(y_data)
		plt.clf()
		print("-------------------------"+title+"-------------------------")
		plt.plot(x_data,y_data,color='green',linewidth=2.0,linestyle='-')
		ax.set_xlabel('Node Num')
		ax.set_ylabel('Improvement compared to MPI')
		ax.set_title(title)
      	#plt.plot(x_data,y_data2,color='blue',linewidth=3.0,linestyle='-.')
		plt.savefig("./节点数量-性能提升图/"+title+".pdf")