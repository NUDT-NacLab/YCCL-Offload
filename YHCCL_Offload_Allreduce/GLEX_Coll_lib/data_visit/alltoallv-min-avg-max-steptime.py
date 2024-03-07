# coding=UTF-8
import os
import sys
FileName="alltoallv_performance.txt"
#draw the step time statistics
title=sys.argv[2]
#print(str(sys.argv))
f=open(sys.argv[1]+FileName,"r")
Tmp=((f.readline()).strip()).split()
statistics_count=int(Tmp[0])
#print(statistics_count)
Data_vec = []
for i in range(0,statistics_count):
    vec = ((f.readline()).strip()).split()
    vec_double  = []
    for d in vec:
        vec_double.append(float(d))
    Data_vec.append(vec_double)
    #print(i,vec_double)
f.close()

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# 处理乱码

x=range(0,statistics_count)
y = Data_vec
# "r" 表示红色，ms用来设置*的大小
plt.plot(x, y, marker='*', ms=10, label="a")
# plt.plot([1, 2, 3, 4], [20, 30, 80, 40], label="b")
plt.xticks(rotation=45)
plt.xlabel("step")
plt.title("step overhead min avg max us")
# upper left 将图例a显示到左上角
plt.savefig('./'+title+'.pdf')
#plt.show()
