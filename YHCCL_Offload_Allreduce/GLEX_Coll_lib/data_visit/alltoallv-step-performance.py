# coding=UTF-8
import os
import sys
FileName="alltoallv_all_statistics.txt"
title=sys.argv[2]
#print(str(sys.argv))
f=open(sys.argv[1]+FileName,"r")

Tmp=((f.readline()).strip()).split()
procn=int(Tmp[0])
statistics_count=int(Tmp[1])
Data_vec = []
for i in range(0,procn):
    vec = ((f.readline()).strip()).split()
    vec_double  = []
    for d in vec:
        vec_double.append(float(d))
    # for j in range(1,len(vec)):
    #     vec_double[j]+=vec_double[j-1]
    ##print(vec_double)
    Data_vec.append(vec_double)
f.close()

# import seaborn as sns
# import matplotlib.pyplot as plt
# from pandas.core.frame import DataFrame
# sns.set_theme(style="whitegrid")

# x=[0,1,2,3]
# c={"a" : a,
# print(data1)

# # Initialize the matplotlib figure
# f, ax = plt.subplots(figsize=(6, 15))

# # Load the example car crash dataset
# crashes = sns.load_dataset('car_crashes',data_home='./',cache=True).sort_values("total", ascending=False)
# print(type(crashes))


# sns.set_color_codes("pastel")
# sns.barplot(x="a", y="x", data=data1,
#             label="Total", color="b")
# sns.set_color_codes("muted")
# sns.barplot(x="b", y="x", data=data1,
#             label="Alcohol-involved", color="b")

# # sns.set_color_codes("pastel")
# # sns.barplot(x="total", y="abbrev", data=crashes,
# #             label="Total", color="b")

# # # Plot the crashes where alcohol was involved
# # sns.set_color_codes("muted")
# # sns.barplot(x="alcohol", y="abbrev", data=crashes,
# #             label="Alcohol-involved", color="b")


# # Add a legend and informative axis label
# ax.legend(ncol=2, loc="lower right", frameon=True)
# ax.set(xlim=(0, 10), ylabel="",
#        xlabel="Automobile collisions per billion miles")
# sns.despine(left=True, bottom=True)


# # import numpy as np
# # import matplotlib as mpl
# # import matplotlib.pyplot as plt
# # import seaborn as sns
# # np.random.seed(sum(map(ord, "aesthetics")))
# # def sinplot(flip=1):
# #     x = np.linspace(0, 14, 100)
# #     for i in range(1, 7):
# #         plt.plot(x, np.sin(x + i * .5) * (7 - i) * flip)
# # sinplot()        
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
 
# Data
x=range(0,procn)
Data_vec_numpy=np.array(Data_vec)
y=Data_vec_numpy.T
#print(y.shape)
labela = []
for v in x:
    labela.append(str(x))
#y=[ [1,4,6,8,9], [2,2,7,10,12], [2,8,5,10,6] ]
 
# Plot
plt.stackplot(x,y,labels=labela)



plt.savefig('./'+title+'.pdf')
plt.show( )