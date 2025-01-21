import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
data = pd.read_excel('dataflow_20231020.xlsx')
parameter = data.iloc[:,0].values

x=data.iloc[:,1].values

y1=data.iloc[:,6].values
y2=data.iloc[:,7].values
y3=data.iloc[:,8].values
y4=data.iloc[:,9].values
x = list(map(int, x))

custom_x = []
x1 = x[0]
xn = x[len(x)-1]
block = int(((xn-x1)/20))
for i in range(21):
    custom_x.append(i*block)

custom_y = []
y = 0
yn = 3
block = float((yn-y)/6)
for i in range(7):
    custom_y.append(i*block)
fig, ax = plt.subplots(figsize=(10, 3))
ax.plot(x, y1, label="inner-product")
ax.plot(x, y2, label="outer-product")
ax.plot(x, y3, label="row-based dataflow 1")
ax.plot(x, y4, label="row-based dataflow 2")
ax.set_xticks(custom_x)
ax.set_yticks(custom_y)
ax.set_xticklabels(custom_x, rotation=45, fontsize='large')
ax.set_yticklabels(custom_y, fontsize='medium')
font_dict = {'family': 'serif', 'color':  '0.0', 'weight': 'normal', 'size': 15}
ax.set_xlabel("Datasize (KB)",fontdict=font_dict)
ax.set_ylabel("Execution time (s)",fontdict=font_dict)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_label_coords(0.5, -0.26)
ax.yaxis.set_label_coords(-0.05, 0.5)
ax.set_ylim(0,3)
ax.set_xlim(0,1000)
plt.subplots_adjust(bottom=0.25, left=0.08, right=0.98)

ax.legend()
fig.savefig("dataflow_data.png")
