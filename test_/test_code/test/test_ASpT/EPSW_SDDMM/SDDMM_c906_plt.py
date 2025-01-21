import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
data = pd.read_excel('SDDMM_r.xlsx')
parameter = data.iloc[:,0].values
block_size = math.floor(len(parameter)/3)
base = [1] * len(parameter)

x=data.iloc[:,0].values
y1=data.iloc[:,1].values
y2=data.iloc[:,2].values
y3=data.iloc[:,3].values
y4=data.iloc[:,4].values
y5=data.iloc[:,5].values
y6=data.iloc[:,6].values
y7=data.iloc[:,7].values
y8=data.iloc[:,8].values
y9=data.iloc[:,9].values
y10=data.iloc[:,10].values
y11=data.iloc[:,11].values

fig, ax = plt.subplots(figsize=(20, 4))
ax.plot(x, y1, label="ASpT_vec")
ax.plot(x, y2, label="ASpT_scalar")
# ax.plot(x, y3, label="CSR_vec_m8")
# ax.plot(x, y4, label="EPSW_0_m8")
# ax.plot(x, y5, label="EPSW_1_m8")
ax.plot(x, y6, label="EPSW_2_m8")
ax.plot(x, y7, label="CSR_vec_m2")
# ax.plot(x, y8, label="EPSW_0_m2")
# ax.plot(x, y9, label="EPSW_1_m2")
# ax.plot(x, y10, label="EPSW_2_m2")
ax.plot(x, y11, label="base")
ax.legend()
# plt.subplots_adjust(hspace=0.5, wspace=0.2)
# fig = plt.figure()
# ax.xlabel('10000B')
# ax.ylabel('ms')

# 显示图形
plt.show()
