import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
data = pd.read_excel('SDDMM_xt.xlsx')
parameter = data.iloc[:,11].values
data1 = data.iloc[0:len(parameter),:].values
color = ["#ccece6", "#99d8c9", "#66c2a4", "#41ae76", "#238b45", "#006d2c", "#00441b"]

x=[]
for i in range(len(parameter)):
    x.append(i)
kb = data1[:len(parameter), 13]
kb = list(map(int, kb))
ASpT_vec_1 = data1[:, 1]
ASpT_sca_1 = data1[:, 2]
CSR_vec_m8_1 = data1[:,3]
CSR_vec_m2_1 = data1[:,7]
EPSW_0_m8_1 = data1[:,4]
EPSW_1_m8_1 = data1[:,5]
EPSW_2_m8_1 = data1[:,6]
EPSW_0_m2_1 = data1[:,8]
EPSW_1_m2_1 = data1[:,9]
EPSW_2_m2_1 = data1[:,10]

ASpT_vec_1_arr = np.array(ASpT_vec_1)
ASpT_sca_1_arr = np.array(ASpT_sca_1)
CSR_vec_m8_1_arr = np.array(CSR_vec_m8_1)
CSR_vec_m2_1_arr = np.array(CSR_vec_m2_1)
EPSW_0_m8_1_arr = np.array(EPSW_0_m8_1)
EPSW_1_m8_1_arr = np.array(EPSW_1_m8_1)
EPSW_2_m8_1_arr = np.array(EPSW_2_m8_1)
EPSW_0_m2_1_arr = np.array(EPSW_0_m2_1)
EPSW_1_m2_1_arr = np.array(EPSW_1_m2_1)
EPSW_2_m2_1_arr = np.array(EPSW_2_m2_1)

bar_width = 0.12
labelx = 0
x = np.array(x)

plt.figure(figsize=(13, 4))
plt.bar(x, ASpT_vec_1_arr, bar_width, color=color[0], label="ASpT_vec")
plt.bar(x+1*bar_width, ASpT_sca_1_arr, bar_width, color=color[1], label="ASpT_scalar")
plt.bar(x+2*bar_width, CSR_vec_m2_1_arr, bar_width, color=color[2], label="CSR_vec_m2")
plt.bar(x+3*bar_width, EPSW_1_m2_1_arr, bar_width, color=color[3], label="EPSW_1_m2")
plt.bar(x+4*bar_width, EPSW_2_m2_1_arr, bar_width, color=color[4], label="EPSW_2_m2")
plt.bar(x+5*bar_width, EPSW_1_m8_1_arr, bar_width, color=color[5], label="EPSW_1_m8")
plt.bar(x+6*bar_width, EPSW_2_m8_1_arr, bar_width, color=color[6], label="EPSW_2_m8")
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xticks(x)
ax.set_xticklabels(kb)
plt.xticks(x, kb, rotation=30, fontsize=12)
plt.yticks(fontsize=12)
ax.yaxis.set_label_coords(-0.04, 1)
ax.xaxis.set_label_coords(1, -0.27)
font_dict = {'family': 'serif', 'color':  '0.0', 'weight': 'normal', 'size': 15}
ax.set_xlabel("The number of non-zero values (KB)", loc="center", alpha = 1, fontdict=font_dict, labelpad=9)
ax.set_ylabel("Execution time compared to CSR", loc="center", alpha = 1, fontdict=font_dict, labelpad=9)
# ax.set_xlim(-1, len(parameter)-1)
ax.legend(bbox_to_anchor=(0.3, 1), ncol=3)

plt.tight_layout()
plt.subplots_adjust(hspace=0.5, wspace=0.2, bottom=0.21)

# 显示图形
plt.savefig("SDDMM_c910.png")
