import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as bar
data = pd.read_excel('SDDMM_r.xlsx')
parameter = data.iloc[:,12].values
block_size = math.floor(len(parameter)/3)

data1 = data.iloc[0:block_size,:].values
data2 = data.iloc[block_size:2*block_size, :].values
data3 = data.iloc[2*block_size:len(parameter), :].values
base = [1] * len(parameter)

x=[]
for i in range(len(parameter)):
    x.append(i+1)
x1 = x[0:block_size]
x2 = x[block_size:2*block_size]
x3 = x[2*block_size:len(parameter)]

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
base_1 = base[0:block_size]

ASpT_vec_2 = data2[:, 1]
ASpT_sca_2 = data2[:, 2]
CSR_vec_m8_2 = data2[:,3]
CSR_vec_m2_2 = data2[:,7]
EPSW_0_m8_2 = data2[:,4]
EPSW_1_m8_2 = data2[:,5]
EPSW_2_m8_2 = data2[:,6]
EPSW_0_m2_2 = data2[:,8]
EPSW_1_m2_2 = data2[:,9]
EPSW_2_m2_2 = data2[:,10]
base_2 = base[block_size:2*block_size]

ASpT_vec_3 = data3[:, 1]
ASpT_sca_3 = data3[:, 2]
CSR_vec_m8_3 = data3[:,3]
CSR_vec_m2_3 = data3[:,7]
EPSW_0_m8_3 = data3[:,4]
EPSW_1_m8_3 = data3[:,5]
EPSW_2_m8_3 = data3[:,6]
EPSW_0_m2_3 = data3[:,8]
EPSW_1_m2_3 = data3[:,9]
EPSW_2_m2_3 = data3[:,10]
base_3 = base[2*block_size:len(parameter)]

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
base_1 = base[0:block_size]

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
base_1_arr = np.array(base_1)

ASpT_vec_2_arr = np.array(ASpT_vec_2)
ASpT_sca_2_arr = np.array(ASpT_sca_2)
CSR_vec_m8_2_arr = np.array(CSR_vec_m8_2)
CSR_vec_m2_2_arr = np.array(CSR_vec_m2_2)
EPSW_0_m8_2_arr = np.array(EPSW_0_m8_2)
EPSW_1_m8_2_arr = np.array(EPSW_1_m8_2)
EPSW_2_m8_2_arr = np.array(EPSW_2_m8_2)
EPSW_0_m2_2_arr = np.array(EPSW_0_m2_2)
EPSW_1_m2_2_arr = np.array(EPSW_1_m2_2)
EPSW_2_m2_2_arr = np.array(EPSW_2_m2_2)
base_2_arr = np.array(base_2)

ASpT_vec_3_arr = np.array(ASpT_vec_3)
ASpT_sca_3_arr = np.array(ASpT_sca_3)
CSR_vec_m8_3_arr = np.array(CSR_vec_m8_3)
CSR_vec_m2_3_arr = np.array(CSR_vec_m2_3)
EPSW_0_m8_3_arr = np.array(EPSW_0_m8_3)
EPSW_1_m8_3_arr = np.array(EPSW_1_m8_3)
EPSW_2_m8_3_arr = np.array(EPSW_2_m8_3)
EPSW_0_m2_3_arr = np.array(EPSW_0_m2_3)
EPSW_1_m2_3_arr = np.array(EPSW_1_m2_3)
EPSW_2_m2_3_arr = np.array(EPSW_2_m2_3)
base_3_arr = np.array(base_3)

parameter1 = parameter[0:block_size]
parameter2 = parameter[block_size:2*block_size]
parameter3 = parameter[2*block_size:len(parameter)]

bar_width = 0.1
labelx = 0

fig, ax= bar.subplots(nrows=3, ncols=1)
ax[0].bar(np.array(x1), base_1_arr, bar_width, label="base")
ax[0].bar(np.array(x1)+bar_width, ASpT_vec_1_arr, bar_width, label="ASpT_vec")
ax[0].bar(np.array(x1)+2*bar_width, ASpT_sca_1_arr, bar_width, label="ASpT_scalar")
ax[0].bar(np.array(x1)+3*bar_width, CSR_vec_m2_1_arr, bar_width, label="CSR_vec_m2")
# ax[0].bar(np.array(x1)+4*bar_width, EPSW_0_m2_1_arr, bar_width, label="EPSW_0_m2")
ax[0].bar(np.array(x1)+4*bar_width, EPSW_1_m2_1_arr, bar_width, label="EPSW_1_m2")
ax[0].bar(np.array(x1)+5*bar_width, EPSW_2_m2_1_arr, bar_width, label="EPSW_2_m2")
# ax[0].bar(np.array(x1)+7*bar_width, CSR_vec_m8_1_arr, bar_width, label="CSR_vec_m8")
# ax[0].bar(np.array(x1)+4*bar_width, EPSW_0_m8_1_arr, bar_width, label="EPSW_0_m8")
ax[0].bar(np.array(x1)+6*bar_width, EPSW_1_m8_1_arr, bar_width, label="EPSW_1_m8")
ax[0].bar(np.array(x1)+7*bar_width, EPSW_2_m8_1_arr, bar_width, label="EPSW_2_m8")
# ax.set_title('SpMM performance comparation on C906')
ax[0].spines['right'].set_visible(False)
ax[0].spines['top'].set_visible(False)
ax[0].set_xticks(x1)
ax[0].set_xticklabels(parameter1, rotation=30, fontsize='small')
ax[0].yaxis.set_label_coords(0, 1)
ax[0].xaxis.set_label_coords(0.98, -0.05)
ax[0].set_xlabel("nnz", ha='left')
ax[0].set_ylabel("%", ha='left')
# ax.set_xlim(-1, len(parameter)-1)
ax[0].legend(bbox_to_anchor=(0.3, 1), ncol=3)


ax[1].bar(np.array(x2), base_2_arr, bar_width, label="base")
ax[1].bar(np.array(x2)+bar_width, ASpT_vec_2_arr, bar_width, label="ASpT_vec")
ax[1].bar(np.array(x2)+2*bar_width, ASpT_sca_2_arr, bar_width, label="ASpT_scalar")
ax[1].bar(np.array(x2)+3*bar_width, CSR_vec_m2_2_arr, bar_width, label="CSR_vec_m2")
# ax[1].bar(np.array(x2)+4*bar_width, EPSW_0_m2_2_arr, bar_width, label="EPSW_0_m2")
ax[1].bar(np.array(x2)+4*bar_width, EPSW_1_m2_2_arr, bar_width, label="EPSW_1_m2")
ax[1].bar(np.array(x2)+5*bar_width, EPSW_2_m2_2_arr, bar_width, label="EPSW_2_m2")
# ax[1].bar(np.array(x2)+7*bar_width, CSR_vec_m8_2_arr, bar_width, label="CSR_vec_m8")
# ax[1].bar(np.array(x2)+4*bar_width, EPSW_0_m8_2_arr, bar_width, label="EPSW_0_m8")
ax[1].bar(np.array(x2)+6*bar_width, EPSW_1_m8_2_arr, bar_width, label="EPSW_1_m8")
ax[1].bar(np.array(x2)+7*bar_width, EPSW_2_m8_2_arr, bar_width, label="EPSW_2_m8")
ax[1].spines['right'].set_visible(False)
ax[1].spines['top'].set_visible(False)
ax[1].set_xticks(x2)
ax[1].set_xticklabels(parameter2, rotation=30, fontsize='small')
ax[1].yaxis.set_label_coords(0, 1)
ax[1].xaxis.set_label_coords(0.98, -0.05)
ax[1].set_xlabel("nnz", ha='left')
ax[1].set_ylabel("%", ha='left')
# ax.set_xlim(-1, len(parameter)-1)
# ax[1].legend()

ax[2].bar(np.array(x3), base_3_arr, bar_width, label="base")
ax[2].bar(np.array(x3)+bar_width, ASpT_vec_3_arr, bar_width, label="ASpT_vec")
ax[2].bar(np.array(x3)+2*bar_width, ASpT_sca_3_arr, bar_width, label="ASpT_scalar")
ax[2].bar(np.array(x3)+3*bar_width, CSR_vec_m2_3_arr, bar_width, label="CSR_vec_m2")
# ax[2].bar(np.array(x3)+4*bar_width, EPSW_0_m2_3_arr, bar_width, label="EPSW_0_m2")
ax[2].bar(np.array(x3)+4*bar_width, EPSW_1_m2_3_arr, bar_width, label="EPSW_1_m2")
ax[2].bar(np.array(x3)+5*bar_width, EPSW_2_m2_3_arr, bar_width, label="EPSW_2_m2")
# ax[2].bar(np.array(x3)+7*bar_width, CSR_vec_m8_3_arr, bar_width, label="CSR_vec_m8")
# ax[2].bar(np.array(x3)+4*bar_width, EPSW_0_m8_3_arr, bar_width, label="EPSW_0_m8")
ax[2].bar(np.array(x3)+6*bar_width, EPSW_1_m8_3_arr, bar_width, label="EPSW_1_m8")
ax[2].bar(np.array(x3)+7*bar_width, EPSW_2_m8_3_arr, bar_width, label="EPSW_2_m8")


# ax.set_title('SpMM performance comparation on C906')
ax[2].spines['right'].set_visible(False)
ax[2].spines['top'].set_visible(False)
ax[2].set_xticks(x3)
ax[2].set_xticklabels(parameter3, rotation=30, fontsize='small')
ax[2].yaxis.set_label_coords(0, 1)
ax[2].xaxis.set_label_coords(0.98, -0.05)
ax[2].set_xlabel("nnz", ha='left')
ax[2].set_ylabel("%", ha='left')
# ax.set_xlim(-1, len(parameter)-1)
# ax[2].legend()

bar.tight_layout()
bar.subplots_adjust(hspace=0.5, wspace=0.2)
# fig = bar.figure()
fig.set_size_inches(20, 6)

# ax.xlabel('10000B')
# ax.ylabel('ms')

# 显示图形
bar.show()
