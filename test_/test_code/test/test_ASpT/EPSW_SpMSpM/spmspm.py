import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
data = pd.read_excel('spmspm_12.1_csr_m8.xlsx')
color = ["#014636", "#016c59", "#02818a", "#3690c0", "#67a9cf", "#a6bddb"]

parameter = len(data.iloc[:,0])  # all data number
kb = data.iloc[0:parameter,3].values
kb = list(map(int, kb))
sparsity = data.iloc[0:parameter,4].values
epsw = data.iloc[0:parameter,5].values #csr_sca
csr_spmspm = data.iloc[0:parameter,6].values #csr_m2
csr_spmm = data.iloc[0:parameter,7].values #epsw0_m2

x=[]
for i in range (parameter):
    x.append(i+1)

kb = np.array(kb)
sparsity = np.array(sparsity)
epsw = np.array(epsw)
csr_spmspm = np.array(csr_spmspm)
csr_spmm = np.array(csr_spmm)


plt.figure(figsize=(7, 3))
# fig, ax= plt.subplots()
plt.plot(np.array(x),csr_spmm,color=color[1], label="CSR_SpMM",marker='D')
plt.plot(np.array(x),csr_spmspm,color=color[3], label="CSR_SpMSpM",marker='s')
plt.plot(np.array(x),epsw, color=color[5], label="EPSW_SpMSpM", marker='o')
ax = plt.gca()
# plt.set_title('SpMM performance comparation on C906')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xticks(x)
ax.set_xticklabels(kb, rotation=30)
plt.xticks(x, kb, fontsize=12)
plt.yticks(fontsize=12)
# plt.ylim(0, 2)
# current_xlim = plt.xlim()
# plt.xlim(current_xlim[0], current_xlim[1])

# plt.set_xlim(-1, len(parameter)-1)
plt.legend(bbox_to_anchor=(0.6, 1), ncol=1, prop = {'size':12})

# plt.tight_layout()
# plt.ypltis.set_label_coords(-0.04, 1)
# plt.xpltis.set_label_coords(1, -0.11)
font_dict = {'family': 'serif', 'color':  '0.0', 'weight': 'normal', 'size': 15}
plt.ylabel("Execution time (ms)", loc="center", alpha = 1, fontdict=font_dict, labelpad=9)
plt.xlabel("Used memory (KB)", loc="center", alpha = 1, fontdict=font_dict, labelpad=3)
plt.subplots_adjust(bottom=0.25, left=0.12, right=0.98)

                                  

# fig.set_size_inches(18, 5)

# 显示图形
plt.savefig("spmspm.png")
