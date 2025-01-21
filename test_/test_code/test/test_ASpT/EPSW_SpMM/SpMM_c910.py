import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
data = pd.read_excel('c910_spmm_11.28.xlsx')
color = ["#014636", "#016c59", "#02818a", "#3690c0", "#67a9cf", "#a6bddb"]

parameter = len(data.iloc[:,0]) - 2  # all data number
bytes = data.iloc[0:parameter,0].values
kb = data.iloc[0:parameter,1].values
kb = list(map(int, kb))
base = data.iloc[0:parameter,2].values #csr_sca
aspt = data.iloc[0:parameter,3].values #aspt_sca
csr = data.iloc[0:parameter,6].values #csr_m2
ell = data.iloc[0:parameter,7].values #epsw0_m8
epsw_1 = data.iloc[0:parameter,9].values #epsw1_m8
epsw_2 = data.iloc[0:parameter,11].values #epsw2_m8

x=[]
for i in range (parameter):
    x.append(i+1)

bytes = np.array(bytes)
kb = np.array(kb)
base = np.array(base)
aspt = np.array(aspt)
csr = np.array(csr)
ell = np.array(ell)
epsw_1 = np.array(epsw_1)
epsw_2 = np.array(epsw_2)

bar_width = 0.16
plt.figure(figsize=(15, 4))
# fig, ax= plt.subplots()
plt.bar(np.array(x), base, bar_width, color=color[0], label="Base", )
plt.bar(np.array(x)+bar_width, aspt, bar_width, color=color[1], label="ASpT")
plt.bar(np.array(x)+2*bar_width, csr, bar_width, color=color[2], label="CSR")
plt.bar(np.array(x)+3*bar_width, ell, bar_width, color=color[3], label="ELLPACK")
plt.bar(np.array(x)+4*bar_width, epsw_1, bar_width, color=color[4], label="EPSW_s")
plt.bar(np.array(x)+5*bar_width, epsw_2, bar_width, color=color[5], label="EPSW_ada")
ax = plt.gca()
# plt.set_title('SpMM performance comparation on C906')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xticks(x)
ax.set_xticklabels(kb, rotation=30)
plt.xticks(x, kb, fontsize=12)
plt.yticks(fontsize=12)
plt.ylim(0, 2)
# current_xlim = plt.xlim()
# plt.xlim(current_xlim[0], current_xlim[1])

# plt.set_xlim(-1, len(parameter)-1)
plt.legend(bbox_to_anchor=(0.7, 0.9), ncol=3, prop = {'size':12})

# plt.tight_layout()
# plt.ypltis.set_label_coords(-0.04, 1)
# plt.xpltis.set_label_coords(1, -0.11)
font_dict = {'family': 'serif', 'color':  '0.0', 'weight': 'normal', 'size': 15}
plt.ylabel("Execution time compared to CSR", loc="center", alpha = 1, fontdict=font_dict, labelpad=9)
plt.xlabel("Used memory (KB)", loc="center", alpha = 1, fontdict=font_dict, labelpad=3)
plt.subplots_adjust(bottom=0.19, left=0.07, right=0.98)

                                  

# fig.set_size_inches(18, 5)

# 显示图形
plt.savefig("spmm_c910.png")
