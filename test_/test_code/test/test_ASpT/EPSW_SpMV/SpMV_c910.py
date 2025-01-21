import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
data = pd.read_excel('c910_spmv_old.xlsx')
color = ["#014636", "#016c59", "#02818a", "#3690c0", "#67a9cf", "#a6bddb"]

parameter = len(data.iloc[:,0]) - 12  # all data number
kb = data.iloc[0:parameter,0].values
kb = list(map(int, kb))
base = data.iloc[0:parameter,1].values #csr_sca
aspt = data.iloc[0:parameter,11].values #aspt_sca
csr = data.iloc[0:parameter,12].values #csr_m2
epsw_2 = data.iloc[0:parameter,13].values #epsw1_m2
epsw_8 = data.iloc[0:parameter,14].values #epsw2_m8

x=[]
for i in range (parameter):
    x.append(i+1)

bytes = np.array(bytes)
kb = np.array(kb)
base = np.array(base)
aspt = np.array(aspt)
csr = np.array(csr)
epsw_2 = np.array(epsw_2)
epsw_8 = np.array(epsw_8)

bar_width = 0.16
plt.figure(figsize=(15, 8))
# fig, ax= plt.subplots()
plt.bar(np.array(x), base, bar_width, color=color[0], label="Base", )
plt.bar(np.array(x)+bar_width, aspt, bar_width, color=color[1], label="ASpT")
plt.bar(np.array(x)+2*bar_width, csr, bar_width, color=color[2], label="CSR")
plt.bar(np.array(x)+3*bar_width, epsw_2, bar_width, color=color[4], label="EPSW_m2")
plt.bar(np.array(x)+4*bar_width, epsw_8, bar_width, color=color[5], label="EPSW_m8")
ax = plt.gca()
# plt.set_title('SpMM performance comparation on C906')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xticks(x)
ax.set_xticklabels(kb, rotation=30)
plt.xticks(x, kb, fontsize=18)
plt.yticks(fontsize=18)
plt.ylim(0, 2)
# current_xlim = plt.xlim()
# plt.xlim(current_xlim[0], current_xlim[1])

# plt.set_xlim(-1, len(parameter)-1)
legend = ax.legend(bbox_to_anchor=(0.3, 1.1), ncol=2, prop = {'size':22})
# legend.set_size(2)

# plt.tight_layout()
# plt.ypltis.set_label_coords(-0.04, 1)
# plt.xpltis.set_label_coords(1, -0.11)
font_dict = {'family': 'serif', 'color':  '0.0', 'weight': 'normal', 'size': 22}
plt.ylabel("Execution time compared to CSR", loc="center", alpha = 1, fontdict=font_dict, labelpad=9)
plt.xlabel("Used memory (KB)", loc="center", alpha = 1, fontdict=font_dict, labelpad=3)
plt.subplots_adjust(bottom=0.12, left=0.08, right=0.99)

                                  

# fig.set_size_inches(18, 5)

# 显示图形
plt.savefig("spmv_c910.png")
