import matplotlib.pyplot as plt
import numpy as np

# 准备数据
color = ["#c7eae5", "#dfc27d", "#80cdc1", "#35978f"]
x_labels = [32, 64, 128, 256, 512, 1024, 4096]  # x轴数据
x = [1,2,3,4,5,6,7]
E32M8_vf = [8, 20, 134, 183, 231, 357, 374]  # y轴数据
E32M4_vf = [16, 28, 150, 206, 258, 377, 412]
E32M2_vf = [17, 36, 99, 136, 169, 205, 272]
E32M1_vf = [30, 62, 207, 271, 341, 409, 544]
E32M8_vf_arr = np.array(E32M8_vf)
E32M4_vf_arr = np.array(E32M4_vf)
E32M2_vf_arr = np.array(E32M2_vf)
E32M1_vf_arr = np.array(E32M1_vf)

E32M8_vf10 = [5,  11, 16,  22,  24,  29,  39]  # y轴数据
E32M4_vf10 = [31, 87, 156, 177, 157, 123, 224]
E32M2_vf10 = [24, 67, 83,  100, 128, 137, 188]
E32M1_vf10 = [22, 68, 113, 97, 129, 153, 185]
E32M8_vf_arr10 = np.array(E32M8_vf10)
E32M4_vf_arr10 = np.array(E32M4_vf10)
E32M2_vf_arr10 = np.array(E32M2_vf10)
E32M1_vf_arr10 = np.array(E32M1_vf10)

fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10, 5))
# custom_x_ticks = [32768, 65536, 98304, 131072, 163840, 196608, 229376, 262144, 294912, 327680]  
custom_x_labels = x_labels
custom_x_ticks = x 
# custom_x_labels = x
bar_width = 0.2
labelx = -0.3
axes[0].bar(np.array(x), E32M8_vf_arr, bar_width, color=color[0], label="E32M8" )
axes[0].bar(np.array(x)+bar_width, E32M4_vf_arr, bar_width, color=color[1],label="E32M4"  )
axes[0].bar(np.array(x)+2*bar_width, E32M2_vf_arr, bar_width, color=color[2],label="E32M2"  )
axes[0].bar(np.array(x)+3*bar_width, E32M1_vf_arr, bar_width, color=color[3],label="E32M1"  )
# axes[0].set_title('Vector-Float Operation on C906')
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)
axes[0].set_xticks([])
axes[0].tick_params(axis='y', labelsize=12)
# axes[0].set_xticklabels(custom_x_labels, rotation=60, fontsize='small')
axes[0].legend()
axes[1].bar(np.array(x), E32M8_vf_arr10, bar_width, color=color[0],label="E32M8"  )
axes[1].bar(np.array(x)+bar_width, E32M4_vf_arr10, bar_width, color=color[1],label="E32M4"  )
axes[1].bar(np.array(x)+2*bar_width, E32M2_vf_arr10, bar_width, color=color[2],label="E32M2"  )
axes[1].bar(np.array(x)+3*bar_width, E32M1_vf_arr10, bar_width, color=color[3],label="E32M1"  )
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)
# axes[1].set_title('Vector-Float Operation on C910')
axes[1].tick_params(axis='y', labelsize=12)
axes[1].set_xticks(custom_x_ticks)
axes[1].set_xticklabels(custom_x_labels, fontsize='large')
axes[1].yaxis.set_label_coords(-0.01, 1.05)
axes[1].xaxis.set_label_coords(1, -0.25)

font_dict = {'family': 'serif', 'color':  '0.0', 'weight': 'normal', 'size': 14}
font_dict_t = {'family': 'serif', 'color':  '0.0', 'weight': 'normal', 'size': 17}
# axes[0].set_xlabel("C906", loc="center", alpha = 1, fontdict=font_dict, labelpad=9)
axes[1].set_xlabel("Used memory (KB)", loc="center", alpha = 1, fontdict=font_dict)
plt.subplots_adjust(hspace=0.02, left=0.07, right=1, bottom=0.16, top=1)
fig.text(0.005, 0.62, 'Execution time (ms)', va='center', rotation='vertical', fontdict=font_dict)
fig.text(0.3, 0.52, 'Vector-Float Operations on C910', va='center', rotation='horizontal', fontdict=font_dict_t)
fig.text(0.3, 0.92, 'Vector-Float Operations on C906', va='center', rotation='horizontal', fontdict=font_dict_t)

# 显示图形
plt.savefig("test_vf.png")
