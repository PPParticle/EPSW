import matplotlib.pyplot as plt
import numpy as np

# 准备数据

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

# E32M8_vi = [8, 22, 93, 216, 229, 279, 411, 375, 421, 469]
# E32M4_vi = [12, 28, 156, 251, 256, 308, 367, 413, 463, 521]
# E32M2_vi = [16, 33, 183, 135, 166, 200, 233, 267, 300, 333]
# E32M1_vi = [29, 61, 203, 270, 337, 405, 473, 540, 608, 676]
# E32M8_vi_arr = np.array(E32M8_vi)
# E32M4_vi_arr = np.array(E32M4_vi)
# E32M2_vi_arr = np.array(E32M2_vi)
# E32M1_vi_arr = np.array(E32M1_vi)

# E32M8_vl = [8, 25, 138, 182, 298, 393, 489, 374, 424, 469]
# E32M4_vl = [9, 23, 101, 133, 184, 201, 236, 272, 306, 339]
# E32M2_vl = [12, 28, 109, 127, 229, 190, 223, 257, 285, 319]
# E32M1_vl = [23, 48, 186, 250, 316, 372, 434, 497, 558, 620]
# E32M8_vl_arr = np.array(E32M8_vl)
# E32M4_vl_arr = np.array(E32M4_vl)
# E32M2_vl_arr = np.array(E32M2_vl)
# E32M1_vl_arr = np.array(E32M1_vl)

fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10, 6))
# custom_x_ticks = [32768, 65536, 98304, 131072, 163840, 196608, 229376, 262144, 294912, 327680]  
custom_x_labels = x_labels
custom_x_ticks = x 
# custom_x_labels = x
bar_width = 0.2
labelx = -0.3


axes[0].bar(np.array(x), E32M8_vf_arr, bar_width, label="E32M8" )
axes[0].bar(np.array(x)+bar_width, E32M4_vf_arr, bar_width, label="E32M4"  )
axes[0].bar(np.array(x)+2*bar_width, E32M2_vf_arr, bar_width, label="E32M2"  )
axes[0].bar(np.array(x)+3*bar_width, E32M1_vf_arr, bar_width, label="E32M1"  )
# axes[0].set_title('float computation')
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)

axes[0].yaxis.set_label_coords(-0.01, 1.05)
axes[0].xaxis.set_label_coords(0.98, -0.05)
axes[0].legend()
# axes[0].set_xticklabels(custom_x_labels)

axes[1].bar(np.array(x), E32M8_vf_arr10, bar_width, label="E32M8"  )
axes[1].bar(np.array(x)+bar_width, E32M4_vf_arr10, bar_width, label="E32M4"  )
axes[1].bar(np.array(x)+2*bar_width, E32M2_vf_arr10, bar_width, label="E32M2"  )
axes[1].bar(np.array(x)+3*bar_width, E32M1_vf_arr10, bar_width, label="E32M1"  )
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)
axes[1].set_xticks(custom_x_ticks)
axes[1].set_xticklabels(custom_x_labels, rotation=60, fontsize='small')
axes[1].yaxis.set_label_coords(-0.01, 1.05)
axes[1].xaxis.set_label_coords(0.98, -0.05)
axes[1].legend()
# axes[1].set_xticklabels(custom_x_labels)

plt.tight_layout()
plt.subplots_adjust(wspace=0.1)

# 显示图形
plt.savefig("test_vf.png")
