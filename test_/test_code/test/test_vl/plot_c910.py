import matplotlib.pyplot as plt
import numpy as np

# 准备数据
x_labels = [32768, 65536, 98304, 131072, 163840, 196608, 229376, 262144, 294912, 327680]  # x轴数据
x = [1,2,3,4,5,6,7,8,9,10]
e32m8_vf = [5,  11, 16,  22,  24,  29,  35,  39,  45,  50]  # y轴数据
e32m4_vf = [31, 87, 156, 177, 157, 123, 196, 224, 227, 232]
e32m2_vf = [24, 67, 83,  100, 128, 137, 183, 188, 185, 176]
e32m1_vf = [22, 68, 113, 97, 129, 153, 182, 185, 173, 185]
e32m8_vf_arr = np.array(e32m8_vf)
e32m4_vf_arr = np.array(e32m4_vf)
e32m2_vf_arr = np.array(e32m2_vf)
e32m1_vf_arr = np.array(e32m1_vf)

e32m8_vi = [2,  11, 18,  25,  37,  41,  47,  54,  51,  52]
e32m4_vi = [15, 31, 50,  67,  92,  117, 141, 138, 173, 176]
e32m2_vi = [27, 55, 68,  112, 124, 145, 202, 170, 205, 222]
e32m1_vi = [33, 69, 109, 122, 173, 200, 197, 209, 229, 218]
e32m8_vi_arr = np.array(e32m8_vi)
e32m4_vi_arr = np.array(e32m4_vi)
e32m2_vi_arr = np.array(e32m2_vi)
e32m1_vi_arr = np.array(e32m1_vi)

e32m8_vl = [4, 7, 9, 20, 15, 38, 22, 47, 56, 45]
e32m4_vl = [12, 24, 40, 55, 68, 84, 86, 115, 129, 118]
e32m2_vl = [16, 37, 58, 85, 122, 129, 139, 156, 175, 190]
e32m1_vl = [17, 39, 64, 36, 111, 166, 115, 204, 190, 184]
e32m8_vl_arr = np.array(e32m8_vl)
e32m4_vl_arr = np.array(e32m4_vl)
e32m2_vl_arr = np.array(e32m2_vl)
e32m1_vl_arr = np.array(e32m1_vl)

fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(20, 3))
# custom_x_ticks = [32768, 65536, 98304, 131072, 163840, 196608, 229376, 262144, 294912, 327680]  
custom_x_labels = x_labels
custom_x_ticks = x 
# custom_x_labels = x
bar_width = 0.2
labelx = -0.3


axes[0].bar(np.array(x), e32m8_vf_arr, bar_width, label="e32m8" )
axes[0].bar(np.array(x)+bar_width, e32m4_vf_arr, bar_width, label="e32m4"  )
axes[0].bar(np.array(x)+2*bar_width, e32m2_vf_arr, bar_width, label="e32m2"  )
axes[0].bar(np.array(x)+3*bar_width, e32m1_vf_arr, bar_width, label="e32m1"  )
axes[0].set_title('float computation')
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)
axes[0].set_xticks(custom_x_ticks)
axes[0].set_xticklabels(custom_x_labels, rotation=60, fontsize='small')
axes[0].yaxis.set_label_coords(-0.01, 1.05)
axes[0].xaxis.set_label_coords(0.98, -0.05)
axes[0].legend()
axes[0].set_xlabel("Byte", ha='left')
axes[0].set_ylabel("ms", ha='left')
# axes[0].set_xticklabels(custom_x_labels)


axes[1].bar(np.array(x), e32m8_vi_arr, bar_width, label="e32m8"  )
axes[1].bar(np.array(x)+bar_width, e32m4_vi_arr, bar_width, label="e32m4"  )
axes[1].bar(np.array(x)+2*bar_width, e32m2_vi_arr, bar_width, label="e32m2"  )
axes[1].bar(np.array(x)+3*bar_width, e32m1_vi_arr, bar_width, label="e32m1"  )
axes[1].set_title('integer computation')
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)
axes[1].set_xticks(custom_x_ticks)
axes[1].set_xticklabels(custom_x_labels, rotation=60, fontsize='small')
axes[1].yaxis.set_label_coords(-0.01, 1.05)
axes[1].xaxis.set_label_coords(0.98, -0.05)
axes[1].legend()
axes[1].set_xlabel("Byte", ha='left')
axes[1].set_ylabel("ms", ha='left')
# axes[1].set_xticklabels(custom_x_labels)


axes[2].bar(np.array(x), e32m8_vl_arr, bar_width, label="e32m8"  )
axes[2].bar(np.array(x)+bar_width, e32m4_vl_arr, bar_width, label="e32m4"  )
axes[2].bar(np.array(x)+2*bar_width, e32m2_vl_arr, bar_width, label="e32m2"  )
axes[2].bar(np.array(x)+3*bar_width, e32m1_vl_arr, bar_width, label="e32m1"  )
axes[2].set_title('memory access')
axes[2].spines['top'].set_visible(False)
axes[2].spines['right'].set_visible(False)
axes[2].set_xticks(custom_x_ticks)
axes[2].set_xticklabels(custom_x_labels, rotation=60, fontsize='small')
axes[2].yaxis.set_label_coords(-0.01, 1.05)
axes[2].xaxis.set_label_coords(0.98, -0.05)
axes[2].legend()
axes[2].set_xlabel("Byte", ha='left')
axes[2].set_ylabel("ms", ha='left')

# for j in range(2):
#     axes[j].yaxis.set_label_coords(labelx, 0)
#     axes[j].xaxis.set_label_coords(labelx, 0)

# 添加标题和坐标轴标签
plt.tight_layout()
plt.subplots_adjust(wspace=0.1)
# plt.xlabel('10000B')
# plt.ylabel('ms')

# 显示图形
plt.show()
