import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re


data = pd.read_excel('parallelism_0.7_v1_初步整理.xlsx', sheet_name='Sheet1')
# print(data.columns.tolist) ['n', 'storage', 'deviation', 'MIN', 'MAX', 'in_row', 'cross_line',...)
labels_x = []
data_list = []
in_row_list = []
cross_line_list = []
for index, row in data.iterrows():
    # print(index)
    text_n = row['n']
    text_s = row['storage']
    text_min = row['MIN']
    text_max = row['MAX']
    text_inrow = row['in_row']
    text_crossline = row['cross_line']
    text_data = row.drop(['n', 'storage', 'deviation', 'MIN', 'MAX', 'in_row', 'cross_line'])
    text_data_list = text_data.tolist()
    text_data_list = text_data_list[:int(text_n)]
    # data_list.append(text_data_list)
    # in_row_list.append(text_inrow)
    # cross_line_list.append(text_crossline)
    if text_n == 16 :
        data_list.append(text_data_list)
        in_row_list.append(text_inrow)
        cross_line_list.append(text_crossline)

    # if(index==0):
    #     print(text_data_list)


    # if text_s :
    #     labels_x.append(text_s) #labels_x ['1KB', '10KB', '256KB', '1024KB', '4096KB', '16384KB']
    # print(text_min)
inrow_list_s = [x/100000 for x in in_row_list] #0.1s
crossline_list_s = [x/100000 for x in cross_line_list]

# x_tick_labels = set(labels_x)
# x_tick_labels = labels_x
# # tmp = np.array(x_tick_labels)
# # tmp.reshape((1,120))
# # x_tick_labels = tmp
# print(x_tick_labels.shape)
# 绘制箱图
fig, ax1 = plt.subplots(figsize=(20,6))

ax1.boxplot(data_list, showfliers=False)
# ax1.set_xticks(range(len(x_tick_labels)))
# ax1.set_xticklabels(x_tick_labels)

ax2 = ax1.twinx()
print(inrow_list_s)
ax2.plot([inrow_list_s], linestyle='-', color='yellow')
ax2.plot([crossline_list_s], linestyle='-', color='blue')
ax2.xaxis.set_visible(False)

# 添加标题和轴标签
ax1.set_title('Boxplot and Line Plot')
ax1.set_xlabel('Storage')
ax1.set_ylabel('Sparstity')
ax2.set_ylabel('time')

# 显示图形
plt.show()


# # 绘制图表
# fig, ax = plt.subplots()
# ax.bar(x, data, align='center', alpha=0.5, ecolor='black', capsize=10)

# # 配置图表
# ax.set_ylabel('Value')
# ax.yaxis.grid(True)

# # 显示图表
# plt.tight_layout()
# plt.show()