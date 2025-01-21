import matplotlib.pyplot as plt
import numpy as np

# 生成一些示例数据
x = [1,2,3,4,5,6,7,8]
y = [111,12,455,567, 1000, 2000, 2123, 4500]

# 创建画布和子图
fig, ax = fig, ax = plt.subplots(figsize=(8, 4)) 

# 设置中线
midline = 5  # 中线位置
ax.axhline(midline, color='gray', linestyle='--')

# 绘制上半部分的折线图，刻度以1000为单位
ax.plot(x, y)
ax.set_yticks(np.arange(0, 11, 1000))
ax.set_ylim(0, 10000)

# 创建第二个y轴用于下半部分，刻度以500为单位
# ax2 = ax.twinx()
# ax2.set_yticks(np.arange(0, 11, 500))
# ax2.set_ylim(0, 5000)

# 显示图形
plt.show()
