import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

file_path = 'data.txt'

# Initialize lists to hold coordinates and values
x_coords = []
y_coords = []
z_values_all = []

# Read the data from the file
print('Reading data from file...')
with open(file_path, 'r') as file:
    lines = file.readlines()

    for i in range(0, len(lines), 4):
        # Parse x, y coordinates
        x, y = map(float, lines[i].split())
        x_coords.append(x)
        y_coords.append(y)

        # Parse and store z values for each state
        z_values = [float(val) for val in lines[i + 1].split()]
        z_values_all.append(z_values)
print('Data reading complete.')

# Convert z_values_all to a numpy array for efficient indexing
z_values_all = np.array(z_values_all)

# 准备3D绘图
fig = plt.figure(dpi=300)
ax = fig.add_subplot(111, projection='3d')

# 使用颜色映射
cmap = plt.get_cmap('viridis')
num_eigenstates = 100

# 绘制每个本征态的结果
for j in range(num_eigenstates):
    color = cmap(float(j) / num_eigenstates)  # 根据本征态编号计算颜色
    ax.plot(x_coords, y_coords, z_values_all[:, j], color=color, zorder=num_eigenstates-j)
    print(f'Plotting eigenstate {j + 1} of {num_eigenstates}...')

# 设置标签和标题
ax.set_xlabel('X Coordinate')
ax.set_ylabel('Y Coordinate')
ax.set_zlabel('Z Value')
ax.set_title('All Eigenstates Visualization')

# 保存图像
plt.savefig('3dplot_all_eigenstates.png')
plt.close(fig)  # 关闭图形以释放内存
