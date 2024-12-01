from main import pre_process, init_element_list, fill_element_list
from element import element
from typing import List
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

def plot_element(e_list:List[element]):
  fig, ax = plt.subplots()
  color = ["red", "green", "black", "purple"]
  #j = 0
  for i in range(216):
    element = e_list[i]
    poly = Polygon(element.node_coordinate, closed=True, facecolor="none",edgecolor=color[i%4])
    #x = 0.5*(element.node_coordinate[0, 0] + element.node_coordinate[2, 0])
    #y = 0.5*(element.node_coordinate[0, 0] + element.node_coordinate[1, 0])
    #ax.text(x, y, j, fontsize=8)
    ax.add_patch(poly)
    #j = j + 1
  ax.set_xlim([-60.0, 60.0])
  ax.set_ylim([-60.0, 60.0])
  plt.axis("equal")
  plt.show()


# theta > 0 && theta < pi/18 && r >40 && r <43.615
# 对模型进行前处理
data = pre_process()
# 初始化四边形单元列表
element_list = init_element_list()
print(len(element_list))
# 根据原子的极坐标将其逐一填充进element_list中
fill_element_list(data, element_list)
# 绘制模型
plot_element(element_list)