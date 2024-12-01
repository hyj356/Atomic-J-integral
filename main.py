from ovito.io import import_file, export_file
from ovito.data import DataCollection
from ovito.modifiers import ComputePropertyModifier as CPM, ExpressionSelectionModifier as ESM, InvertSelectionModifier as ISM, DeleteSelectedModifier as DM, AffineTransformationModifier as AFM, CalculateDisplacementsModifier as CDM
from math import pi, cos, sin, sqrt
from element import element
import numpy as np
import pickle
from typing import List

# 裂纹尖端的xy坐标, 这里我们需要将整个模型平移以使得裂纹尖端与坐标原点对齐, 便于计算和公式推导
origin_x, origin_y = 178.5, 72.5
dr = 3.615        # 约等于一个晶格常数
n_dr = 6         # 在r方向上划分了多少份
r0 = 40           # 小圆的半径
r1 = r0 + n_dr*dr    # 大圆的半径
depth = 10.5961   # 盒子的厚度
n_dθ = 36      # 在theta方向上划分了多少份
dθ = 2*pi / n_dθ      # dθ大约为10°
num_of_element = n_dr * n_dθ   # 通过上述划分, 可以获得216个四边形单元

def pre_process():
  pipeline = import_file("./Data/Deform_*.gz", input_format="lammps/dump",
  sort_particles=True)
  # 将原子进行平移, 这样可以把极坐标原点与直角坐标原点对齐
  pipeline.modifiers.append(AFM(transformation=[[1, 0, 0, -origin_x], [0, 1, 0, -origin_y], [0, 0, 1, 0]]))
  # 计算平移之后的极坐标r
  pipeline.modifiers.append(CPM(expressions="sqrt( Position.X^2+Position.Y^2+Position.Z^2)",
                                output_property="r"))
  # 计算极坐标系下的theta, 这时候我们希望theta=[0,2PI]
  pipeline.modifiers.append(ESM(expression="Position.X > 0 && Position.Y > 0"))
  pipeline.modifiers.append(CPM(expressions="atan(Position.Y/Position.X)", 
                                output_property="theta", only_selected=True))
  pipeline.modifiers.append(ESM(expression="Position.X > 0 && Position.Y < 0"))
  pipeline.modifiers.append(CPM(expressions="atan(Position.Y/Position.X)+2*pi", 
                                output_property="theta", only_selected=True))
  pipeline.modifiers.append(ESM(expression="Position.X < 0"))
  pipeline.modifiers.append(CPM(expressions="atan(Position.Y/Position.X)+pi", 
                                output_property="theta", only_selected=True))
  # 计算平均化之后的应力张量S11, S22和S12, 其中S12 = S21
  pipeline.modifiers.append(CPM(cutoff_radius=5.0, expressions="v_sy/(NumNeighbors+1)",
                                neighbor_expressions="v_sy/(NumNeighbors+1)",
                                output_property="ave_syy"))
  pipeline.modifiers.append(CPM(cutoff_radius=5.0, expressions="v_sx/(NumNeighbors+1)",
                                neighbor_expressions="v_sx/(NumNeighbors+1)",
                                output_property="ave_sxx"))
  pipeline.modifiers.append(CPM(cutoff_radius=5.0, expressions="v_sxy/(NumNeighbors+1)",
                                neighbor_expressions="v_sxy/(NumNeighbors+1)",
                                output_property="ave_sxy"))
  pipeline.modifiers.append(ESM(expression=f'r > {r0} && r < {r1}'))
  pipeline.modifiers.append(ISM())
  pipeline.modifiers.append(DM())
  data = pipeline.compute(frame=1)

  # 导出文件
  export_file(data=data, file="pre_process.gz", format='lammps/dump',columns =
    ['Particle Identifier', 'Particle Type', "r", "theta", 'Position.X', 'Position.Y', 'Position.Z',
    "ave_syy", "ave_sxx", "ave_sxy", "Mass", "c_pot"])
  # 返回计算完成的数据集
  return data

def fill_element_list(data:DataCollection, e_list:List[element]):
  '''
  此函数将对处理之后的轨迹文件进行遍历处理, 将相关信息逐一填充进去
  '''
  # 遍历所有原子的r
  for i, r in enumerate(data.particles["r"]):
    theta = data.particles["theta"][i]
    int_r, int_θ = round((r-r0)//dr), round(theta//dθ)
    # 计算第i个原子应该位于第几个element
    index = int_r*n_dθ + int_θ
    # 将第i个原子的xy坐标填充到对应的列表中
    e_list[index].cart_x.append(data.particles['Position.X'][i])
    e_list[index].cart_y.append(data.particles['Position.Y'][i])
    # 将第i个原子的极坐标填充到对应的列表中
    e_list[index].polar_r.append(data.particles['r'][i])
    e_list[index].polar_θ.append(data.particles['theta'][i])
    # 将第i个原子的原子id写入到对应的list中
    e_list[index].id_list.append(data.particles['Particle Identifier'][i])
    # 增加单元内部原子的数量
    e_list[index].N += 1

def init_element_list():
  '''
  此函数用于初始化单元的4个节点的全局xy坐标
  '''
  element_list = [element(dr=dr, dθ=dθ, depth=depth) for _ in range(num_of_element)]
  for i in range(n_dr):
    for j in range(n_dθ):
      # 计算极坐标
      r, theta = r0 + i * dr, j * dθ
      # 按照逆时针方向计算第一个点的xy坐标
      element_list[i*n_dθ + j].node_coordinate[0, 0] = r * cos(theta)
      element_list[i*n_dθ + j].node_coordinate[0, 1] = r * sin(theta)
      # 计算第二个点的xy坐标
      element_list[i*n_dθ + j].node_coordinate[1, 0] = (r + dr) * cos(theta)
      element_list[i*n_dθ + j].node_coordinate[1, 1] = (r + dr) * sin(theta)
      # 计算第三个点的xy坐标
      element_list[i*n_dθ + j].node_coordinate[2, 0] = (r + dr) * cos(theta + dθ)
      element_list[i*n_dθ + j].node_coordinate[2, 1] = (r + dr) * sin(theta + dθ)
      # 计算第四个点的xy坐标
      element_list[i*n_dθ + j].node_coordinate[3, 0] = r * cos(theta + dθ)
      element_list[i*n_dθ + j].node_coordinate[3, 1] = r * sin(theta + dθ)
      # 存储单元的极坐标r
      element_list[i*n_dθ + j].r = i*n_dr
  return element_list

def calculate_strain_energy(e_list:List[element], filename:str, start_frame:int, end_frame:int):
    # 导入文件, 注意sort_particles必须设置为True
    pipeline = import_file(filename, sort_particles=True)
    initial_pe = pipeline.compute(frame=start_frame).particles["c_pot"]
    after_pe = pipeline.compute(frame=end_frame).particles["c_pot"]
    # 遍历每个单元中的原子, 计算势能差
    for element in e_list:
      for id in element.id_list:
        element.delta_pe += after_pe[id-1] - initial_pe[id-1]
      # 计算单元的应变能密度
      element.w = element.delta_pe / element.volume

def calculate_stress_tensor(e_list:List[element], filename:str, frame:int):
  # 导入文件并通过计算获取3个应力张量的数据
  pipeline = import_file(filename, sort_particles=True)
  # 计算平均化之后的应力张量S11, S22和S12, 其中S12 = S21
  pipeline.modifiers.append(CPM(cutoff_radius=5.0, expressions="v_sy/(NumNeighbors+1)",
                                neighbor_expressions="v_sy/(NumNeighbors+1)",
                                output_property="ave_syy"))
  pipeline.modifiers.append(CPM(cutoff_radius=5.0, expressions="v_sx/(NumNeighbors+1)",
                                neighbor_expressions="v_sx/(NumNeighbors+1)",
                                output_property="ave_sxx"))
  pipeline.modifiers.append(CPM(cutoff_radius=5.0, expressions="v_sxy/(NumNeighbors+1)",
                                neighbor_expressions="v_sxy/(NumNeighbors+1)",
                                output_property="ave_sxy"))
  data = pipeline.compute(frame=frame)
  s11, s22, s12 = data.particles["ave_sxx"], data.particles["ave_syy"], data.particles["ave_sxy"]
  # 遍历每个单元, 计算其中的应力张量
  for element in e_list:
    if element.N == 0:
      continue
    else:
      for id in element.id_list:
        element.s11 += s11[id-1]
        element.s22 += s22[id-1]
        element.s12 += s12[id-1]
    # 计算单元的平均应力分量
    element.s11 = element.s11 / element.N
    element.s22 = element.s22 / element.N
    element.s12 = element.s12 / element.N
    # 注意每个单元中的S21 = S12
    element.s21 = element.s12

def calculate_volume_area(e_list:List[element]):
  for element in e_list:
    element.cal_area()
    element.cal_volume()

def write_VTK(e_list:List[element], filename:str):
  fdata = open(file=filename, mode='w')
  fdata.write("# vtk DataFile Version 3.0\n")
  fdata.write("# Visuliaze Stress tensor\n")
  fdata.write("ASCII\nDATASET UNSTRUCTURED_GRID\n\n")
  fdata.write(f"POINTS  {n_dθ * (n_dr+1)} float\n")
  # 写入所有节点坐标
  for i in range(n_dr+1):
    r1 = r0 + i*dr    # 大圆的半径
    t = np.linspace(0, 2 * np.pi, 37)
    x = r1 * np.cos(t)
    y = r1 * np.sin(t)
    for j in range(n_dθ):
      fdata.write(f"{x[j]:.2f} {y[j]:.2f} 0.00\n")
  fdata.write('\n')
  # 写入单元序号
  fdata.write(f"CELLS {n_dr*n_dθ} {5*n_dr*n_dθ}\n")
  for i in range(n_dr):
    for j in range(n_dθ):
      if (i*n_dθ + j + 1) % 36 ==0:
        fdata.write(f"4 {i*n_dθ + j} {i*n_dθ + j + n_dθ} {i*n_dθ + j + 1} {i*n_dθ + j - 35}\n")
      else:
        fdata.write(f"4 {i*n_dθ + j} {i*n_dθ + j + n_dθ} {i*n_dθ + j + n_dθ + 1} {i*n_dθ + j + 1}\n")
  fdata.write('\n')

  # 指明单元类型
  fdata.write(f"CELL_TYPES {n_dr*n_dθ}\n")
  for _ in range(n_dθ*n_dr):
    fdata.write('9\n')

  # 将应力数据写入到VTK文件中
  fdata.write(f"CELL_DATA {n_dθ*n_dr}\n")

  # 写入应力分量s11
  fdata.write(f"SCALARS s11 float 1\nLOOKUP_TABLE default\n")
  for element in e_list:
    fdata.write(f"{element.s11:.7f}\n")
  fdata.write('\n')
  # 写入应力分量s22
  fdata.write(f"SCALARS s22 float 1\nLOOKUP_TABLE default\n")
  for element in e_list:
    fdata.write(f"{element.s22:.7f}\n")
  fdata.write('\n')
  # 写入应力分量s12
  fdata.write(f"SCALARS s12 float 1\nLOOKUP_TABLE default\n")
  for element in e_list:
    fdata.write(f"{element.s12:.7f}\n")
  fdata.write('\n')
  # 写入应力分量s21
  fdata.write(f"SCALARS s21 float 1\nLOOKUP_TABLE default\n")
  for element in e_list:
    fdata.write(f"{element.s21:.7f}\n")

  fdata.close()

def print_element(e_list:List[element]):
  fdata = open("element_list.txt", mode='w')
  for i, element in enumerate(e_list):
    fdata.write(f"{i}, {element.N}, {element.s11}, {element.s22}, {element.s21}\n")
  fdata.close()

def write_element(e_list:List[element], filename:str):
  '''
  此函数将每个单元类中的数据以二进制写入到文件中
  e_list: 含有所有单元的相关数据的列表
  filename: 二进制文件的名称
  返回值: 无
  '''
  with open(file=filename, mode='wb') as f:
    pickle.dump(e_list, f)

def read_element(filename:str) -> List[element]:
  '''
  此函数从二进制文件中读取序列化的单元类列表
  filename: 二进制文件的名称
  返回值: 无
  '''
  with open(file=filename, mode='wb') as f:
    my_restored_list = pickle.load(f)
  return my_restored_list

def calculate_displacement(e_list:List[element], filename:str, start_frame:int, end_frame:int):
  '''
  该函数耗费较长时间, 除非用C/C++, Fortran这类高性能语言重写, 否则没有太好的改善方法
  '''
  # 计算每个原子的位移
  pipeline = import_file(filename, sort_particles=True)
  # 将模型进行平移, 以确保对准原点
  pipeline.modifiers.append(AFM(transformation=[[1, 0, 0, -origin_x], [0, 1, 0, -origin_y], [0, 0, 1, 0]]))
  # 参考第0帧计算位移
  pipeline.modifiers.append(CDM(reference_frame=start_frame))
  # 计算获取位移
  data = pipeline.compute(frame=end_frame)
  mass,u1,u2 = data.particles["Mass"],data.particles["Displacement.X"],data.particles["Displacement.Y"]
  # 构建邻居列表, 原文献使用0.3nm, 也就是3埃米, 这里偷懒直接照抄, 其实不建议照抄
  # 接下来遍历每个单元
  for element in e_list:
    # 遍历每个单元的4个角点
    for i, node in enumerate(element.node_coordinate):
      # 调用expression select获取距离角度在3埃米的截断半径以内的原子, 注意我们假设是二维模型, 所以实际上
      # 选择区域是一个圆柱形的立体区域
      pipeline.modifiers.append(ESM(expression=f"(Position.X-{node[0]})^2 + (Position.Y - {node[1]})^2 < 9"))
      data = pipeline.compute(frame=end_frame)
      mask = (data.particles["Selection"] == 1) # 获取掩码数组
      # 根据掩码数组提取原子在x, y方向上的位移和质量掩码
      _u1, _u2, _mass= u1[mask], u2[mask], mass[mask]
      # 考虑到裂纹附近存在大量空腔, 因此极有可能在指定的截断半径下选不到任何原子, 因此需要分类讨论处理
      if(np.sum(_mass) == 0):
        element.node_u1[i], element.node_u2[i] = 0.0, 0.0
        continue
      else:
        # 根据文献提供的加权平均公式(7)计算单元四个角点的x和y方向位移
        element.node_u1[i] = np.sum(_u1*_mass) / np.sum(_mass)
        element.node_u2[i] = np.sum(_u2*_mass) / np.sum(_mass)

def calculate_J_integral(e_list:List[element]):
  '''
  此函数用于计算样品变形之后的J积分
  '''
  # 初始化J积分的值
  J_integral = 0.0
  # 将ev/Å²转换为J/m²
  unit = 16.02176634  
  # 遍历每个单元
  for element in e_list:
    if element.N <= 6:   # 如果某个有限元单元内部原子数量小于10, 认为它对J积分不起贡献
      continue
    else:
      w = element.w       # 应变能密度w
      s11 = element.s11   # σ11
      partial_u1_x1 = element.cal_partial_u1_x1() # ∂u1/∂x1
      s12 = element.s12   # σ12
      partial_u2_x1 = element.cal_partial_u2_x1() # ∂u2/∂x1
      partial_g_x1 = element.cal_partial_g_x1(r0=r0, r1=r1)   # ∂g/∂x1
      s21 = element.s21   # σ21
      s22 = element.s22   # σ22
      partial_g_x2 = element.cal_partial_g_x2(r0=r0, r1=r1)   # ∂g/∂x2
      dAm = element.dA  # dAm
      J_integral += ((w - s11*partial_u1_x1-s12*partial_u2_x1)*partial_g_x1 + \
                    (s21*partial_u1_x1+s22*partial_u2_x1)*partial_g_x2) * dAm

  return J_integral * unit

if __name__ == "__main__":
  # 对模型进行前处理
  data = pre_process()
  # 初始化四边形单元列表
  element_list = init_element_list()
  print("有限元单元列表初始化完成")
  # 根据原子的极坐标将其逐一填充进element_list中
  fill_element_list(data, element_list)
  print("已将所有原子填充至对应单元内部")
  # 计算每个单元4个节点的位移
  calculate_displacement(element_list, "./Data/Deform_*.gz", start_frame=0,end_frame=1)
  print("单元节点位移计算完成")
  # 计算单元的面积和体积
  calculate_volume_area(e_list=element_list)
  print("单元面积与体积计算完成")
  # 计算变形前后的应变能密度
  calculate_strain_energy(element_list, "./Data/Deform_*.gz", start_frame=0,end_frame=1)
  print("单元应变能密度计算完成")
  # 计算每个单元内部的平均应力分量
  calculate_stress_tensor(element_list, "./Data/Deform_40000.gz", frame=0)
  print("单元应力分量计算完成")
  # 将结果输出到VTK文件中进行可视化
  #write_VTK(e_list=element_list, filename="./VTK/stress_tensor.vtk")
  #print_element(e_list=element_list)
  # 计算J积分
  J_integral = calculate_J_integral(e_list=element_list)
  print(f"The vaule of J integral is {J_integral} J/m²")