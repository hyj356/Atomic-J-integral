from ovito.io import import_file
from ovito.modifiers import ExpressionSelectionModifier as ESM, CalculateDisplacementsModifier as CDM,AffineTransformationModifier as AFM
import numpy as np
# from ovito.data import CutoffNeighborFinder as CNF

pipeline = import_file('./Data/Deform_*.gz', sort_particles=True)
pipeline.modifiers.append(AFM(transformation=[[1, 0, 0, -178], [0, 1, 0, -72.5], [0, 0, 1, 0]]))
pipeline.modifiers.append(CDM(reference_frame=0))
pipeline.modifiers.append(ESM(expression="(Position.X-50)^2 + (Position.Y - 0.0)^2 < 9"))
data = pipeline.compute(frame=1)

mass,ux,uy = data.particles["Mass"],data.particles["Displacement.X"],data.particles["Displacement.Y"]
mask = (data.particles["Selection"] == 1)
_ux, _uy, _mass= ux[mask], uy[mask], mass[mask]
print(len(ux[mask]))
print(np.sum(_ux*_mass) / np.sum(_mass))
print(np.sum(_uy*_mass) / np.sum(_mass))
# # 以0.3nm为截断半径构建邻居列表
# neighbor = CNF(cutoff=3, data_collection=data)
# # 找出在距离指定坐标的截断半径内部的邻居原子ID
# for neigh in neighbor.find_at([180.0,72.0,0.0]):
#   print(neigh.index+1)
# import vtk
# reader = vtk.vtkUnstructuredGridReader()
# reader.SetFileName("./VTK/stress_tensor.vtk")
# try:
#     reader.Update()
#     print("VTK file format appears to be correct.")
# except:
#     print("There may be an issue with the VTK file format.")
# import numpy as np
# from math import pi
# origin_x, origin_y = 178, 72.5
# dr = 3.615        # 约等于一个晶格常数
# n_dr = 6         # 在r方向上划分了多少份
# r0 = 40           # 小圆的半径
# r1 = r0 + n_dr*dr    # 大圆的半径
# depth = 10.5961   # 盒子的厚度
# n_dθ = 36      # 在theta方向上划分了多少份
# dθ = 2*pi / n_dθ      # dθ大约为10°
# num_of_element = 6 * 36   # 通过上述划分, 可以获得216个四边形单元
# fdata = open("coordinate.txt", mode='w')
# for i in range(n_dr+1):
#     r1 = r0 + i*dr    # 大圆的半径
#     t = np.linspace(0, 2 * np.pi, 37)
#     x = r1 * np.cos(t)
#     y = r1 * np.sin(t)
#     for j in range(n_dθ):
#       fdata.write(f"{i*n_dθ+j} {x[j]:.2f} {y[j]:.2f} 0.00\n")
# fdata.close()