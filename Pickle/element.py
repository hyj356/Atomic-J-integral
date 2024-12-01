import numpy as np

class element():
  def __init__(self, dr, dθ, depth):
    self.dr = dr
    self.dθ = dθ
    self.depth = depth
    self.w = 0.0    # 单元的应变能密度
    self.N = 0      # 单元内部原子的个数
    self.pe = 0.0   # 单元内的总势能
    self.dA = 0.0   # 单元的面积
    self.volume = 0.0   # 单元的体积
    self.r = 0.0    # 记录单元的极坐标r
    self.id_list =[]  # 记录单元内部的原子ID
    self.cart_x = [] # 记录单元内部的所有原子的x坐标
    self.cart_y = [] # 记录单元内部的所有原子的y坐标
    self.polar_r = []  # 记录单元内部的所有原子的极坐标r
    self.polar_θ = []  # 记录单元内部的所有原子的极坐标θ
    self.mass = []    # 原子质量
    self.delta_pe = 0.0      # 原子变形前后势能差
    self.s11 = 0.0    # 单元内部的应力分量s11的总和
    self.s22 = 0.0    # 单元内部的应力分量s22的总和
    self.s12 = self.s21 = 0.0 # 单元内部的应力分量s21和s12的总和, 在MD中, S12 = S21
    self.ux = []      # 单元中原子在x方向上的变形量
    self.uy = []      # 单元中原子在y方向上的变形量
    self.node_coordinate = np.zeros((4, 2)) # 记录四个节点的xy坐标

  def _cal_volume(self):
    '''
    计算单元的体积Vm
    '''
    self.volume = self.cal_area() * self.depth

  def cal_area(self):
    '''
    计算单元的面积dAm
    '''
    return 0.5 * ((self.r + self.dr)**2 -self.r**2) * self.dθ