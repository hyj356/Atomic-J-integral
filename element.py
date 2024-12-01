import numpy as np
from math import sqrt

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
    self.delta_pe = 0.0      # 原子变形前后势能差
    self.s11 = 0.0    # 单元内部的应力分量s11的总和
    self.s22 = 0.0    # 单元内部的应力分量s22的总和
    self.s12 = self.s21 = 0.0 # 单元内部的应力分量s21和s12的总和, 在MD中, S12 = S21
    self.node_coordinate = np.zeros((4, 2)) # 记录四个节点的xy坐标
    self.node_u1 = np.zeros(4)          # 记录四个节点在加载前后的x方向上的位移
    self.node_u2 = np.zeros(4)          # 记录四个节点在加载前后的y方向上的位移

  def cal_volume(self):
    '''
    计算单元的体积Vm
    '''
    self.volume = self.dA * self.depth

  def cal_area(self):
    '''
    计算单元的面积dAm
    '''
    self.dA = 0.5 * ((self.r + self.dr)**2 -self.r**2) * self.dθ
  def __inv(self, matrix:np.ndarray):
    '''
    此函数用于计算2x2矩阵的逆矩阵
    '''
    a,b,c,d = matrix[0, 0],matrix[0, 1],matrix[1, 0],matrix[1, 1]
    factor = a * d - b * c
    inverse_matrix = np.array([[d / factor, -b / factor],
                               [-c / factor, a / factor]])
    return inverse_matrix
  
  def cal_jacobi_matrix(self, xi:float, eta:float):
    '''
    此函数用于计算单元的雅可比矩阵\n
    xi: 希腊字母ε, 自然坐标之一, 通常表示水平方向, 取值范围在-1到1之间\n
    eta: 希腊字母η, 自然坐标之一, 通常表示竖直方向, 取值范围在-1到1之间
    '''
    jacobi = np.zeros((2, 2))
    # 计算∂x/∂ε
    jacobi[0, 0] = -0.25*self.node_coordinate[0, 0] * (1-eta) + \
                   0.25*self.node_coordinate[1, 0] * (1-eta)  + \
                   0.25*self.node_coordinate[2, 0] * (1+eta)    \
                   -0.25*self.node_coordinate[3, 0] * (1+eta)
    # 计算∂y/∂ε
    jacobi[0, 1] = -0.25*self.node_coordinate[0, 1] * (1-eta) + \
                   0.25*self.node_coordinate[1, 1] * (1-eta)  + \
                   0.25*self.node_coordinate[2, 1] * (1+eta)    \
                   -0.25*self.node_coordinate[3, 1] * (1+eta) 
    # 计算∂x/∂η
    jacobi[1, 0] = -0.25*self.node_coordinate[0, 0] * (1-xi)    \
                   -0.25*self.node_coordinate[1, 0] * (1+xi)  + \
                   0.25*self.node_coordinate[2, 0] * (1+xi)   + \
                   0.25*self.node_coordinate[3, 0] * (1-xi)
    # 计算∂y/∂η
    jacobi[1, 1] = -0.25*self.node_coordinate[0, 1] * (1-xi)    \
                   -0.25*self.node_coordinate[1, 1] * (1+xi)  + \
                   0.25*self.node_coordinate[2, 1] * (1+xi)   + \
                   0.25*self.node_coordinate[3, 1] * (1-xi)
    return jacobi
  
  def cal_jacobi_det(self, xi:float, eta:float):
    '''
    此函数用于计算单元的雅可比行列式的值\n
    xi: 希腊字母ε, 自然坐标之一, 通常表示水平方向, 取值范围在-1到1之间\n
    eta:希腊字母η, 自然坐标之一, 通常表示竖直方向, 取值范围在-1到1之间
    '''
    jacobi = np.zeros((2, 2))
    # 计算∂x/∂ε
    jacobi[0, 0] = -0.25*self.node_coordinate[0, 0] * (1-eta) + \
                   0.25*self.node_coordinate[1, 0] * (1-eta)  + \
                   0.25*self.node_coordinate[2, 0] * (1+eta)    \
                   -0.25*self.node_coordinate[3, 0] * (1+eta)
    # 计算∂y/∂ε
    jacobi[0, 1] = -0.25*self.node_coordinate[0, 1] * (1-eta) + \
                   0.25*self.node_coordinate[1, 1] * (1-eta)  + \
                   0.25*self.node_coordinate[2, 1] * (1+eta)    \
                   -0.25*self.node_coordinate[3, 1] * (1+eta) 
    # 计算∂x/∂η
    jacobi[1, 0] = -0.25*self.node_coordinate[0, 0] * (1-xi)    \
                   -0.25*self.node_coordinate[1, 0] * (1+xi)  + \
                   0.25*self.node_coordinate[2, 0] * (1+xi)   + \
                   0.25*self.node_coordinate[3, 0] * (1-xi)
    # 计算∂y/∂η
    jacobi[1, 1] = -0.25*self.node_coordinate[0, 1] * (1-xi)    \
                   -0.25*self.node_coordinate[1, 1] * (1+xi)  + \
                   0.25*self.node_coordinate[2, 1] * (1+xi)   + \
                   0.25*self.node_coordinate[3, 1] * (1-xi)
    # 行列式的值等于ad-bc
    return jacobi[0, 0]*jacobi[1, 1] - jacobi[0, 1]*jacobi[1, 0]

  
  def __cal_partial_u1_x1(self, xi:float, eta:float):
    '''
    此函数用于计算单元单元在1方向上的应变, 即ε1(ξ, eta) \n
    xi: 希腊字母ε, 自然坐标之一, 通常表示水平方向, 取值范围在-1到1之间\n
    eta:希腊字母η, 自然坐标之一, 通常表示竖直方向, 取值范围在-1到1之间
    '''
    # 首先计算雅可比矩阵
    jacobi = self.cal_jacobi_matrix(xi, eta)
    # 接着计算逆矩阵
    inv_jacobi = self.__inv(matrix=jacobi)
    # 计算形函数对自然坐标的导数, 即[∂Ni/∂ε; ∂Ni/∂η], 其中下标i = 1,2,3,4
    partial_N1_natrue = np.array([-0.25 * (1 - eta), -0.25 * (1 - xi)]) # [∂N1/∂ε; ∂N1/∂η]
    partial_N2_natrue = np.array([ 0.25 * (1 - eta), -0.25 * (1 + xi)]) # [∂N2/∂ε; ∂N2/∂η]
    partial_N3_natrue = np.array([ 0.25 * (1 + eta),  0.25 * (1 + xi)]) # [∂N3/∂ε; ∂N3/∂η]
    partial_N4_natrue = np.array([-0.25 * (1 + eta),  0.25 * (1 - xi)]) # [∂N4/∂ε; ∂N4/∂η]
    # 计算形函数对直角坐标的导数, 即[∂Ni/∂x; ∂Ni/∂y], 其中下标i = 1,2,3,4
    partial_N1_cart = np.matmul(inv_jacobi, partial_N1_natrue)
    partial_N2_cart = np.matmul(inv_jacobi, partial_N2_natrue)
    partial_N3_cart = np.matmul(inv_jacobi, partial_N3_natrue)
    partial_N4_cart = np.matmul(inv_jacobi, partial_N4_natrue)
    # 乘以4个节点位移, 然后累加求和, 获得单元在1方向上的应变
    element_strain1 = partial_N1_cart[0]*self.node_u1[0] + partial_N2_cart[0]*self.node_u1[1] + \
                      partial_N3_cart[0]*self.node_u1[2] + partial_N4_cart[0]*self.node_u1[3]
    return element_strain1
  
  def cal_partial_u1_x1(self):
    '''
    此函数用于计算单元的∂u1/∂x1的值 \n
    xi: 希腊字母ε, 自然坐标之一, 通常表示水平方向, 取值范围在-1到1之间\n
    eta:希腊字母η, 自然坐标之一, 通常表示竖直方向, 取值范围在-1到1之间
    '''
    # 高斯积分点的坐标位置
    W1 = [-1/sqrt(3), -1/sqrt(3), 1/sqrt(3), 1/sqrt(3)]
    W2 = [-1/sqrt(3), 1/sqrt(3), -1/sqrt(3), 1/sqrt(3)]
    Am = self.cal_Am()
    ε11 = 0.0   # ε11 = ∂u1/∂x1
    for ξ, η in zip(W1, W2):
      jacobi = self.cal_jacobi_matrix(xi=ξ, eta=η)
      ε1 = self.__cal_partial_u1_x1(xi=ξ, eta=η)
      ε11 += ε1 * np.linalg.det(jacobi)
    return ε11 / Am

  def __cal_partial_u2_x1(self, xi:float, eta:float):
    '''
    此函数用于计算单元单元在1方向上的应变, 即ε2\n
    xi: 希腊字母ε, 自然坐标之一, 通常表示水平方向, 取值范围在-1到1之间\n
    eta:希腊字母η, 自然坐标之一, 通常表示竖直方向, 取值范围在-1到1之间
    '''
    # 首先计算雅可比矩阵
    jacobi = self.cal_jacobi_matrix(xi, eta)
    # 接着计算逆矩阵
    inv_jacobi = self.__inv(matrix=jacobi)
    # 计算形函数对自然坐标的导数, 即[∂Ni/∂ε; ∂Ni/∂η], 其中下标i = 1,2,3,4
    partial_N1_natrue = np.array([-0.25 * (1 - eta), -0.25 * (1 - xi)]) # [∂N1/∂ε; ∂N1/∂η]
    partial_N2_natrue = np.array([ 0.25 * (1 - eta), -0.25 * (1 + xi)]) # [∂N2/∂ε; ∂N2/∂η]
    partial_N3_natrue = np.array([ 0.25 * (1 + eta),  0.25 * (1 + xi)]) # [∂N3/∂ε; ∂N3/∂η]
    partial_N4_natrue = np.array([-0.25 * (1 + eta),  0.25 * (1 - xi)]) # [∂N4/∂ε; ∂N4/∂η]
    # 计算形函数对直角坐标的导数, 即[∂Ni/∂x; ∂Ni/∂y], 其中下标i = 1,2,3,4
    partial_N1_cart = np.matmul(inv_jacobi, partial_N1_natrue)
    partial_N2_cart = np.matmul(inv_jacobi, partial_N2_natrue)
    partial_N3_cart = np.matmul(inv_jacobi, partial_N3_natrue)
    partial_N4_cart = np.matmul(inv_jacobi, partial_N4_natrue)
    # 乘以4个节点位移, 然后累加求和, 获得单元在2方向上的应变
    element_strain2 = partial_N1_cart[1]*self.node_u2[0] + partial_N2_cart[1]*self.node_u2[1] + \
                      partial_N3_cart[1]*self.node_u2[2] + partial_N4_cart[1]*self.node_u2[3]
    return element_strain2
  
  def cal_partial_u2_x1(self):
    '''
    此函数用于计算单元单元在2方向上的应变, 即∂u2/∂x1\n
    xi: 希腊字母ε, 自然坐标之一, 通常表示水平方向, 取值范围在-1到1之间\n
    eta:希腊字母η, 自然坐标之一, 通常表示竖直方向, 取值范围在-1到1之间
    '''
    # 高斯积分点的坐标位置
    W1 = [-1/sqrt(3), -1/sqrt(3), 1/sqrt(3), 1/sqrt(3)]
    W2 = [-1/sqrt(3), 1/sqrt(3), -1/sqrt(3), 1/sqrt(3)]
    Am = self.cal_Am()
    ε22 = 0.0   # ε11 = ∂u1/∂x1
    for ξ, η in zip(W1, W2):
      jacobi = self.cal_jacobi_matrix(xi=ξ, eta=η)
      ε2 = self.__cal_partial_u2_x1(xi=ξ, eta=η)
      ε22 += ε2 * np.linalg.det(jacobi)
    return ε22 / Am
  
  def g_func(self, r:float, r0:float, r1:float):
    '''
    辅助函数g(r) = (r - r0) / (r1 - r0), 其中r为极坐标下原子到原点的距离
    '''
    return (r - r0) / (r1 - r0)
  
  def __cal_partial_g_x1(self, xi:float, eta:float, r0:float, r1:float):
    '''
    此函数用于计算辅助函数g对x1的偏导, 即g1 = ∂Ni/∂x1, i = 1, 2, 3, 4\n
    xi: 希腊字母ε, 自然坐标之一, 通常表示水平方向, 取值范围在-1到1之间\n
    eta:希腊字母η, 自然坐标之一, 通常表示竖直方向, 取值范围在-1到1之间
    '''
    # 首先计算雅可比矩阵
    jacobi = self.cal_jacobi_matrix(xi, eta)
    # 接着计算逆矩阵
    inv_jacobi = self.__inv(matrix=jacobi)
    # 计算形函数对自然坐标的导数, 即[∂Ni/∂ε; ∂Ni/∂η], 其中下标i = 1,2,3,4
    partial_N1_natrue = np.array([-0.25 * (1 - eta), -0.25 * (1 - xi)]) # [∂N1/∂ε; ∂N1/∂η]
    partial_N2_natrue = np.array([ 0.25 * (1 - eta), -0.25 * (1 + xi)]) # [∂N2/∂ε; ∂N2/∂η]
    partial_N3_natrue = np.array([ 0.25 * (1 + eta),  0.25 * (1 + xi)]) # [∂N3/∂ε; ∂N3/∂η]
    partial_N4_natrue = np.array([-0.25 * (1 + eta),  0.25 * (1 - xi)]) # [∂N4/∂ε; ∂N4/∂η]
    # 计算形函数对直角坐标的导数, 即[∂Ni/∂x; ∂Ni/∂y], 其中下标i = 1,2,3,4
    partial_N1_cart = np.matmul(inv_jacobi, partial_N1_natrue)
    partial_N2_cart = np.matmul(inv_jacobi, partial_N2_natrue)
    partial_N3_cart = np.matmul(inv_jacobi, partial_N3_natrue)
    partial_N4_cart = np.matmul(inv_jacobi, partial_N4_natrue)
    # 乘以4个节点的g(r)值, 然后求和
    r_1, r_2, r_3, r_4 = np.linalg.norm(self.node_coordinate[0, :], ord=2), \
                         np.linalg.norm(self.node_coordinate[1, :], ord=2), \
                         np.linalg.norm(self.node_coordinate[2, :], ord=2), \
                         np.linalg.norm(self.node_coordinate[3, :], ord=2)
    partial_g_x1 = partial_N1_cart[0]*self.g_func(r_1, r0, r1) + partial_N2_cart[0]*self.g_func(r_2, r0, r1)\
                 + partial_N3_cart[0]*self.g_func(r_3, r0, r1) + partial_N4_cart[0]*self.g_func(r_4, r0, r1)
    
    return partial_g_x1
  
  def cal_partial_g_x1(self, r0:float, r1:float):
    # 高斯积分点的坐标位置
    W1 = [-1/sqrt(3), -1/sqrt(3), 1/sqrt(3), 1/sqrt(3)]
    W2 = [-1/sqrt(3), 1/sqrt(3), -1/sqrt(3), 1/sqrt(3)]
    Am = self.cal_Am()
    partial_g_x1 = 0.0    # g1 = ∂g/∂x1
    for ξ, η in zip(W1, W2):
      jacobi = self.cal_jacobi_matrix(xi=ξ, eta=η)
      g1 = self.__cal_partial_g_x1(xi=ξ, eta=η, r0=r0, r1=r1)
      partial_g_x1 += g1 * np.linalg.det(jacobi)
    return partial_g_x1 / Am

  def __cal_partial_g_x2(self, xi:float, eta:float, r0:float, r1:float):
    '''
    此函数用于计算辅助函数g对x2的偏导, 即g1 = ∂Ni/∂x2, i = 1, 2, 3, 4\n
    xi: 希腊字母ε, 自然坐标之一, 通常表示水平方向, 取值范围在-1到1之间\n
    eta:希腊字母η, 自然坐标之一, 通常表示竖直方向, 取值范围在-1到1之间
    '''
    # 首先计算雅可比矩阵
    jacobi = self.cal_jacobi_matrix(xi, eta)
    # 接着计算逆矩阵
    inv_jacobi = self.__inv(matrix=jacobi)
    # 计算形函数对自然坐标的导数, 即[∂Ni/∂ε; ∂Ni/∂η], 其中下标i = 1,2,3,4
    partial_N1_natrue = np.array([-0.25 * (1 - eta), -0.25 * (1 - xi)]) # [∂N1/∂ε; ∂N1/∂η]
    partial_N2_natrue = np.array([ 0.25 * (1 - eta), -0.25 * (1 + xi)]) # [∂N2/∂ε; ∂N2/∂η]
    partial_N3_natrue = np.array([ 0.25 * (1 + eta),  0.25 * (1 + xi)]) # [∂N3/∂ε; ∂N3/∂η]
    partial_N4_natrue = np.array([-0.25 * (1 + eta),  0.25 * (1 - xi)]) # [∂N4/∂ε; ∂N4/∂η]
    # 计算形函数对直角坐标的导数, 即[∂Ni/∂x; ∂Ni/∂y], 其中下标i = 1,2,3,4
    partial_N1_cart = np.matmul(inv_jacobi, partial_N1_natrue)
    partial_N2_cart = np.matmul(inv_jacobi, partial_N2_natrue)
    partial_N3_cart = np.matmul(inv_jacobi, partial_N3_natrue)
    partial_N4_cart = np.matmul(inv_jacobi, partial_N4_natrue)
    # 乘以4个节点的g(r)值, 然后求和
    r_1, r_2, r_3, r_4 = np.linalg.norm(self.node_coordinate[0, :], ord=2), \
                         np.linalg.norm(self.node_coordinate[1, :], ord=2), \
                         np.linalg.norm(self.node_coordinate[2, :], ord=2), \
                         np.linalg.norm(self.node_coordinate[3, :], ord=2)
    partial_g_x2 = partial_N1_cart[1]*self.g_func(r_1, r0, r1) + partial_N2_cart[1]*self.g_func(r_2, r0, r1)\
                 + partial_N3_cart[1]*self.g_func(r_3, r0, r1) + partial_N4_cart[1]*self.g_func(r_4, r0, r1)
    
    return partial_g_x2
  
  def cal_partial_g_x2(self, r0:float, r1:float):
    # 高斯积分点的坐标位置
    W1 = [-1/sqrt(3), -1/sqrt(3), 1/sqrt(3), 1/sqrt(3)]
    W2 = [-1/sqrt(3), 1/sqrt(3), -1/sqrt(3), 1/sqrt(3)]
    Am = self.cal_Am()
    partial_g_x2 = 0.0    # g2 = ∂g/∂x2
    for ξ, η in zip(W1, W2):
      jacobi = self.cal_jacobi_matrix(xi=ξ, eta=η)
      g2 = self.__cal_partial_g_x2(xi=ξ, eta=η, r0=r0, r1=r1)
      partial_g_x2 += g2 * np.linalg.det(jacobi)
    return partial_g_x2 / Am
  
  def cal_Am(self):
    '''
    此函数用于计算参考文献中的公式(11)中的Am
    '''
    Am = 0.0
    W1 = [-1/sqrt(3), -1/sqrt(3), 1/sqrt(3), 1/sqrt(3)]
    W2 = [-1/sqrt(3), 1/sqrt(3), -1/sqrt(3), 1/sqrt(3)]
    for w1, w2 in zip(W1, W2):
      Am += self.cal_jacobi_det(w1, w2)
    return Am